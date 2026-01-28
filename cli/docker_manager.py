"""Docker Compose management."""
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rich.console import Console

console = Console()


@dataclass
class ContainerStatus:
    """Status of a Docker container."""
    name: str
    running: bool
    container_id: str | None
    healthy: bool
    port: int | None


class DockerManager:
    """Manages services via docker compose."""
    
    def __init__(self, compose_file: Path | None = None):
        """Initialize Docker manager.
        
        Args:
            compose_file: Path to docker-compose.yml. Defaults to current directory.
        """
        self.compose_file = compose_file or Path("docker-compose.yml")
        self.compose_cmd = self._detect_compose_command()
    
    def _detect_compose_command(self) -> list[str]:
        """Detect whether to use 'docker-compose' or 'docker compose'."""
        try:
            subprocess.run(
                ["docker", "compose", "version"],
                capture_output=True,
                check=True,
            )
            return ["docker", "compose"]
        except (subprocess.CalledProcessError, FileNotFoundError):
            try:
                subprocess.run(
                    ["docker-compose", "version"],
                    capture_output=True,
                    check=True,
                )
                return ["docker-compose"]
            except (subprocess.CalledProcessError, FileNotFoundError):
                console.print("[yellow]⚠ Docker compose not found[/]")
                return ["docker", "compose"]  # Default, will fail gracefully
    
    def is_available(self) -> bool:
        """Check if Docker is available and running."""
        try:
            result = subprocess.run(
                ["docker", "info"],
                capture_output=True,
                timeout=5,
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def get_running_services(self) -> list[str]:
        """Get list of services currently running in Docker.
        
        Returns:
            List of service names that are currently running.
        """
        if not self.is_available():
            return []
        
        try:
            result = subprocess.run(
                [*self.compose_cmd, "ps", "--format", "json"],
                capture_output=True,
                text=True,
                cwd=self.compose_file.parent if self.compose_file.parent.exists() else None,
            )
            if result.returncode != 0:
                return []
            
            services = []
            output = result.stdout.strip()
            
            # Handle both single JSON object and newline-separated JSON objects
            if not output:
                return []
            
            for line in output.split("\n"):
                if line.strip():
                    try:
                        container = json.loads(line)
                        state = container.get("State", "").lower()
                        if state == "running":
                            name = container.get("Service", "")
                            if name and name not in services:
                                services.append(name)
                    except json.JSONDecodeError:
                        continue
            
            return services
        except Exception:
            return []
    
    def start(
        self, 
        services: list[str] | None = None, 
        scale: dict[str, int] | None = None,
        build: bool = False,
        no_deps: bool = False,
    ) -> bool:
        """Start docker services.
        
        Args:
            services: Specific services to start. None for all.
            scale: Dict of service name to replica count for scaling.
            build: Whether to build images before starting.
            no_deps: Whether to skip starting dependencies.
        
        Returns:
            True if successful, False otherwise.
        """
        cmd = [*self.compose_cmd, "up", "-d"]
        
        if build:
            cmd.append("--build")

        if no_deps:
            cmd.append("--no-deps")
        
        if scale:
            for svc, count in scale.items():
                cmd.extend(["--scale", f"{svc}={count}"])
        
        if services:
            cmd.extend(services)
        
        console.print(f"[blue]Running: {' '.join(cmd)}[/]")
        result = subprocess.run(cmd)
        
        if result.returncode == 0:
            console.print("[green]✓ Docker services started[/]")
            return True
        else:
            console.print("[red]✗ Failed to start Docker services[/]")
            return False
    
    def stop(self, services: list[str] | None = None) -> bool:
        """Stop docker services.
        
        Args:
            services: Specific services to stop. None for all (uses 'down').
        
        Returns:
            True if successful, False otherwise.
        """
        if services:
            cmd = [*self.compose_cmd, "stop", *services]
        else:
            cmd = [*self.compose_cmd, "down"]
        
        console.print(f"[blue]Running: {' '.join(cmd)}[/]")
        result = subprocess.run(cmd)
        
        if result.returncode == 0:
            console.print("[green]✓ Docker services stopped[/]")
            return True
        else:
            console.print("[red]✗ Failed to stop Docker services[/]")
            return False
    
    def restart(self, services: list[str] | None = None) -> bool:
        """Restart docker services.
        
        Args:
            services: Specific services to restart. None for all.
        
        Returns:
            True if successful, False otherwise.
        """
        cmd = [*self.compose_cmd, "restart"]
        
        if services:
            cmd.extend(services)
        
        console.print(f"[blue]Running: {' '.join(cmd)}[/]")
        result = subprocess.run(cmd)
        
        return result.returncode == 0
    
    def get_status(self) -> list[ContainerStatus]:
        """Get status of all docker services.
        
        Returns:
            List of ContainerStatus objects for each service.
        """
        if not self.is_available():
            return []
        
        try:
            result = subprocess.run(
                [*self.compose_cmd, "ps", "--format", "json"],
                capture_output=True,
                text=True,
            )
            
            if result.returncode != 0:
                return []
            
            statuses = []
            output = result.stdout.strip()
            
            if not output:
                return []
            
            for line in output.split("\n"):
                if line.strip():
                    try:
                        container = json.loads(line)
                        state = container.get("State", "").lower()
                        statuses.append(ContainerStatus(
                            name=container.get("Service", "unknown"),
                            running=state == "running",
                            container_id=container.get("ID", "")[:12] if container.get("ID") else None,
                            healthy=container.get("Health", "").lower() == "healthy",
                            port=self._extract_port(container.get("Ports", "")),
                        ))
                    except json.JSONDecodeError:
                        continue
            
            return statuses
        except Exception:
            return []
    
    def _extract_port(self, ports_str: str) -> int | None:
        """Extract first exposed port from ports string.
        
        Format examples:
            "0.0.0.0:8000->8000/tcp"
            "8000/tcp"
        """
        if not ports_str:
            return None
        
        if "->" in ports_str:
            try:
                # Format: "0.0.0.0:8000->8000/tcp"
                host_port = ports_str.split(":")[1].split("->")[0]
                return int(host_port)
            except (IndexError, ValueError):
                pass
        
        return None
    
    def logs(
        self, 
        service: str, 
        follow: bool = False, 
        lines: int = 100,
    ) -> None:
        """Show logs for a docker service.
        
        Args:
            service: Service name to show logs for.
            follow: Whether to follow log output.
            lines: Number of lines to show.
        """
        cmd = [*self.compose_cmd, "logs"]
        
        if follow:
            cmd.append("--follow")
        
        cmd.extend(["--tail", str(lines), service])
        
        console.print(f"[blue]═══ Docker logs: {service} ═══[/]")
        subprocess.run(cmd)
