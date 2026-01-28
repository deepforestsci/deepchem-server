"""Local process management with PID and log tracking."""
import os
import socket
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import httpx
import psutil
from rich.console import Console

from .config import Config, HealthCheck

console = Console()


@dataclass
class ServiceStatus:
    """Status of a local service."""
    name: str
    running: bool
    pid: int | None
    healthy: bool
    port: int | None


class LocalProcessManager:
    """Manages local service processes with PID and log tracking."""
    
    def __init__(self, config: Config):
        """Initialize the local process manager.
        
        Args:
            config: CLI configuration with service definitions.
        """
        self.config = config
        self.pid_dir = Path(config.settings.pid_dir)
        self.log_dir = Path(config.settings.log_dir)
        self._ensure_dirs()
    
    def _ensure_dirs(self) -> None:
        """Create required directories."""
        for dir_path in self.config.settings.create_dirs:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def _expand_env_value(self, value: str) -> str:
        """Expand environment variable syntax in value.
        
        Handles ${VAR:-default} syntax.
        """
        if not value.startswith("${"):
            return value
        
        # Parse ${VAR:-default} or ${VAR}
        inner = value[2:-1]  # Remove ${ and }
        
        if ":-" in inner:
            var_name, default = inner.split(":-", 1)
            return os.environ.get(var_name, default)
        else:
            return os.environ.get(inner, "")
    
    def _build_env(self, env_vars: dict[str, str]) -> dict[str, str]:
        """Build environment dict with expanded variables."""
        full_env = os.environ.copy()
        
        for key, value in env_vars.items():
            full_env[key] = self._expand_env_value(value)
        
        return full_env
    
    def start_dependency(self, name: str) -> int | None:
        """Start an infrastructure dependency (e.g., redis).
        
        Args:
            name: Name of the dependency to start.
        
        Returns:
            PID of the started process, or None if failed.
        """
        dep = self.config.dependencies.get(name)
        if not dep:
            console.print(f"[red]✗ Unknown dependency: {name}[/]")
            return None
        
        return self._start_process(
            name=name,
            command=dep.local.command,
            env=dep.local.env,
            health_check=dep.local.health_check,
        )
    
    def start_service(self, name: str, replicas: int = 1) -> list[int]:
        """Start a service with optional replicas.
        
        Args:
            name: Name of the service to start.
            replicas: Number of replicas to start (for scalable services).
        
        Returns:
            List of PIDs for started processes.
        """
        svc = self.config.services.get(name)
        if not svc:
            console.print(f"[red]✗ Unknown service: {name}[/]")
            return []
        
        pids = []
        for i in range(replicas):
            suffix = f".{i+1}" if replicas > 1 else ""
            pid = self._start_process(
                name=f"{name}{suffix}",
                command=svc.local.command,
                env=svc.local.env,
                health_check=svc.local.health_check,
            )
            if pid:
                pids.append(pid)
            
            # Small delay between replicas
            if replicas > 1 and i < replicas - 1:
                time.sleep(0.5)
        
        return pids
    
    def _start_process(
        self, 
        name: str, 
        command: list[str], 
        env: dict[str, str],
        health_check: HealthCheck,
    ) -> int | None:
        """Start a process with logging.
        
        Args:
            name: Process name (used for PID and log files).
            command: Command to execute.
            env: Environment variables.
            health_check: Health check configuration.
        
        Returns:
            PID of started process, or None if failed.
        """
        # Check if already running
        if self.is_running(name):
            console.print(f"[yellow]  ⚠ {name} is already running[/]")
            return self._get_pid(name)
        
        # Build environment
        full_env = self._build_env(env)
        
        # Open log file
        log_file = self.log_dir / f"{name}.log"
        
        try:
            log_fd = open(log_file, "a")
            
            # Write startup marker
            log_fd.write(f"\n{'='*60}\n")
            log_fd.write(f"Starting {name} at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_fd.write(f"Command: {' '.join(command)}\n")
            log_fd.write(f"{'='*60}\n")
            log_fd.flush()
            
            # Start process
            proc = subprocess.Popen(
                command,
                stdout=log_fd,
                stderr=subprocess.STDOUT,
                env=full_env,
                start_new_session=True,  # Detach from terminal
            )
            
            # Save PID
            pid_file = self.pid_dir / f"{name}.pid"
            pid_file.write_text(str(proc.pid))
            
            # Wait for health check
            if self._wait_for_health(name, health_check, proc.pid):
                console.print(f"[green]  ✓ {name} started (PID: {proc.pid})[/]")
                return proc.pid
            else:
                console.print(f"[red]  ✗ {name} failed health check[/]")
                self._kill_process(proc.pid)
                pid_file.unlink(missing_ok=True)
                return None
                
        except FileNotFoundError:
            console.print(f"[red]  ✗ Command not found: {command[0]}[/]")
            console.print(f"[dim]    Make sure {command[0]} is installed and in PATH[/]")
            return None
        except Exception as e:
            console.print(f"[red]  ✗ Failed to start {name}: {e}[/]")
            return None
    
    def _wait_for_health(
        self, 
        name: str, 
        health_check: HealthCheck, 
        pid: int,
        timeout: int | None = None,
    ) -> bool:
        """Wait for service to become healthy.
        
        Args:
            name: Service name.
            health_check: Health check configuration.
            pid: Process ID to verify is still running.
            timeout: Override timeout from health check config.
        
        Returns:
            True if healthy, False if timeout or failed.
        """
        check_timeout = timeout or health_check.timeout
        check_type = health_check.type
        
        # For process-only checks, just verify it started
        if check_type == "process":
            time.sleep(0.5)  # Brief pause to let process crash if it will
            return psutil.pid_exists(pid)
        
        # Poll for health
        start_time = time.time()
        while time.time() - start_time < check_timeout:
            # First check process is still alive
            if not psutil.pid_exists(pid):
                return False
            
            if check_type == "http":
                try:
                    resp = httpx.get(health_check.url, timeout=2)
                    if resp.status_code == 200:
                        return True
                except Exception:
                    pass
                    
            elif check_type == "command":
                try:
                    result = subprocess.run(
                        health_check.command,
                        capture_output=True,
                        text=True,
                        timeout=5,
                    )
                    expected = health_check.expected or ""
                    if expected in result.stdout:
                        return True
                except Exception:
                    pass
            
            time.sleep(0.5)
        
        return False
    
    def stop_services(self, services: list[str] | None = None) -> None:
        """Stop services (all if none specified).
        
        Args:
            services: List of service names to stop. None for all.
        """
        # Find all PID files
        pid_files = list(self.pid_dir.glob("*.pid"))
        
        if not pid_files:
            console.print("[yellow]No running services found[/]")
            return
        
        stopped_count = 0
        for pid_file in pid_files:
            name = pid_file.stem
            
            # Filter if specific services requested
            if services:
                base_name = name.split(".")[0]  # Handle worker.1, worker.2
                if base_name not in services:
                    continue
            
            try:
                pid = int(pid_file.read_text().strip())
                self._kill_process(pid)
                pid_file.unlink(missing_ok=True)
                console.print(f"[green]  ✓ Stopped {name} (PID: {pid})[/]")
                stopped_count += 1
            except ValueError:
                console.print(f"[yellow]  ⚠ Invalid PID file for {name}[/]")
                pid_file.unlink(missing_ok=True)
            except Exception as e:
                console.print(f"[yellow]  ⚠ Could not stop {name}: {e}[/]")
        
        if stopped_count == 0 and services:
            console.print(f"[yellow]No matching services found for: {', '.join(services)}[/]")
    
    def _kill_process(self, pid: int) -> None:
        """Gracefully kill a process.
        
        Sends SIGTERM first, then SIGKILL if process doesn't terminate.
        """
        try:
            proc = psutil.Process(pid)
            
            # Try graceful termination first
            proc.terminate()
            
            try:
                proc.wait(timeout=5)
            except psutil.TimeoutExpired:
                # Force kill if graceful shutdown failed
                proc.kill()
                proc.wait(timeout=2)
                
        except psutil.NoSuchProcess:
            pass  # Already dead
        except psutil.AccessDenied:
            console.print(f"[yellow]  ⚠ Permission denied killing PID {pid}[/]")
    
    def is_running(self, name: str) -> bool:
        """Check if a service is running.
        
        Args:
            name: Service name to check.
        
        Returns:
            True if running, False otherwise.
        """
        pid = self._get_pid(name)
        if pid is not None:
            return psutil.pid_exists(pid)

        dep = self.config.dependencies.get(name)
        if dep:
            if dep.local.health_check.type != "process":
                if self._check_health_quick(dep.local.health_check):
                    return True
            if dep.port and self._is_port_open(dep.port):
                return True
            return False

        svc = self.config.services.get(name)
        if svc:
            if svc.local.health_check.type != "process":
                if self._check_health_quick(svc.local.health_check):
                    return True
            if svc.port and self._is_port_open(svc.port):
                return True

        return False
    
    def _get_pid(self, name: str) -> int | None:
        """Get PID for a service from its PID file.
        
        Args:
            name: Service name.
        
        Returns:
            PID if found and valid, None otherwise.
        """
        pid_file = self.pid_dir / f"{name}.pid"
        if not pid_file.exists():
            return None
        
        try:
            return int(pid_file.read_text().strip())
        except (ValueError, OSError):
            return None
    
    def get_all_status(self) -> list[ServiceStatus]:
        """Get status of all services.
        
        Returns:
            List of ServiceStatus objects for each dependency and service.
        """
        statuses = []
        
        # Check dependencies
        for name, dep in self.config.dependencies.items():
            running = self.is_running(name)
            pid = self._get_pid(name) if running else None
            healthy = False
            
            if running:
                healthy = self._check_health_quick(dep.local.health_check)
            
            statuses.append(ServiceStatus(
                name=name,
                running=running,
                pid=pid,
                healthy=healthy,
                port=dep.port,
            ))
        
        # Check services (including replicas)
        for name, svc in self.config.services.items():
            # Check for replicas (worker.1, worker.2, etc.)
            replica_patterns = [f"{name}.*.pid", f"{name}.pid"]
            replica_pids = []
            
            for pattern in replica_patterns:
                replica_pids.extend(self.pid_dir.glob(pattern))
            
            # Deduplicate
            replica_pids = list(set(replica_pids))
            
            if not replica_pids:
                # Service not running
                statuses.append(ServiceStatus(
                    name=name,
                    running=False,
                    pid=None,
                    healthy=False,
                    port=svc.port,
                ))
            else:
                for pid_file in sorted(replica_pids):
                    replica_name = pid_file.stem
                    pid = self._get_pid(replica_name)
                    running = pid is not None and psutil.pid_exists(pid)
                    
                    healthy = False
                    if running:
                        healthy = self._check_health_quick(svc.local.health_check)
                    
                    statuses.append(ServiceStatus(
                        name=replica_name,
                        running=running,
                        pid=pid if running else None,
                        healthy=healthy,
                        port=svc.port,
                    ))
        
        return statuses
    
    def _check_health_quick(self, health_check: HealthCheck) -> bool:
        """Quick health check (non-blocking).
        
        Args:
            health_check: Health check configuration.
        
        Returns:
            True if healthy, False otherwise.
        """
        check_type = health_check.type
        
        if check_type == "http":
            try:
                resp = httpx.get(health_check.url, timeout=2)
                return resp.status_code == 200
            except Exception:
                return False
                
        elif check_type == "command":
            try:
                result = subprocess.run(
                    health_check.command,
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                expected = health_check.expected or ""
                return expected in result.stdout
            except Exception:
                return False
        
        # For process type, just return True (already checked running)
        return True

    def _is_port_open(self, port: int, host: str = "127.0.0.1") -> bool:
        """Check if a local TCP port is accepting connections."""
        try:
            with socket.create_connection((host, port), timeout=0.3):
                return True
        except OSError:
            return False
    
    def show_logs(
        self, 
        service: str, 
        follow: bool = False, 
        lines: int = 100,
    ) -> None:
        """Show logs for a service.
        
        Args:
            service: Service name to show logs for.
            follow: Whether to follow (tail -f) the logs.
            lines: Number of lines to show.
        """
        # Find matching log files
        log_patterns = [f"{service}*.log", f"{service}.log"]
        log_files = []
        
        for pattern in log_patterns:
            log_files.extend(self.log_dir.glob(pattern))
        
        # Deduplicate and sort
        log_files = sorted(set(log_files))
        
        if not log_files:
            console.print(f"[yellow]No logs found for {service}[/]")
            console.print(f"[dim]  Log directory: {self.log_dir}[/]")
            return
        
        for log_file in log_files:
            console.print(f"[blue]═══ {log_file.name} ═══[/]")
            
            if follow:
                # Use tail -f for following
                try:
                    subprocess.run(["tail", "-f", str(log_file)])
                except KeyboardInterrupt:
                    console.print("\n[dim]Stopped following logs[/]")
            else:
                try:
                    subprocess.run(["tail", f"-n{lines}", str(log_file)])
                except FileNotFoundError:
                    # Fallback if tail not available
                    content = log_file.read_text()
                    log_lines = content.split("\n")
                    for line in log_lines[-lines:]:
                        console.print(line)
