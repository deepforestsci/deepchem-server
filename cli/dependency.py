"""Dependency resolution for service startup."""
from collections import defaultdict

from rich.console import Console

from .config import Config
from .docker_manager import DockerManager
from .process_manager import LocalProcessManager


console = Console()


class DependencyResolver:
    """Resolves and manages service dependencies for startup ordering."""

    def __init__(self, config: Config, mode: str):
        """Initialize dependency resolver.
        
        Args:
            config: CLI configuration with service definitions.
            mode: Execution mode ('local' or 'docker').
        """
        self.config = config
        self.mode = mode
        self.docker_mgr = DockerManager()
        self.local_mgr = LocalProcessManager(config)

    def resolve_order(self, services: list[str] | None) -> tuple[list[str], list[str]]:
        """Resolve startup order based on dependencies.
        
        Args:
            services: Specific services to start, or None for all.
        
        Returns:
            Tuple of (dependencies to start, services in dependency order).
        """
        # If no services specified, start all
        if not services:
            services = list(self.config.services.keys())

        # Collect all required infrastructure dependencies
        required_deps: set[str] = set()
        for svc_name in services:
            svc = self.config.services.get(svc_name)
            if svc:
                for dep in svc.depends_on:
                    if dep in self.config.dependencies:
                        required_deps.add(dep)

        # Topological sort for services
        ordered_services = self._topo_sort(services)

        return list(required_deps), ordered_services

    def start_services(
        self,
        services: list[str] | None = None,
        num_workers: int = 2,
        skip_deps: bool = False,
    ) -> bool:
        """Start services in dependency order.
        
        Args:
            services: Specific services to start, or None for all.
            num_workers: Number of worker replicas.
            skip_deps: Skip starting infrastructure dependencies.
        
        Returns:
            True if all services started successfully.
        """
        deps, ordered_services = self.resolve_order(services)

        # In docker mode, do a single compose up (optionally scaling),
        # and build images to avoid stale-code mismatches.
        if self.mode == "docker":
            # Determine services to start
            services_to_start: list[str] = []
            if not skip_deps:
                services_to_start.extend(deps)
            services_to_start.extend(ordered_services)

            # Deduplicate while preserving order
            seen: set[str] = set()
            services_to_start = [s for s in services_to_start if not (s in seen or seen.add(s))]  # type: ignore

            scale: dict[str, int] | None = None
            if "worker" in services_to_start and num_workers > 1:
                scale = {"worker": num_workers}

            console.print("[bold]Starting services (Docker Compose)...[/]")
            ok = self.docker_mgr.start(
                services=services_to_start if services_to_start else None,
                scale=scale,
                no_deps=skip_deps,
            )
            if not ok:
                return False

            console.print()
            console.print("[bold green]✓ All services started successfully![/]")
            self._print_urls()
            return True

        # Check what's already running in Docker (if in local mode)
        docker_running: set[str] = set()
        if self.mode == "local":
            docker_running = set(self.docker_mgr.get_running_services())
            if docker_running:
                console.print(f"[yellow]⚠ Services already running in Docker: "
                              f"{', '.join(sorted(docker_running))}[/]")
                console.print("[yellow]  These will be skipped in local mode.[/]")
                console.print()

        # 1. Start dependencies first
        if not skip_deps and deps:
            console.print("[bold]Starting dependencies...[/]")
            for dep_name in deps:
                if dep_name in docker_running:
                    console.print(f"[yellow]  >> {dep_name}: Running in Docker, skipping[/]")
                    continue

                if self._is_running(dep_name):
                    console.print(f"[green]  ✓ {dep_name}: Already running[/]")
                    continue

                console.print(f"[blue]  ▶ Starting {dep_name}...[/]")
                if not self._start_dependency(dep_name):
                    console.print(f"[red]✗ Failed to start dependency: {dep_name}[/]")
                    return False
            console.print()

        # 2. Start services in order
        console.print("[bold]Starting services...[/]")
        for svc_name in ordered_services:
            if svc_name in docker_running:
                console.print(f"[yellow]  >> {svc_name}: Running in Docker, skipping[/]")
                continue

            # Check dependencies are running
            svc = self.config.services.get(svc_name)
            if svc:
                for dep in svc.depends_on:
                    if not self._is_running(dep) and dep not in docker_running:
                        console.print(f"[red]✗ Cannot start {svc_name}: "
                                      f"dependency '{dep}' is not running[/]")
                        console.print(f"[dim]  Try: deepchem-server-cli start {dep}[/]")
                        return False

            console.print(f"[blue]  ▶ Starting {svc_name}...[/]")

            if svc and svc.scalable and num_workers > 1:
                if not self._start_service(svc_name, replicas=num_workers):
                    return False
            else:
                if not self._start_service(svc_name):
                    return False

        console.print()
        console.print("[bold green]✓ All services started successfully![/]")
        self._print_urls()

        return True

    def _print_urls(self) -> None:
        """Print service URLs after successful startup."""
        console.print()
        console.print("[bold]Service URLs:[/]")

        for name, svc in self.config.services.items():
            if svc.port:
                console.print(f"  • {svc.name}: [cyan]http://localhost:{svc.port}[/]")

        for name, dep in self.config.dependencies.items():
            if dep.port:
                console.print(f"  • {dep.name}: [cyan]localhost:{dep.port}[/]")

    def _is_running(self, name: str) -> bool:
        """Check if a service/dependency is running.
        
        Checks both local processes and Docker containers.
        """
        if self.mode == "docker":
            return name in self.docker_mgr.get_running_services()
        else:
            # Check local first
            if self.local_mgr.is_running(name):
                return True
            # Also check Docker (for hybrid scenarios)
            return name in self.docker_mgr.get_running_services()

    def _start_dependency(self, name: str) -> bool:
        """Start an infrastructure dependency."""
        if self.mode == "docker":
            return self.docker_mgr.start([name])
        else:
            pid = self.local_mgr.start_dependency(name)
            return pid is not None

    def _start_service(self, name: str, replicas: int = 1) -> bool:
        """Start an application service."""
        if self.mode == "docker":
            scale = {name: replicas} if replicas > 1 else None
            return self.docker_mgr.start([name], scale=scale)
        else:
            pids = self.local_mgr.start_service(name, replicas=replicas)
            return len(pids) == replicas

    def _topo_sort(self, services: list[str]) -> list[str]:
        """Topological sort of services by dependencies.
        
        Only considers service-to-service dependencies (not infrastructure deps).
        """
        # Build dependency graph
        graph: dict[str, list[str]] = defaultdict(list)
        in_degree: dict[str, int] = defaultdict(int)

        # Initialize all requested services
        for svc_name in services:
            in_degree[svc_name] = in_degree.get(svc_name, 0)

        for svc_name in services:
            svc = self.config.services.get(svc_name)
            if svc:
                for dep in svc.depends_on:
                    # Only consider service-to-service deps, not infra deps
                    if dep in services and dep in self.config.services:
                        graph[dep].append(svc_name)
                        in_degree[svc_name] += 1

        # Kahn's algorithm
        queue = [s for s in services if in_degree[s] == 0]
        result = []

        while queue:
            node = queue.pop(0)
            result.append(node)

            for neighbor in graph[node]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)

        # If we couldn't process all nodes, there's a cycle
        if len(result) != len(services):
            # Return original order as fallback
            console.print("[yellow]⚠ Circular dependency detected, using original order[/]")
            return services

        return result
