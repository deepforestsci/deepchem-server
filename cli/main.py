"""DeepChem Server CLI - Service Manager.

A unified CLI tool to manage all DeepChem Server services,
supporting both local development and Docker deployment modes.
"""
import os
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .config import load_config
from .dependency import DependencyResolver
from .docker_manager import DockerManager
from .process_manager import LocalProcessManager


app = typer.Typer(
    name="deepchem-server-cli",
    help="Manage DeepChem Server services (local and Docker modes)",
    add_completion=False,
    no_args_is_help=True,
)
console = Console()


def get_mode(mode_override: str | None = None) -> str:
    """Get execution mode from override, environment, or default to local.
    
    Priority:
    1. Explicit --mode flag
    2. DEEPCHEM_MODE environment variable
    3. Default: 'local'
    """
    if mode_override:
        return mode_override
    return os.environ.get("DEEPCHEM_MODE", "local")


def version_callback(value: bool) -> None:
    """Show version and exit."""
    if value:
        console.print(f"deepchem-server-cli version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            "-v",
            callback=version_callback,
            is_eager=True,
            help="Show version and exit",
        ),
    ] = False,
) -> None:
    """DeepChem Server CLI - Manage services locally or with Docker."""
    pass


@app.command()
def start(
    services: Annotated[
        Optional[list[str]],
        typer.Argument(help="Services to start (default: all). Options: datastore, gateway, worker"),
    ] = None,
    mode: Annotated[
        Optional[str],
        typer.Option(
            "--mode",
            "-m",
            help="Execution mode: local or docker (default: local, or DEEPCHEM_MODE env)",
        ),
    ] = None,
    workers: Annotated[
        int,
        typer.Option(
            "--workers",
            "-w",
            help="Number of worker replicas",
        ),
    ] = 2,
    skip_deps: Annotated[
        bool,
        typer.Option(
            "--skip-deps",
            help="Skip starting dependencies (redis)",
        ),
    ] = False,
) -> None:
    """Start DeepChem services.
    
    Examples:
    
        deepchem-server-cli start                    # Start all services locally
        
        deepchem-server-cli start gateway worker     # Start specific services
        
        deepchem-server-cli start --mode docker      # Start with Docker
        
        deepchem-server-cli start --workers 4        # Start with 4 workers
    """
    effective_mode = get_mode(mode)

    console.print()
    console.print(f"[bold]DeepChem Server - Starting services ({effective_mode} mode)[/]")
    console.print("=" * 50)
    console.print()

    try:
        config = load_config()
    except FileNotFoundError as e:
        console.print(f"[red]✗ Configuration error: {e}[/]")
        raise typer.Exit(1)

    resolver = DependencyResolver(config, effective_mode)

    success = resolver.start_services(
        services=services,
        num_workers=workers,
        skip_deps=skip_deps,
    )

    if not success:
        raise typer.Exit(1)


@app.command()
def stop(
    services: Annotated[
        Optional[list[str]],
        typer.Argument(help="Services to stop (default: all)"),
    ] = None,
    mode: Annotated[
        Optional[str],
        typer.Option(
            "--mode",
            "-m",
            help="Execution mode: local or docker",
        ),
    ] = None,
) -> None:
    """Stop running services.
    
    Examples:
    
        deepchem-server-cli stop                     # Stop all local services
        
        deepchem-server-cli stop gateway             # Stop specific service
        
        deepchem-server-cli stop --mode docker       # Stop Docker services
    """
    effective_mode = get_mode(mode)

    console.print()
    console.print(f"[bold]DeepChem Server - Stopping services ({effective_mode} mode)[/]")
    console.print("=" * 50)
    console.print()

    try:
        config = load_config()
    except FileNotFoundError as e:
        console.print(f"[red]✗ Configuration error: {e}[/]")
        raise typer.Exit(1)

    if effective_mode == "docker":
        docker_mgr = DockerManager()
        docker_mgr.stop(services)
    else:
        local_mgr = LocalProcessManager(config)
        local_mgr.stop_services(services)

    console.print()
    console.print("[green]✓ Stop complete[/]")


@app.command()
def restart(
    services: Annotated[
        Optional[list[str]],
        typer.Argument(help="Services to restart (default: all)"),
    ] = None,
    mode: Annotated[
        Optional[str],
        typer.Option("--mode", "-m", help="Execution mode: local or docker"),
    ] = None,
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", help="Number of worker replicas"),
    ] = 2,
) -> None:
    """Restart services (stop then start).
    
    Examples:
    
        deepchem-server-cli restart                  # Restart all services
        
        deepchem-server-cli restart worker           # Restart workers only
    """
    effective_mode = get_mode(mode)

    console.print()
    console.print(f"[bold]DeepChem Server - Restarting services ({effective_mode} mode)[/]")
    console.print("=" * 50)
    console.print()

    try:
        config = load_config()
    except FileNotFoundError as e:
        console.print(f"[red]✗ Configuration error: {e}[/]")
        raise typer.Exit(1)

    # Stop first
    console.print("[bold]Stopping services...[/]")
    if effective_mode == "docker":
        docker_mgr = DockerManager()
        docker_mgr.stop(services)
    else:
        local_mgr = LocalProcessManager(config)
        local_mgr.stop_services(services)

    console.print()

    # Then start
    console.print("[bold]Starting services...[/]")
    resolver = DependencyResolver(config, effective_mode)

    success = resolver.start_services(
        services=services,
        num_workers=workers,
        skip_deps=False,
    )

    if not success:
        raise typer.Exit(1)


@app.command()
def status(
    mode: Annotated[
        Optional[str],
        typer.Option(
            "--mode",
            "-m",
            help="Check status for mode: local or docker",
        ),
    ] = None,
) -> None:
    """Show status of all services.
    
    Examples:
    
        deepchem-server-cli status                   # Show local service status
        
        deepchem-server-cli status --mode docker     # Show Docker service status
    """
    effective_mode = get_mode(mode)

    try:
        config = load_config()
    except FileNotFoundError as e:
        console.print(f"[red]✗ Configuration error: {e}[/]")
        raise typer.Exit(1)

    console.print()

    # Build status table
    table = Table(
        title=f"DeepChem Server Status ({effective_mode.title()} Mode)",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Service", style="white")
    table.add_column("Status", justify="center")
    table.add_column("PID/Container", justify="right")
    table.add_column("Health", justify="center")
    table.add_column("Port", justify="right")

    if effective_mode == "docker":
        docker_mgr = DockerManager()
        statuses = docker_mgr.get_status()

        if not statuses:
            console.print("[yellow]No Docker services found[/]")
            console.print("[dim]  Run 'deepchem-server-cli start --mode docker' to start services[/]")
            return

        for svc in statuses:
            status_str = "[green]● Running[/]" if svc.running else "[red]○ Stopped[/]"
            health_str = "[green]✓ Healthy[/]" if svc.healthy else (
                "[yellow]⚠ Unknown[/]" if svc.running else "[dim]-[/]")

            table.add_row(
                svc.name,
                status_str,
                svc.container_id or "-",
                health_str,
                str(svc.port) if svc.port else "-",
            )
    else:
        local_mgr = LocalProcessManager(config)
        statuses = local_mgr.get_all_status()  # type: ignore

        if not any(s.running for s in statuses):
            console.print("[yellow]No local services running[/]")
            console.print("[dim]  Run 'deepchem-server-cli start' to start services[/]")
            return

        for svc in statuses:
            status_str = "[green]● Running[/]" if svc.running else "[red]○ Stopped[/]"
            health_str = "[green]✓ Healthy[/]" if svc.healthy else (
                "[yellow]⚠ Unhealthy[/]" if svc.running else "[dim]-[/]")

            table.add_row(
                svc.name,
                status_str,
                str(svc.pid) if svc.pid else "-",  # type: ignore
                health_str,
                str(svc.port) if svc.port else "-",
            )

    console.print(table)
    console.print()


@app.command()
def logs(
    service: Annotated[
        str,
        typer.Argument(help="Service to show logs for (e.g., gateway, worker, datastore)"),
    ],
    follow: Annotated[
        bool,
        typer.Option("--follow", "-f", help="Follow log output (like tail -f)"),
    ] = False,
    lines: Annotated[
        int,
        typer.Option("--lines", "-n", help="Number of lines to show"),
    ] = 100,
    mode: Annotated[
        Optional[str],
        typer.Option("--mode", "-m", help="Execution mode: local or docker"),
    ] = None,
) -> None:
    """Show service logs.
    
    Examples:
    
        deepchem-server-cli logs gateway             # Show gateway logs
        
        deepchem-server-cli logs worker --follow     # Follow worker logs
        
        deepchem-server-cli logs datastore -n 50     # Show last 50 lines
    """
    effective_mode = get_mode(mode)

    try:
        config = load_config()
    except FileNotFoundError as e:
        console.print(f"[red]✗ Configuration error: {e}[/]")
        raise typer.Exit(1)

    console.print()

    if effective_mode == "docker":
        docker_mgr = DockerManager()
        docker_mgr.logs(service, follow=follow, lines=lines)
    else:
        local_mgr = LocalProcessManager(config)
        local_mgr.show_logs(service, follow=follow, lines=lines)


if __name__ == "__main__":
    app()
