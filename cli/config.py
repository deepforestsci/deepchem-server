"""Configuration loader for services.yaml."""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class HealthCheck:
    """Health check configuration."""
    type: str  # "http", "command", or "process"
    url: str | None = None
    command: list[str] | None = None
    expected: str | None = None
    timeout: int = 30


@dataclass
class LocalConfig:
    """Local execution configuration."""
    command: list[str]
    env: dict[str, str] = field(default_factory=dict)
    health_check: HealthCheck = field(default_factory=lambda: HealthCheck(type="process"))


@dataclass
class DockerConfig:
    """Docker execution configuration."""
    service: str


@dataclass
class DependencyConfig:
    """Infrastructure dependency configuration (e.g., Redis)."""
    name: str
    description: str
    local: LocalConfig
    docker: DockerConfig
    port: int | None = None


@dataclass
class ServiceConfig:
    """Application service configuration."""
    name: str
    description: str
    depends_on: list[str]
    local: LocalConfig
    docker: DockerConfig
    port: int | None = None
    scalable: bool = False
    default_replicas: int = 1


@dataclass
class Settings:
    """Global CLI settings."""
    data_dir: str
    pid_dir: str
    log_dir: str
    create_dirs: list[str]


@dataclass
class Config:
    """Complete CLI configuration."""
    dependencies: dict[str, DependencyConfig]
    services: dict[str, ServiceConfig]
    settings: Settings


def _parse_health_check(data: dict[str, Any]) -> HealthCheck:
    """Parse health check configuration from YAML data."""
    if not data:
        return HealthCheck(type="process")
    
    return HealthCheck(
        type=data.get("type", "process"),
        url=data.get("url"),
        command=data.get("command"),
        expected=data.get("expected"),
        timeout=data.get("timeout", 30),
    )


def _parse_local_config(data: dict[str, Any]) -> LocalConfig:
    """Parse local execution configuration from YAML data."""
    return LocalConfig(
        command=data.get("command", []),
        env=data.get("env", {}),
        health_check=_parse_health_check(data.get("health_check", {})),
    )


def load_config(config_path: Path | None = None) -> Config:
    """Load configuration from services.yaml.
    
    Args:
        config_path: Optional path to config file. Defaults to services.yaml
                    in the CLI package directory.
    
    Returns:
        Parsed Config object with dependencies, services, and settings.
    """
    if config_path is None:
        # Look for config in cli package directory
        config_path = Path(__file__).parent / "services.yaml"
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path) as f:
        raw = yaml.safe_load(f)
    
    # Parse dependencies
    dependencies: dict[str, DependencyConfig] = {}
    for name, dep_data in raw.get("dependencies", {}).items():
        dependencies[name] = DependencyConfig(
            name=dep_data.get("name", name),
            description=dep_data.get("description", ""),
            local=_parse_local_config(dep_data.get("local", {})),
            docker=DockerConfig(service=dep_data.get("docker", {}).get("service", name)),
            port=dep_data.get("port"),
        )
    
    # Parse services
    services: dict[str, ServiceConfig] = {}
    for name, svc_data in raw.get("services", {}).items():
        services[name] = ServiceConfig(
            name=svc_data.get("name", name),
            description=svc_data.get("description", ""),
            depends_on=svc_data.get("depends_on", []),
            local=_parse_local_config(svc_data.get("local", {})),
            docker=DockerConfig(service=svc_data.get("docker", {}).get("service", name)),
            port=svc_data.get("port"),
            scalable=svc_data.get("scalable", False),
            default_replicas=svc_data.get("default_replicas", 1),
        )
    
    # Parse settings
    settings_data = raw.get("settings", {})
    settings = Settings(
        data_dir=settings_data.get("data_dir", "./data"),
        pid_dir=settings_data.get("pid_dir", ".deepchem/pids"),
        log_dir=settings_data.get("log_dir", ".deepchem/logs"),
        create_dirs=settings_data.get("create_dirs", []),
    )
    
    return Config(
        dependencies=dependencies,
        services=services,
        settings=settings,
    )
