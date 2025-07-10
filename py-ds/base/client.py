"""
Base client class with common functionality for all DeepChem API clients.
"""

import requests
from typing import Optional, Dict, Any
from ..settings import Settings


class BaseClient:
    """
    Base client class containing common functionality for all DeepChem API clients.

    This class provides shared methods for HTTP requests, configuration validation,
    and other common operations that are used by specific client implementations.
    """

    def __init__(self, settings: Optional[Settings] = None, base_url: Optional[str] = None):
        """
        Initialize BaseClient.

        Args:
            settings: Settings instance for configuration
            base_url: Base URL for the API (overrides settings if provided)
        """
        if settings is None:
            settings = Settings()

        self.settings = settings
        self.base_url = base_url or settings.get_base_url()
        self.session = requests.Session()

    def _check_configuration(self) -> None:
        """
        Check if profile and project are configured.

        Raises:
            ValueError: If profile or project is not configured
        """
        if not self.settings.is_configured():
            missing = []
            if not self.settings.get_profile():
                missing.append("profile")
            if not self.settings.get_project():
                missing.append("project")
            raise ValueError(
                f"Missing required settings: {', '.join(missing)}. "
                f"Please configure using settings.set_profile() and settings.set_project()"
            )

    def _make_request(self, method: str, endpoint: str, **kwargs) -> Dict[str, Any]:
        """
        Make HTTP request to the API.

        Args:
            method: HTTP method ('GET', 'POST', etc.)
            endpoint: API endpoint (without base URL)
            **kwargs: Additional arguments for requests

        Returns:
            JSON response as dictionary

        Raises:
            requests.exceptions.RequestException: If request fails
        """
        url = f"{self.base_url}{endpoint}"

        try:
            response = self.session.request(method, url, **kwargs)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            raise requests.exceptions.RequestException(f"API request failed: {e}")

    def _get_profile_and_project(
        self, profile_name: Optional[str] = None, project_name: Optional[str] = None
    ) -> tuple[str, str]:
        """
        Get profile and project names, either from parameters or settings.

        Args:
            profile_name: Profile name (uses settings if not provided)
            project_name: Project name (uses settings if not provided)

        Returns:
            Tuple of (profile, project) names

        Raises:
            ValueError: If required settings are missing
        """
        profile = profile_name or self.settings.get_profile()
        project = project_name or self.settings.get_project()

        if not profile or not project:
            self._check_configuration()

        # After _check_configuration(), we know both values are not None
        assert profile is not None and project is not None
        return profile, project

    def healthcheck(self) -> Dict[str, Any]:
        """
        Check server health status.

        Returns:
            Health status response

        Raises:
            requests.exceptions.RequestException: If API request fails
        """
        return self._make_request("GET", "/healthcheck")

    def get_settings(self) -> Settings:
        """
        Get the current settings instance.

        Returns:
            Settings instance
        """
        return self.settings

    def get_base_url(self) -> str:
        """
        Get the current base URL.

        Returns:
            Base URL for the API
        """
        return self.base_url

    def set_base_url(self, url: str) -> None:
        """
        Set a new base URL for this client instance.

        Args:
            url: New base URL for the API
        """
        self.base_url = url.rstrip("/")

    def close(self) -> None:
        """
        Close the HTTP session.
        """
        self.session.close()
