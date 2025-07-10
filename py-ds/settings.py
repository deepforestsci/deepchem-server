"""
Settings module for managing profile and project configurations.
"""

import json
import os
from typing import Optional, Dict, Any
from pathlib import Path


class Settings:
    """
    Settings class for managing profile and project configurations.

    This class handles storing and loading user settings including profile
    and project information, and persists them to a settings.json file.
    """

    def __init__(self, settings_file: str = "settings.json"):
        """
        Initialize Settings instance.

        Args:
            settings_file: Path to the settings JSON file (default: "settings.json")
        """
        self.settings_file = Path(settings_file)
        self.profile: Optional[str] = None
        self.project: Optional[str] = None
        self.base_url: str = "http://localhost:8000"
        self._additional_settings: Dict[str, Any] = {}

        self.load()

    def set_profile(self, profile_name: str) -> None:
        """
        Set the profile name.

        Args:
            profile_name: Name of the profile
        """
        self.profile = profile_name
        self.save()

    def set_project(self, project_name: str) -> None:
        """
        Set the project name.

        Args:
            project_name: Name of the project
        """
        self.project = project_name
        self.save()

    def set_base_url(self, url: str) -> None:
        """
        Set the base URL for the API.

        Args:
            url: Base URL for the DeepChem server
        """
        self.base_url = url.rstrip("/")
        self.save()

    def get_profile(self) -> Optional[str]:
        """
        Get the current profile name.

        Returns:
            Current profile name or None if not set
        """
        return self.profile

    def get_project(self) -> Optional[str]:
        """
        Get the current project name.

        Returns:
            Current project name or None if not set
        """
        return self.project

    def get_base_url(self) -> str:
        """
        Get the base URL for the API.

        Returns:
            Base URL for the DeepChem server
        """
        return self.base_url

    def set_setting(self, key: str, value: Any) -> None:
        """
        Set a custom setting.

        Args:
            key: Setting key
            value: Setting value
        """
        self._additional_settings[key] = value
        self.save()

    def get_setting(self, key: str, default: Any = None) -> Any:
        """
        Get a custom setting.

        Args:
            key: Setting key
            default: Default value if key not found

        Returns:
            Setting value or default
        """
        return self._additional_settings.get(key, default)

    def save(self) -> None:
        """
        Save current settings to the JSON file.
        """
        settings_data = {
            "profile": self.profile,
            "project": self.project,
            "base_url": self.base_url,
            "additional_settings": self._additional_settings,
        }

        try:
            with open(self.settings_file, "w") as f:
                json.dump(settings_data, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save settings to {self.settings_file}: {e}")

    def load(self) -> None:
        """
        Load settings from the JSON file if it exists.
        """
        if not self.settings_file.exists():
            return

        try:
            with open(self.settings_file, "r") as f:
                settings_data = json.load(f)

            self.profile = settings_data.get("profile")
            self.project = settings_data.get("project")
            self.base_url = settings_data.get("base_url", "http://localhost:8000")
            self._additional_settings = settings_data.get("additional_settings", {})

        except Exception as e:
            print(f"Warning: Could not load settings from {self.settings_file}: {e}")

    def reset(self) -> None:
        """
        Reset all settings to default values and remove the settings file.
        """
        self.profile = None
        self.project = None
        self.base_url = "http://localhost:8000"
        self._additional_settings = {}

        if self.settings_file.exists():
            try:
                os.remove(self.settings_file)
            except Exception as e:
                print(f"Warning: Could not remove settings file {self.settings_file}: {e}")

    def is_configured(self) -> bool:
        """
        Check if both profile and project are configured.

        Returns:
            True if both profile and project are set, False otherwise
        """
        return self.profile is not None and self.project is not None

    def __str__(self) -> str:
        """
        String representation of current settings.

        Returns:
            String representation of settings
        """
        return f"Settings(profile='{self.profile}', project='{self.project}', base_url='{self.base_url}')"

    def __repr__(self) -> str:
        """
        Detailed representation of current settings.

        Returns:
            Detailed representation of settings
        """
        return (
            f"Settings(profile='{self.profile}', project='{self.project}', "
            f"base_url='{self.base_url}', settings_file='{self.settings_file}')"
        )
