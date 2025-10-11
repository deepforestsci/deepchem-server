"""
Unit tests for Settings class.
"""

import json
import os
from pathlib import Path
import tempfile
from typing import Any
from unittest.mock import mock_open, patch

from pyds.settings import Settings


class TestSettings:
    """Unit tests for Settings class."""

    def test_init_default(self, temp_settings_file: str) -> None:
        """Test Settings initialization with defaults."""
        settings = Settings(settings_file=temp_settings_file)

        assert settings.profile is None
        assert settings.project is None
        assert settings.base_url == "http://localhost:8000"
        assert settings._additional_settings == {}

    def test_init_with_parameters(self, temp_settings_file: str) -> None:
        """Test Settings initialization with parameters."""
        additional = {"key": "value"}
        settings = Settings(
            settings_file=temp_settings_file,
            profile="test_profile",
            project="test_project",
            base_url="http://custom:9000",
            additional_settings=additional,
        )

        assert settings.profile == "test_profile"
        assert settings.project == "test_project"
        assert settings.base_url == "http://custom:9000"
        assert settings._additional_settings == additional

    def test_set_profile(self, test_settings: Settings) -> None:
        """Test set_profile method."""
        test_settings.set_profile("new_profile")

        assert test_settings.profile == "new_profile"

    def test_set_project(self, test_settings: Settings) -> None:
        """Test set_project method."""
        test_settings.set_project("new_project")

        assert test_settings.project == "new_project"

    def test_set_base_url(self, test_settings: Settings) -> None:
        """Test set_base_url method."""
        test_settings.set_base_url("http://newserver:8080/")

        # Should strip trailing slash
        assert test_settings.base_url == "http://newserver:8080"

    def test_get_profile(self, test_settings: Settings) -> None:
        """Test get_profile method."""
        assert test_settings.get_profile() == "test_profile"

    def test_get_project(self, test_settings: Settings) -> None:
        """Test get_project method."""
        assert test_settings.get_project() == "test_project"

    def test_get_base_url(self, test_settings: Settings) -> None:
        """Test get_base_url method."""
        assert test_settings.get_base_url() == "http://localhost:8000"

    def test_set_get_setting(self, test_settings: Settings) -> None:
        """Test set_setting and get_setting methods."""
        test_settings.set_setting("custom_key", "custom_value")

        assert test_settings.get_setting("custom_key") == "custom_value"
        assert test_settings.get_setting("nonexistent", "default") == "default"

    def test_touch_creates_file(self, temp_settings_file: str) -> None:
        """Test touch method creates settings file."""
        # Remove the file first
        if os.path.exists(temp_settings_file):
            os.remove(temp_settings_file)

        Settings(settings_file=temp_settings_file)

        assert os.path.exists(temp_settings_file)

        # Check file content
        with open(temp_settings_file, "r") as f:
            data = json.load(f)

        expected: dict[str, Any] = {"profile": None, "project": None, "base_url": None, "additional_settings": {}}
        assert data == expected

    def test_save_and_load(self, temp_settings_file: str) -> None:
        """Test save and load methods."""
        settings = Settings(
            settings_file=temp_settings_file,
            profile="saved_profile",
            project="saved_project",
            base_url="http://saved:8000",
        )
        settings.set_setting("saved_key", "saved_value")

        # Create new instance to test loading
        new_settings = Settings(settings_file=temp_settings_file)

        assert new_settings.profile == "saved_profile"
        assert new_settings.project == "saved_project"
        assert new_settings.base_url == "http://saved:8000"
        assert new_settings.get_setting("saved_key") == "saved_value"

    def test_load_nonexistent_file(self, temp_settings_file: str) -> None:
        """Test load method with nonexistent file."""
        os.remove(temp_settings_file)

        settings = Settings()
        settings.settings_file = Path(temp_settings_file)
        settings.load()

        # Should not raise exception and maintain defaults
        assert settings.profile is None
        assert settings.project is None

    def test_load_invalid_json(self, temp_settings_file: str) -> None:
        """Test load method with invalid JSON."""
        with open(temp_settings_file, "w") as f:
            f.write("invalid json content")

        with patch("builtins.print") as mock_print:
            settings = Settings(settings_file=temp_settings_file)
            mock_print.assert_called()

        # Should maintain defaults
        assert settings.profile is None
        assert settings.project is None

    def test_save_permission_error(self, test_settings: Settings) -> None:
        """Test save method with permission error."""
        with patch("builtins.open", mock_open()) as mock_file:
            mock_file.side_effect = PermissionError("Permission denied")

            with patch("builtins.print") as mock_print:
                test_settings.save()
                mock_print.assert_called()

    def test_reset(self, test_settings: Settings) -> None:
        """Test reset method."""
        settings_file = test_settings.settings_file

        test_settings.reset()

        assert test_settings.profile is None
        assert test_settings.project is None
        assert test_settings.base_url == "http://localhost:8000"
        assert test_settings._additional_settings == {}
        assert not settings_file.exists()

    def test_reset_file_removal_error(self, test_settings: Settings) -> None:
        """Test reset method with file removal error."""
        with patch("os.remove", side_effect=PermissionError("Permission denied")):
            with patch("builtins.print") as mock_print:
                test_settings.reset()
                mock_print.assert_called()

    def test_is_configured_true(self, test_settings: Settings) -> None:
        """Test is_configured returns True when both profile and project are set."""
        assert test_settings.is_configured() is True

    def test_is_configured_false_no_profile(self, temp_settings_file: str) -> None:
        """Test is_configured returns False when profile is missing."""
        settings = Settings(settings_file=temp_settings_file, project="test_project")
        assert settings.is_configured() is False

    def test_is_configured_false_no_project(self, temp_settings_file: str) -> None:
        """Test is_configured returns False when project is missing."""
        settings = Settings(settings_file=temp_settings_file, profile="test_profile")
        assert settings.is_configured() is False

    def test_is_configured_false_both_missing(self, temp_settings_file: str) -> None:
        """Test is_configured returns False when both are missing."""
        settings = Settings(settings_file=temp_settings_file)
        assert settings.is_configured() is False

    def test_str_representation(self, test_settings: Settings) -> None:
        """Test string representation."""
        expected = "Settings(profile='test_profile', project='test_project', base_url='http://localhost:8000')"
        assert str(test_settings) == expected

    def test_repr_representation(self, test_settings: Settings) -> None:
        """Test detailed representation."""
        settings_file = test_settings.settings_file
        expected = (f"Settings(profile='test_profile', project='test_project', "
                    f"base_url='http://localhost:8000', settings_file='{settings_file}')")
        assert repr(test_settings) == expected

    def test_init_file_creation_when_not_exists(self) -> None:
        """Test Settings creates file when it doesn't exist."""
        with tempfile.NamedTemporaryFile(delete=True) as f:
            temp_file = f.name

        # File should not exist now
        assert not os.path.exists(temp_file)

        Settings(settings_file=temp_file)

        # File should be created
        assert os.path.exists(temp_file)

        # Clean up
        os.remove(temp_file)

    def test_load_with_missing_keys(self, temp_settings_file: str) -> None:
        """Test load method handles missing keys gracefully."""
        # Create file with partial data
        partial_data = {"profile": "partial_profile"}
        with open(temp_settings_file, "w") as f:
            json.dump(partial_data, f)

        settings = Settings(settings_file=temp_settings_file)

        assert settings.profile == "partial_profile"
        assert settings.project is None
        assert settings.base_url == "http://localhost:8000"
        assert settings._additional_settings == {}

    def test_base_url_default_when_none_in_file(self, temp_settings_file: str) -> None:
        """Test base_url defaults to localhost when None in saved file."""
        data: dict[str, Any] = {"profile": "test", "project": "test", "base_url": None, "additional_settings": {}}
        with open(temp_settings_file, "w") as f:
            json.dump(data, f)

        settings = Settings(settings_file=temp_settings_file)

        assert settings.base_url == "http://localhost:8000"
