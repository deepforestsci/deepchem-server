"""
Tests for the Settings module.
"""

import json
import os
import tempfile

from pyds.settings import Settings


class TestSettings:
    """Test cases for Settings class."""

    def test_init_with_new_file(self, temp_settings_file):
        """Test Settings initialization with a new file."""
        settings = Settings(
            settings_file=temp_settings_file,
            profile="test-profile",
            project="test-project",
            base_url="http://localhost:8000",
        )

        assert settings.get_profile() == "test-profile"
        assert settings.get_project() == "test-project"
        assert settings.get_base_url() == "http://localhost:8000"
        assert os.path.exists(temp_settings_file)

    def test_init_with_existing_file(self, temp_settings_file):
        """Test Settings initialization with existing file."""
        # Create initial settings
        initial_data = {
            "profile": "initial-profile",
            "project": "initial-project",
            "base_url": "http://initial:8000",
        }
        with open(temp_settings_file, "w") as f:
            json.dump(initial_data, f)

        # Load settings without overriding
        settings = Settings(settings_file=temp_settings_file)

        assert settings.get_profile() == "initial-profile"
        assert settings.get_project() == "initial-project"
        assert settings.get_base_url() == "http://initial:8000"

    def test_init_with_override(self, temp_settings_file):
        """Test Settings initialization with override parameters."""
        # Create initial settings
        initial_data = {
            "profile": "initial-profile",
            "project": "initial-project",
            "base_url": "http://initial:8000",
        }
        with open(temp_settings_file, "w") as f:
            json.dump(initial_data, f)

        # Load settings with overrides
        settings = Settings(settings_file=temp_settings_file,
                            profile="new-profile",
                            project="new-project")

        assert settings.get_profile() == "new-profile"
        assert settings.get_project() == "new-project"
        assert settings.get_base_url(
        ) == "http://initial:8000"  # Not overridden

    def test_set_and_get_profile(self, settings_instance):
        """Test setting and getting profile."""
        settings_instance.set_profile("new-profile")
        assert settings_instance.get_profile() == "new-profile"

    def test_set_and_get_project(self, settings_instance):
        """Test setting and getting project."""
        settings_instance.set_project("new-project")
        assert settings_instance.get_project() == "new-project"

    def test_set_and_get_base_url(self, settings_instance):
        """Test setting and getting base URL."""
        settings_instance.set_base_url("http://newhost:9000")
        assert settings_instance.get_base_url() == "http://newhost:9000"

    def test_set_and_get_custom_setting(self, settings_instance):
        """Test setting and getting custom settings."""
        settings_instance.set_setting("custom_key", "custom_value")
        assert settings_instance.get_setting("custom_key") == "custom_value"

    def test_get_nonexistent_setting(self, settings_instance):
        """Test getting a non-existent setting."""
        assert settings_instance.get_setting("nonexistent") is None

    def test_get_nonexistent_setting_with_default(self, settings_instance):
        """Test getting a non-existent setting with default value."""
        default_value = "default"
        assert settings_instance.get_setting("nonexistent",
                                             default_value) == default_value

    def test_is_configured_true(self, settings_instance):
        """Test is_configured returns True when profile and project are set."""
        assert settings_instance.is_configured() is True

    def test_is_configured_false_no_profile(self, temp_settings_file):
        """Test is_configured returns False when profile is missing."""
        settings = Settings(settings_file=temp_settings_file,
                            project="test-project")
        assert settings.is_configured() is False

    def test_is_configured_false_no_project(self, temp_settings_file):
        """Test is_configured returns False when project is missing."""
        settings = Settings(settings_file=temp_settings_file,
                            profile="test-profile")
        assert settings.is_configured() is False

    def test_save_and_load(self, temp_settings_file):
        """Test saving and loading settings."""
        settings = Settings(settings_file=temp_settings_file,
                            profile="test-profile",
                            project="test-project")
        settings.set_setting("custom_key", "custom_value")

        # Create new instance to test loading
        new_settings = Settings(settings_file=temp_settings_file)

        assert new_settings.get_profile() == "test-profile"
        assert new_settings.get_project() == "test-project"
        assert new_settings.get_setting("custom_key") == "custom_value"

    def test_touch_creates_file(self):
        """Test that touch method creates a new file."""
        with tempfile.NamedTemporaryFile(delete=True) as f:
            temp_file = f.name

        # File should not exist
        assert not os.path.exists(temp_file)

        settings = Settings(settings_file=temp_file)

        # File should now exist
        assert os.path.exists(temp_file)

        # Cleanup
        os.unlink(temp_file)

    def test_get_all_settings(self, settings_instance):
        """Test getting all settings."""
        settings_instance.set_setting("key1", "value1")
        settings_instance.set_setting("key2", "value2")

        # Test individual getters since get_all() doesn't exist
        assert settings_instance.get_profile() == "test-profile"
        assert settings_instance.get_project() == "test-project"
        assert settings_instance.get_base_url() == "http://localhost:8000"
        assert settings_instance.get_setting("key1") == "value1"
        assert settings_instance.get_setting("key2") == "value2"

    def test_delete_setting(self, settings_instance):
        """Test deleting a setting by resetting."""
        settings_instance.set_setting("temp_key", "temp_value")
        assert settings_instance.get_setting("temp_key") == "temp_value"

        # Since there's no delete method, test reset functionality
        settings_instance.reset()
        assert settings_instance.get_setting("temp_key") is None

    def test_reset_settings(self, settings_instance):
        """Test resetting all settings."""
        settings_instance.set_setting("key1", "value1")

        # Reset should not raise any exception
        settings_instance.reset()
        assert settings_instance.get_setting("key1") is None

    def test_clear_settings(self, settings_instance):
        """Test clearing all settings using reset."""
        settings_instance.set_setting("key1", "value1")
        settings_instance.set_setting("key2", "value2")

        settings_instance.reset()

        assert settings_instance.get_setting("key1") is None
        assert settings_instance.get_setting("key2") is None
        assert settings_instance.get_profile() is None
        assert settings_instance.get_project() is None

    def test_persistence_across_instances(self, temp_settings_file):
        """Test that settings persist across different instances."""
        # Create first instance and set values
        settings1 = Settings(
            settings_file=temp_settings_file,
            profile="persistent-profile",
            project="persistent-project",
        )
        settings1.set_setting("persistent_key", "persistent_value")

        # Create second instance and verify values persist
        settings2 = Settings(settings_file=temp_settings_file)

        assert settings2.get_profile() == "persistent-profile"
        assert settings2.get_project() == "persistent-project"
        assert settings2.get_setting("persistent_key") == "persistent_value"

    def test_additional_settings_init(self, temp_settings_file):
        """Test initialization with additional settings."""
        additional = {"custom1": "value1", "custom2": "value2"}

        settings = Settings(settings_file=temp_settings_file,
                            additional_settings=additional)

        assert settings.get_setting("custom1") == "value1"
        assert settings.get_setting("custom2") == "value2"
