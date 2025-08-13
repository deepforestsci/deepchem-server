"""
Tests for the pyds package __init__.py module.
"""

import pytest

import pyds
from pyds import BaseClient, Data, Primitives, Settings

class TestPackageInit:
    """Test cases for package initialization."""

    def test_version_exists(self):
        """Test that __version__ is defined."""
        assert hasattr(pyds, "__version__")
        assert isinstance(pyds.__version__, str)
        assert pyds.__version__ == "0.1.0"

    def test_all_exports(self):
        """Test that __all__ contains expected exports."""
        expected_exports = ["Settings", "Primitives", "Data", "BaseClient"]
        assert hasattr(pyds, "__all__")
        assert pyds.__all__ == expected_exports

    def test_settings_import(self):
        """Test that Settings can be imported from package."""
        assert hasattr(pyds, "Settings")
        assert pyds.Settings is Settings

    def test_primitives_import(self):
        """Test that Primitives can be imported from package."""
        assert hasattr(pyds, "Primitives")
        assert pyds.Primitives is Primitives

    def test_data_import(self):
        """Test that Data can be imported from package."""
        assert hasattr(pyds, "Data")
        assert pyds.Data is Data

    def test_base_client_import(self):
        """Test that BaseClient can be imported from package."""
        assert hasattr(pyds, "BaseClient")
        assert pyds.BaseClient is BaseClient

    def test_direct_imports(self):
        """Test that direct imports work correctly."""
        from pyds import BaseClient, Data, Primitives, Settings

        # Test that we can instantiate classes
        settings = Settings(settings_file="test.json")
        assert isinstance(settings, Settings)

        data_client = Data(settings=settings)
        assert isinstance(data_client, Data)

        primitives_client = Primitives(settings=settings)
        assert isinstance(primitives_client, Primitives)

        base_client = BaseClient(settings=settings)
        assert isinstance(base_client, BaseClient)

        # Cleanup
        import os

        if os.path.exists("test.json"):
            os.unlink("test.json")

    def test_star_import(self):
        """Test that star import works correctly."""
        # Test star import functionality by importing all items individually
        # to avoid using wildcard import in class
        import pyds

        # Check that all items in __all__ can be accessed
        for item in pyds.__all__:
            assert hasattr(pyds, item)

    def test_module_attributes(self):
        """Test that module has expected attributes."""
        # Check that the module has the expected structure
        assert hasattr(pyds, "settings")
        assert hasattr(pyds, "primitives")
        assert hasattr(pyds, "data")
        assert hasattr(pyds, "base")

    def test_submodule_imports(self):
        """Test that submodules can be imported."""
        from pyds.base import BaseClient as BaseClientClass
        from pyds.data import Data as DataClass
        from pyds.primitives import Primitives as PrimitivesClass
        from pyds.settings import Settings as SettingsClass

        assert SettingsClass is Settings
        assert PrimitivesClass is Primitives
        assert DataClass is Data
        assert BaseClientClass is BaseClient

    def test_import_error_handling(self):
        """Test that import errors are handled gracefully."""
        # Test importing non-existent attribute
        with pytest.raises(AttributeError):
            import pyds

            _ = getattr(pyds, "NonExistentClass")

    def test_package_structure(self):
        """Test that package structure is as expected."""
        import pkgutil

        import pyds

        # Get all modules in the package
        package_modules = [name for _, name, _ in pkgutil.iter_modules(pyds.__path__)]

        expected_modules = ["settings", "primitives", "data", "base"]
        for module in expected_modules:
            assert module in package_modules

    def test_circular_import_prevention(self):
        """Test that there are no circular import issues."""
        # Try importing all modules in different orders
        # Reload the main package
        import importlib

        import pyds.primitives

        importlib.reload(pyds)

        # All imports should still work
        assert hasattr(pyds, "Settings")
        assert hasattr(pyds, "Primitives")
        assert hasattr(pyds, "Data")
        assert hasattr(pyds, "BaseClient")

    def test_lazy_import_behavior(self):
        """Test that imports work properly on first access."""
        # This tests that the imports in __init__.py work correctly
        # when accessed for the first time

        # Access each class to ensure they're properly imported
        settings_class = pyds.Settings
        primitives_class = pyds.Primitives
        data_class = pyds.Data
        base_class = pyds.BaseClient

        # Verify they are the correct classes
        assert settings_class.__name__ == "Settings"
        assert primitives_class.__name__ == "Primitives"
        assert data_class.__name__ == "Data"
        assert base_class.__name__ == "BaseClient"

    def test_compatibility_imports(self):
        """Test that imports maintain compatibility."""
        # Test that different import styles work

        # Style 1: Import package and access attributes
        import pyds

        Settings1 = pyds.Settings

        # Style 2: Direct import
        from pyds import Settings as Settings2

        # Style 3: Import from submodule
        from pyds.settings import Settings as Settings3

        # All should be the same class
        assert Settings1 is Settings2 is Settings3
