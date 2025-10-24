"""
Unit tests to verify that pyds imports work correctly after refactoring.
"""


class TestImports:
    """Test suite for verifying pyds package imports."""

    def test_basic_imports(self):
        """Test basic imports from pyds."""
        from pyds import Settings, Data, BaseClient

        assert Settings is not None
        assert Data is not None
        assert BaseClient is not None

    def test_primitive_imports(self):
        """Test importing primitive classes."""
        from pyds import Primitive, Featurize, Train, Evaluate, Infer, TVTSplit

        assert Primitive is not None
        assert Featurize is not None
        assert Train is not None
        assert Evaluate is not None
        assert Infer is not None
        assert TVTSplit is not None

    def test_all_exports_accessible(self):
        """Test that all items in __all__ are accessible."""
        import pyds

        expected_exports = [
            "Settings",
            "Data",
            "BaseClient",
            "Primitive",
            "Featurize",
            "Train",
            "Evaluate",
            "Infer",
            "TVTSplit",
        ]

        for export in expected_exports:
            assert export in pyds.__all__, f"{export} not in __all__"
            assert hasattr(pyds, export), f"{export} not accessible from pyds"

    def test_package_structure(self):
        """Test that the package structure is as expected."""
        import pyds
        import pyds.base
        import pyds.base.client
        import pyds.data
        import pyds.settings
        import pyds.primitives
        import pyds.primitives.base
        import pyds.primitives.featurize
        import pyds.primitives.train
        import pyds.primitives.evaluate
        import pyds.primitives.infer
        import pyds.primitives.splitter

        assert True
