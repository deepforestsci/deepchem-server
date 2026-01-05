"""Allow running the CLI as a module: python -m cli"""
from .main import app

if __name__ == "__main__":
    app()
