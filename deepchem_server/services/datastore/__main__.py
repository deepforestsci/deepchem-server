"""
Entry point for running DeepchemDatastore service standalone.

Usage:
    python -m deepchem_server.services.datastore
    
    # Or with custom settings:
    DATASTORE_BASE_DIR=/data DATASTORE_PORT=8081 DATASTORE_API_KEY=secret \
        python -m deepchem_server.services.datastore
"""
import argparse
import logging
import os

import uvicorn

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run the datastore service."""
    parser = argparse.ArgumentParser(description="DeepchemDatastore Service")
    parser.add_argument(
        "--host",
        default=os.getenv("DATASTORE_HOST", "0.0.0.0"),
        help="Host to bind to (default: 0.0.0.0)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.getenv("DATASTORE_PORT", "8081")),
        help="Port to bind to (default: 8081)",
    )
    parser.add_argument(
        "--base-dir",
        default=os.getenv("DATASTORE_BASE_DIR", "./datastore_data"),
        help="Base directory for storage (default: ./datastore_data)",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload for development",
    )

    args = parser.parse_args()

    # Set environment variables for the service
    os.environ["DATASTORE_BASE_DIR"] = args.base_dir
    os.environ["DATASTORE_PORT"] = str(args.port)

    logger.info("Starting DeepchemDatastore Service")
    logger.info(f"  Host: {args.host}")
    logger.info(f"  Port: {args.port}")
    logger.info(f"  Base Dir: {args.base_dir}")

    uvicorn.run(
        "deepchem_server.services.datastore.server.api:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
