"""
Entry point for running RQ workers.

Usage:
    python -m deepchem_server.services.jobs.run_worker
    
Environment variables:
    REDIS_URL: Redis connection URL
    JOBS_DATABASE_URL: Database connection URL
    DATADIR: Data directory for DiskDataStore
    OBJC_DISABLE_INITIALIZE_FORK_SAFETY: Set to YES on macOS to avoid fork crashes
"""
import argparse
import logging
import os
import platform

# Fix for macOS fork() crash with Objective-C libraries (DeepChem/RDKit)
# This MUST be set before importing anything that uses Objective-C
if platform.system() == "Darwin":
    os.environ["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

from redis import Redis
from rq import Queue, SimpleWorker, Worker

from deepchem_server.services.jobs.models import init_db
from deepchem_server.services.jobs.queue import JobQueue

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run an RQ worker."""
    parser = argparse.ArgumentParser(description="DeepChem Job Worker")
    parser.add_argument(
        "--queue",
        default="deepchem_jobs",
        help="Queue name to process (default: deepchem_jobs)",
    )
    parser.add_argument(
        "--redis-url",
        default=os.getenv("REDIS_URL", "redis://localhost:6379/0"),
        help="Redis URL",
    )
    parser.add_argument(
        "--recover",
        action="store_true",
        help="Recover stale jobs on startup",
    )
    parser.add_argument(
        "--simple",
        action="store_true",
        help="Use SimpleWorker (no forking, runs in main process)",
    )

    args = parser.parse_args()

    # Initialize database
    init_db()

    # Recover stale jobs if requested
    if args.recover:
        logger.info("Recovering stale jobs...")
        job_queue = JobQueue(redis_url=args.redis_url)
        recovered = job_queue.recover_stale_jobs()
        logger.info(f"Recovered {recovered} stale jobs")

    # Connect to Redis
    redis_conn = Redis.from_url(args.redis_url)
    queue = Queue(args.queue, connection=redis_conn)

    logger.info(f"Starting worker for queue: {args.queue}")
    logger.info(f"Redis: {args.redis_url}")

    # Use SimpleWorker on macOS to avoid fork issues with DeepChem/RDKit
    # SimpleWorker runs jobs in the main process without forking
    worker: Worker | SimpleWorker
    if args.simple or platform.system() == "Darwin":
        logger.info("Using SimpleWorker (no forking) for macOS compatibility")
        worker = SimpleWorker([queue], connection=redis_conn)
    else:
        worker = Worker([queue], connection=redis_conn)

    worker.work()


if __name__ == "__main__":
    main()
