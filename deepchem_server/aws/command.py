import argparse
import json
import logging
import sys

from deepchem_server.core import config
from deepchem_server.core.datastore import S3DataStore
from deepchem_server.core.primitives.compute import ComputeWorkflow

logging.basicConfig(stream=sys.stdout, level=logging.INFO, force=True)
logger = logging.getLogger("deepchem-batch")


def main(args: argparse.Namespace) -> None:
    try:
        datastore = S3DataStore(
            profile_name=args.profile,
            project_name=args.project,
            bucket_name=args.bucket,
        )
        config.set_datastore(datastore)
        try:
            program = json.loads(args.program)
        except json.JSONDecodeError as exc:
            logger.error(f"Invalid JSON in --program: {exc}")
            raise
        logger.info(f"Running program: {program.get('program_name', 'unknown')}")
        ComputeWorkflow(program).execute()
        logger.info("Done")
    except Exception:
        logger.exception("Job failed")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="deepchem-server Batch container entrypoint")
    parser.add_argument("--program", required=True, help="JSON-encoded program dict")
    parser.add_argument("--profile", required=True, help="Profile name")
    parser.add_argument("--project", required=True, help="Project name")
    parser.add_argument("--bucket",  required=True, help="S3 bucket name")
    main(parser.parse_args())
