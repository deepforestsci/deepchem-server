"""Large library screening script: upload -> partition -> featurize each part -> infer each part."""

import sys
import time
import uuid
from pathlib import Path

from pyds import BaseClient, Data, Featurize, Infer, Partition, Settings  # type: ignore


repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root / "pyds"))

BASE_URL = "http://deepchem-server"
PROFILE = "test_profile"
PROJECT = "test_project"
DATASET_PATH = repo_root / "deepchem_server" / "core" / "tests" / "assets" / "zinc250k.csv"

SCREENING_MODEL = "deepchem://test_profile/test_project/rf_logp_5a4c5f9f"

N_PARTITIONS = 25

if not DATASET_PATH.exists():
    raise FileNotFoundError(f"Dataset not found: {DATASET_PATH}")

run_id = uuid.uuid4().hex[:8]
print(f"Case Study 2 - Run ID: {run_id} - {time.strftime('%Y-%m-%d %H:%M:%S')}")

settings = Settings(profile=PROFILE, project=PROJECT, base_url=BASE_URL)
base_client = BaseClient(settings=settings)
data_client = Data(settings=settings)
featurize_client = Featurize(settings=settings)
partition_client = Partition(settings=settings)
infer_client = Infer(settings=settings)

print(f"healthcheck: {base_client.healthcheck()}")

# 1) Upload dataset
upload_result = data_client.upload_data(
    file_path=DATASET_PATH,
    filename=f"library_{run_id}.csv",
    description=f"Screening library upload for run {run_id}",
)
library_address = upload_result["dataset_address"]
print(f"Uploaded library to: {library_address}")

# 2) Partition the CSV dataset into N smaller parts
partition_result = partition_client.run(
    dataset_address=library_address,
    n_partition=N_PARTITIONS,
    shuffle=False,
)
partitions = partition_result["partitioned_dataset_addresses"]
print(f"Partitioned into {len(partitions)} parts:")
for i, part_addr in enumerate(partitions):
    print(f"  partition[{i}]: {part_addr}")

# 3) Featurize each partition and run inference on it
prediction_files = []
for i, part_addr in enumerate(partitions):
    featurize_result = featurize_client.run(
        dataset_address=part_addr,
        featurizer="ecfp",
        output=f"library_ecfp_{run_id}_part{i}",
        dataset_column="smiles",
        label_column="logp",
        feat_kwargs={
            "radius": 2,
            "size": 1024
        },
    )
    featurized_part_address = featurize_result["featurized_file_address"]
    print(f"  featurized partition[{i}]: {featurized_part_address}")

    infer_result = infer_client.run(
        model_address=SCREENING_MODEL,
        data_address=part_addr,
        dataset_column="smiles",
        output=f"screening_{run_id}_part{i}",
    )
    pred_address = infer_result["inference_results_address"]
    prediction_files.append(pred_address)
    print(f"  predictions for partition[{i}]: {pred_address}")

print(f"Screening complete. {len(prediction_files)} prediction files produced.")
