"""Basic Zinc5k pipeline script: featurize -> train -> evaluate -> infer."""

import sys
import time
import uuid
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root / "pyds"))

repo_root = Path(__file__).resolve().parents[1]

from pyds import (
    BaseClient,
    Data,
    Evaluate,
    Featurize,
    Infer,
    Settings,
    TVTSplit,
    Train,
)

# Global config
BASE_URL = "http://deepchem-server"
PROFILE = "test_profile"
PROJECT = "test_project"
DATASET_PATH = repo_root / "deepchem_server" / "core" / "tests" / "assets" / "zinc5k.csv"

if not DATASET_PATH.exists():
    raise FileNotFoundError(f"Dataset not found: {DATASET_PATH}")

run_id = uuid.uuid4().hex[:8]
print(f"Case Study 1 - Run ID: {run_id} - {time.strftime('%Y-%m-%d %H:%M:%S')}")

settings = Settings(profile=PROFILE, project=PROJECT, base_url=BASE_URL)
base_client = BaseClient(settings=settings)
data_client = Data(settings=settings)
featurize_client = Featurize(settings=settings)
split_client = TVTSplit(settings=settings)
train_client = Train(settings=settings)
evaluate_client = Evaluate(settings=settings)
infer_client = Infer(settings=settings)

# Healthcheck
health = base_client.healthcheck()
print(f"healthcheck: {health}")

# 1) Upload
upload_result = data_client.upload_data(
    file_path=DATASET_PATH,
    filename=f"zinc5k_{run_id}.csv",
    description=f"Zinc5k upload for run {run_id}",
)
dataset_address = upload_result["dataset_address"]
print(f"Uploaded dataset to: {dataset_address}")

# 2) Featurize
featurize_result = featurize_client.run(
    dataset_address=dataset_address,
    featurizer="ecfp",
    output=f"zinc5k_ecfp_{run_id}",
    dataset_column="smiles",
    label_column="logp",
    feat_kwargs={
        "radius": 2,
        "size": 1024,
    },
)
featurized_address = featurize_result["featurized_file_address"]
print(f"Featurized dataset address: {featurized_address}")

# 3) Split
split_result = split_client.run(
    splitter_type="random",
    dataset_address=featurized_address,
    frac_train=0.8,
    frac_valid=0.1,
    frac_test=0.1,
)
train_addr, valid_addr, test_addr = split_result["train_valid_test_split_results_address"]
print(f"Train dataset address: {train_addr}")
print(f"Valid dataset address: {valid_addr}")
print(f"Test dataset address: {test_addr}")

# 4) Train
train_result = train_client.run(
    dataset_address=train_addr,
    model_type="random_forest_regressor",
    model_name=f"rf_logp_{run_id}",
    init_kwargs={"n_estimators": 100, "random_state": 42},
    train_kwargs={},
)
model_address = train_result["trained_model_address"]
print(f"Trained model address: {model_address}")

# 5) Evaluate
evaluate_result = evaluate_client.run(
    dataset_addresses=[test_addr],
    model_address=model_address,
    metrics=["pearson_r2_score", "rms_score", "mae_error"],
    output_key=f"eval_logp_{run_id}",
    is_metric_plots=False,
)
evaluation_address = evaluate_result["evaluation_result_address"]
print(f"Evaluation results address: {evaluation_address}")

# 6) Infer
inference_output = f"infer_logp_{run_id}"

infer_result = infer_client.run(
    model_address=model_address,
    data_address=test_addr,
    output=inference_output,
    dataset_column="smiles",
)

inference_address = infer_result["inference_results_address"]
print(f"Inference results address: {inference_address}")
