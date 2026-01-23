import multiprocessing as mp

import os
import socket
import subprocess
import sys
import time
from pathlib import Path
import shutil
from urllib.error import URLError
from urllib.request import urlopen

import pandas as pd
import pytest

from deepchem_server.core.datastore import DiskDataStore
from deepchem_server.services.datastore.client import DatastoreClient


@pytest.fixture(scope="session", autouse=True)
def set_multiprocessing_start_method():
    """Set multiprocessing start method to 'fork' for the test session.

    On macOS, the default is 'spawn' which creates fresh Python processes
    that don't inherit global state. This causes 'Datastore not set' errors
    in multicore featurization tests. Using 'fork' ensures child processes
    inherit the global datastore configuration.
    """
    try:
        mp.set_start_method('fork', force=True)
    except RuntimeError:
        # Already set, ignore
        pass


@pytest.fixture(scope="session")
def datastore_api_key() -> str:
    # The datastore server requires DATASTORE_API_KEY to be set. For tests we
    # default to a known dev key, but allow overrides via environment.
    return os.getenv("DATASTORE_API_KEY", "dev-api-key")


@pytest.fixture(scope="session")
def datastore_service(tmp_path_factory: pytest.TempPathFactory, datastore_api_key: str):
    """Start a demo datastore service for the pytest session.

    This avoids reusing a persisted external datastore between test runs.

    Override:
      - Set DEEPCHEM_TEST_DATASTORE_URL to use an already-running service and
        skip starting a subprocess.
    """
    override_url = os.getenv("DEEPCHEM_TEST_DATASTORE_URL")
    if override_url:
        yield (override_url.rstrip("/"), None)
        return

    host = "127.0.0.1"
    base_dir = Path(tmp_path_factory.mktemp("datastore_service"))

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((host, 0))
        port = s.getsockname()[1]

    # Keep logs outside the datastore base dir so we can safely wipe the base dir
    # between tests to avoid FileExistsError collisions.
    log_dir = Path(tmp_path_factory.mktemp("datastore_service_logs"))
    log_path = log_dir / "datastore.log"
    log_file = open(log_path, "w", encoding="utf-8")

    env = os.environ.copy()
    env["DATASTORE_API_KEY"] = datastore_api_key
    env["PYTHONUNBUFFERED"] = "1"

    cmd = [
        sys.executable,
        "-m",
        "deepchem_server.services.datastore",
        "--host",
        host,
        "--port",
        str(port),
        "--base-dir",
        str(base_dir),
    ]

    proc = subprocess.Popen(
        cmd,
        stdout=log_file,
        stderr=subprocess.STDOUT,
        env=env,
        cwd=os.getcwd(),
    )

    health_url = f"http://{host}:{port}/api/v1/healthcheck"
    deadline = time.time() + 30.0
    while time.time() < deadline:
        if proc.poll() is not None:
            log_file.flush()
            log_file.close()
            logs = ""
            try:
                logs = log_path.read_text(encoding="utf-8", errors="replace")
            except Exception:
                pass
            raise RuntimeError(
                "Datastore service exited before becoming healthy.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Logs:\n{logs}"
            )
        try:
            with urlopen(health_url, timeout=1.0) as resp:  # nosec B310
                if resp.status == 200:
                    break
        except URLError:
            time.sleep(0.1)
            continue
    else:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            proc.kill()
        log_file.flush()
        log_file.close()
        logs = ""
        try:
            logs = log_path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            pass
        raise TimeoutError(
            "Timed out waiting for datastore service to become healthy.\n"
            f"Health URL: {health_url}\n"
            f"Command: {' '.join(cmd)}\n"
            f"Logs:\n{logs}"
        )

    url = f"http://{host}:{port}"
    try:
        yield (url, base_dir)
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            proc.kill()
        log_file.flush()
        log_file.close()


@pytest.fixture(scope="session")
def datastore_url(datastore_service) -> str:
    url, _ = datastore_service
    return url


@pytest.fixture(scope="session")
def datastore_base_dir(datastore_service) -> Path | None:
    _, base_dir = datastore_service
    return base_dir


@pytest.fixture(autouse=True)
def clean_datastore_storage(datastore_base_dir: Path | None) -> None:
    """Ensure each test starts with an empty datastore.

    This prevents collisions from reused constant output/model keys like
    'feat_test', 'feat', 'gcn_feat', 'infer.csv', 'test_model', etc.

    If DEEPCHEM_TEST_DATASTORE_URL is set (external datastore), we skip wiping.
    """
    if datastore_base_dir is None:
        return

    for child in datastore_base_dir.iterdir():
        try:
            if child.is_dir():
                shutil.rmtree(child, ignore_errors=True)
            else:
                child.unlink(missing_ok=True)
        except Exception:
            # Best-effort cleanup: ignore transient filesystem errors.
            pass


@pytest.fixture(scope="session")
def datastore_client(datastore_url: str, datastore_api_key: str):
    return DatastoreClient(url=datastore_url, api_key=datastore_api_key)


@pytest.fixture
def disk_datastore(tmp_path, datastore_url: str, datastore_api_key: str):
    client = DatastoreClient(url=datastore_url, api_key=datastore_api_key)
    fds = DiskDataStore(
        client=client, profile_name="test", project_name="user", basedir=str(tmp_path)
    )
    return fds


@pytest.fixture
def alternate_disk_datastore(tmp_path, datastore_url: str, datastore_api_key: str):
    client = DatastoreClient(url=datastore_url, api_key=datastore_api_key)
    fds = DiskDataStore(
        client=client,
        profile_name="alternate-test",
        project_name="alternate-user",
        basedir=str(tmp_path),
    )
    return fds


@pytest.fixture
def tmp_csv(tmp_path):
    path = tmp_path / "temp.csv"
    col = [1, 2, 3, 4, 5]
    df = pd.DataFrame({'col': col})
    df.to_csv(path, index=False)
    return str(path)
