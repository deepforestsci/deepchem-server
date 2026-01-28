"""
Resource models for pyds client.

Contains Job and DeepchemData classes that wrap API responses with methods.
"""

from __future__ import annotations

import enum
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

if TYPE_CHECKING:
    from .base import BaseClient


class JobStatus(enum.Enum):
    """Job status enum matching server-side statuses."""
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class Job:
    """
    Represents a submitted job on the DeepChem server.

    Provides methods to check status, wait for completion, cancel, and retrieve results.

    Attributes
    ----------
    id : str
        Unique job identifier.
    status : JobStatus
        Current status of the job.
    program_name : str
        Name of the program/primitive being executed.
    result : Optional[str]
        Result data (typically a datastore address) when job completes.
    error : Optional[str]
        Error message if the job failed.
    created_at : Optional[datetime]
        Timestamp when the job was created.
    updated_at : Optional[datetime]
        Timestamp when the job was last updated.
    started_at : Optional[datetime]
        Timestamp when job execution started.
    completed_at : Optional[datetime]
        Timestamp when job execution completed.
    """

    id: str
    _client: "BaseClient" = field(repr=False)
    status: JobStatus = JobStatus.QUEUED
    program_name: str = ""
    result: Optional[str] = None
    error: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    profile: str = ""
    project: str = ""
    retry_count: int = 0
    max_retries: int = 3

    @classmethod
    def from_dict(cls, data: Dict[str, Any], client: "BaseClient") -> "Job":
        """
        Create a Job from API response dictionary.

        Parameters
        ----------
        data : dict
            API response containing job details.
        client : BaseClient
            Client instance for making subsequent API calls.

        Returns
        -------
        Job
            Initialized Job instance.
        """
        status_str = data.get("status", "queued")
        try:
            status = JobStatus(status_str)
        except ValueError:
            status = JobStatus.QUEUED

        def parse_datetime(val: Optional[str]) -> Optional[datetime]:
            if val:
                try:
                    return datetime.fromisoformat(val.replace("Z", "+00:00"))
                except (ValueError, AttributeError):
                    pass
            return None

        return cls(
            id=data.get("id") or data.get("job_id") or "",
            _client=client,
            status=status,
            program_name=data.get("program_name", ""),
            result=data.get("result"),
            error=data.get("error"),
            created_at=parse_datetime(data.get("created_at")),
            updated_at=parse_datetime(data.get("updated_at")),
            started_at=parse_datetime(data.get("started_at")),
            completed_at=parse_datetime(data.get("completed_at")),
            profile=data.get("profile", ""),
            project=data.get("project", ""),
            retry_count=data.get("retry_count", 0),
            max_retries=data.get("max_retries", 3),
        )

    def refresh(self) -> "Job":
        """
        Fetch the latest job status from the server and update this instance.

        Returns
        -------
        Job
            This Job instance with updated attributes.

        Raises
        ------
        Exception
            If the API request fails.
        """
        response = self._client._get(f"/v1/jobs/{self.id}")
        data = self._client._validate_response(response)
        updated = Job.from_dict(data, self._client)

        # Update self with new values
        self.status = updated.status
        self.result = updated.result
        self.error = updated.error
        self.updated_at = updated.updated_at
        self.started_at = updated.started_at
        self.completed_at = updated.completed_at
        self.program_name = updated.program_name or self.program_name
        self.profile = updated.profile or self.profile
        self.project = updated.project or self.project
        self.retry_count = updated.retry_count
        self.max_retries = updated.max_retries

        return self

    def wait(self, timeout: Optional[float] = None, interval: float = 2.0) -> "Job":
        """
        Block until the job reaches a terminal state (COMPLETED or FAILED).

        Parameters
        ----------
        timeout : float, optional
            Maximum time to wait in seconds. None means wait indefinitely.
        interval : float
            Polling interval in seconds (default: 2.0).

        Returns
        -------
        Job
            This Job instance with final status.

        Raises
        ------
        TimeoutError
            If the timeout is reached before the job completes.
        Exception
            If the job failed.
        """
        start_time = time.time()

        while not self.is_terminal:
            if timeout is not None:
                elapsed = time.time() - start_time
                if elapsed >= timeout:
                    raise TimeoutError(f"Job {self.id} did not complete within {timeout} seconds")

            time.sleep(interval)
            self.refresh()

        # Check for failure - either explicit FAILED status or has error
        if self.status == JobStatus.FAILED or self.error:
            raise Exception(f"Job {self.id} failed: {self.error}")

        return self

    def cancel(self) -> bool:
        """
        Cancel the job if it's pending or queued.

        Returns
        -------
        bool
            True if cancellation was successful, False otherwise.

        Raises
        ------
        Exception
            If the API request fails.
        """
        response = self._client._delete(f"/v1/jobs/{self.id}")
        try:
            self._client._validate_response(response)
            self.status = JobStatus.FAILED
            self.error = "Cancelled by user"
            return True
        except Exception:
            return False

    def result_as_data(self) -> Optional["DeepchemData"]:
        """
        Convert the job result to a DeepchemData object if applicable.

        Returns
        -------
        DeepchemData or None
            DeepchemData object if the result is a valid datastore address,
            None otherwise.
        """
        if not self.result:
            return None

        # Server may return result with surrounding quotes, strip them
        result = self.result.strip('"').strip("'")

        if result.startswith("deepchem://"):
            return DeepchemData.from_address(result, self._client)
        return None

    @property
    def is_terminal(self) -> bool:
        """Check if the job is in a terminal state.
        
        A job is terminal if:
        - Status is COMPLETED or FAILED
        - OR job has an error (even if retries remain, treat as terminal for client)
        """
        if self.status in {JobStatus.COMPLETED, JobStatus.FAILED}:
            return True
        # Also terminal if there's an error - server doesn't actually retry
        if self.error:
            return True
        return False

    @property
    def is_running(self) -> bool:
        """Check if the job is currently running."""
        return self.status == JobStatus.RUNNING

    @property
    def is_successful(self) -> bool:
        """Check if the job completed successfully."""
        return self.status == JobStatus.COMPLETED


@dataclass
class DeepchemData:
    """
    Represents a data object stored in the DeepChem datastore.

    Provides methods to download data and access metadata.

    Attributes
    ----------
    address : str
        Full deepchem:// address of the data.
    profile : str
        Profile name from the address.
    project : str
        Project name from the address.
    key : str
        Data key/path within the project.
    """

    address: str
    _client: "BaseClient" = field(repr=False)
    profile: str = ""
    project: str = ""
    key: str = ""

    @classmethod
    def from_address(cls, address: str, client: "BaseClient") -> "DeepchemData":
        """
        Create a DeepchemData from a deepchem:// address string.

        Parameters
        ----------
        address : str
            Full deepchem:// address (e.g., "deepchem://profile/project/path/to/data").
        client : BaseClient
            Client instance for making subsequent API calls.

        Returns
        -------
        DeepchemData
            Initialized DeepchemData instance.
        """
        profile = ""
        project = ""
        key = ""

        if address.startswith("deepchem://"):
            parts = address[len("deepchem://"):].split("/", 2)
            if len(parts) >= 1:
                profile = parts[0]
            if len(parts) >= 2:
                project = parts[1]
            if len(parts) >= 3:
                key = parts[2]

        return cls(
            address=address,
            _client=client,
            profile=profile,
            project=project,
            key=key,
        )

    @classmethod
    def from_dict(cls, data: Dict[str, Any], client: "BaseClient") -> "DeepchemData":
        """
        Create a DeepchemData from API response dictionary.

        Parameters
        ----------
        data : dict
            API response containing dataset_address key.
        client : BaseClient
            Client instance for making subsequent API calls.

        Returns
        -------
        DeepchemData
            Initialized DeepchemData instance.
        """
        address = data.get("dataset_address", "")
        return cls.from_address(address, client)

    def download(
        self,
        destination: Union[str, Path],
        format: str = "auto",
    ) -> Path:
        """
        Download the data to a local file.

        Parameters
        ----------
        destination : str or Path
            Local path where the data should be saved.
        format : str
            Download format: 'zip' for directories, 'auto' for automatic detection.

        Returns
        -------
        Path
            Path to the downloaded file.

        Raises
        ------
        Exception
            If the download fails.
        """
        destination = Path(destination)

        params = {
            "profile_name": self.profile,
            "project_name": self.project,
            "address": self.address,
            "format": format,
        }

        response = self._client._get("/v1/data/download", params=params, stream=True)

        if response.status_code >= 400:
            raise Exception(f"Download failed: HTTP {response.status_code}")

        # Ensure parent directory exists
        destination.parent.mkdir(parents=True, exist_ok=True)

        with open(destination, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        return destination

    def __str__(self) -> str:
        return self.address
