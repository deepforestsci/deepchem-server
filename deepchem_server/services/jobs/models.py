"""
Job models for persistent job state tracking.

Uses SQLAlchemy with MySQL backend (SQLite for testing).
"""
import enum
import os
from datetime import datetime, timedelta
from typing import Any, Dict, Optional

from sqlalchemy import Column, DateTime, Enum, Integer, String, Text, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


DATABASE_URL = os.getenv("JOBS_DATABASE_URL")

# SQLAlchemy setup
Base = declarative_base()
engine = None
SessionLocal = None


def init_db(database_url: Optional[str] = None) -> None:
    """Initialize database connection and create tables.
    
    Parameters
    ----------
    database_url : str, optional
        Database connection URL. If not provided, uses JOBS_DATABASE_URL env var.
        
    Raises
    ------
    ValueError
        If no database URL is provided and JOBS_DATABASE_URL is not set.
    """
    global engine, SessionLocal
    url = database_url or DATABASE_URL

    if url is None:
        raise ValueError("JOBS_DATABASE_URL environment variable is required. "
                         "Set it to a valid database URL, e.g.: "
                         "mysql+pymysql://user:pass@localhost:3306/deepchem_jobs")

    # Configure engine based on database type
    if url.startswith("sqlite"):
        # SQLite doesn't support connection pooling options
        engine = create_engine(url, pool_pre_ping=True)
    else:
        # MySQL/other databases: use connection pooling
        engine = create_engine(
            url,
            pool_pre_ping=True,
            pool_size=5,
            max_overflow=10,
            pool_recycle=3600,  # Recycle connections after 1 hour
        )

    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    Base.metadata.create_all(bind=engine)


def get_session():
    """Get a database session."""
    if SessionLocal is None:
        init_db()
    return SessionLocal()


class JobStatus(enum.Enum):
    """Job status enum."""
    PENDING = "pending"  # Created but not yet queued
    QUEUED = "queued"  # In Redis queue waiting
    RUNNING = "running"  # Worker is processing
    COMPLETED = "completed"  # Successfully finished
    FAILED = "failed"  # Failed after retries


class Job(Base):  # type: ignore
    """Job model for persistent storage."""

    __tablename__ = "jobs"

    id = Column(String(36), primary_key=True)
    status = Column(Enum(JobStatus), default=JobStatus.PENDING, nullable=False)  # type: ignore
    program_name = Column(String(100), nullable=False)
    program_json = Column(Text, nullable=False)
    profile = Column(String(100), nullable=False)
    project = Column(String(100), nullable=False)

    # Result and error tracking
    result = Column(Text, nullable=True)
    error = Column(Text, nullable=True)
    retry_count = Column(Integer, default=0)
    max_retries = Column(Integer, default=3)

    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    expires_at = Column(DateTime, nullable=True)

    def __init__(
        self,
        id: str,
        program_name: str,
        program_json: str,
        profile: str,
        project: str,
        max_retries: int = 3,
        retention_days: int = 7,
    ):
        self.id = id  # type: ignore
        self.program_name = program_name  # type: ignore
        self.program_json = program_json  # type: ignore
        self.profile = profile  # type: ignore
        self.project = project  # type: ignore
        self.max_retries = max_retries  # type: ignore
        self.expires_at = datetime.utcnow() + timedelta(days=retention_days)  # type: ignore

    def to_dict(self) -> Dict[str, Any]:
        """Convert job to dictionary."""
        return {
            "id": self.id,
            "status": self.status.value,
            "program_name": self.program_name,
            "profile": self.profile,
            "project": self.project,
            "result": self.result,
            "error": self.error,
            "retry_count": self.retry_count,
            "max_retries": self.max_retries,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": self.completed_at.isoformat() if self.completed_at else None,
        }
