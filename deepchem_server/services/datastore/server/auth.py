"""
API Key Authentication middleware for DeepchemDatastore service.
"""

import os
from typing import Optional

from fastapi import HTTPException, Security, status
from fastapi.security import APIKeyHeader


API_KEY_ENV = "DATASTORE_API_KEY"
API_KEY_HEADER = "X-API-Key"

api_key_header = APIKeyHeader(name=API_KEY_HEADER, auto_error=False)


def get_api_key() -> str:
    """Get the configured API key from environment or use default."""
    api_key = os.getenv(API_KEY_ENV)
    if api_key is None:
        raise ValueError(f"API key not found in environment variable {API_KEY_ENV}")
    return api_key


async def verify_api_key(api_key: Optional[str] = Security(api_key_header)) -> str:
    """Verify the API key from request header.

    Parameters
    ----------
    api_key : str, optional
        API key from X-API-Key header

    Returns
    -------
    str
        The verified API key

    Raises
    ------
    HTTPException
        If API key is missing or invalid
    """
    expected_key = get_api_key()

    if api_key is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing API key. Provide X-API-Key header.",
        )

    if api_key != expected_key:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Invalid API key.",
        )

    return api_key
