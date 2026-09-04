"""
Unit tests for BaseClient class.
"""

import socket
import threading
import time
from typing import Any, Dict
from unittest.mock import Mock, patch

import pytest
import requests
import responses

from pyds.base.client import BaseClient
from pyds.settings import Settings


def _localhost_listener_that_never_accepts() -> tuple[int, threading.Event, threading.Thread]:
    """
    Bind a TCP socket on 127.0.0.1, listen, and never call accept().

    Clients can complete connect and send a request; the peer will not send an HTTP
    response, so requests read timeouts can be exercised.

    Returns
    -------
    tuple[int, threading.Event, threading.Thread]
        Listening port, event to stop the background thread, and that thread.
    """
    port_ref: list[int] = []
    shutdown = threading.Event()

    def run() -> None:
        srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        srv.bind(("127.0.0.1", 0))
        srv.listen(1)
        port_ref.append(int(srv.getsockname()[1]))
        shutdown.wait(timeout=120.0)
        srv.close()

    thread = threading.Thread(target=run, daemon=True)
    thread.start()
    while not port_ref:
        time.sleep(0.01)
    return port_ref[0], shutdown, thread


class TestBaseClient:
    """Unit tests for BaseClient class."""

    def test_init_with_settings(self, test_settings: Settings) -> None:
        """Test BaseClient initialization with settings."""
        client = BaseClient(settings=test_settings)

        assert client.settings == test_settings
        assert client.base_url == "http://localhost:8000"
        assert isinstance(client.session, requests.Session)

    def test_init_with_base_url_override(self, test_settings: Settings) -> None:
        """Test BaseClient initialization with base URL override."""
        custom_url = "http://custom:9000"
        client = BaseClient(settings=test_settings, base_url=custom_url)

        assert client.base_url == custom_url

    def test_init_without_settings(self) -> None:
        """Test BaseClient initialization without settings creates default."""
        with patch.object(Settings, "__init__", return_value=None) as _:
            mock_settings = Mock()
            mock_settings.get_base_url.return_value = "http://localhost:8000"

            with patch("pyds.base.client.Settings", return_value=mock_settings):
                client = BaseClient()

                assert client.settings == mock_settings

    def test_check_configuration_success(self, base_client: BaseClient) -> None:
        """Test _check_configuration with valid settings."""
        # Should not raise any exception
        base_client._check_configuration()

    def test_check_configuration_missing_profile(self, test_settings_not_configured: Settings) -> None:
        """Test _check_configuration with missing profile."""
        client = BaseClient(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings: profile, project"):
            client._check_configuration()

    def test_check_configuration_missing_project(self, temp_settings_file: str) -> None:
        """Test _check_configuration with missing project."""
        settings = Settings(
            settings_file=temp_settings_file,
            profile="test_profile",
            base_url="http://localhost:8000",
        )
        client = BaseClient(settings=settings)

        with pytest.raises(ValueError, match="Missing required settings: project"):
            client._check_configuration()

    @responses.activate
    def test_make_request_success(self, base_client: BaseClient) -> None:
        """Test successful HTTP request."""
        responses.add(responses.GET, "http://localhost:8000/test", json={"status": "ok"}, status=200)

        response = base_client._make_request("GET", "/test")

        assert response.status_code == 200
        assert response.json() == {"status": "ok"}

    def test_make_request_failure(self, base_client: BaseClient) -> None:
        """Test failed HTTP request."""
        with patch.object(
                base_client.session,
                "request",
                side_effect=requests.RequestException("Connection error"),
        ):
            with pytest.raises(Exception, match="API request failed: Connection error"):
                base_client._make_request("GET", "/test")

    @responses.activate
    def test_get_request(self, base_client: BaseClient) -> None:
        """Test GET request method."""
        responses.add(responses.GET, "http://localhost:8000/test", json={"method": "GET"}, status=200)

        response = base_client._get("/test")

        assert response.status_code == 200
        assert response.json() == {"method": "GET"}

    @responses.activate
    def test_post_request(self, base_client: BaseClient) -> None:
        """Test POST request method."""
        responses.add(responses.POST, "http://localhost:8000/test", json={"method": "POST"}, status=200)

        response = base_client._post("/test", json={"data": "test"})

        assert response.status_code == 200
        assert response.json() == {"method": "POST"}

    @responses.activate
    def test_put_request(self, base_client: BaseClient) -> None:
        """Test PUT request method."""
        responses.add(responses.PUT, "http://localhost:8000/test", json={"method": "PUT"}, status=200)

        response = base_client._put("/test")

        assert response.status_code == 200
        assert response.json() == {"method": "PUT"}

    @responses.activate
    def test_delete_request(self, base_client: BaseClient) -> None:
        """Test DELETE request method."""
        responses.add(responses.DELETE, "http://localhost:8000/test", json={"method": "DELETE"}, status=200)

        response = base_client._delete("/test")

        assert response.status_code == 200
        assert response.json() == {"method": "DELETE"}

    def test_validate_response_success(self, base_client: BaseClient) -> None:
        """Test successful response validation."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"result": "success"}

        result = base_client._validate_response(mock_response)

        assert result == {"result": "success"}

    def test_validate_response_http_error_with_detail(self, base_client: BaseClient) -> None:
        """Test response validation with HTTP error and detail."""
        mock_response = Mock()
        mock_response.status_code = 400
        mock_response.json.return_value = {"detail": "Bad request"}

        with pytest.raises(Exception, match="Bad request"):
            base_client._validate_response(mock_response)

    def test_validate_response_http_error_without_detail(self, base_client: BaseClient) -> None:
        """Test response validation with HTTP error but no detail."""
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.json.side_effect = Exception("Invalid JSON")

        with pytest.raises(Exception, match="HTTP 500"):
            base_client._validate_response(mock_response)

    def test_get_profile_project_with_params(self, base_client: BaseClient) -> None:
        """Test _get_profile_project with explicit parameters."""
        profile, project = base_client._get_profile_project("custom_profile", "custom_project")

        assert profile == "custom_profile"
        assert project == "custom_project"

    def test_get_profile_project_from_settings(self, base_client: BaseClient) -> None:
        """Test _get_profile_project using settings."""
        profile, project = base_client._get_profile_project()

        assert profile == "test_profile"
        assert project == "test_project"

    def test_get_profile_project_missing_settings(self, test_settings_not_configured: Settings) -> None:
        """Test _get_profile_project with missing settings."""
        client = BaseClient(settings=test_settings_not_configured)

        with pytest.raises(ValueError, match="Missing required settings"):
            client._get_profile_project()

    def test_get_profile_and_project_with_params(self, base_client: BaseClient) -> None:
        """Test _get_profile_and_project with explicit parameters."""
        profile, project = base_client._get_profile_and_project("custom_profile", "custom_project")

        assert profile == "custom_profile"
        assert project == "custom_project"

    def test_get_profile_and_project_from_settings(self, base_client: BaseClient) -> None:
        """Test _get_profile_and_project using settings."""
        profile, project = base_client._get_profile_and_project()

        assert profile == "test_profile"
        assert project == "test_project"

    @responses.activate
    def test_healthcheck_success(self, base_client: BaseClient, sample_healthcheck_response: Dict[str, Any]) -> None:
        """Test successful healthcheck."""
        responses.add(
            responses.GET,
            "http://localhost:8000/healthcheck",
            json=sample_healthcheck_response,
            status=200,
        )

        result = base_client.healthcheck()

        assert result == sample_healthcheck_response

    @responses.activate
    def test_healthcheck_failure(self, base_client: BaseClient) -> None:
        """Test failed healthcheck."""
        responses.add(
            responses.GET,
            "http://localhost:8000/healthcheck",
            json={"detail": "Service unavailable"},
            status=503,
        )

        with pytest.raises(Exception, match="Service unavailable"):
            base_client.healthcheck()

    def test_get_settings(self, base_client: BaseClient) -> None:
        """Test get_settings method."""
        settings = base_client.get_settings()
        assert settings == base_client.settings

    def test_get_base_url(self, base_client: BaseClient) -> None:
        """Test get_base_url method."""
        url = base_client.get_base_url()
        assert url == "http://localhost:8000"

    def test_set_base_url(self, base_client: BaseClient) -> None:
        """Test set_base_url method."""
        new_url = "http://newserver:9000/"
        base_client.set_base_url(new_url)

        # Should strip trailing slash
        assert base_client.base_url == "http://newserver:9000"

    def test_close_session(self, base_client: BaseClient) -> None:
        """Test close method."""
        mock_session = Mock()
        base_client.session = mock_session

        base_client.close()

        mock_session.close.assert_called_once()


class TestBaseClientLiveTransport:
    """
    Tests for live server transport behavior.
    """

    def test_make_request_connection_refused_live(self, test_settings: Settings) -> None:
        """Test connection refused live."""
        client = BaseClient(settings=test_settings, base_url="http://127.0.0.1:1")
        with pytest.raises(Exception) as exc_info:
            client._make_request("GET", "/healthcheck", timeout=5.0)
        msg = str(exc_info.value).lower()
        assert "api request failed" in msg
        assert "refused" in msg or "failed to establish" in msg

    def test_healthcheck_connection_refused_live(self, test_settings: Settings) -> None:
        """Test healthcheck connection refused live."""
        client = BaseClient(settings=test_settings, base_url="http://127.0.0.1:1")
        with pytest.raises(Exception) as exc_info:
            client.healthcheck()
        assert "api request failed" in str(exc_info.value).lower()

    def test_make_request_read_timeout_live(self, test_settings: Settings) -> None:
        """Test make request read timeout live."""
        port, shutdown, thread = _localhost_listener_that_never_accepts()
        try:
            client = BaseClient(settings=test_settings, base_url=f"http://127.0.0.1:{port}")
            with pytest.raises(Exception) as exc_info:
                client._make_request("GET", "/", timeout=2.0)
        finally:
            shutdown.set()
            thread.join(timeout=2.0)
        msg = str(exc_info.value).lower()
        assert "api request failed" in msg
        assert "read" in msg and "timed out" in msg

    def test_healthcheck_read_timeout_live(self, test_settings: Settings) -> None:
        """Test healthcheck read timeout live."""
        port, shutdown, thread = _localhost_listener_that_never_accepts()
        try:
            client = BaseClient(settings=test_settings, base_url=f"http://127.0.0.1:{port}")
            with pytest.raises(Exception) as exc_info:
                client.healthcheck(timeout=2.0)
        finally:
            shutdown.set()
            thread.join(timeout=2.0)
        msg = str(exc_info.value).lower()
        assert "api request failed" in msg
        assert "read" in msg and "timed out" in msg

    def test_make_request_connect_timeout_live(self, test_settings: Settings) -> None:
        """Test make request connect timeout live."""
        client = BaseClient(settings=test_settings, base_url="http://192.0.2.1")
        with pytest.raises(Exception) as exc_info:
            client._make_request("GET", "/", timeout=1.5)
        msg = str(exc_info.value).lower()
        assert "api request failed" in msg
        if "timed out" not in msg and "timeout" not in msg:
            pytest.skip("connect did not time out as expected for 192.0.2.1; "
                        f"got {exc_info.value!r}",)
