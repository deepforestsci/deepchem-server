API Reference
=============

Complete reference documentation for the DeepChem Server API, including REST endpoints and core modules.

.. toctree::
   :maxdepth: 2

   rest_api
   core_modules
   utils

Overview
--------

DeepChem Server provides a comprehensive API for molecular machine learning workflows. The API is organized into several layers:

* **REST API**: HTTP endpoints for data upload, featurization, and retrieval
* **Core Modules**: Python modules that implement the server functionality  
* **Utilities**: Helper functions and utilities used throughout the system

For interactive API testing and the most up-to-date endpoint documentation, visit http://localhost:8000/docs when your server is running.

Authentication
--------------

Currently, DeepChem Server does not require authentication. All endpoints are publicly accessible when the server is running.

.. note::
   In production deployments, you should implement appropriate authentication and authorization mechanisms.

Error Handling
--------------

All API endpoints return consistent error responses with appropriate HTTP status codes. Common status codes include:

* **200 OK**: Request succeeded
* **400 Bad Request**: Invalid request parameters
* **404 Not Found**: Resource not found
* **422 Unprocessable Entity**: Validation error
* **500 Internal Server Error**: Server error

For detailed error response formats and examples, visit http://localhost:8000/docs. 