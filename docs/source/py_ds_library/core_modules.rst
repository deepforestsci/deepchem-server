Core Modules
============

The pyds library consists of three core modules that provide the foundation for all client operations.

Settings Class
--------------

The Settings class manages persistent configuration for profile, project, and server settings.

.. autoclass:: pyds.settings.Settings
   :members:
   :undoc-members:

BaseClient Class
----------------

The BaseClient class provides common functionality for all API clients, including HTTP requests, configuration validation, and health checks.

.. autoclass:: pyds.base.client.BaseClient
   :members:
   :undoc-members:

Data Class
----------

The Data class handles data management operations, primarily file uploads to the datastore.

.. autoclass:: pyds.data.Data
   :members:
   :undoc-members:
