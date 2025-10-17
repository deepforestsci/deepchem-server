Configuration
=============

The pyds library uses a centralized configuration system through the Settings class, which persists configuration data to a JSON file.

Settings File Format
--------------------

The Settings class persists configuration to a JSON file (``.pyds.settings.json``) with the following structure:

.. code-block:: json

   {
     "profile": "my_profile",
     "project": "my_project",
     "base_url": "http://localhost:8000",
     "additional_settings": {
       "custom_key": "custom_value"
     }
   }

Configuration Methods
---------------------

Settings can be configured programmatically or through the JSON file:

* **Programmatic**: Use ``set_profile()``, ``set_project()``, ``set_base_url()``
* **File-based**: Edit ``.pyds.settings.json`` directly
* **Mixed**: Load from file and override programmatically

Error Handling
--------------

The pyds library includes comprehensive error handling:

* **Configuration Errors**: Raised when required settings (profile/project) are missing
* **File Errors**: Raised when upload files don't exist or are inaccessible
* **API Errors**: Raised when HTTP requests fail or return error responses
* **Validation Errors**: Raised when parameters don't meet requirements

All errors include descriptive messages to help with debugging and troubleshooting.

Integration with DeepChem Server
--------------------------------

The pyds library interfaces with the following DeepChem Server endpoints:

**Data Endpoints:**

* ``POST /data/uploaddata``: Upload data to datastore

**Primitive Endpoints:**

* ``POST /primitive/featurize``: Submit featurization jobs
* ``POST /primitive/train``: Submit training jobs
* ``POST /primitive/evaluate``: Submit evaluation jobs
* ``POST /primitive/infer``: Submit inference jobs
* ``POST /primitive/train-valid-test-split``: Submit splitting jobs

**System Endpoints:**

* ``GET /healthcheck``: Check server health status

For interactive API testing and the most up-to-date endpoint documentation, visit http://localhost:8000/docs when your DeepChem Server is running.
