REST API
========

The DeepChem Server REST API provides HTTP endpoints for all major operations. For interactive documentation, examples, and testing, visit http://localhost:8000/docs when your server is running.

Base URL
--------

The default base URL is:

.. code-block:: text

   http://localhost:8000

All endpoints are relative to this base URL.

Interactive Documentation
-------------------------

**Swagger UI**: http://localhost:8000/docs

This provides:

* Complete endpoint documentation with interactive testing
* Request/response examples
* Parameter descriptions and validation
* Authentication details
* Schema definitions

**ReDoc**: http://localhost:8000/redoc

Alternative documentation interface with:

* Clean, readable documentation
* Detailed schema information
* Comprehensive API overview

Router Module Documentation
---------------------------

For detailed implementation information, the router endpoints are documented using auto-generated docstrings:

.. automodule:: deepchem_server.routers.data
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: deepchem_server.routers.primitives
   :members:
   :undoc-members:
   :show-inheritance:

Available Featurizers
---------------------

DeepChem Server supports the following molecular featurizers:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Featurizer
     - DeepChem Class
     - Description
   * - ``ecfp``
     - ``CircularFingerprint``
     - Extended Connectivity Fingerprints (Morgan/ECFP)
   * - ``graphconv``
     - ``ConvMolFeaturizer``
     - Graph convolution molecular featurizer
   * - ``weave``
     - ``WeaveFeaturizer``
     - Weave molecular featurizer for graph networks
   * - ``molgraphconv``
     - ``MolGraphConvFeaturizer``
     - Molecular graph convolution featurizer

For detailed information about these featurizers, including parameters, usage examples, and implementation details, please refer to the `DeepChem Featurizers documentation <https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html>`_.

Each featurizer accepts different parameters that can be passed in the ``feat_kwargs`` field of your featurization requests. Consult the DeepChem documentation for the specific parameters supported by each featurizer.

Response Formats
----------------

All API responses follow consistent JSON formats as documented in the interactive API documentation at http://localhost:8000/docs.

Error Handling
--------------

Common HTTP status codes:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Code
     - Status
     - Description
   * - 200
     - OK
     - Request succeeded
   * - 400
     - Bad Request
     - Invalid request parameters
   * - 404
     - Not Found
     - Resource not found
   * - 422
     - Unprocessable Entity
     - Validation error
   * - 500
     - Internal Server Error
     - Server error occurred

For detailed error response formats and examples, visit http://localhost:8000/docs.

FastAPI Configuration
---------------------

DeepChem Server is built with FastAPI. For server configuration, deployment, and advanced settings, refer to the `FastAPI documentation <https://fastapi.tiangolo.com/>`_.

Best Practices
--------------

1. **Use the Interactive Documentation**: Visit http://localhost:8000/docs for the most up-to-date API information
2. **Test Endpoints Interactively**: Use the "Try it out" feature in the Swagger UI
3. **Check Response Schemas**: Review the response examples in the documentation
4. **Handle Errors Gracefully**: Implement proper error handling in your client code
5. **Consult DeepChem Documentation**: For detailed featurizer information, refer to the `DeepChem Featurizers documentation <https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html>`_

For detailed examples, parameter descriptions, and interactive testing, always refer to http://localhost:8000/docs. 