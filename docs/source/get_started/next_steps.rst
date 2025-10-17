Next Steps
----------

1. **Installation**: Set up DeepChem Server on your system
2. **Setup Python Client**: Install the pyds library (see :doc:`PyDS Library Docs </py_ds_library/getting_started>`)
3. **Quick Start**: Run your first molecular featurization job

**Key Endpoints**

Once running, you'll primarily work with the following endpoints:

* ``POST /data/uploaddata``: Upload molecular datasets
* ``POST /primitive/featurize``: Apply molecular featurization
* ``POST /primitive/train``: Train machine learning models
* ``POST /primitive/evaluate``: Evaluate model performance
* ``POST /primitive/infer``: Run inference on new data
* ``GET /healthcheck``: Check server status