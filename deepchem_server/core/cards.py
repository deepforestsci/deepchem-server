import ast
import json
import datetime
from typing import Optional, Dict


class Card:
    """Baseclass for cards"""

    def __init__(self):
        self.last_updated_time = datetime.datetime.now().strftime(
            '%d-%B-%Y %H:%M:%S')

    def __bytes__(self):
        return bytes(self.to_json(), encoding='utf8')

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__)

    def update_card(self, key, value):
        setattr(self, key, value)


class DataCard(Card):
    """Class for storing data card attributes

    Parameters
    ----------
    address: str
      Address of the reference object in the datastore.
    file_type: str
      The file extension - ex. csv filetype, .json file type etc
    data_type: str
      The type of object stored at the location pointef by filename - ex: pd.DataFrame, dast.dataframe.DataFrame
    description: str
      A description about the datastore
    featurizer: str
      The featurizer used in the dataset
    intended_use: str
      Notes on dataset - the intended use of the dataset
    caveats: str
      Notes on dataset - the caveats in using the dataset.
    feat_kwargs: dict
      Keyword arguments for featurizer (used when featurizer is not None)

    Notes
    -----
    1. Difference between data_type and file_type:
      An example can illustrate this better. A csv file (file_type) can either be
      a pandas.DataFrame or dask.dataframe.DataFrame or just a csv file. The file_type
      holds the file extension ('csv') while data_type refers to the data
      object (pandas.DataFrame, dask.dataframe.DataFrame, etc)
    """
    SUPPORTED_DATA_TYPES = [
        'pandas.DataFrame',
        'dc.data.NumpyDataset',
        'dc.data.DiskDataset',
        'json',
        'text/plain',
        'png',
        'binary',
    ]
    SUPPORTED_FILE_TYPES = [
        'csv', 'dir', 'json', 'pdb', 'fasta', 'fastq', 'png', 'sdf', 'dcd',
        'txt', 'xml', 'py', 'pdbqt', 'zip', 'smi', 'smiles', 'bz2', 'cxsmiles',
        'onnx', 'hdf5', 'log'
    ]

    def __init__(self,
                 address: str,
                 file_type: str,
                 data_type: str,
                 shape=None,
                 description: Optional[str] = None,
                 featurizer: Optional[str] = None,
                 intended_use: Optional[str] = None,
                 caveats: Optional[str] = None,
                 feat_kwargs: Optional[Dict] = None,
                 **kwargs):
        super().__init__()
        data_type = self.validate_datatype(data_type)
        if not isinstance(address, str):
            raise TypeError("address must be a string")
        if not isinstance(file_type, str):
            raise TypeError("file_type must by a string")
        if not isinstance(data_type, str):
            raise TypeError("data_type must be a string")
        assert file_type in self.SUPPORTED_FILE_TYPES, 'Filetype {} is not supported. Supported data type are {}'.format(
            file_type, self.SUPPORTED_FILE_TYPES)
        self.shape = shape
        self.address = address
        self.file_type = file_type
        self.data_type = data_type
        self.description = description
        self.featurizer = featurizer
        self.intended_use = intended_use
        self.caveats = caveats
        self.feat_kwargs = feat_kwargs
        # FIXME we should not depend on kwargs as internal parameters. These are features
        # for users to store additional details.
        for key, value in kwargs.items():
            setattr(self, key, value)

    def validate_datatype(self, data_type):
        # expand data type name
        if data_type == 'DataFrame':
            data_type = 'pandas.DataFrame'
        elif data_type == 'DiskDataset':
            data_type = 'dc.data.DiskDataset'
        elif data_type == 'NumpyDataset':
            data_type = 'dc.data.NumpyDataset'
        assert data_type in self.SUPPORTED_DATA_TYPES, 'Datatype {} is not supported. Supported file type are {}'.format(
            data_type, self.SUPPORTED_DATA_TYPES)
        return data_type

    @classmethod
    def from_json(cls, json_data: str):
        args = json.loads(json_data)
        return cls(**args)
        # Note: The above method may fail if `DataCard` contains nested objects.
        # return json.loads(json_data, object_hook=lambda d: DataCard(**d))

    @classmethod
    def from_bytes(cls, card_bytes: bytes):
        return DataCard.from_json(card_bytes.decode('utf8'))

    def get_n_samples(self):
        if self.shape is None:
            raise ValueError("the dataset does not have shape")
        if self.data_type == 'pandas.DataFrame':
            return self.shape[0]  # shape of dataframe is (n_rows x n_cols)
        elif self.data_type in ['dc.data.DiskDataset', 'dc.data.NumpyDataset']:
            # A deepchem dataset has shape (X_shape, y_shape, w_shape, ids_shape)
            x_shape = self.shape[0]
            return x_shape[0]
        else:
            return self.shape

    def to_json(self):
        json_str = json.dumps(self, default=lambda o: o.__dict__)
        return json_str.replace('_shape', 'shape')

    @property
    def shape(self):
        # FIXME This might pose security risk if user
        # arbitarily sets card shape
        return ast.literal_eval(self._shape)

    @shape.setter
    def shape(self, value):
        self._shape = str(value)


class ModelCard(Card):
    pass
