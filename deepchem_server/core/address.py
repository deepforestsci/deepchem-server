import os
from deepchem_server.core import config
from typing import Optional

DEEPCHEM_ADDRESS_PREFIX = 'deepchem://'


class DeepchemAddress(object):
    """A uniform representation to refer deepchem Objects.

    DeepchemAddress provides access to storage location of the
    object by inferring it from the DeepchemAddress provided.

    Parameters
    ----------
    address : str
        The address of the object.
    kind : str, optional
        The kind of object. Can be 'data' or 'model', by default 'data'.

    Examples
    --------
    >>> address = 'deepchem://profile/project/zinc.csv'
    >>> deepchem_address = DeepchemAddress(address)
    """
    address_prefix: str = 'deepchem://'

    def __init__(self, address: str, kind: Optional[str] = "data") -> None:
        """Initialize DeepchemAddress.

        Parameters
        ----------
        address : str
            The address of the object.
        kind : str, optional
            The kind of object, by default 'data'.
        """
        if address.startswith(DeepchemAddress.address_prefix):
            self.address = address
        else:
            self.address = DeepchemAddress.address_prefix + address
        self.kind = kind
        parsed_address = DeepchemAddress.parse_address(self.address)
        self.profile = parsed_address['profile']
        self.project = parsed_address['project']
        self.key = parsed_address['key']

    @classmethod
    def make_deepchem_address_from_filename(cls, end: str) -> str:
        """Return a deepchem address string from a filename.

        Parameters
        ----------
        end : str
            The filename whose DeepchemAddress we are creating.

        Returns
        -------
        str
            The DeepchemAddress of the file in the format
            deepchem://<storage_loc>/<end>.

        Raises
        ------
        ValueError
            If no datastore is configured.

        Examples
        --------
        >>> DeepchemAddress.make_deepchem_address_from_filename('temp.txt')
        deepchem://test_company/test_user/working_dir/temp.txt
        """
        datastore = config.get_datastore()
        if datastore is None:
            raise ValueError("No datastore configured")
        return DeepchemAddress.address_prefix + os.path.join(
            datastore.storage_loc, end)

    @classmethod
    def get_key(cls, address: str) -> str:
        """Return the key from an address.

        A key is used to refer to one of DeepChem's dataset or model.

        Parameters
        ----------
        address : str
            The address string whose key we are extracting.

        Returns
        -------
        str
            The extracted key from the address.

        Examples
        --------
        The following are all examples for different formats of the same address

        Example 1:
        ----------
        >>> dataset_address = 'deepchem://deepchem/data/delaney'
        >>> key = DeepchemAddress.get_key(dataset_address)
        >>> key
        delaney

        Example 2:
        ----------
        >>> dataset_address = 'deepchem/data/delaney'
        >>> key = DeepchemAddress.get_key(dataset_address)
        >>> key
        deepchem/data/delaney

        Example 3:
        ----------
        >>> dataset_address = 'delaney'
        >>> key = DeepchemAddress.get_key(dataset_address)
        delaney
        """
        if address.startswith(DeepchemAddress.address_prefix):
            address = address[len(DeepchemAddress.address_prefix):]
            return '/'.join(address.split('/')[2:])
        return address

    @classmethod
    def parse_address(cls, address: str) -> dict:
        """Return different components of the address.

        Parameters
        ----------
        address : str
            The deepchem address of the object.

        Returns
        -------
        dict
            Dictionary containing 'profile', 'project', and 'key' components.

        Raises
        ------
        ValueError
            If the address format is invalid.

        Examples
        --------
        >>> address = 'deepchem://user/test/file'
        >>> parsed_address = DeepchemAddress.parse_address(address)
        >>> parsed_address
        {'profile': 'user', 'project': 'test', 'key': 'file'}
        """
        if address.startswith(DeepchemAddress.address_prefix):
            address = address[len(DeepchemAddress.address_prefix):]
        tokens = address.split('/')
        if len(tokens) < 3:
            raise ValueError('Invalid deepchem address')
        parsed_address = dict()
        parsed_address['profile'] = tokens[0]
        parsed_address['project'] = tokens[1]
        parsed_address['key'] = '/'.join(tokens[2:])
        return parsed_address

    @classmethod
    def get_path(cls,
                 storage_loc: str,
                 address: str,
                 format: Optional[str] = 's3',
                 base_dir: Optional[str] = None) -> str:
        """Return the path of the object in the storage from the address.

        When the format is ``local``, the ``base_dir`` is used as the base
        directory and ensures that the path returned matches the OS path format.

        Parameters
        ----------
        storage_loc : str
            The storage location of the object (used in case the address is not
            in default deepchem address format).
        address : str
            The deepchem address of the object.
        format : {'s3', 'local'}, optional
            The format of the path to be returned, by default 's3'.
        base_dir : str, optional
            The base directory to be used in case of 'local' format.

        Returns
        -------
        str
            The path of the object in the specified format.

        Raises
        ------
        ValueError
            If the format is not 's3' or 'local'.

        Examples
        --------
        All the following examples return the same path - profile/project/key

        Example 1:
        ----------
        >>> address = 'deepchem://profile/project/key'
        >>> storage_loc = 'profile/project'
        >>> path = DeepchemAddress.get_path(storage_loc, address)
        >>> path
        profile/project/key

        Example 2:
        ----------
        >>> address = 'profile/project/key'
        >>> storage_loc = 'profile/project'
        >>> path = DeepchemAddress.get_path(storage_loc, address)
        >>> path
        profile/project/key

        Example 3:
        ----------
        >>> address = 'key'
        >>> storage_loc = 'profile/project'
        >>> path = DeepchemAddress.get_path(storage_loc, address)
        >>> path
        profile/project/key
        """
        try:
            # Address is of the form deepchem://profile/project/key
            address_parsed = DeepchemAddress.parse_address(address)
            profile = address_parsed["profile"]
            project = address_parsed["project"]
            key = address_parsed["key"]
            if format == 's3':
                return profile + '/' + project + '/' + key
            elif format == 'local':
                key = key.replace('/', os.sep)
                if base_dir is not None:
                    return os.path.join(base_dir, profile, project, key)
                else:
                    return os.path.join(profile, project, key)
        except ValueError:
            # Address is not in the form deepchem://profile/project/key
            if not address.startswith(storage_loc):
                address = DEEPCHEM_ADDRESS_PREFIX + storage_loc.strip(
                    '/') + '/' + address
            address_key = DeepchemAddress.get_key(address)
            if format == 's3':
                return storage_loc + address_key
            elif format == 'local':
                address_key = address_key.replace('/', os.sep)
                return os.path.join(storage_loc, address_key)
        # if the format is neither s3 nor local
        raise ValueError(f"Unknown format: {format}")

    @classmethod
    def get_parent_key(cls, address: str) -> str:
        """Return the parent key of the object.

        Parameters
        ----------
        address : str
            The deepchem address of the object or the key of the object.

        Returns
        -------
        str
            The parent key path.

        Examples
        --------
        >>> address = 'deepchem://profile/project/parent1/parent2/key'
        >>> parent_key = DeepchemAddress.get_parent_key(address)
        >>> parent_key
        parent1/parent2

        >>> address = 'profile/project/parent1/parent2/key'
        >>> parent_key = DeepchemAddress.get_parent_key(address)
        >>> parent_key
        parent1/parent2
        """
        object_key = DeepchemAddress.get_key(address)
        if '/' not in object_key:
            return ''
        return '/'.join(object_key.split('/')[:-1]) + '/'

    @classmethod
    def get_object_name(cls, address: str) -> str:
        """Return the name of the object.

        Parameters
        ----------
        address : str
            The deepchem address of the object or the key of the object.

        Returns
        -------
        str
            The object name.

        Examples
        --------
        >>> address = 'deepchem://profile/project/parent1/parent2/key'
        >>> object_name = DeepchemAddress.get_object_name(address)
        >>> object_name
        key

        >>> address = 'profile/project/parent1/parent2/key'
        >>> object_name = DeepchemAddress.get_object_name(address)
        >>> object_name
        key
        """
        object_key = DeepchemAddress.get_key(address)
        if '/' not in object_key:
            return object_key
        return object_key.split('/')[-1]

    def __str__(self) -> str:
        """Return string representation of the address.

        Returns
        -------
        str
            The address string.
        """
        return self.address

    def __repr__(self) -> str:
        """Return string representation of the address.

        Returns
        -------
        str
            The address string.
        """
        return self.address
