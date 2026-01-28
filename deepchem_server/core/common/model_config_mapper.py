import logging
from typing import Any, Dict, List, Literal, Optional

from _collections_abc import dict_keys, dict_values


logger = logging.getLogger(__name__)


class DeepChemModelConfigMapper:
    """Mappings between models and their configuration in Deepchem.

    This class contains mappings between the models and their configuration in
    Deepchem. It is used to generate the model cards while uploading
    models.

    The main purpose of this class is to validate and parse the model parameters
    from the config.yaml file and generate the model cards.

    **The config.yaml file contains the following parameters:**

    1. model_class (required): The model class in Deepchem.

    2. init_args (optional): The init arguments for the model.

    3. train_args (optional): The train arguments for the model.

    4. description (optional): The description of the model (will be stored in the model card).

    5. featurizer (optional): The featurizer for the model (will be stored in the model card).

    Sample config.yaml file:
    ------------------------
    .. code-block:: yaml

        model_class: GCNModel
        init_args:
          n_tasks: 1
          mode: classification
          batch: 2
          learning_rate: 0.0003
        train_args:
          nb_epoch: 1
        description: Description of the model (will be stored in the model card)

    Parameters
    ----------
    model_class : Any
        The model class in Deepchem.
    model_class_name : str, optional
        The name of the model class. If not provided, will be inferred.
    required_init_params : list, optional
        A list of required init parameters.
    optional_init_params : list, optional
        A list of optional init parameters.
    required_train_params : list, optional
        A list of required train parameters.
    optional_train_params : list, optional
        A list of optional train parameters.
    tasks : dict, optional
        A Dictionary of tasks mapped to their respective parameter name supported by the model.

    Examples
    --------
    >>> from deepchem_server.core.model_config_mapper import DeepChemModelConfigMapper
    >>> from deepchem.models import GCNModel
    >>> model = DeepChemModelConfigMapper(
    ... model_class=GCNModel,
    ... required_init_params=["init_param"],
    ... optional_init_params=["init_param1", "init_param2"],
    ... required_train_params=["train_param"],
    ... optional_train_params=["train_param1", "train_param2"])
    >>> model.get_model_class_name()
    'gcn'
    >>> model.get_model_class()
    <class 'deepchem.models.torch_models.gcn.GCNModel'>
    >>> model
    <class 'deepchem.models.torch_models.gcn.GCNModel'>
    >>> model.add_init_params(["test_required_init_param"])
    >>> model.get_init_params("required")
    ['init_param', 'test_required_init_param']
    >>> model.add_init_params(["test_optional_init_param"], "optional")
    >>> model.get_init_params("optional")
    ['init_param1', 'init_param2', 'test_optional_init_param']
    >>> model.get_init_params()
    {'required': ['init_param', 'test_required_init_param'],
    'optional': ['init_param1', 'init_param2', 'test_optional_init_param']}
    >>> model.add_tasks({"task1": "task", "task2": "mode"})
    >>> model.get_tasks()
    {'task1': 'task', 'task2': 'mode'}

    In the above example, the model tasks are mapped to their respective parameter name supported by the model.
    For example, the task "task1" is mapped to parameter "task" and the task so, during model initialization,
    if "task1" is provided as a task, then the parameter "task" will be used to initialize the model. Similarly,
    if "task2" is provided as a task, then the parameter "mode" will be used to initialize the model.
    """

    @staticmethod
    def parse_params(required_params: Optional[List], optional_params: Optional[List]) -> Dict:
        """Parse the required and optional parameters of the model.

        Returns a dictionary with the required and optional parameters.

        Parameters
        ----------
        required_params : list, optional
            A list of required parameters.
        optional_params : list, optional
            A list of optional parameters.

        Returns
        -------
        dict
            A dictionary containing the required and optional parameters.
        """
        if required_params is None:
            required_params = []
        if optional_params is None:
            optional_params = []
        return {
            'required': required_params,
            'optional': optional_params,
        }

    @staticmethod
    def get_class_name(model_class: Any) -> str:
        """Try to detect the model name for the model.

        Parameters
        ----------
        model_class : Any
            The model class.

        Returns
        -------
        str
            The model class name.
        """
        try:
            model_class_name = model_class.__name__
        except AttributeError:
            if model_class.__class__.__name__ == 'SklearnModel':
                model_class_name = model_class.model.__class__.__name__
            else:
                model_class_name = str(model_class.__class__)
                logger.error(f"{model_class.__class__}: Model class name not provided and could not be inferred.")
        return model_class_name

    def __init__(
        self,
        model_class: Any,
        model_class_name: Optional[str] = None,
        required_init_params: Optional[List] = None,
        optional_init_params: Optional[List] = None,
        required_train_params: Optional[List] = None,
        optional_train_params: Optional[List] = None,
        tasks: Optional[Dict] = None,
    ) -> None:
        """Initialize DeepChemModelConfigMapper."""
        if model_class_name is None:
            model_class_name = self.get_class_name(model_class)

        if tasks is None:
            tasks = {}

        self.model_config_mapping = {
            "model_class": model_class,
            "model_class_name": model_class_name,
            "init_params": self.parse_params(required_init_params, optional_init_params),
            "train_params": self.parse_params(required_train_params, optional_train_params),
            "tasks": tasks,
        }

    def add_init_params(self, init_params: List, kind: Literal["required", "optional"] = "required") -> None:
        """Add the init parameters to the model config mapping.

        Parameters
        ----------
        init_params : list
            A list of init parameters.
        kind : {'required', 'optional'}, optional
            Whether the init parameters are required or optional, by default 'required'.

        Returns
        -------
        None
        """
        self.model_config_mapping['init_params'][kind].extend(init_params)

    def add_train_params(self, train_params: List, kind: Literal["required", "optional"] = "required") -> None:
        """Add the train parameters to the model config mapping.

        Parameters
        ----------
        train_params : list
            A list of train parameters.
        kind : {'required', 'optional'}, optional
            Whether the train parameters are required or optional, by default 'required'.

        Returns
        -------
        None
        """
        self.model_config_mapping['train_params'][kind].extend(train_params)

    def add_tasks(self, tasks: Dict) -> None:
        """Add the tasks to the model config mapping.

        Parameters
        ----------
        tasks : dict
            A dictionary of tasks mapped to their respective parameter name supported by the model.

        Returns
        -------
        None
        """
        self.model_config_mapping['tasks'].update(tasks)

    def get_model_class(self) -> Any:
        """Return the model class for the model.

        Returns
        -------
        Any
            The model class for the model.
        """
        return self.model_config_mapping['model_class']

    def get_model_class_name(self) -> str:
        """Return the model class name for the model.

        Returns
        -------
        str
            The model class name for the model.
        """
        return self.model_config_mapping['model_class_name']

    def get_init_params(self, kind: Literal["required", "optional", None] = None) -> Dict:
        """Return the initialization parameters for the model.

        Parameters
        ----------
        kind : {'required', 'optional', None}, optional
            If kind is None, then the function returns all the init parameters
            for the model. If kind is "required", then the function returns
            only the required init parameters. If kind is "optional", then the
            function returns only the optional init parameters.

        Returns
        -------
        dict
            Returns a dictionary containing the init parameters for the model.
        """
        if kind is not None:
            return self.model_config_mapping['init_params'][kind]
        return self.model_config_mapping['init_params']

    def get_train_params(self, kind: Literal["required", "optional", None] = None) -> Dict:
        """Return the train parameters for the model.

        Parameters
        ----------
        kind : {'required', 'optional', None}, optional
            If kind is None, then the function returns all the train parameters
            for the model. If kind is "required", then the function returns
            only the required train parameters. If kind is "optional", then the
            function returns only the optional train parameters.

        Returns
        -------
        dict
            Returns a dictionary containing the train parameters for the model.
        """
        if kind is not None:
            return self.model_config_mapping['train_params'][kind]
        return self.model_config_mapping['train_params']

    def get_tasks(self) -> Dict:
        """Return the tasks for the model.

        Returns
        -------
        dict
            Returns a Dictionary containing the tasks mapped to their respective parameter name of the model.
        """
        return self.model_config_mapping['tasks']

    def __getitem__(self, item: str) -> Any:
        """Return the mentioned item from the model config mapping.

        Parameters
        ----------
        item : str
            The item to be returned from the model config mapping.

        Returns
        -------
        Any
            The item from the model config mapping.
        """
        return self.model_config_mapping[item]

    def __str__(self) -> str:
        """Return the model class name for the model.

        Returns
        -------
        str
            The model class name.

        Examples
        --------
        >>> from deepchem_server.core.model_config_mapper import DeepChemModelConfigMapper
        >>> from deepchem.models import GCNModel
        >>> model = DeepChemModelConfigMapper(
        ... model_class=GCNModel,
        ... required_init_params=["init_param"],
        ... optional_init_params=["init_param1", "init_param2"],
        ... required_train_params=["train_param"],
        ... optional_train_params=["train_param1", "train_param2"])
        >>> str(model)
        'GCNModel'
        """
        return self.model_config_mapping['model_class_name']

    def __repr__(self) -> Any:
        """Return the model class for the model.

        Returns
        -------
        Any
            The model class.

        Examples
        --------
        >>> from deepchem_server.core.model_config_mapper import DeepChemModelConfigMapper
        >>> from deepchem.models import GCNModel
        >>> model = DeepChemModelConfigMapper(
        ... model_class=GCNModel,
        ... required_init_params=["init_param"],
        ... optional_init_params=["init_param1", "init_param2"],
        ... required_train_params=["train_param"],
        ... optional_train_params=["train_param1", "train_param2"])
        >>> model
        <class 'deepchem.models.torch_models.gcn.GCNModel'>
        """
        return repr(self.get_model_class())


class ModelAddressWrapper(dict):
    """Wrapper for deepchem-server model name and deepchem model config.

    This class is used to wrap the deepchem-server model name and deepchem model config.
    It is used as a custom dictionary to map the deepchem-server model name to the
    deepchem model config.

    Examples
    --------
    >>> from deepchem_server.core.model_config_mapper import ModelAddressWrapper, DeepChemModelConfigMapper
    >>> from deepchem.models import GCNModel
    >>> model = DeepChemModelConfigMapper(
    ... model_class=GCNModel,
    ... required_init_params=["init_param"],
    ... optional_init_params=["init_param1", "init_param2"],
    ... required_train_params=["train_param"],
    ... optional_train_params=["train_param1", "train_param2"])
    >>> model_address_map = ModelAddressWrapper({"gcn": model})
    >>> model_address_map
    {'gcn': <class 'deepchem.models.torch_models.gcn.GCNModel'>}
    >>> model_address_map['gcn']
    <class 'deepchem.models.torch_models.gcn.GCNModel'>
    >>> model_address_map.get_model_class_name('gcn')
    'gcn'

    >>> # using key value pairs
    >>> from sklearn.linear_model import LinearRegression
    >>> from deepchem.models import SklearnModel
    >>> model = DeepChemModelConfigMapper(
    ... model_class=SklearnModel,
    ... required_init_params=None,
    ... optional_init_params=["fit_intercept", "copy_X", "n_jobs", "positive"],
    ... required_train_params=None,
    ... optional_train_params=None)
    >>> model_address_map['linear_regression'] = model
    >>> model_address_map['linear_regression']
    <class 'deepchem.models.sklearn_models.SklearnModel'>
    """

    def __init__(self, *args, **kwargs) -> None:
        """Initialize ModelAddressWrapper.

        Parameters
        ----------
        *args
            Variable length argument list. Expected dict as first argument.
        **kwargs
            Arbitrary keyword arguments for model mappings.

        Raises
        ------
        TypeError
            If more than 1 positional argument is provided or if the first
            argument is not a dict.
        """
        super().__init__()
        if args:
            if len(args) > 1:
                raise TypeError(f"ModelAddressWrapper expected at most 1 arguments, got {len(args)}")
            arg = args[0]  # type: ignore[index]
            if isinstance(arg, dict):
                for key, value in arg.items():
                    self.__setitem__(key, value)
            else:
                raise TypeError(f"ModelAddressWrapper expected dict, got {type(arg)}")
        if kwargs:
            for key, value in kwargs.items():
                self.__setitem__(key, value)

    def get_model_config(
            self,
            key: str,
            kind: Literal["model_name", "class_name"] = "model_name") -> Optional[DeepChemModelConfigMapper]:
        """Return the model config map given the model key.

        Parameters
        ----------
        key : str
            The name/key of the model.
        kind : {'model_name', 'class_name'}, optional
            Whether the key is the model name or the model class, by default 'model_name'.

        Returns
        -------
        DeepChemModelConfigMapper or None
            The model config map for the model, or None if not found.
        """
        if kind == "model_name":
            return self.__dict__[key]
        elif kind == "class_name":
            for model_name, model_config_map in self.__dict__.items():
                if model_config_map.get_model_class_name() == key:
                    return model_config_map
        return None

    def get_model_name_from_class_name(self, model_class_name: str) -> Optional[str]:
        """Return the model name for the model class name.

        The class will be used when parsing the config.yaml file,
        since we don't have the model name in the config.yaml file.

        Parameters
        ----------
        model_class_name : str
            The model class name for the model.

        Returns
        -------
        str or None
            The model name for the model, or None if not found.
        """
        for model_name, model_config_map in self.__dict__.items():
            if model_config_map.get_model_class_name() == model_class_name:
                return model_name
        return None

    def get_model_class_name(self, key: str) -> str:
        """Return the model class name for the model key.

        Since using a key to access the ModelAddressWrapper returns the model
        class, this function reduces the code complexity.

        The below code snippets are equivalent:

        >>> from deepchem.models import GCNModel
        >>> model_address_map = ModelAddressWrapper({"gcn": DeepChemModelConfigMapper(model_class=GCNModel)})
        >>> model_address_map.get_model_class_name("gcn")
        'GCNModel'

        >>> model_address_map.get_model_config("gcn").get_model_class_name()
        'GCNModel'

        Parameters
        ----------
        key : str
            The name/key of the model.

        Returns
        -------
        str
            The model class name for the model.
        """
        return self.__dict__[key].get_model_class_name()

    def get_model_class_names(self) -> List[str]:
        """Return the model class names for the models.

        Returns
        -------
        list of str
            The model class names for the models.
        """
        return [self.__dict__[key].get_model_class_name() for key in self.__dict__.keys()]

    def __setitem__(self, key: str, value: DeepChemModelConfigMapper) -> None:
        """Set item in the wrapper.

        Parameters
        ----------
        key : str
            The model name key.
        value : DeepChemModelConfigMapper
            The model config mapper to store.

        Returns
        -------
        None
        """
        self.__dict__[key] = value

    def __getitem__(self, key: str) -> Any:
        """Get item from the wrapper.

        Parameters
        ----------
        key : str
            The model name key.

        Returns
        -------
        Any
            The model config mapper.
        """
        return self.__dict__[key].get_model_class()

    def __contains__(self, key) -> bool:
        """Check if the key is in the wrapper.

        Parameters
        ----------
        key : str
            The model name key.

        Returns
        -------
        bool
            True if the key is in the wrapper, False otherwise.
        """
        return key in self.__dict__.keys()

    def keys(self) -> dict_keys:
        """Return the keys of the wrapper.

        Returns
        -------
        dict_keys
            The keys of the wrapper.
        """
        return self.__dict__.keys()

    def values(self) -> dict_values:
        """Return the values of the wrapper.

        Returns
        -------
        dict_values
            The values of the wrapper.
        """
        return self.__dict__.values()
