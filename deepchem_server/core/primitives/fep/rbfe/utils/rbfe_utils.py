import json
from typing import Type

from openfe import SolventComponent  # type: ignore
from openfe.protocols.openmm_rfe.equil_rfe_settings import RelativeHybridTopologyProtocolSettings  # type: ignore
import pint  # type: ignore

from deepchem_server.core.primitives.fep.rbfe.system_setup import get_default_RBFE_simulation_settings


PINT_QUANTITY_TYPE = type(pint.Quantity("1.0 kcal/mol"))  # type: ignore


class RBFESettingsUtils:
    """RBFESettingsUtility class containing utility methods for RBFESettings objects.

    The RBFESettingsUtils class provides several utility methods for the ease of handling of RBFESettings objects.
    For example, this includes methods for serializing and deserializing RBFESettings objects to and from JSON strings,
    which are primarily used to pass RBFESettings objects between the frontend and backend.
    """

    class _RBFESettingsEncoder(json.JSONEncoder):
        """Custom JSON encoder for RBFESettings objects.
        """

        def default(self, obj: Type[RelativeHybridTopologyProtocolSettings]):
            """Overridden default method of JSONEncoder.

            Due to non-serializable nature of pint.Quantity objects, this method is overridden to serialize pint.Quantity objects to their string representations.
            """
            if isinstance(obj, pint.Quantity):  # type: ignore
                return obj.__str__()
            else:
                try:
                    return json.JSONEncoder.default(self, obj)
                except TypeError:
                    return obj.__dict__

    class _RBFESettingsDecoder(json.JSONDecoder):
        """Custom JSON decoder for RBFESettings objects.
        """

        def decode(self, settings_json: str):  # type: ignore
            """Overridden decode method of JSONDecoder.

            Due to non-serializable nature of pint.Quantity objects, this method is overridden to deserialize pint.Quantity objects from their string representations.
            """
            settings_dict = json.loads(settings_json)

            rbfe_settings = get_default_RBFE_simulation_settings()

            for setting in settings_dict.keys():
                for setting_parameter in settings_dict[setting].keys():

                    setting_attribute = getattr(getattr(rbfe_settings, setting), setting_parameter)

                    if type(setting_attribute) is PINT_QUANTITY_TYPE:
                        setattr(
                            (getattr(rbfe_settings, setting)),
                            setting_parameter,
                            pint.Quantity(  # type: ignore
                                settings_dict[setting][setting_parameter]))
                    else:
                        setattr((getattr(rbfe_settings, setting)), setting_parameter,
                                settings_dict[setting][setting_parameter])

            return rbfe_settings

    @staticmethod
    def loads(settings_json: str) -> Type[RelativeHybridTopologyProtocolSettings]:
        """Loads a RBFESettings object from a JSON string.

        For a comprehensive list of parameters that can be loaded, please refer to the OpenFE API reference https://docs.openfree.energy/en/stable/reference/api/openmm_rfe.html#protocol-settings

        Parameters
        ----------
        settings_json : str
            JSON string containing the RBFESettings object.

            The JSON is not required to encompass all fields within the RBFESettings object. However,
            it should encompass those specific fields of the RBFESettings object which deviate from the default configurations.

        Returns
        -------
        rbfe_settings : Type[RelativeHybridTopologyProtocolSettings]
            The deserialized RBFESettings object loaded from the JSON string.

        Example
        -------
        >>> # An rbfe_settings JSON string with two overridden parameters.
        >>> overridden_rbfe_settings = '{"simulation_settings" : {"equilibration_length" : "2 ns", "production_length" : "20 ns"}}'
        >>> rbfe_settings = RBFESettingsUtils.loads(overridden_rbfe_settings)
        """
        rbfe_settings = json.loads(settings_json, cls=RBFESettingsUtils._RBFESettingsDecoder)
        return rbfe_settings

    @staticmethod
    def dumps(rbfe_settings: Type[RelativeHybridTopologyProtocolSettings]) -> str:
        """Serializes a RBFESettings object to a JSON string.

        Parameters
        ----------
        rbfe_settings : Type[RelativeHybridTopologyProtocolSettings]
            The RBFESettings object to be serialized.

        Returns
        -------
        settings_json : str
            The serialized JSON representation of the RelativeHybridTopologyProtocolSettings object.

        Example
        -------
        >>> from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
        >>> rbfe_settings = RelativeHybridTopologyProtocol.default_settings()
        >>> rbfe_settings_json = RBFESettingsUtils.dumps(rbfe_settings)
        """
        if not isinstance(rbfe_settings, RelativeHybridTopologyProtocolSettings):
            raise TypeError("rbfe_settings must be a RelativeHybridTopologyProtocolSettings object")
        settings_json = json.dumps(rbfe_settings, cls=RBFESettingsUtils._RBFESettingsEncoder)
        return settings_json


class SolventComponentUtils:
    """SolventComponentUtils class containing utility methods for SolventComponent objects.

    The SolventComponentUtils class provides several utility methods for the ease of handling of SolventComponent objects.
    For example, this includes methods for serializing and deserializing SolventComponent objects to and from JSON strings,
    which are primarily used to pass SolventComponent objects between the frontend and backend.
    """

    @staticmethod
    def loads(solvent_json: str) -> Type[SolventComponent]:
        """Loads a SolventComponent object from a JSON string.

        The JSON string must contain the following keys: "positive_ion", "negative_ion", "neutralize", "ion_concentration".
        The expected values for the keys are as follows:

        positive_ion : str
            The positive ion of the solvent.
        negative_ion : str
            The negative ion of the solvent.
        neutralize : bool
            Whether to neutralize the solvent.
        ion_concentration : str
            The ion concentration of the solvent specified as a parsable pint.Quantity.

        Parameters
        ----------
        solvent_json : str
            JSON string containing the SolventComponent object.

        Returns
        -------
        solvent : Type[SolventComponent]
            The deserialized SolventComponent object loaded from the JSON string.

        Example
        -------
        >>> # Please note the keys in the JSON string.
        >>> solvent_json = '{"positive_ion": "Na", "negative_ion": "Cl", "neutralize": true, "ion_concentration": "0.15 mol/l"}'
        >>> solvent = SolventComponentUtils.loads(valid_json)
        """
        solvent_dict = json.loads(solvent_json)

        solvent = SolventComponent(
            positive_ion=solvent_dict['positive_ion'],
            negative_ion=solvent_dict['negative_ion'],
            ion_concentration=pint.Quantity(  # type: ignore
                solvent_dict['ion_concentration']),
            neutralize=solvent_dict['neutralize'])

        return solvent

    @staticmethod
    def dumps(solvent: Type[SolventComponent]) -> str:
        """Serializes a SolventComponent object to a JSON string.

        Parameters
        ----------
        solvent : Type[SolventComponent]
            The SolventComponent object to be serialized.

        Returns
        -------
        solvent_json
            The serialized JSON representation of the SolventComponent object.

        Example
        -------
        >>> solvent = SolventComponent(positive_ion='Na',
                               negative_ion='Cl',
                               neutralize=True,
                               ion_concentration=0.15 * unit.mole / unit.litre)
        >>> solvent_json = SolventComponentUtils.dumps(solvent)
        """
        if not isinstance(solvent, SolventComponent):
            raise TypeError("solvent must be a SolventComponent object")
        solvent_dict = {
            'positive_ion': solvent.positive_ion,
            'negative_ion': solvent.negative_ion,
            'ion_concentration': solvent.ion_concentration.__str__(),
            'neutralize': solvent.neutralize
        }

        solvent_json = json.dumps(solvent_dict)
        return solvent_json
