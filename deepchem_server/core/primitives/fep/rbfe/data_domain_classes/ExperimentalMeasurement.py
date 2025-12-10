import logging
from typing import Union

import pint


class ExperimentalMeasurement:
    """DataClass to store ExperimentalMeasurement Values

    Parameters
    ----------

    value: pint.Quantity
        The measured experimental value

    uncertainty: pint.Quantity
        The uncertainty in the measured experimental quantity
    """

    def __init__(self, value: Union[str, pint.Quantity], uncertainty: Union[str, pint.Quantity]):
        try:
            self.value = get_validated_input(value)
        except Exception as e:
            logging.error(f"Value {value} not a valid pint Quantity!")
            raise e

        try:
            self.uncertainty = get_validated_input(uncertainty)
        except Exception as e:
            logging.error(f"Uncertainty{uncertainty} not a valid pint Quantity!")
            raise e


def get_validated_input(quantity: Union[str, pint.Quantity]):  # type: ignore
    """Validates the input Quantity String

    If the input is a pint.Quantity, it returns it as is, otherwise it converts it to a pint.Quantity object and then returns it.

    Parameters
    ----------
    quantity_string : Union[str, pint.Quantity]
        The input quantity object.

    Returns
    -------
    quantity
        The validated Quantity

    Raises
    ------
    TypeError
        If the input is invalid
    """
    if isinstance(quantity, type(pint.Quantity('1 kcal/mol'))):  # type: ignore
        quantity = quantity
    elif isinstance(quantity, str):
        quantity = pint.Quantity(quantity)  # type: ignore
    else:
        raise TypeError()

    return quantity
