from typing import Type, Union

import pint  # type: ignore


class EdgeSimulationResult:
    """Dataclass to store the result of a simulation for a single edge.

    Parameters
    ----------
    componentA_name : str
        The name of the first SmallMoleculeComponent in the edge.
    componentB_name : str
        The name of the second SmallMoleculeComponent in the edge.
    is_dry_run : bool
        Whether the simulation was run in dry run mode.
    dry_run_status : str
        The status of the dry run upon completion (Success or Failure).
    complex_dG : pint.Quantity
        The change in free energy of the protein-ligand-complex when it is formed from its constituent parts in the gas phase
    complex_dG_uncertainty : pint.Quantity
        The calculation uncertainty in complex_dG.
    solvent_dG : pint.Quantity
        The change in free energy of the protein-ligand-complex when it is solvated in a solvent.
    solvent_dG_uncertainty : pint.Quantity
        The calculation uncertainty in solvent_dG.
    edge_ddG : pint.Quantity
        The difference between the complex_dG and the solvent_dG.


    References
    ----------
        1. Relative Binding Free Energy Calculations in Drug Discovery: Recent Advances and Practical Considerations - https://doi.org/10.1021/acs.jcim.7b00564
        2. Best Practices for Alchemical Free Energy Calculations - https://livecomsjournal.org/index.php/livecoms/article/view/v2i1e18378
    """

    def __init__(
            self,
            componentA_name: str,
            componentB_name: str,
            is_dry_run: bool,
            dry_run_status: Union[str, None] = None,
            complex_dG: Type[pint.Quantity] = None,  # type: ignore
            complex_dG_uncertainty: Type[pint.Quantity] = None,  # type: ignore
            solvent_dG: Type[pint.Quantity] = None,  # type: ignore
            solvent_dG_uncertainty: Type[pint.Quantity] = None,  # type: ignore
            edge_ddG: Type[pint.Quantity] = None,  # type: ignore
            edge_ddG_uncertainty: Type[pint.Quantity] = None  # type: ignore
    ):
        self.componentA_name = componentA_name
        self.componentB_name = componentB_name
        self.is_dry_run = is_dry_run
        self.dry_run_status = dry_run_status
        self.complex_dG = (complex_dG or '').__str__()  # type: ignore
        self.complex_dG_uncertainty = (complex_dG_uncertainty or '').__str__()  # type: ignore
        self.solvent_dG = (solvent_dG or '').__str__()  # type: ignore
        self.solvent_dG_uncertainty = (solvent_dG_uncertainty or '').__str__()  # type: ignore
        self.edge_ddG = (edge_ddG or '').__str__()  # type: ignore
        self.edge_ddG_uncertainty = (edge_ddG_uncertainty or '').__str__()  # type: ignore
