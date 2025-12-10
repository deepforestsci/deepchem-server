from typing import Type

from gufe.protocols.protocoldag import ProtocolDAG  # type: ignore
from openfe import SmallMoleculeComponent  # type: ignore


class RunnableEdge:
    """A RunnableEdge containing information needed to run a single-edge RBFE simulation.

    For each edge in the generated perturbation network, we need to create 4 systems.
        1. Component A complex
        2. Component A solvent
        3. Component B complex
        4. Component B solvent

    Using these 4 systems, two alchemical transformations (represented as ProtocolDAG Objects) are created.
        1. Complex DAG: Component A complex -> Component B complex
        2. Solvent DAG: Component A solvent -> Component B solvent

    The RunnableEdge object encapsulates all of this information for a single edge in the perturbation network.
    Using the RelativeHybridTopologyProtocol object, the two DAGs can then be used to run the RBFE simulations.

    References
    ----------
        https://docs.openfree.energy/en/stable/reference/index.html
    """

    def __init__(self, componentA: Type[SmallMoleculeComponent], componentB: Type[SmallMoleculeComponent],
                 solvent_dag: Type[ProtocolDAG], complex_dag: type[ProtocolDAG]):
        """
        Parameters
        ----------
        componentA : Type[SmallMoleculeComponent]
            The first component of the edge.
        componentB : Type[SmallMoleculeComponent]
            The second component of the edge.
        solvent_dag : Type[ProtocolDAG]
            The solvent dag for the edge generated from the rbfe_transform.create() method
        complex_dag : Type[ProtocolDAG]
            The complex dag generated for the edge from the rbfe_transform.create() method
        """
        self.componentA = componentA
        self.componentB = componentB
        self.solvent_dag = solvent_dag
        self.complex_dag = complex_dag
