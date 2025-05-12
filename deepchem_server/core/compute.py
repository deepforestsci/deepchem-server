import deepchem_server.core as core

program_map = {
    'featurize': core.featurize
}

class ComputeWorkflow:
    """A Compute Workflow is a workflow that runs on Deepchem Server.

    The compute workflow encapsulates the concept of a
    workflow on Deepchem Server. A workflow is a series of 
    Deepchem Server primitive commands that update the state of the
    underlying datastore. Users can specify workflows
    using a simple programming language and can execute
    workflows if their personal budgets permit.
    """

    def __init__(self, program):
        """Initialize Worfklow."""
        self.program = program

    def execute(self):
        """Runs the set of commands in this workflow."""
        if 'program_name' not in self.program:
            raise ValueError("program_name not found in program")
        program_name = self.program['program_name']
        params = {key: value for key, value in self.program.items() if key != 'program_name'}
        if program_name not in program_map:
            raise ValueError(f"{program_name} not in available programs")

        output = program_map[program_name](**params)
        return output
