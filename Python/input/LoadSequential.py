from input.LoadModel import LoadModel
from input.LoadState import LoadState
def LoadSequential(orb):
    """
    This function defines the mission's sequence by loading the model,
    initial state, and defining the mission's phases.

    The sequence types can be among the following:
    - "Propag": Free propagation
    - "DVPropag": Propagation after an impulsive boost
    - "TBPOptim": Optimize initial state for a 3-body problem orbit
    - "Lambert": Solve and propagate using Lambert's problem
    - "LambertOptim": Find the best Lambert transfer between two orbits
    """
    # Load the model and initial state
    orb = LoadModel("Ref", orb)  #
    orb = LoadState("RefCapstone", orb)  #

    # Initialize 'seq' key in 'orb' dictionary
    orb['seq'] = {}  # Add this line to ensure 'seq' exists

    # Set the initial time for the sequence
    orb['seq']['Time'] = orb['sat']['t0']

    # Define sequence 'a'
    orb['seq']['a'] = {
        'type': "Propag",
        'span': 2 * 86400  # Propagate for 2 days
    }

    # Define sequence 'b'
    orb['seq']['b'] = {
        'type': "Lambert",
        'stop': "ELFO",     # Target state from LoadState
        'span': 2 * 3600    # Propagate for 2 hours
    }

    # Define sequence 'c'
    orb['seq']['c'] = {
        'type': "DVPropag",
        'Orbi': "ELFO",     # Target orbit to match velocity with
        'span': 3 * 3600    # Propagate for 3 hours
    }

    return orb