import itertools
import numpy as np


def neighbour_count(state, periodic=False):
    """
    Count the number of live neighbours of a regular n-d array.
    Parameters
    ----------

    state : numpy.ndarray of bools
       The Game of Life state.
    periodic : bool
       Whether boundaries are periodic.
    Returns
    ------
    ndarray of ints
        Number of living neighbours for each cell.


    The algorithm looks in each direction, and increases the count in a cell
    if that neighbour is alive. This version is generic and will perform
    the count for any system dimension.

    """

    dim = state.shape
    ncount = np.zeros(dim)

    # itertools.product([A,B,C, ...]) with A, B, (e.g.) lists gives
    # an iterator spanning every possible combination of tuples
    # taking the first element from A, the second from B etc.
    #
    # We use it here to automate looking one place back, level, and
    # one place forward in every dimension of the array
    #
    # This also uses two other Python syntax generalizations.
    # For functions f(*[a, b, c]) == f(a, b, c) and
    # 3*[(-1, 1)] = [(-1, 1), (-1, 1), (-1, 1)]
    combos = list(itertools.product(*(state.ndim*[(-1, 0, 1)])))

    # Because it include every combination, it also has the "do nothing"
    # case to count the cell itself. We don't want that, so we remove
    # the all zeros entry from the list.
    combos.remove(state.ndim*(0, ))

    if periodic:
        # In the periodic case, we can use the np.roll function to shift
        # things.
        for combo in combos:
            ncount[...] += 1*np.roll(state, combo, range(state.ndim))
    else:
        # In the non periodic case, we loop over the combinations and
        # deal only with the slices which are actually relevant
        # e.g. in 1D we want to do
        # neighbour_count[0:dim-1] += X[1:dim]
        # neighbour_count[1:dim] += X[0:dim-1]

        _slices = (slice(None, -1), slice(None, None), slice(1, None))

        def lhs_slice(combo):
            """Return the slice of the neighbour_count mesh to increment."""
            return tuple(_slices[c+1] for c in combo)

        def rhs_slice(combo):
            """Return the slice of the X mesh we're testing."""
            return tuple(_slices[1-c] for c in combo)

        for combo in combos:
            ncount[lhs_slice(combo)] += 1*state[rhs_slice(combo)]

    return ncount


def life(initial_state, nt, periodic=False):
    """
    Perform iterations of Conway's Game of Life.

    Parameters
    ----------
    initial_state : array_like or list of lists
         Initial 2d state of grid in an array of booleans.
    nt : int
         Number of steps of Life to perform.
    periodic : bool
         If true, then grid is assumed periodic.

    Returns
    -------
    array_like
         Final state of grid in array of booleans
    """

    state = np.array(initial_state, bool)

    # main loop
    for _ in range(nt):
        ncount = neighbour_count(state,
                                 periodic)
        state = state*(ncount == 2)+(ncount == 3)

    return state
