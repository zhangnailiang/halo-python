from prop.prophpop import prophpop
from prop.lambert import lambert
from input.LoadState import LoadState
def LambertOptim(orb):
    """
    Optimizes the Lambert transfer between two orbits.

    Parameters:
    - orb: Orbital data structure.
    """
    import numpy as np
    from scipy.optimize import fmin

    def DVTOT(T):
        """
        Computes the total Delta-V for the Lambert transfer, including a penalty term.

        Parameters:
        - T: Array containing [t1, t2, span].
        """
        X1 = Pos(T[0], orb['sat']['X0iner'], orb)

        target = orb['seq']['a']['target']
        orb2 = LoadState(target, orb)
        X2 = Pos(T[1], orb2['sat']['X0iner'], orb)

        V1, V2, extremal_distances, exitflag = lambert(
            X1[:3], X2[:3], T[2] / 86400, 0, orb['centralPlanet']['GM']
        )

        # Assessment of Lambert solution
        X2t = Pos(T[2], np.concatenate((X1[:3], V1)), orb)
        pun = np.linalg.norm(X2t[:3] - X2[:3]) * 1 / 2e3  # Punishment term

        DV1 = np.linalg.norm(X1[3:6] - V1)
        DV2 = np.linalg.norm(X2[3:6] - V2)
        DV = DV1 + DV2 + pun
        return DV

    T0 = [orb['seq']['a']['t1'], orb['seq']['a']['t2'], orb['seq']['a']['span']]
    T = fmin(DVTOT, T0, disp=True)
    t1, t2, span = T

    target = orb['seq']['a']['target']
    orb2 = LoadState(target, orb)
    X2 = Pos(t2, orb2['sat']['X0iner'], orb)
    return t1, t2, X2, span

def Pos(t, X0, orb):
    """
    Propagates an initial state X0 for a time t.

    Parameters:
    - t: Time to propagate.
    - X0: Initial state vector.
    - orb: Orbital data structure.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    options = {'rtol': 1e-7, 'atol': 1e-12}
    sol = solve_ivp(
        prophpop, [0, t], X0, method='DOP853',
        rtol=options['rtol'], atol=options['atol'], args=(orb,)
    )
    X = sol.y[:, -1]
    return X