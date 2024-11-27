def TBPOptim(Xi, orb):
    """
    Optimize the initial state for a three-body problem orbit to achieve
    periodicity in the Earth-Moon rotational frame.

    Parameters:
    - Xi: Initial state vector.
    - orb: Orbital data structure.
    """
    import numpy as np
    from scipy.optimize import fsolve

    def Propag(X):
        import numpy as np
        from scipy.integrate import solve_ivp
        from prop.prophpop import prophpop

        n = 3  # Number of periods to consider
        options = {'rtol': 1e-7, 'atol': 1e-14}

        # Time spans for each period
        span = [orb['seq']['Time'] + i * orb['seq']['a']['T'] for i in range(n + 1)]

        Ytraj = []
        Xr_initial = MatRI(span[0], X)[2]

        dY = np.zeros(6 * n)
        for i in range(n):
            sol = solve_ivp(
                prophpop, [span[i], span[i + 1]], X, method='DOP853',
                rtol=options['rtol'], atol=options['atol'], args=(orb,)
            )
            X = sol.y[:, -1]
            _, _, Ytr = MatRI(span[i + 1], X)
            dY[6 * i:6 * i + 3] = Ytr[:3] - Xr_initial[:3]
            dY[6 * i + 3:6 * i + 6] = Ytr[3:6] - Xr_initial[3:6]
        return dY

    X, infodict, ier, mesg = fsolve(Propag, Xi, full_output=True)
    if ier != 1:
        print(f"Optimization did not converge: {mesg}")
    return X

def MatRI(t, X):
    """
    Computes the rotation matrix from rotational to inertial frame at time t,
    its inverse, and the state vector X in the rotational frame.
    """
    import numpy as np
    import spiceypy as spice

    SE = spice.spkezr('EARTH', t, 'J2000', 'NONE', 'MOON')[0]
    dt = 1
    SEm = spice.spkezr('EARTH', t - dt, 'J2000', 'NONE', 'MOON')[0]

    uR = -SE[:3] / np.linalg.norm(SE[:3])
    uRm = -SEm[:3] / np.linalg.norm(SEm[:3])
    uTh = -SE[3:6] / np.linalg.norm(SE[3:6])
    uZ = np.cross(uR, uTh)

    Mri = np.column_stack((uR, uTh, uZ))
    Mir = np.linalg.inv(Mri)
    Omg = np.linalg.norm(uR - uRm) / dt / np.linalg.norm(uTh)

    Xr = np.zeros(6)
    Xr[:3] = Mir.dot(X[:3])
    Xr[3:6] = Mir.dot(X[3:6] - np.cross(Omg * uZ, X[:3]))
    return Mri, Mir, Xr