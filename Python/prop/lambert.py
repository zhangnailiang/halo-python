def lambert(r1vec, r2vec, tf, m, muC):
    """
    Solves the Lambert problem using an adjusted Lancaster-Blanchard method.

    Parameters:
    - r1vec: ndarray, initial position vector.
    - r2vec: ndarray, final position vector.
    - tf: float, flight time (in days). Positive indicates a short path, negative indicates a long path.
    - m: int, number of revolutions. Positive for right branch, negative for left branch.
    - muC: float, gravitational parameter.

    Returns:
    - V1: ndarray, initial velocity vector.
    - V2: ndarray, final velocity vector.
    - extremal_distances: list, [minimum distance, maximum distance].
    - exitflag: int, status flag (1 for success, negative values indicate errors).
    """
    import numpy as np

    # Tolerance and constants
    tol = 1e-14
    bad = False
    days = 86400  # Seconds in a day

    # Non-dimensionalization
    r1 = np.linalg.norm(r1vec)
    r1_unit = r1vec / r1
    V = np.sqrt(muC / r1)
    r2vec = r2vec / r1
    T = r1 / V
    tf = tf * days / T  # Convert tf to non-dimensional time

    # Calculate transfer angle
    mr2vec = np.linalg.norm(r2vec)
    cos_dth = np.dot(r1_unit, r2vec) / mr2vec
    cos_dth = np.clip(cos_dth, -1, 1)  # Ensure within [-1, 1] range
    dth = np.arccos(cos_dth)

    # Determine long or short path
    longway = np.sign(tf)
    tf = abs(tf)
    if longway < 0:
        dth = 2 * np.pi - dth

    # Determine left or right branch
    leftbranch = np.sign(m)
    m = abs(m)

    # Calculate relevant parameters
    c = np.sqrt(1 + mr2vec ** 2 - 2 * mr2vec * np.cos(dth))  # Chord length
    s = (1 + mr2vec + c) / 2  # Semi-perimeter
    a_min = s / 2  # Semi-major axis of the minimum energy orbit
    Lambda = np.sqrt(mr2vec) * np.cos(dth / 2) / s  # Lambda parameter
    crossprd = np.cross(r1_unit, r2vec)
    mcr = np.linalg.norm(crossprd)
    if mcr == 0:
        # Handle zero vector case
        nrmunit = np.array([0, 0, 0])
    else:
        nrmunit = crossprd / mcr  # Normal unit vector

    # Initial guess
    if m == 0:
        inn1 = -0.5233
        inn2 = 0.5233
        x1 = np.log(1 + inn1)
        x2 = np.log(1 + inn2)
    else:
        if leftbranch < 0:
            inn1 = -0.5234
            inn2 = -0.2234
        else:
            inn1 = 0.7234
            inn2 = 0.5234
        x1 = np.tan(inn1 * np.pi / 2)
        x2 = np.tan(inn2 * np.pi / 2)

    xx = np.array([inn1, inn2])
    aa = a_min / (1 - xx ** 2)
    # Calculate beta and alfa
    term = np.sqrt((s - c) / (2 * aa))
    term = np.clip(term, -1, 1)  # Ensure within [-1, 1] range
    bbeta = longway * 2 * np.arcsin(term)
    aalfa = 2 * np.arccos(np.clip(xx, -1, 1))

    # Initial estimation of flight time
    y12 = aa * np.sqrt(aa) * ((aalfa - np.sin(aalfa)) - (bbeta - np.sin(bbeta)) + 2 * np.pi * m)

    # Initial error estimation
    if m == 0:
        logt = np.log(tf)
        y1 = np.log(y12[0]) - logt
        y2 = np.log(y12[1]) - logt
    else:
        y1 = y12[0] - tf
        y2 = y12[1] - tf

    # Newton-Raphson iteration to solve for x
    err = np.inf
    iterations = 0
    while err > tol:
        iterations += 1
        xnew = (x1 * y2 - y1 * x2) / (y2 - y1)
        if m == 0:
            x = np.exp(xnew) - 1
        else:
            x = np.tan(xnew * np.pi / 2)
        a = a_min / (1 - x ** 2)
        if x < 1:  # Elliptical orbit
            beta = longway * 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
            alfa = 2 * np.arccos(np.clip(x, -1, 1))
            psi = (alfa - beta) / 2
            eta2 = 2 * a * np.sin(psi) ** 2 / s
            eta = np.sqrt(eta2)
        else:  # Hyperbolic orbit
            alfa = 2 * np.arccosh(x)
            beta = longway * 2 * np.arcsinh(np.sqrt((c - s) / (2 * a)))
            psi = (alfa - beta) / 2
            eta2 = -2 * a * np.sinh(psi) ** 2 / s
            eta = np.sqrt(eta2)
        # Calculate flight time
        if a > 0:
            tof = a * np.sqrt(a) * ((alfa - np.sin(alfa)) - (beta - np.sin(beta)) + 2 * np.pi * m)
        else:
            tof = -a * np.sqrt(-a) * ((np.sinh(alfa) - alfa) - (np.sinh(beta) - beta))
        # Update error estimation
        if m == 0:
            ynew = np.log(tof) - np.log(tf)
        else:
            ynew = tof - tf
        x1, x2 = x2, xnew
        y1, y2 = y2, ynew
        err = abs(x2 - xnew)
        if iterations > 15:
            bad = True
            break

    # Check if iteration failed
    if bad:
        exitflag = -1
        V1 = np.full(3, np.nan)
        V2 = np.full(3, np.nan)
        extremal_distances = [np.nan, np.nan]
        return V1, V2, extremal_distances, exitflag

    # Calculate final x value
    if m == 0:
        x = np.exp(xnew) - 1
    else:
        x = np.tan(xnew * np.pi / 2)

    # Calculate semi-major axis a
    a = a_min / (1 - x ** 2)

    # Recalculate psi, eta, etc.
    if x < 1:
        beta = longway * 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
        alfa = 2 * np.arccos(np.clip(x, -1, 1))
        psi = (alfa - beta) / 2
        eta2 = 2 * a * np.sin(psi) ** 2 / s
        eta = np.sqrt(eta2)
    else:
        beta = longway * 2 * np.arcsinh(np.sqrt((c - s) / (2 * a)))
        alfa = 2 * np.arccosh(x)
        psi = (alfa - beta) / 2
        eta2 = -2 * a * np.sinh(psi) ** 2 / s
        eta = np.sqrt(eta2)

    # Calculate unit normal vector
    ih = longway * nrmunit

    # Normalize r2 vector
    r2_unit = r2vec / mr2vec

    # Calculate velocity direction
    crsprd1 = np.cross(ih, r1_unit)
    crsprd2 = np.cross(ih, r2_unit)

    # Calculate radial and transverse components of departure and arrival velocities
    Vr1 = 1 / (eta * np.sqrt(a_min)) * (2 * Lambda * a_min - Lambda - x * eta)
    Vt1 = np.sqrt(mr2vec / (a_min * eta2) * np.sin(dth / 2) ** 2)

    Vt2 = Vt1 / mr2vec
    Vr2 = (Vt1 - Vt2) / np.tan(dth / 2) - Vr1

    # Calculate departure and arrival velocity vectors
    V1 = (Vr1 * r1_unit + Vt1 * crsprd1) * V
    V2 = (Vr2 * r2_unit + Vt2 * crsprd2) * V

    # Set success flag
    exitflag = 1

    # Calculate minimum and maximum distances
    extremal_distances = minmax_distances(
        r1_unit * r1, r1, r2_unit * r1, mr2vec * r1, dth, a * r1, V1, V2, m, muC
    )

    return V1, V2, extremal_distances, exitflag


def minmax_distances(r1vec, r1, r2vec, r2, dth, a, V1, V2, m, muC):
    """
    Compute minimum and maximum distances to the central body during the transfer.

    Parameters:
    - r1vec: ndarray, initial position vector.
    - r1: float, magnitude of the initial position.
    - r2vec: ndarray, final position vector.
    - r2: float, magnitude of the final position.
    - dth: float, transfer angle.
    - a: float, semi-major axis.
    - V1: ndarray, initial velocity vector.
    - V2: ndarray, final velocity vector.
    - m: int, number of revolutions.
    - muC: float, gravitational parameter.

    Returns:
    - extremal_distances: list, [minimum distance, maximum distance].
    """
    import numpy as np

    # By default, minimum and maximum distances are the minimum and maximum of r1 and r2
    minimum_distance = min(r1, r2)
    maximum_distance = max(r1, r2)

    # Check if the long path was used
    longway = abs(dth) > np.pi

    # Calculate eccentricity vector
    h_vec = np.cross(r1vec, V1)
    h = np.linalg.norm(h_vec)
    e_vec = (np.cross(V1, h_vec) / muC) - (r1vec / r1)
    e = np.linalg.norm(e_vec)

    # Calculate pericenter and apocenter
    pericenter = a * (1 - e)
    if e < 1:
        apocenter = a * (1 + e)
    else:
        apocenter = np.inf

    if m > 0:
        # Multiple revolutions - orbit passes both pericenter and apocenter
        minimum_distance = pericenter
        maximum_distance = apocenter
    else:
        # Calculate theta1 and theta2
        pm1 = np.sign(np.dot(np.cross(e_vec, r1vec), h_vec))
        pm2 = np.sign(np.dot(np.cross(e_vec, r2vec), h_vec))
        theta1 = pm1 * np.arccos(np.clip(np.dot(e_vec / e, r1vec / r1), -1, 1))
        theta2 = pm2 * np.arccos(np.clip(np.dot(e_vec / e, r2vec / r2), -1, 1))

        # Check if pericenter or apocenter was passed
        if theta1 * theta2 < 0:
            if abs(abs(theta1) + abs(theta2) - dth) < 5 * np.finfo(float).eps * abs(dth):
                minimum_distance = pericenter
            else:
                maximum_distance = apocenter
        elif longway:
            minimum_distance = pericenter
            if e < 1:
                maximum_distance = apocenter

    extremal_distances = [minimum_distance, maximum_distance]
    return extremal_distances
