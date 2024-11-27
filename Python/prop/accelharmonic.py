import numpy as np

def accelharmonic(r, E, n_max, m_max, Cnm, Snm, GM, R_ref):
    """
    Calculate the inertial acceleration due to the gravity field of a central planet
    considering its spherical harmonics.

    Parameters:
    r : np.ndarray
        Position vector of the satellite in an inertial frame [km].
    E : np.ndarray
        Rotation matrix from inertial frame to body-fixed frame [3x3].
    n_max : int
        Maximum degree of harmonics coefficients.
    m_max : int
        Maximum order of harmonics coefficients.
    Cnm : np.ndarray
        Normalized harmonics coefficients (cosine terms) [n x n].
    Snm : np.ndarray
        Normalized harmonics coefficients (sine terms) [n x n].
    GM : float
        Gravitational parameter of the central planet [km^3/s^2].
    R_ref : float
        Reference radius of the central planet [km].

    Returns:
    a : np.ndarray
        Acceleration vector of the satellite in the inertial frame [km/s^2].
    """
    # Body-fixed position
    r_bf = E @ r

    # Auxiliary quantities
    d = np.linalg.norm(r_bf)
    latgc = np.arcsin(r_bf[2] / d)
    lon = np.arctan2(r_bf[1] / (d * np.cos(latgc)), r_bf[0] / (d * np.cos(latgc)))

    # Calculate Legendre polynomials and their derivatives
    pnm, dpnm = legendre(n_max, m_max, latgc)

    dUdr = 0
    dUdlatgc = 0
    dUdlon = 0
    q3 = 0
    q2 = q3
    q1 = q2

    for n in range(n_max + 1):
        b1 = (-GM / d**2) * (R_ref / d)**n * (n + 1)
        b2 = (GM / d) * (R_ref / d)**n
        b3 = (GM / d) * (R_ref / d)**n
        for m in range(m_max + 1):
            cos_m_lon = np.cos(m * lon)
            sin_m_lon = np.sin(m * lon)

            q1 += pnm[n, m] * (Cnm[n, m] * cos_m_lon + Snm[n, m] * sin_m_lon)
            q2 += dpnm[n, m] * (Cnm[n, m] * cos_m_lon + Snm[n, m] * sin_m_lon)
            q3 += m * pnm[n, m] * (Snm[n, m] * cos_m_lon - Cnm[n, m] * sin_m_lon)

        dUdr += q1 * b1
        dUdlatgc += q2 * b2
        dUdlon += q3 * b3
        q3 = 0
        q2 = q3
        q1 = q2
    # Body-fixed acceleration
    r2xy = r_bf[0]**2 + r_bf[1]**2

    ax = (1 / d * dUdr - r_bf[2] / (d**2 * np.sqrt(r2xy)) * dUdlatgc) * r_bf[0] - (1 / r2xy * dUdlon) * r_bf[1]
    ay = (1 / d * dUdr - r_bf[2] / (d**2 * np.sqrt(r2xy)) * dUdlatgc) * r_bf[1] + (1 / r2xy * dUdlon) * r_bf[0]
    az = 1 / d * dUdr * r_bf[2] + np.sqrt(r2xy) / d**2 * dUdlatgc

    # a_bf = np.array([ax, ay, az])
    a_bf = np.array([ax, ay, az]).T

    # Inertial acceleration
    a = E.T @ a_bf
    return a


def legendre(n, m, fi):
    """
    Compute the normalized associated Legendre polynomials and their derivatives.

    Parameters:
    n : int
        Maximum degree.
    m : int
        Maximum order.
    fi : float
        Geocentric latitude in radians.

    Returns:
    pnm : np.ndarray
        Associated Legendre polynomials.
    dpnm : np.ndarray
        Derivatives of the associated Legendre polynomials.
    """
    pnm = np.zeros((n + 1, m + 1))
    dpnm = np.zeros((n + 1, m + 1))

    # Initial conditions
    pnm[0, 0] = 1
    dpnm[0, 0] = 0
    pnm[1, 1] = np.sqrt(3) * np.cos(fi)
    dpnm[1, 1] = -np.sqrt(3) * np.sin(fi)

    # Diagonal coefficients
    for i in range(1, n):
        pnm[i + 1, i + 1] = np.sqrt((2 * i + 1) / (2 * i)) * np.cos(fi) * pnm[i, i]
        dpnm[i + 1, i + 1] = np.sqrt((2 * i + 1) / (2 * i)) * (np.cos(fi) * dpnm[i, i] - np.sin(fi) * pnm[i, i])

    # Horizontal first step coefficients
    for i in range(n):
        pnm[i + 1, i] = np.sqrt(2 * i + 1) * np.sin(fi) * pnm[i, i]
        dpnm[i + 1, i] = np.sqrt(2 * i + 1) * (np.cos(fi) * pnm[i, i] + np.sin(fi) * dpnm[i, i])

    # Horizontal second step coefficients
    j = 0
    k = 2
    while j <= m:
        for i in range(k, n + 1):
            pnm[i, j] = np.sqrt((2 * i + 1) / ((i - j) * (i + j))) * (
                np.sqrt(2 * i - 1) * np.sin(fi) * pnm[i - 1, j] - np.sqrt(((i + j - 1) * (i - j - 1)) / (2 * i - 3)) * pnm[i - 2, j]
            )
        j += 1
        k += 1

    j = 0
    k = 2
    while j <= m:
        for i in range(k, n + 1):
            dpnm[i, j] = np.sqrt((2 * i + 1) / ((i - j) * (i + j))) * (
                np.sqrt(2 * i - 1) * (np.sin(fi) * dpnm[i - 1, j] + np.cos(fi) * pnm[i - 1, j])
                - np.sqrt(((i + j - 1) * (i - j - 1)) / (2 * i - 3)) * dpnm[i - 2, j]
            )
        j += 1
        k += 1

    return pnm, dpnm
