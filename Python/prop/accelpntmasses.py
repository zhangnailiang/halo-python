def accelpntmasses(r, pointMasses, GM, t, frame, centralPlanet, model):
    """
    Compute the acceleration due to point masses (e.g., Sun, Jupiter).

    Parameters:
    - r: ndarray, spacecraft position vector [km]
    - pointMasses: list of strings, names of point masses (e.g., ['SUN', 'JUPITER BARYCENTER'])
    - GM: list or ndarray, gravitational parameters of the point masses [km^3/s^2]
    - t: float, time [seconds past J2000]
    - frame: str, reference frame (e.g., 'J2000')
    - centralPlanet: str, name of the central planet (e.g., 'MOON')
    - model: dict, model parameters (not used in this function)

    Returns:
    - acc_pointMasses: ndarray, total acceleration due to point masses [km/s^2]
    """
    import numpy as np
    import spiceypy as spice

    acc_pointMasses = np.zeros(3)
    for j in range(len(pointMasses)):
        # Get the state vector of the point mass relative to the central planet
        temp = spice.spkezr(pointMasses[j], t, frame, 'NONE', centralPlanet)[0]
        r0j = temp[:3]  # Position vector of point mass j

        # Compute the vectors
        rjS = r - r0j  # Vector from point mass to spacecraft

        # Compute distances cubed
        r0j_norm3 = np.linalg.norm(r0j) ** 3
        rjS_norm3 = np.linalg.norm(rjS) ** 3

        # Update the acceleration
        acc_pointMasses += GM[j] * (-r0j / r0j_norm3 - rjS / rjS_norm3)

    return acc_pointMasses
