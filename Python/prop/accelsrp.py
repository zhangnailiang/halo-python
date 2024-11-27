def accelsrp(X, srp, const, stringBody, t, stringFrame, stringCentrBody, model):
    """
    Compute the acceleration due to solar radiation pressure.

    Parameters:
    - X: ndarray, spacecraft state vector [km]
    - srp: dict, solar radiation pressure parameters
        - srp['A']: cross-sectional area of the spacecraft [m^2]
        - srp['m']: mass of the spacecraft [kg]
        - srp['CR']: coefficient of reflectivity
    - const: dict, constants
        - const['c']: speed of light [km/s]
        - const['Ls']: solar luminosity [W]
    - stringBody: str, name of the celestial body (e.g., 'SUN')
    - t: float, time [seconds past J2000]
    - stringFrame: str, reference frame (e.g., 'J2000')
    - stringCentrBody: str, name of the central body (e.g., 'MOON')
    - model: dict, model parameters
        - model['centralPlanet']['RE']: radius of the central planet [km]

    Returns:
    - accsrp: ndarray, acceleration due to solar radiation pressure [km/s^2]
    """
    import numpy as np
    import spiceypy as spice

    # Get position of the body (e.g., Sun) relative to the central body
    XB = spice.spkezr(stringBody, t, stringFrame, 'NONE', stringCentrBody)[0]

    # Unit vector from spacecraft to the body
    u = XB[:3] - X[:3]
    u_norm = np.linalg.norm(u)
    if u_norm == 0:
        u_norm = 1e-10  # Avoid division by zero
    u = u / u_norm

    # Parameters
    A = srp['A']  # Cross-sectional area [m^2]
    m = srp['m']  # Mass [kg]
    CR = srp['CR']  # Coefficient of reflectivity
    c = const['c'] * 1e3  # Speed of light [m/s]
    Ls = const['Ls']  # Solar luminosity [W]

    # Distances
    r = np.linalg.norm(XB[:3])  # Distance from central body to the body [km]
    rL = np.linalg.norm(X[:3])  # Distance from central body to spacecraft [km]
    d = np.sqrt(r ** 2 - model['centralPlanet']['RE'] ** 2)  # Distance from body's surface to central body [km]
    dL = np.linalg.norm(X[:3] - XB[:3])  # Distance from spacecraft to the body [km]

    # Angles
    cos_eta = d / r
    cos_etaL = (r ** 2 + dL ** 2 - rL ** 2) / (2 * r * dL)

    # Check if spacecraft is in shadow
    if (cos_etaL > cos_eta) and (dL > d):
        accsrp = np.zeros(3)  # No acceleration due to SRP
    else:
        # Compute pressure
        P = Ls / (4 * np.pi * c * dL ** 2 * 1e6)  # Solar pressure at distance dL [N/m^2]
        accsrp = -(P * A * CR * 1e-3 / m) * u  # Acceleration [km/s^2]

    return accsrp

