def accelalb(X, sat, const, stringBody, t, stringFrame, stringCentrBody, model):
    """
    Compute the acceleration due to Earth's albedo.

    Parameters:
    - X: ndarray, position vector of the spacecraft [km]
    - sat: dict, satellite parameters
        - sat['srp']['A']: satellite cross-sectional area [m^2]
        - sat['srp']['m']: satellite mass [kg]
        - sat['srp']['CR']: satellite's coefficient of reflectivity
        - sat['alb']['CR']: albedo coefficient
    - const: dict, constants
        - const['c']: speed of light [km/s]
        - const['Ls']: solar luminosity [W]
    - stringBody: str, name of the celestial body (e.g., 'EARTH')
    - t: float, time [seconds past J2000]
    - stringFrame: str, reference frame (e.g., 'J2000')
    - stringCentrBody: str, name of the central body (e.g., 'MOON')
    - model: dict, model parameters
        - model['centralPlanet']['RE']: radius of the central planet [km]

    Returns:
    - accalb: ndarray, acceleration due to Earth's albedo [km/s^2]
    """
    import numpy as np
    import spiceypy as spice

    # Get position of the body (e.g., Earth) relative to the central body
    XB = spice.spkezr(stringBody, t, stringFrame, 'NONE', stringCentrBody)[0]
    # Get position of the Sun relative to the central body
    XS = spice.spkezr('SUN', t, stringFrame, 'NONE', stringCentrBody)[0]

    # Unit vector from spacecraft to the body
    u = XB[:3] - X[:3]
    u = u / np.linalg.norm(u)

    # Satellite parameters
    A = sat['srp']['A']  # Cross-sectional area [m^2]
    m = sat['srp']['m']  # Mass [kg]
    CR = sat['srp']['CR']  # Coefficient of reflectivity
    Calb = sat['alb']['CR']  # Albedo coefficient

    # Constants
    c = const['c'] * 1e3  # Convert to m/s
    Ls = const['Ls']  # Solar luminosity [W]

    # Compute distances
    rS = np.linalg.norm(XS[:3] - XB[:3])  # Distance from Sun to body [km]
    r = np.linalg.norm(XB[:3])  # Distance from body to central body [km]
    rL = np.linalg.norm(X[:3])  # Distance from spacecraft to central body [km]
    d = np.sqrt(r**2 - model['centralPlanet']['RE']**2)  # Distance from body's surface to central body [km]
    dL = np.linalg.norm(X[:3] - XB[:3])  # Distance from spacecraft to body [km]

    # Angles
    cos_eta = d / r
    cos_etaL = (r**2 + dL**2 - rL**2) / (2 * r * dL)

    # Check if spacecraft is illuminated by Earth's albedo
    if (cos_etaL > cos_eta) and (dL > d):
        accalb = np.zeros(3)  # No acceleration due to albedo
    else:
        # Compute pressure
        PE = Ls / (4 * np.pi * c * rS**2 * 1e6)  # Solar flux at body [W/m^2]
        PE_mean = PE / 4  # Average solar flux on Earth's surface
        P = PE_mean * Calb * (6371 / dL)**2  # Pressure due to albedo [N/m^2]
        accalb = -(P * A * CR / m) * u * 1e-3  # Acceleration [km/s^2]

    return accalb