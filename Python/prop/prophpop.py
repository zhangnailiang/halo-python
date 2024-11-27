import numpy as np
from scipy.integrate import odeint
from spiceypy import spkezr, pxfrm2
from prop.accelalb import accelalb
from prop.accelsrp import accelsrp
from prop.accelpntmasses import accelpntmasses
from prop.accelharmonic import accelharmonic

import numpy as np
import spiceypy as spice

def prophpop(t, X0, model):
    """
    Computes the derivative of the state vector for propagation.

    Parameters:
    - t: Time.
    - X0: State vector at time t.
    - model: Model containing necessary parameters and functions.
    """
    import numpy as np
    import spiceypy as spice

    # Extract position and velocity
    xJ2000 = X0[:3]
    vJ2000 = X0[3:6]

    # Acceleration due to point masses (e.g., Sun)
    acc_pointMasses = accelpntmasses(
        xJ2000,
        model['pointMasses']['stringName'],
        model['pointMasses']['GM'],
        t,
        model['frame']['integr'],
        model['centralPlanet']['stringName'],
        model  # Added this parameter
    )

    # Acceleration due to central planet gravity field
    Rgeog_iner = spice.pxfrm2(model['frame']['from'], model['frame']['to'], t, t)
    acc_centralPlanet = accelharmonic(
        xJ2000,
        Rgeog_iner.T,
        model['prop']['harmonics']['degree'],
        model['prop']['harmonics']['order'],
        model['prop']['harmonics']['Cnm'],
        model['prop']['harmonics']['Snm'],
        model['centralPlanet']['GM'],
        model['centralPlanet']['RE']
    )

    # Acceleration due to Earth gravity field
    REarth_iner = spice.pxfrm2(model['frame']['fromE'], model['frame']['to'], t, t)
    R_Moon_Earth = spice.spkezr('EARTH', t, model['frame']['to'], 'NONE', 'MOON')[0][:3]
    acc_earth = (
        accelharmonic(
            xJ2000 - R_Moon_Earth,
            REarth_iner.T,
            model['prop']['harmonics']['degreeE'],
            model['prop']['harmonics']['orderE'],
            model['prop']['harmonics']['ECnm'],
            model['prop']['harmonics']['ESnm'],
            model['Earth']['GM'],
            model['Earth']['RE']
        ) - accelharmonic(
            -R_Moon_Earth,
            REarth_iner.T,
            model['prop']['harmonics']['degreeE'],
            model['prop']['harmonics']['orderE'],
            model['prop']['harmonics']['ECnm'],
            model['prop']['harmonics']['ESnm'],
            model['Earth']['GM'],
            model['Earth']['RE']
        )
    )

    # Acceleration due to general relativity
    gamma = 1
    beta = 1
    c = model['const']['c']
    acc_genRel = (
        model['sat']['rel'] * model['centralPlanet']['GM'] / (c ** 2 * np.linalg.norm(xJ2000) ** 3)
        * (
            (2 * (beta + gamma) * model['centralPlanet']['GM'] / np.linalg.norm(xJ2000) - gamma * np.linalg.norm(vJ2000) ** 2) * xJ2000
            + 2 * (1 + gamma) * np.dot(xJ2000, vJ2000) * vJ2000
        )
    )

    # Acceleration due to solar radiation pressure from the Sun
    acc_SRPSun = accelsrp(
        xJ2000,
        model['sat']['srp'],
        model['const'],
        'SUN',
        t,
        model['frame']['integr'],
        model['centralPlanet']['stringName'],
        model
    )

    # Acceleration due to Earth albedo
    acc_alb = accelalb(
        xJ2000,
        model['sat'],
        model['const'],
        'EARTH',
        t,
        model['frame']['integr'],
        model['centralPlanet']['stringName'],
        model
    )

    # Total acceleration
    xdot = vJ2000
    vdot = acc_centralPlanet + acc_earth + acc_pointMasses + acc_SRPSun + acc_alb + acc_genRel

    # Return derivative of state vector
    return np.concatenate((xdot, vdot))
