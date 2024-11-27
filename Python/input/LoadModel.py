from prop.normalizedharmonics import normalizedharmonics
def LoadModel(modeldef, orb):
    """
    Load the model used by the propagator.

    Parameters:
    - modeldef: str, model definition ('Ref' or 'Fast')
    - orb: dict, orbital data structure to update

    Returns:
    - orb: dict, updated orbital data structure
    """
    import os

    if modeldef == "Ref":
        HarmD = 150  # Maximum degree of the harmonics
        HarmO = 150  # Maximum order of the harmonics
        isE = 1      # Consider the Earth as a perturbation
        isS = 1      # Consider the Sun as a perturbation
        RC = 1.3     # Reflection coefficient
        isGR = 1     # Use general relativity
        isJ = 1      # Consider Jupiter as a perturbation
        HarmDE = 3   # Maximum degree of the Earth harmonics
        HarmOE = 3   # Maximum order of the Earth harmonics
        AC = 0.3     # Albedo coefficient

        Asat = 2     # Satellite area [m^2]
        msat = 1E3   # Satellite mass [kg]

    elif modeldef == "Fast":
        HarmD = 8
        HarmO = 8
        isE = 1
        isS = 0
        RC = 1.3
        isGR = 0
        isJ = 0
        HarmDE = 0
        HarmOE = 0
        AC = 0

        Asat = 2
        msat = 1E3

    else:
        raise ValueError(f"Unknown model definition: {modeldef}")

    # FORCE MODEL
    orb['frame'] = {}
    orb['frame']['integr'] = 'J2000'   # Integration frame (inertial)
    orb['frame']['from'] = 'MOON_PA'   # Frame where lunar gravity potential is defined
    orb['frame']['fromE'] = 'ITRF93'   # Frame where Earth's gravity potential is defined
    orb['frame']['to'] = orb['frame']['integr']

    orb['centralPlanet'] = {}
    orb['Earth'] = {}
    orb['pointMasses'] = {}
    orb['const'] = {}
    orb['prop'] = {}
    orb['sat'] = {}

    orb['centralPlanet']['stringName'] = 'Moon'
    orb['Earth']['stringName'] = 'Earth'
    orb['pointMasses']['stringName'] = ['Sun', 'JUPITER BARYCENTER']

    # Physical Constants
    orb['const']['G'] = 6.67428e-20   # Universal gravitational constant [km^3 kg^-1 s^-2]
    # Moon
    orb['centralPlanet']['RE'] = 1738                       # Radius of the Moon [km]
    orb['centralPlanet']['GM'] = 4902.7999671               # Gravitational parameter of the Moon [km^3 s^-2]
    # Earth
    orb['Earth']['RE'] = 6378.1363                          # Radius of the Earth [km]
    orb['Earth']['GM'] = 398600.4415 * isE                  # Gravitational parameter of the Earth [km^3 s^-2]
    # Point Masses: Sun and Jupiter
    orb['pointMasses']['M'] = [1.9884e30 * isS, 1.89813e27 * isJ]  # Masses [kg]
    orb['pointMasses']['GM'] = [orb['const']['G'] * M for M in orb['pointMasses']['M']]  # Gravitational parameters [km^3 s^-2]
    orb['pointMasses']['numb'] = len(orb['pointMasses']['M'])

    # Gravity models
    orb['prop']['harmonics'] = {}
    orb['prop']['harmonics']['degree'] = HarmD
    orb['prop']['harmonics']['order'] = HarmO
    orb['prop']['harmonics']['degreeE'] = HarmDE
    orb['prop']['harmonics']['orderE'] = HarmOE

    # File paths for gravity models
    current_dir = os.getcwd()
    # orb['prop']['harmonics']['filepath'] = os.path.join(current_dir, 'input',  'Moon_AIUB-GRL350B.txt')
    orb['prop']['harmonics']['filepath'] = os.path.join(current_dir, 'input/gravity_models/Moon_AIUB-GRL350B.txt')
    orb['prop']['harmonics']['filepathE'] = os.path.join(current_dir, 'input/gravity_models/Earth_EGM2008.txt')
    # orb['prop']['harmonics']['filepathE'] = os.path.join(current_dir, 'input',  'Earth_EGM2008.txt')

    # Harmonics coefficients
    orb['prop']['harmonics']['Cnm'], orb['prop']['harmonics']['Snm'] = normalizedharmonics(
        orb['prop']['harmonics']['filepath'], orb['prop']['harmonics']['degree']
    )
    orb['prop']['harmonics']['ECnm'], orb['prop']['harmonics']['ESnm'] = normalizedharmonics(
        orb['prop']['harmonics']['filepathE'], orb['prop']['harmonics']['degreeE']
    )

    # Solar Radiation Pressure and Earth Albedo
    orb['const']['c'] = 299792.458  # Speed of light [km/s]
    orb['const']['Ls'] = 3.839e26   # Solar luminosity [W]
    orb['sat']['srp'] = {}
    orb['sat']['alb'] = {}
    orb['sat']['srp']['A'] = Asat   # Satellite area [m^2]
    orb['sat']['srp']['m'] = msat   # Satellite mass [kg]
    orb['sat']['srp']['CR'] = RC    # Reflection coefficient
    orb['sat']['alb']['CR'] = AC    # Albedo coefficient

    # General Relativity
    orb['sat']['rel'] = isGR

    return orb