import spiceypy as spice
import numpy as np

def LoadState(orbit, orb):
    # LoadState: This function gathers the states of some selected orbits of interest.
    # The Orbits called "RefSpacecraft" are the orbit of interest used in the paper.
    # The NRHO and DRO14 are CR3BP orbits (close to the Capstone and Orion one) fitted to ephemeris model.
    # ELFO is an Elliptical Lunar Frozen Orbit.
    # Two examples/templates for Keplerian coordinates, Kepl and Kepl2, were also added.
    # The user can add new ones or modify the "Perso" Orbit in last position.

    if orbit == "NRHO":
        ref_s = 'J2000'  # Reference System (J2000/MOON_ME)
        coord_t = 'Cartesian'  # Coordinates Type (Cartesian/Keplerian)

        X, Y, Z = -14007.49388066819, 37926.40572749509, -59472.47612911355  # Coordinates
        VX, VY, VZ = 0.03815140145811324, 0.05511836319581054, 0.028138675143097314  # Velocities
        time = '2024 Jan 01 00:00:00.000'  # start time [yyyy month dd hh:mm:ss.---]

    elif orbit == "DRO14":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = 63650.58084169839, -24686.81610635431, -15448.941380173477
        VX, VY, VZ = -0.14532575232216205, -0.290549113768819, -0.1502956531877957
        time = '2024 Jan 01 00:00:00.000'

    elif orbit == "ELFO":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = 12.931, 235.558, 2388.38
        VX, VY, VZ = 1.80519, -0.0991, -3.83e-07
        time = '2024 Jan 01 00:00:00.000'

    elif orbit == "RefClementine":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = -2098.4745703999997, 1298.5243538999998, -3494.2586871999997
        VX, VY, VZ = 0.8948704814799997, 0.17693000307000004, -0.15601097269
        time = '1994 Apr 15 15:00:00.000'

    elif orbit == "RefCapstone":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = -16983.14075642353, 21213.5584242304, -58035.6304537942
        VX, VY, VZ = -0.04457645856905286, -0.05137482728591616, 0.1416421540429468
        time = '2022 Nov 25 00:00:00.000'

    elif orbit == "RefOrion":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = 41570.98496247831, -50015.47166299814, -28760.46002667214
        VX, VY, VZ = -0.2465381485425481, -0.182296424457602, -0.07642694222326024
        time = '2022 Nov 29 16:00:00.000'

    elif orbit == "RefLRO":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = 820.1514640535461, 1091.2367697664858, 1262.7894519284164
        VX, VY, VZ = 0.22458861495091112, 1.1177228435015545, -1.1391038535297986
        time = '2020 Feb 01 00:00:00.000'

    elif orbit == "Kepl":
        ref_s = 'MOON_ME'
        coord_t = 'Keplerian'

        SMA, ECC, INC, LAN, AOP, MA = 1738 + 100, 0, 90, 0, 45, 0
        time = '2020 Feb 01 00:00:00.000'

    elif orbit == "Kepl2":
        ref_s = 'MOON_ME'
        coord_t = 'Keplerian2'

        PER, APO, INC, LAN, AOP, MA = 1738 + 100, 1738 + 1000, 90, 0, 45, 0
        time = '2020 Feb 01 00:00:00.000'

    elif orbit == "Perso":
        ref_s = 'J2000'
        coord_t = 'Cartesian'

        X, Y, Z = -175.001861, -1175.17121, 1410.69401
        VX, VY, VZ = 0.791300523, 1.03420255, 0.971046194
        time = '2022 Nov 26 00:00:00.000'

    # Initialization
    orb['frame']['initstate'] = ref_s
    et = spice.str2et(time)
    orb['sat']['t0'] = et

    if coord_t == 'Cartesian':
        orb['sat']['X0'] = np.array([X, Y, Z, VX, VY, VZ])
        if ref_s not in ['ICRF', 'J2000']:
            Rgeog_iner = spice.sxform(ref_s, orb['frame']['to'], et)
            orb['sat']['X0iner'] = np.dot(Rgeog_iner, orb['sat']['X0'])
        else:
            orb['sat']['X0iner'] = orb['sat']['X0']
    else:
        if coord_t == 'Keplerian':
            orb['sat']['keplstate'] = {
                'SMA': SMA,
                'ECC': ECC,
                'INC': np.radians(INC),
                'LAN': np.radians(LAN),
                'AOP': np.radians(AOP),
                'MA': np.radians(MA),
            }
        else:
            SMA = (APO + PER) / 2
            ECC = 1 - PER / SMA
            orb['sat']['keplstate'] = {
                'SMA': SMA,
                'ECC': ECC,
                'INC': np.radians(INC),
                'LAN': np.radians(LAN),
                'AOP': np.radians(AOP),
                'MA': np.radians(MA),
            }
        conics_params = [
            orb['sat']['keplstate']['SMA'] * (1 - orb['sat']['keplstate']['ECC']),
            orb['sat']['keplstate']['ECC'],
            orb['sat']['keplstate']['INC'],
            orb['sat']['keplstate']['LAN'],
            orb['sat']['keplstate']['AOP'],
            orb['sat']['keplstate']['MA'],
            et,
            orb['centralPlanet']['GM'],
        ]
        orb['sat']['X0'] = spice.conics(conics_params, et)
        if ref_s not in ['ICRF', 'J2000']:
            Rgeog_iner = spice.pxform(ref_s, orb['frame']['to'], et)
            orb['sat']['X0iner'] = np.concatenate((
                np.dot(Rgeog_iner, orb['sat']['X0'][:3]),
                np.dot(Rgeog_iner, orb['sat']['X0'][3:])
            ))
        else:
            orb['sat']['X0iner'] = orb['sat']['X0']

    return orb