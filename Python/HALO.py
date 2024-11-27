import time
import numpy as np
import json
import os

import spiceypy as spice
from scipy.integrate import solve_ivp
import pickle
from input.LoadState import LoadState
from input.LoadModel import LoadModel
from metakernelcheck import metakernelcheck
from prop.TBPOptim import TBPOptim
from prop.accelalb import accelalb
from prop.accelpntmasses import accelpntmasses
from prop.accelharmonic import accelharmonic
from prop.accelsrp import accelsrp
from prop.lambert import lambert
from prop.LambertOptim import LambertOptim
from input.LoadSequential import LoadSequential
from prop.prophpop import prophpop


#start time for runtime calculation
start_time = time.time()
# Start time
ts = time.process_time()

# SPICE LIBRARIES
spice.kclear()
metakernelcheck()
spice.furnsh('metakernel.tm')


# Initialization
orb = {}
orb = LoadSequential(orb)  # Ensure to provide this function
Time = orb['seq']['Time']




# Running the sequence
SeqNames = list(orb['seq'].keys())
options = {'rtol': 1e-7, 'atol': 1e-14}
TStep = 30

for i in range(1, len(SeqNames)):
    Seq = orb['seq'][SeqNames[i]]

    if Seq['type'] == "Propag":
        
        span = Seq['span']
        Tspan = (Time, Time + span)
        t_eval = np.arange(Time, Time + span + TStep, TStep)
        X0 = orb['sat']['X0iner']

        block_start_time = time.time()


        
        # Use solve_ivp with method similar to ode113
        sol = solve_ivp(
            prophpop, Tspan, X0, method='DOP853', t_eval=t_eval,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )
        orb['seq'][SeqNames[i]]['t'] = sol.t
        orb['seq'][SeqNames[i]]['XJ2000'] = sol.y.T
        block_end_time = time.time()
        block_runtime = block_end_time - block_start_time
        print("block_runtime:", block_runtime)
        

    elif Seq['type'] == "DVPropag":
        X = orb['sat']['X0iner']
        orb = LoadState(Seq['Orbi'], orb)  # Ensure to provide this function
        newX = orb['sat']['X0iner']


        print(f"Diff. pos. target = {np.linalg.norm(newX[:3] - X[:3])} km")
        print(f"DV2 = {np.linalg.norm(newX[3:6] - X[3:6]) * 1e3} m/s")
        span = Seq['span']
        Tspan = (Time, Time + span)
        t_eval = np.arange(Time, Time + span + TStep, TStep)
        X0 = np.concatenate((X[:3], newX[3:6]))

        sol = solve_ivp(
            prophpop, Tspan, X0, method='DOP853', t_eval=t_eval,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )
        orb['seq'][SeqNames[i]]['t'] = sol.t
        orb['seq'][SeqNames[i]]['XJ2000'] = sol.y.T

    elif Seq['type'] == "TBPOptim":
        Period = Seq['T']
        Tspan = (Time, Time + 4 * Period)
        t_eval = np.arange(Time, Time + 4 * Period + TStep, TStep)
        X0 = orb['sat']['X0iner']
        X = TBPOptim(X0, orb)  # Ensure to provide this function

        sol1 = solve_ivp(
            prophpop, Tspan, X0, method='DOP853', t_eval=t_eval,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )
        sol2 = solve_ivp(
            prophpop, Tspan, X, method='DOP853', t_eval=t_eval,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )

        orb['seq'][SeqNames[i]]['t'] = sol1.t
        orb['seq'][SeqNames[i]]['XJ2000'] = sol1.y.T
        orb['seq'][SeqNames[i]]['XC'] = sol2.y.T
        print('The best solution was saved in output')

    elif Seq['type'] == "Lambert":
        span = Seq['span']
        X0 = orb['sat']['X0iner']
        r1 = X0[:3]

        orb = LoadState(Seq['stop'], orb)
        r2 = orb['sat']['X0iner'][:3]

        V1, V2, extremal_distances, exitflag = lambert(
            r1, r2, span / 86400, 0, orb['centralPlanet']['GM']
        )  # Ensure to provide this function
        print(f"DV1 = {np.linalg.norm(V1 - X0[3:6]) * 1e3} m/s")

        Tspan = (Time, Time + span)
        t_eval = np.arange(Time, Time + span + TStep, TStep)
        X0_new = np.concatenate((r1, V1))

        sol = solve_ivp(
            prophpop, Tspan, X0_new, method='DOP853', t_eval=t_eval,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )
        orb['seq'][SeqNames[i]]['t'] = sol.t
        orb['seq'][SeqNames[i]]['XJ2000'] = sol.y.T

    elif Seq['type'] == "LambertOptim":
        t1, t2, X2, span = LambertOptim(orb)  # Ensure to provide this function

        # First propagation
        Tspan1 = (Time, Time + t1)
        t_eval1 = np.arange(Time, Time + t1 + TStep, TStep)
        Time += t1
        X0 = orb['sat']['X0iner']

        sol1 = solve_ivp(
            prophpop, Tspan1, X0, method='DOP853', t_eval=t_eval1,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )
        X1 = sol1.y[:, -1]

        # Lambert computation
        V1, V2, extremal_distances, exitflag = lambert(
            X1[:3], X2[:3], span / 86400, 0, orb['centralPlanet']['GM']
        )
        print(f"DV1 = {np.linalg.norm(V1 - X1[3:6]) * 1e3} m/s")

        # Second propagation
        Tspan2 = (Time, Time + span)
        t_eval2 = np.arange(Time, Time + span + TStep, TStep)
        Time += span
        X0_new = np.concatenate((X1[:3], V1))

        sol2 = solve_ivp(
            prophpop, Tspan2, X0_new, method='DOP853', t_eval=t_eval2,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )
        X2t = sol2.y[:, -1]

        print(f"Diff. pos. target = {np.linalg.norm(X2t[:3] - X2[:3])} km")
        print(f"DV2 = {np.linalg.norm(V2 - X2[3:6]) * 1e3} m/s")

        # Third propagation
        Tspan3 = (Time, Time + t1)
        t_eval3 = np.arange(Time, Time + t1 + TStep, TStep)
        X0_final = np.concatenate((X2t[:3], X2[3:6]))

        sol3 = solve_ivp(
            prophpop, Tspan3, X0_final, method='DOP853', t_eval=t_eval3,
            rtol=options['rtol'], atol=options['atol'], args=(orb,)
        )

        orb['seq'][SeqNames[i]]['t'] = np.concatenate((sol1.t, sol2.t, sol3.t))
        orb['seq'][SeqNames[i]]['XJ2000'] = np.vstack((sol1.y.T, sol2.y.T, sol3.y.T))
        span = t1  # To agree with next lines

    # Update the initial state for the next iteration
    if 'XJ2000' in orb['seq'][SeqNames[i]]:
        orb['sat']['X0iner'] = orb['seq'][SeqNames[i]]['XJ2000'][-1]
    else:
        raise KeyError(f"'XJ2000' not found in sequence {SeqNames[i]}, check the propagation step for issues.")
    Time += span
    orb['seq']['Time'] = Time

def convert_ndarray_to_list(data):
    """ Recursively convert all ndarrays in the dictionary to lists for JSON serialization. """
    if isinstance(data, dict):
        return {key: convert_ndarray_to_list(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_ndarray_to_list(element) for element in data]
    elif isinstance(data, np.ndarray):
        return data.tolist()
    else:
        return data

# Convert all ndarrays in `orb` to lists
orb_serializable = convert_ndarray_to_list(orb)

# Save to JSON
output_path_json = os.path.join('output', 'ORBdata.json')
with open(output_path_json, 'w') as json_file:
    json.dump(orb_serializable, json_file, indent=4)
print(f"Orbital data saved to {output_path_json}")


end_time = time.time()
runtime = end_time - start_time
print("Run Time:", runtime)
