# High-precision Analyser of Lunar Orbits

# Description

HALO is a mission design tool build on a high accuracy propagator for lunar orbits. It intends to help to study any orbits in the vicinity of the Moon. HALO comes with a scientific paper that precisely define all the models and algorithms used and go through a sensitivity analysis on lunar orbits of interest. This paper High-precision_Analyser_of_Lunar_Orbit.pdf can be found at the root of the GitHub, it is highly recommended to get a good look at the paper when using HALO.

# Initialisation

After copying the project on your local machine (git clone [url]), you have to:

1. install the dependencies listed in requirements.txt by running:
```language
# pip install -r requirements.txt
```
2. Download the de430.bsp ephemeris at [NAIF de430](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/).   
   Place it in MATLAB/HALO/ker and in Python/input.


### List of all Python modules
- matplotlib
- spiceypy
- numpy
- tqdm
- pandas
- pymatreader
- pyparsing
- python-dateutil
- pytz
- scipy

# Usage

The project's code is developed in python.

The main usage of HALO is to construct a mission by changing the input files like LoadSequential.py on PYTHON, then run the mission by running the main script HALO.py, and finally use python scripts like Visualisation.py to process the output file. Because of relative path, HALO.py as to be run from the python folder.
Read the following parts to learn what can be done with HALO.

How to run the project:
```language
# python HALO.py
```

### Architechture of the Python
- "HALO.py" is the main file that allows to run the tool.  
- "metakernelcheck.py" creates a metakernel.tm file that specifies essential SPICE kernel files for loading planetary, time, and orientation data.  
- "Visualisation.py" allows to plot a propagation after processing, the user can choose the frame, or choose to add the converged trajectory if there is one, same with a CR3BP trajectory.  
- "TrajectoryErrors.py" allows to compare the propagation of a true ephemeris state of a satellite with its real evolution.   
- "CR3BPfitting.py" converts the initial state from a CR3BP orbit chosen on https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/intro in the Earth-Moon rotationnal and barycentered frame into the ephemeris Moon centered inertial frame in order to use it as an initial state for the propagator.   
- "process.py" gathers all the processing functions used in the other python files.

- Directories:
  - input: This directory contains scripts for loading and processing input data for the simulation, such as initial conditions or configuration files.
  - ker: This may contain SPICE kernel files, which store data related to planetary positions, instrument information, or constants.
  - prop: This directory contain various trajectory propagation and optimization algorithms, such as TBPOptim, accelpntmasses, accelsrp, etc., for different modeling needs.
  - output: The output folder contains all the information on the processed mission.

### Process an ephemeris 3-body-problem orbit
You can get a CR3BP orbit from the [JPL website](https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/intro). This orbit is a solution in a simplified model and in the Earth-Moon rotational frame.   
CR3BPfitting.py will fit the simplified model on an ephemeris one at a specific epoch (chosen at 1 Jan 2024 in the example) and convert the initial position to be used for propagation.   
Then a sequence of type "TBPOptim" allows to compute an ephemeris orbit from this fitted initial guess.    
The following sequential gives an example for a NRHO close to the Capstone one:  

   orb = LoadState("NRHO", orb)  
   orb.seq.Time = sp.utc2et("2024-01-01T00:00:00")  # Convert UTC to ET  
   orb.seq.a.type = "TBPOptim"  
   orb.seq.a.T = 5.7444e+05  # Orbital period in seconds  

Finally, Visualisation.py allows to visualise the initial guess and the converged solution by putting the Converged option to 1. The RotationalF option can also be put to 1 to observe the closed orbit in the rotational frame.   

The same can be done for DRO14, with a period of 14 days, but the convergence is a bit harder so that we have to change n=2 (in the algorithm's file) and use only two periods due to the long period of the DRO.

### Some useful functions
sp.dafec(sp.spklef("input/LRO_ES_90_202003_GRGM900C_L600.bsp"),1000)    :    Allows to assess the comment of a bsp file
sp.datetime2et(sp.datetime(y, m, d, h, m, s))                           :    Converts a UTC time to ET

See also part 4.3.5 of the paper for time and coordinate systems handling in either MATLAB and Python.

# Authors and acknowledgment

Creator of the Github and first contributor to the project: Quentin Granier, Research Associate.  
Supervisor of the Project: Special thanks to Dr. Yang Yang for his continuous support and guidance throughout this project. Contact: yiyinfeixiong@gmail.com.
