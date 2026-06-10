   ____    _____  ______    ____            _                   ____                            __  ___   ______
  / __ \  / ___/ /_  __/   / __ \          (_)   ____          / __ \    ____   ___    ____    /  |/  /  / ____/
 / / / /  \__ \   / /     / /_/ /         / /   / __ \        / / / /   / __ \ / _ \  / __ \  / /|_/ /  / /     
/ /_/ /  ___/ /  / /     / _, _/         / /   / / / /       / /_/ /   / /_/ //  __/ / / / / / /  / /  / /___   
\____/  /____/  /_/     /_/ |_|         /_/   /_/ /_/        \____/   / .___/ \___/ /_/ /_/ /_/  /_/   \____/   
                                                                     /_/                                        

================
|              |
|   READ ME    |
|              |
================

Contained in these Model files are two files, not including this README.

	(1) OSTR_OpenMC_Model.py
			- 	Requires the following python packages:
					-	numpy
					-	os
					-	argparse
					- 	openmc
			- 	Constructed with OpenMC v0.15.0
			- 	The primary script
			- 	Generates the input files for the "openmc" executable to run
			- 	Imports OSTR_Materials.py for the materials definitions
	(2) OSTR_Materials.py
			- 	Requires the following python packages:
					-	os
					- 	openmc
			- 	Contains the materials definitions for use in the primary script


*******************************************************************************************
*  To run the OSTR OpenMC Model, you must place files (1) and (2) in the same directory.  * 
*******************************************************************************************


The model can be run with python using the following command line arguments:

    -n              = (int)         Number of particle histories to track per cycle
    -i              = (int)         Number of inactive cycles to run before tracking statistics
    -a              = (int)         Number of active cycles for which statistics are tracked
    -m              = (bool)        Generate the model with temperatures at full power (1 MW)
    --run           = (bool)        Run the OpenMC model after it is generated
    --plot          = (bool)        Plot any specified plots from the model after it is generated
    --rods          = (float) x4    Four floats for the height of the control rods, i.e. <TR> <SA> <SH> <REG>
    --num_threads   = (int)         Number of OpenMP threads to run OpenMC with
    --dir           = (str)         Path for which to generate AND run OpenMC or OpenMC plotter(if specified)

Integer, Float, and String arguments require values after them. Booleans are True if included, False if not. 

An example command:
 ____________________________________________________________________________________________________________________________
|                                                                                                                            |
| python OSTR_OpenMC_Model.py -n 25000 -i 30 -a 470 -m --run --plot --rods 30 40.1 50.2 60.3 --num_threads 16 --dir ./folder |
|____________________________________________________________________________________________________________________________|

This will generate a model of the OSTR with the Transient Rod 30% withdrawn, Safety Rod 40.1% withdrawn, Shim Rod 50.2%
withdrawn, and the Regulating Rod 60.3% withdrawn with temperatures set to those at 1 MW (full power). This model is
set to track 25000 particles per cycle, with 30 inactive cycles and 470 active cycles (total of 500). The input 
files will be generated in "./folder". OpenMC will then run in that same directory using 16 threads, and plots 
will be generated in that same directory.
