# Model 3: Aqueous Model
# This model contains the complete OpenMC geometry setup for TRIGA reactor
# with aqueous neptunium solution as the target material in the lazy susan

from argparse import ArgumentParser
import numpy as np
import os
import openmc
from openmc import stats
from openmc.model import RightCircularCylinder as RCC
import openmc.deplete

openmc.config['chain_file'] = '/nfs/stak/users/bichselc/Paul TRIGA/aqueous_test/chain_casl_pwr.xml'

# Import the materials module
import OSTR_Materials_Model3 as mats

print("Model 3: Aqueous Model loaded successfully")
print("This model uses aqueous neptunium solution as target material")
print("Target solution: 1.0392 M Np-237 in water")
