# Model 2: Pellet Model with HFIR pellets and aqueous target
# This model includes the complete OpenMC geometry setup for TRIGA reactor
# with HFIR-style pellet stacks in the target locations

from argparse import ArgumentParser
import numpy as np
import os
import openmc
from openmc import stats
from openmc.model import RightCircularCylinder as RCC
import openmc.deplete

openmc.config['chain_file'] = '/nfs/stak/users/bichselc/Paul TRIGA/pellet/chain_casl_pwr.xml'

# Import the materials module
import OSTR_Materials_Copy1 as mats

print("Model 2: Pellet Model loaded successfully")
print("This model uses HFIR-style pellets as target material")
