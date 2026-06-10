import openmc
import os
import numpy as np

# os.environ['OPENMC_CROSS_SECTIONS'] = '/home/chase/openmc/Cross_Section_Libraries/endfb71_hdf5/cross_sections.xml'

cmap = {}

##########################################################################################################
#                                                                                                        #
# --------------------------------------------- MATERIALS ---------------------------------------------- #
#                                                                                                        #
##########################################################################################################

##############
#   TARGET   #
##############

target = openmc.Material()

target.add_nuclide('Np237', 1)
target.add_element('O', 2)

target.set_density('g/cm3', 11.10)
target.name = 'Target'
target.id = 9300
target.depletable = True
target.volume = 599.2388

cmap[target] = (255, 0, 0)

############
#   FUEL   #
############

# 19.75% Enriched, 30 wt% Uranium, 1.07 wt% Erbium, 1.6:1 Hydrogen to Zirconium (mole ratio)

fuel = openmc.Material()

fuel.add_nuclide('U235', 5.925, 'wo')
fuel.add_nuclide('U238', 24.075, 'wo')
fuel.add_element('H', 1.197294, 'wo')
fuel.add_element('Zr', 67.73271, 'wo')
fuel.add_element('Er', 1.1, 'wo')

fuel.add_s_alpha_beta('c_H_in_ZrH')
fuel.add_s_alpha_beta('c_Zr_in_ZrH')

fuel.set_density('g/cm3', 7.1879)
fuel.name = 'Fuel'
fuel.id = 9200
fuel.depletable = False

cmap[fuel] = (170, 175, 170)

################
#   GRAPHITE   #
################

# Nuclear Grade Graphite (10% Porosity)

graphite = openmc.Material()

graphite.add_element('C', 100, 'ao')

graphite.add_s_alpha_beta('c_Graphite')

graphite.set_density('g/cm3', 1.71)
graphite.name = 'Graphite'
graphite.id = 600

cmap[graphite] = (60, 60, 60)


#################
#   ZIRCONIUM   #
#################

# Zirconium Metal

zirconium = openmc.Material()

zirconium.add_element('Zr', 100, 'ao')

zirconium.set_density('g/cm3', 1.75)
zirconium.name = 'Zirconium'
zirconium.id = 4000

cmap[zirconium] = (205, 205, 200)


##################
#   MOLYBDENUM   #
##################

# Molybdenum Metal

molybdenum = openmc.Material()

molybdenum.add_element('Mo', 100, 'ao')

molybdenum.set_density('g/cm3', 10.223)
molybdenum.name = 'Molybdenum'
molybdenum.id = 4200

cmap[molybdenum] = (150, 150, 170)


#####################
#   BORON CARBIDE   #
#####################

# Boron Carbide Sintered Pellets (B4C)

b4c = openmc.Material()

b4c.add_element('B', 4, 'ao')
b4c.add_element('C', 1, 'ao')

b4c.set_density('g/cm3', 2.49)
b4c.name = 'Boron Carbide'
b4c.id = 500

cmap[b4c] = (40, 40, 50)

#############
#   WATER   #
#############

# Water at 31°C and 1.6 atm (at roughly 20 feet deep)

water = openmc.Material()

water.add_element('H', 2, 'ao')
water.add_element('O', 1, 'ao')

water.add_s_alpha_beta('c_H_in_H2O')
water.add_s_alpha_beta('c_D_in_D2O')

water.set_density('g/cm3', 0.995)
water.name = 'Water'
water.id = 100

cmap[water] = (140, 160, 215)


######################
#  6061-T6 Aluminum  #
######################

# From "Compendium of Material Composition Data for Radiation Transport Modeling" - PNNL (2021)

aluminum = openmc.Material()

aluminum.add_element('Al', 97.200, 'wo')
aluminum.add_element('Mg', 1.0000, 'wo')
aluminum.add_element('Si', 0.6000, 'wo')
aluminum.add_element('Ti', 0.0876, 'wo')
aluminum.add_element('Cr', 0.1950, 'wo')
aluminum.add_element('Mn', 0.0876, 'wo')
aluminum.add_element('Fe', 0.4088, 'wo')
aluminum.add_element('Cu', 0.2750, 'wo')
aluminum.add_element('Zn', 0.1460, 'wo')

aluminum.add_s_alpha_beta('c_Al27')

aluminum.set_density('g/cm3', 2.7)
aluminum.name = '6061-T6 Aluminum'
aluminum.id = 1300

cmap[aluminum] = (180, 180, 180)


##############
#  Titanium  #
##############

# Titanium Metal

titanium = openmc.Material()

titanium.add_element('Ti', 100, 'ao')

titanium.set_density('g/cm3', 4.506)
titanium.name = 'Titanium'
titanium.id = 2200

cmap[titanium] = (180, 190, 180)


#############################
#  SAE 304 Stainless Steel  #
#############################

# From "Compendium of Material Composition Data for Radiation Transport Modeling" - PNNL (2021)

steel = openmc.Material()

steel.add_element('Fe', 68.345, 'wo')
steel.add_element('C', 0.0800, 'wo')
steel.add_element('Si', 1.0000, 'wo')
steel.add_element('P', 0.0450, 'wo')
steel.add_element('S', 0.0300, 'wo')
steel.add_element('Cr', 19.000, 'wo')
steel.add_element('Mn', 2.0000, 'wo')
steel.add_element('Ni', 9.5000, 'wo')

steel.add_s_alpha_beta('c_Fe56')

steel.set_density('g/cm3', 8.03)
steel.name = 'Stainless Steel 304'
steel.id = 2600

cmap[steel] = (150, 130, 130)


#################################
#  Steel/Water Mix for Flutes   #
#################################

# From the OSTR MCNP Deck. Unsure of the exact ratios used for this mixture, but it comes from blending
# the geometry of the flutes (made of steel) into the water that passes by them.

flute_mix = openmc.Material()

# Water
flute_mix.add_element('H', 0.48777, 'ao')
flute_mix.add_element('O', 0.04676, 'ao')

# Steel
flute_mix.add_element('Fe', 0.17925, 'ao')
flute_mix.add_element('C', 0.00098, 'ao')
flute_mix.add_element('Si', 0.00521, 'ao')
flute_mix.add_element('P', 0.00021, 'ao')
flute_mix.add_element('S', 0.00014, 'ao')
flute_mix.add_element('Cr', 0.05352, 'ao')
flute_mix.add_element('Mn', 0.00533, 'ao')
flute_mix.add_element('Ni', 0.02371, 'ao')

flute_mix.add_s_alpha_beta('c_H_in_H2O')
flute_mix.add_s_alpha_beta('c_D_in_D2O')
flute_mix.add_s_alpha_beta('c_Fe56')

flute_mix.set_density('atom/b-cm', 0.095865)
flute_mix.name = 'Homogeneous Flute'
flute_mix.id = 26100

cmap[flute_mix] = (200, 195, 190)

#########
#  AIR  #
#########

# Dry Air at STP

air = openmc.Material()

air.add_element('C', 0.0150, 'ao')
air.add_element('N', 78.4429, 'ao')
air.add_element('O', 21.0750, 'ao')
air.add_element('Ar', 0.4671, 'ao')

air.set_density('kg/m3', 1.205)
air.name = 'Air'
air.id = 700

cmap[air] = (220, 220, 220)

#############
#  CADMIUM  #
#############

# Cadmium Metal

cadmium = openmc.Material()

cadmium.add_element('Cd', 100, 'ao')

cadmium.set_density('g/cm3', 8.649)
cadmium.name = 'Cadmium'
cadmium.id = 4800

cmap[cadmium] = (180, 170, 180)

##########
#  LEAD  #
##########

# Lead Metal

lead = openmc.Material()

lead.add_element('Pb', 100, 'ao')

lead.set_density('g/cm3', 10.95)
lead.name = 'Lead'
lead.id = 8200

cmap[lead] = (70, 70, 70)
