
from argparse import ArgumentParser
import numpy as np
import os

import openmc
from openmc import stats
from openmc.model import RightCircularCylinder as RCC
import openmc.deplete

openmc.config['chain_file'] = '/nfs/stak/users/bichselc/Paul TRIGA/pellet/chain_casl_pwr.xml' # This will have to be editted to your files for modifications

def ostr_model(n_particles, n_inactive, n_active, rod_heights, megawatt=False):
    """
    Returns an openmc.Model() of the Oregon State TRIGA Reactor.

    :param n_particles:     (int) number of particles
    :param n_inactive:      (int) number of inactive cycles
    :param n_active:        (int) number of active cycles
    :param rod_heights:     (tuple) a tuple containing the four control rod heights (Tr, Sa, Sh, Reg)
    :param megawatt:        (bool) if True, change temperatures to operating temps

    :return:                (openmc.Model) OpenMC model of the Oregon State TRIGA Reactor

    """

    # Some print statements for important info

    print(f"\n\n   ___    ___   _   _     _____   ___   ___    ___     _       __  __   _           ___   ___ \n"
          "  / _ \\  / __| | | | |   |_   _| | _ \\ |_ _|  / __|   /_\\     |  \\/  | | |__       |_ _| |_ _|\n"
          " | (_) | \\__ \\ | |_| |     | |   |   /  | |  | (_ |  / _ \\    | |\\/| | | / /  _     | |   | | \n"
          "  \\___/  |___/  \\___/      |_|   |_|_\\ |___|  \\___| /_/ \\_\\   |_|  |_| |_\\_\\ (_)   |___| |___|\n\n")

    # Instantiate the model class
    model = openmc.Model()

    print("Rod Heights:\n")
    print(f"Transient: {rod_heights[0]:1.0f}%")
    print(f"Safety: {rod_heights[1]:1.0f}%")
    print(f"Shim: {rod_heights[2]:1.0f}%")
    print(f"Regulating: {rod_heights[3]:1.0f}%\n")

    # These are the materials in a separate file.
    import OSTR_Materials_Copy1 as mats

    print("Materials Generated.\n")

    # Default false, but if true, model will be at 1 MW temps
    if megawatt:
        print("Full Power (1 MW) Temperatures:\n")
        # Temperatures at 1 MW
        bulk_water_temperature = 308.15  # 35°C
        core_water_temperature = 393.15  # 100°C
        fuel_temperature = 623.15  # 350°C
        print(f"Average Fuel Temperature is {(fuel_temperature - 273.15):1.2f}ºC")
        print(f"Average Core Temperature is {(core_water_temperature - 273.15):1.2f}ºC")
        print(f"Average Tank Temperature is {(bulk_water_temperature - 273.15):1.2f}ºC\n")
    else:
        print("Low Power (<10 kW) Temperatures:\n")
        # Temperatures at 100 W (All ~31°C)
        bulk_water_temperature = 304.15
        core_water_temperature = 304.15
        fuel_temperature = 304.15
        print(f"Average Fuel Temperature is {(fuel_temperature - 273.15):1.2f}ºC")
        print(f"Average Core Temperature is {(core_water_temperature - 273.15):1.2f}ºC")
        print(f"Average Tank Temperature is {(bulk_water_temperature - 273.15):1.2f}ºC\n")

    # Set the control rod heights
    trans_pos = rod_heights[0]
    safety_pos = rod_heights[1]
    shim_pos = rod_heights[2]
    reg_pos = rod_heights[3]


    
    def make_hfir_pellet(center, temp):
        x0, y0, z0 = center
    
    # Dimensions in cm
        hp = 0.5 * 2.54          # Total Pellet Height
        ht = 0.05 * 2.54         # Taper Height (for 45 deg)
        ro = 0.165 * 2.54        # Outer Radius
        ri = 0.115 * 2.54        # Inner Radius (at the bore)

    #    Surfaces 
    # Outer Cylinder
        s_outer = openmc.ZCylinder(x0=x0, y0=y0, r=ro)
    
    # Top and Bottom Planes of the single pellet
        p_bot = openmc.ZPlane(z0=z0)
        p_top = openmc.ZPlane(z0=z0 + hp)
    
    # Internal Bore Cylinder
        s_bore = openmc.ZCylinder(x0=x0, y0=y0, r=ri)
    
    # Taper Cones (45 degrees means radius increases 1:1 with height)
    # Bottom Taper: starts at ri at p_bot, opens up to ro
        s_taper_bot = openmc.Cone(x0=x0, y0=y0, z0=z0 - (ro - ri), r2=1.0)
    
    # Top Taper: starts at ri at p_top, opens up to ro
        s_taper_top = openmc.Cone(x0=x0, y0=y0, z0=z0 + hp + (ro - ri), r2=1.0)

    #     Regions 
    # Pellet Meat: Inside outer cyl, outside both cones, and outside the bore
    # Plus, it must be between the top and bottom planes
        meat_region = (-s_outer & +p_bot & -p_top & +s_taper_bot & -s_taper_top & +s_bore)
    
    # Void Region: Everything else inside the outer envelope
        void_region = (-s_outer & +p_bot & -p_top & ~meat_region)

    #     Cells 
        c_meat = openmc.Cell(fill=mats.npo2, region=meat_region)
        c_void = openmc.Cell(fill=mats.air, region=void_region)
    
        return openmc.Universe(cells=[c_meat, c_void])



    
    # Functions to generate core components
    def make_element(position, bounding_cyl, temp, element_type):
        """
        Generates a universe containing the requested element at the position specified by the bounding cylinder.

        :param position:        (int) Grid Position (201 is B1, 302 is C2, etc.) only for naming/ID/debugging
        :param bounding_cyl:    (openmc.model.RightCircularCylinder) cylinder bounding the cell
        :param temp:            (float) temperature in [K]
        :param element_type:    (str) 'fuel', 'refl', 'source', 'ct', 'clicit', 'icit', 'rabbit'

        :return:                (openmc.Universe) universe containing the desired element_type
        """

        # X and Y coordinate
        x0 = bounding_cyl.cyl.x0
        y0 = bounding_cyl.cyl.y0

        # Z coordinate (plane equation is Ax + By + Cz + D = 0, therefore D/C is the z-coordinate.)
        z0 = bounding_cyl.bottom.d / bounding_cyl.bottom.c

        # Calculate the bounding cylinder
        s_bound = RCC(center_base=(x0, y0, z0), height=65.78, radius=1.87325, axis='z')
        c_bound = openmc.Cell(fill=mats.water, region=(+s_bound))

        if element_type == 'fuel':
            '''
            
            From the bottom up, a fuel element is built as follows:
                
                (1)     Bottom flute (modeled as homogeneous mixture)
                (2)     Bottom axial reflector
                (3)     Molybdenum disk (1/32" thick)
                (4a)    A 15" section of TRIGA fuel
                (4b)    A 15" zirconium rod through the TRIGA fuel
                (5)     Top axial reflector (assumed same size as bottom, it is not actually)
                (6)     An expansion gap
                (7)     Top flute (modeled as homogeneous mixture)
                
            The cladding extends from the top of the bottom flute to the bottom of the top flute.
            It is 0.020" in thickness. There is assumed to be no radial gap between the fuel components
            and the inner surface of the cladding.
            
            '''

            # Fuel Component Dimensions
            bottom_flute_height = 4.54
            graph_height = 9.665
            moly_height = 0.079375
            fuel_height = 38.1
            upper_graph_height = 9.585625
            void_gap_height = 1.27
            upper_flute_height = 2.54
            clad_height = graph_height + moly_height + fuel_height + upper_graph_height + void_gap_height

            # Build z-coordinates for bottom of cylinders for each component
            bottom_flute_z = z0
            bottom_graph_z = bottom_flute_z + bottom_flute_height
            moly_z = bottom_graph_z + graph_height
            fuel_meat_z = moly_z + moly_height
            upper_graph_z = fuel_meat_z + fuel_height
            void_gap_z = upper_graph_z + upper_graph_height
            upper_flute_z = void_gap_z + void_gap_height

            # Surfaces for each component
            s_bottom_flute = RCC(center_base=(x0, y0, bottom_flute_z), height=bottom_flute_height, radius=1.87325, axis='z',
                                 name=f"Bottom Flute of Fuel Element in Position {position}")
            s_bottom_refl = RCC(center_base=(x0, y0, bottom_graph_z), height=graph_height, radius=1.82245, axis='z',
                                name=f"Bottom Graphite Reflector of Fuel Element in Position {position}")
            s_moly_disk = RCC(center_base=(x0, y0, moly_z), height=moly_height, radius=1.82245, axis='z',
                              name=f" Moly Disk in of Fuel Element in Position {position}")
            s_fuel_meat = RCC(center_base=(x0, y0, fuel_meat_z), height=fuel_height, radius=1.82245, axis='z',
                              name=f"Fuel Meat of Fuel Element in Position {position}")
            s_zirc_rod = RCC(center_base=(x0, y0, fuel_meat_z), height=fuel_height, radius=0.3175, axis='z',
                             name=f"Zirconium Rod of Fuel Element in Position {position}")
            s_top_refl = RCC(center_base=(x0, y0, upper_graph_z), height=graph_height, radius=1.82245, axis='z',
                             name=f"Top Reflector of Fuel Element in Position {position}")
            s_void_gap = RCC(center_base=(x0, y0, void_gap_z), height=void_gap_height, radius=1.82245, axis='z',
                             name=f"Expansion Gap of Fuel Element in Position {position}")
            s_upper_flute = RCC(center_base=(x0, y0, upper_flute_z), height=upper_flute_height, radius=1.87325, axis='z',
                                name=f"Top Flute of Fuel Element in Position {position}")
            s_cladding = RCC(center_base=(x0, y0, bottom_graph_z), height=clad_height, radius=1.87325, axis='z',
                             name=f"Cladding of Fuel Element in Position {position}")

            # Cells for each component
            c_bottom_flute = openmc.Cell(fill=mats.flute_mix, region=(-s_bottom_flute),
                                         name=f"Bottom Flute of FE in Pos. {position}")
            c_bottom_refl = openmc.Cell(fill=mats.graphite, region=(-s_bottom_refl),
                                        name=f"Bottom Graphite Slug of FE in Pos. {position}")
            c_moly_disk = openmc.Cell(fill=mats.molybdenum, region=(-s_moly_disk),
                                      name=f"Moly Disk of FE in Pos. {position}")
            c_fuel_meat = openmc.Cell(fill=mats.fuel, region=(-s_fuel_meat & +s_zirc_rod),
                                      name=f"Fuel Meat of FE in Pos. {position}")
            c_zirc_rod = openmc.Cell(fill=mats.zirconium, region=(-s_zirc_rod),
                                     name=f"Zirconium Rod of FE in Pos. {position}")
            c_top_refl = openmc.Cell(fill=mats.graphite, region=(-s_top_refl),
                                     name=f"Top Graphite Slug of FE in Pos. {position}")
            c_void_gap = openmc.Cell(fill=mats.air, region=(-s_void_gap),
                                     name=f"Expansion Gap of FE in Pos. {position}")
            c_upper_flute = openmc.Cell(fill=mats.flute_mix, region=(-s_upper_flute),
                                        name=f"Top Flute of FE in Pos. {position}")
            c_cladding = openmc.Cell(fill=mats.steel,
                                     region=(-s_cladding & +s_bottom_refl & +s_fuel_meat & +s_top_refl & +s_void_gap),
                                     name=f"Cladding of FE in Pos. {position}")

            # List of cells in fuel element
            fuel_rod_cells = [c_bottom_flute,
                              c_bottom_refl,
                              c_moly_disk,
                              c_fuel_meat,
                              c_zirc_rod,
                              c_top_refl,
                              c_void_gap,
                              c_upper_flute,
                              c_cladding,
                              c_bound]

            for cell in fuel_rod_cells:
                cell.temperature = temp

            # Make fuel element universe
            fuel_element = openmc.Universe(universe_id=position, name=f"Fuel Element in Position {position}",
                                           cells=fuel_rod_cells)
            return fuel_element

        elif element_type == 'refl':
            '''

            From the bottom up, a reflector element is built as follows:
                (1)     Bottom flute (modeled as homogeneous mixture)
                (2)     A graphite slug
                (6)     An expansion gap
                (7)     Top flute (modeled as homogeneous mixture)

            The cladding extends from the top of the bottom flute to the bottom of the top flute.
            It is 0.020" in thickness. There is assumed to be no radial gap between the components
            and the inner surface of the cladding.

            '''

            # Refl Element Component Dimensions
            bottom_flute_height = 4.54
            graph_height = 9.665 + 38.1 + 9.665
            void_gap_height = 1.27
            upper_flute_height = 2.54
            clad_height = graph_height + void_gap_height

            # Build z-coordinates for bottom of cylinders for each component
            bottom_flute_z = z0
            graph_z = bottom_flute_z + bottom_flute_height
            void_gap_z = graph_z + graph_height
            upper_flute_z = void_gap_z + void_gap_height

            # Surfaces for each component
            s_bottom_flute = RCC(center_base=(x0, y0, bottom_flute_z), height=bottom_flute_height, radius=1.87325,
                                 axis='z',
                                 name=f"Bottom Flute of Reflector Element in Position {position}")
            s_graphite = RCC(center_base=(x0, y0, graph_z), height=graph_height, radius=1.82245, axis='z',
                             name=f"Graphite of Reflector Element in Position {position}")
            s_void_gap = RCC(center_base=(x0, y0, void_gap_z), height=void_gap_height, radius=1.82245, axis='z',
                             name=f"Expansion Gap of Reflector Element in Position {position}")
            s_upper_flute = RCC(center_base=(x0, y0, upper_flute_z), height=upper_flute_height, radius=1.87325, axis='z',
                                name=f"Top Flute of Reflector Element in Position {position}")
            s_cladding = RCC(center_base=(x0, y0, graph_z), height=clad_height, radius=1.87325, axis='z',
                             name=f"Cladding of Reflector Element in Position {position}")

            # Cells for each component
            c_bottom_flute = openmc.Cell(fill=mats.flute_mix, region=(-s_bottom_flute),
                                         name=f"Bottom Flute of RE in Pos. {position}")
            c_graphite = openmc.Cell(fill=mats.graphite, region=(-s_graphite),
                                     name=f"Graphite of RE in Pos. {position}")
            c_void_gap = openmc.Cell(fill=mats.air, region=(-s_void_gap),
                                     name=f"Expansion Gap of RE in Pos. {position}")
            c_upper_flute = openmc.Cell(fill=mats.flute_mix, region=(-s_upper_flute),
                                        name=f"Top Flute of RE in Pos. {position}")
            c_cladding = openmc.Cell(fill=mats.steel,
                                     region=(-s_cladding & +s_graphite & +s_void_gap),
                                     name=f"Cladding of FE in Pos. {position}")

            # List of cells in fuel element
            refl_element_cells = [c_bottom_flute,
                                  c_graphite,
                                  c_void_gap,
                                  c_upper_flute,
                                  c_cladding,
                                  c_bound]

            for cell in refl_element_cells:
                cell.temperature = temp

            # Make fuel element universe
            refl_element = openmc.Universe(universe_id=position, name=f"Reflector Element in Position {position}",
                                           cells=refl_element_cells)

            return refl_element

        elif element_type == 'source':
            '''
            
            Identical to a reflector element, but filled with air and clad in aluminum instead of steel.
            
            '''

            # Source Holder Dimensions
            bottom_flute_height = 4.54
            air_height = 9.665 + 38.1 + 9.665 + 1.27
            upper_flute_height = 2.54
            clad_height = air_height

            # Build z-coordinates for bottom of cylinders for each component
            bottom_flute_z = z0
            air_z = bottom_flute_z + bottom_flute_height
            upper_flute_z = air_z + air_height

            # Surfaces for each component
            s_bottom_flute = RCC(center_base=(x0, y0, bottom_flute_z), height=bottom_flute_height, radius=1.87325,
                                 axis='z', name=f"Bottom Flute of Source Holder in Position {position}")
            s_air = RCC(center_base=(x0, y0, air_z), height=air_height, radius=1.82245, axis='z',
                        name=f"Air in Source Holder in Position {position}")
            s_upper_flute = RCC(center_base=(x0, y0, upper_flute_z), height=upper_flute_height, radius=1.87325, axis='z',
                                name=f"Top Fitting of Source Holder in Position {position}")
            s_cladding = RCC(center_base=(x0, y0, air_z), height=clad_height, radius=1.87325, axis='z',
                             name=f"Cladding of Source Holder in Position {position}")

            # Cells for each component
            c_bottom_flute = openmc.Cell(fill=mats.aluminum, region=(-s_bottom_flute),
                                         name=f"Bottom Fitting of Source in Pos. {position}")
            c_air = openmc.Cell(fill=mats.air, region=(-s_air), name=f"Air in Source in Pos. {position}")
            c_upper_flute = openmc.Cell(fill=mats.aluminum, region=(-s_upper_flute),
                                        name=f"Top Fitting of Source in Pos. {position}")
            c_cladding = openmc.Cell(fill=mats.aluminum, region=(-s_cladding & +s_air),
                                     name=f"Cladding of FE in Pos. {position}")

            # List of cells in fuel element
            source_holder_cells = [c_bottom_flute,
                                   c_air,
                                   c_upper_flute,
                                   c_cladding,
                                   c_bound]

            for cell in source_holder_cells:
                cell.temperature = temp

            # Make fuel element universe
            source_holder = openmc.Universe(universe_id=position, name=f"Source Holder in Position {position}",
                                            cells=source_holder_cells)

            return source_holder

        elif element_type == 'clicit':
            '''
            
            A CLICIT tube is a series of nested tubes. From inside to outside, there is the air
            in the tube, followed by a thin aluminum liner. This inner tube is wrapped in cadmium metal.
            This cadmium metal is encased in an outer aluminum tube.
            
            '''

            # CLICIT (Will Be Translated to Correct Position)
            s_clicit_tube = RCC(center_base=(x0, y0, z0), height=65.78, radius=1.87325, axis='z')
            s_clicit_cad = RCC(center_base=(x0, y0, z0 + 4.54), height=61.24, radius=1.63880, axis='z')
            s_clicit_al = RCC(center_base=(x0, y0, z0), height=85.78, radius=1.58750, axis='z')
            s_clicit_air = RCC(center_base=(x0, y0, z0 + 4.6797), height=81.1003, radius=1.44018, axis='z')

            bound_region = (-s_clicit_tube | -s_clicit_al)

            c_bound = openmc.Cell(fill=mats.water, region=(~bound_region))

            # CLICIT Tube
            c_clicit_air = openmc.Cell(fill=mats.air, region=(-s_clicit_air),
                                       name=f"Air in CLICIT in {position}")
            c_clicit_tube = openmc.Cell(fill=mats.aluminum, region=(-s_clicit_al & +s_clicit_air),
                                        name=f"CLICIT Tube in {position}")
            c_clicit_cad = openmc.Cell(fill=mats.cadmium, region=(-s_clicit_cad & +s_clicit_al),
                                       name=f"CLICIT Cadmium in {position}")
            c_clicit_clad = openmc.Cell(fill=mats.aluminum, region=(-s_clicit_tube & +s_clicit_cad),
                                        name=f"In-Core Portion of CLICIT in {position}")
            c_water_around_clicit = openmc.Cell(fill=mats.water, region=(+s_clicit_tube),
                                                name=f"Water around CLICIT {position}")

            clicit_cells = [c_clicit_air, c_clicit_tube, c_clicit_cad, c_clicit_clad, c_water_around_clicit, c_bound]

            for cell in clicit_cells:
                cell.temperature = temp

            clicit = openmc.Universe(universe_id=position, name=f"CLICIT in {position}", cells=clicit_cells)

            return clicit

        elif element_type == 'ct':
            '''
            
            The central thimble is a tube filled with an aluminum plug in the center of the core with 
            a 1/4" diameter, water-filled hole.
            
            '''

            # Central Thimble (Will Be Translated to Correct Position)
            s_central_thimble = RCC(center_base=(x0, y0, z0), height=85.78, radius=1.87325, axis='z')
            s_ct_plug = RCC(center_base=(x0, y0, z0), height=65.78, radius=1.63880, axis='z')
            s_ct_tube = RCC(center_base=(x0, y0, z0), height=85.78, radius=1.63880, axis='z')
            s_ct_hole = RCC(center_base=(x0, y0, z0 + 4.54), height=61.24, radius=0.3175, axis='z')

            # Central Thimble
            c_central_thimble_tube = openmc.Cell(fill=mats.aluminum, region=(-s_central_thimble & +s_ct_tube),
                                                 name="Central Thimble Tube")
            c_central_thimble_water = openmc.Cell(fill=mats.water, region=(-s_ct_tube & +s_ct_plug.top),
                                                  name="Water In CT Tube")
            c_central_thimble_hole = openmc.Cell(fill=mats.water, region =(-s_ct_hole),
                                                 name="Water in CT Hole")
            c_central_thimble_plug = openmc.Cell(fill=mats.aluminum, region=(-s_ct_plug & +s_ct_hole),
                                                 name="Central Thimble Plug")
            c_ct_bound = openmc.Cell(fill=mats.water, region=(+s_central_thimble))

            u_ct = openmc.Universe(cells=[c_central_thimble_water, c_central_thimble_plug, c_central_thimble_tube,
                                          c_central_thimble_hole, c_ct_bound], name="Central Thimble")
            return u_ct

        elif element_type == 'icit':
            '''
            
            Identical to CLICIT without the cadmium metal.
            
            '''

            # ICIT
            s_in_core_tube = RCC(center_base=(x0, y0, z0), height=85.78, radius=1.87325, axis='z')
            s_in_core_air = RCC(center_base=(x0, y0, z0 + 4.54), height=81.24, radius=1.63880, axis='z')

            # ICIT
            c_icit_tube = openmc.Cell(fill=mats.aluminum, region=(-s_in_core_tube & +s_in_core_air), name="ICIT Tube")
            c_icit_air = openmc.Cell(fill=mats.air, region=(-s_in_core_air), name="ICIT Air")
            c_icit_bound = openmc.Cell(fill=mats.water, region=(+s_in_core_tube), name="ICIT Bound")
            u_ICIT = openmc.Universe(cells=[c_icit_air, c_icit_tube, c_icit_bound], name="ICIT")
            return u_ICIT

        elif element_type == 'rabbit':
            '''
            
            Identical to ICIT, but made of titanium.
            
            '''


            # RABBIT
            s_rabbit_tube = RCC(center_base=(x0, y0, z0), height=85.78, radius=1.87325, axis='z')
            s_rabbit_air = RCC(center_base=(x0, y0, z0+4.54), height=81.24, radius=1.63880, axis='z')

            c_rabbit_tube = openmc.Cell(fill=mats.titanium, region=(-s_rabbit_tube & +s_rabbit_air), name="Rabbit Tube")
            c_rabbit_air = openmc.Cell(fill=mats.air, region=(-s_rabbit_air), name="Rabbit Air")
            c_rabbit_bound = openmc.Cell(fill=mats.water, region=(+s_rabbit_tube), name="Rabbit Bound")
            u_rabbit = openmc.Universe(cells=[c_rabbit_tube, c_rabbit_air, c_rabbit_bound], name="Rabbit")
            return u_rabbit


        
        elif element_type == 'target':

            # 1. Dimensions 
            bottom_flute_height = 4.54
            graph_height = 9.665
            moly_height = 0.079375
            target_zone_height = 38.1
            upper_graph_height = 9.585625
            void_gap_height = 1.27
            upper_flute_height = 2.54

            clad_height = graph_height + moly_height + target_zone_height + upper_graph_height + void_gap_height


            # 2. Axial positions
            bottom_flute_z = z0
            bottom_graph_z = bottom_flute_z + bottom_flute_height
            moly_z = bottom_graph_z + graph_height
            target_meat_z = moly_z + moly_height
            upper_graph_z = target_meat_z + target_zone_height
            void_gap_z = upper_graph_z + upper_graph_height
            upper_flute_z = void_gap_z + void_gap_height

            
            # 3. Surfaces
            s_bottom_flute = RCC((x0, y0, bottom_flute_z), bottom_flute_height, 1.87325)
            s_bottom_refl  = RCC((x0, y0, bottom_graph_z), graph_height, 1.82245)
            s_moly_disk    = RCC((x0, y0, moly_z), moly_height, 1.82245)

            s_target_boundary = RCC((x0, y0, target_meat_z), target_zone_height, 1.82245)
            s_zirc_rod        = RCC((x0, y0, target_meat_z), target_zone_height, 0.2921)

            s_top_refl  = RCC((x0, y0, upper_graph_z), graph_height, 1.82245)
            s_void_gap  = RCC((x0, y0, void_gap_z), void_gap_height, 1.82245)
            s_upper_flute = RCC((x0, y0, upper_flute_z), upper_flute_height, 1.87325)

            s_cladding = RCC((x0, y0, bottom_graph_z), clad_height, 1.87325)

            
            # 4. Standard cells
            c_bottom_flute = openmc.Cell(fill=mats.flute_mix, region=-s_bottom_flute)
            c_bottom_refl  = openmc.Cell(fill=mats.graphite, region=-s_bottom_refl)
            
            # Molybedenum Replaced with Graphite
            c_moly_disk    = openmc.Cell(fill=mats.graphite, region=-s_moly_disk)


            # 5. HFIR pellet stack with surrounding void

            pellet_outer_radius = 0.165 * 2.54   # 0.4191 cm
            pellet_inner_radius = 0.115 * 2.54   # 0.2921 cm
            pellet_height = 0.5 * 2.54           # 1.27 cm

            z_current = target_meat_z

            pellet_cells = []
            void_cells = []

            for i in range(30):

                # axial planes
                p_bot = openmc.ZPlane(z0=z_current)
                p_top = openmc.ZPlane(z0=z_current + pellet_height)

                # pellet cylinders
                s_pellet_outer = openmc.ZCylinder(
                    x0=x0,
                    y0=y0,
                    r=pellet_outer_radius
                )

                s_pellet_inner = openmc.ZCylinder(
                    x0=x0,
                    y0=y0,
                    r=pellet_inner_radius
                )

                # pellet annulus
                pellet_region = (
                    -s_pellet_outer &
                    +s_pellet_inner &
                    +p_bot & -p_top
                )

                pellet_cell = openmc.Cell(
                    fill=mats.npo2,
                    region=pellet_region,
                    name=f"HFIR Pellet {i} in Pos {position}"
                )

                pellet_cells.append(pellet_cell)

                # surrounding void region
                void_region = (
                    -s_target_boundary &
                    +s_pellet_outer &
                    +p_bot & -p_top
                )

                void_cell = openmc.Cell(
                    fill=mats.air,
                    region=void_region,
                    name=f"Void Around Pellet {i} in Pos {position}"
                )

                void_cells.append(void_cell)

                z_current += pellet_height


            # 6. Zirconium rod replaced with Graphite (continuous)
            c_zirc_rod = openmc.Cell(
                fill=mats.graphite,
                region=-s_zirc_rod
    )


            # 7. Upper sections
            c_top_refl  = openmc.Cell(fill=mats.graphite, region=-s_top_refl)
            c_void_gap  = openmc.Cell(fill=mats.air, region=-s_void_gap)
            c_upper_flute = openmc.Cell(fill=mats.flute_mix, region=-s_upper_flute)


            # 8. Cladding (fully enclosing)
            c_cladding = openmc.Cell(
                fill=mats.steel,
                region=(
                    -s_cladding &
                    +s_bottom_refl &
                    +s_moly_disk &
                    +s_target_boundary &
                    +s_top_refl &
                    +s_void_gap
                )
            )


            # 9. Boundary water
            c_bound = openmc.Cell(fill=mats.water, region=(+s_cladding))


            # 10. Assemble
            target_cells = [
                c_bottom_flute,
                c_bottom_refl,
                c_moly_disk,
                *pellet_cells,
                *void_cells,
                c_zirc_rod,
                c_top_refl,
                c_void_gap,
                c_upper_flute,
                c_cladding,
                c_bound
            ]

            for cell in target_cells:
                cell.temperature = temp

            target_element = openmc.Universe(
                universe_id=position,
                name=f"Target Element in Position {position}",
                cells=target_cells
            )

            return target_element


        
        else:
            raise Exception(f"Element Type '{element_type}' not supported")


    def make_control_rod(position, bounding_cyl, core_temp, fuel_temp, percent_withdrawn, cr_type='ffcr'):
        """
        Generates an air-followed or fuel-followed control rod universe.

        :param position:            (float) grid position, i.e. 304 for C4 (for ID/naming/debug purpose only)
        :param bounding_cyl:        (openmc.model.RightCircularCylinder) cylinder bounding the cell
        :param core_temp:           (float) average core temperature in [K]
        :param fuel_temp:           (float) average fuel temperature in [K]
        :param percent_withdrawn:   (float) percent withdrawn (0 - 100)
        :param cr_type:             (str) 'afcr' or 'ffcr' (air-followed or fuel followed control rod)

        :return:                    (openmc.Universe) Universe containing the control rod

        """

        # X and Y coordinate
        x0 = bounding_cyl.cyl.x0
        y0 = bounding_cyl.cyl.y0

        # Z coordinate
        z0 = bounding_cyl.bottom.d / bounding_cyl.bottom.c

        # Total length of control rod travel
        length_of_travel = 38.1 # [cm]

        # Calculate the z-coordinates of each component, then translate based on withdrawal
        bottom_of_rod = z0 + (percent_withdrawn / 100) * length_of_travel
        end_cap = z0 + 1.99 + (percent_withdrawn / 100) * length_of_travel
        bottom_can = z0 + 16.185 + (percent_withdrawn / 100) * length_of_travel
        double_magneform = z0 + 19.445 + (percent_withdrawn / 100) * length_of_travel
        fuel = z0 + 57.545 + (percent_withdrawn / 100) * length_of_travel
        fuel_void = z0 + 58.185 + (percent_withdrawn / 100) * length_of_travel
        magneform = z0 + 59.545 + (percent_withdrawn / 100) * length_of_travel
        poison = z0 + 97.645 + (percent_withdrawn / 100) * length_of_travel
        poison_void = z0 + 100.645 + (percent_withdrawn / 100) * length_of_travel
        magneform2 = z0 + 102.275 + (percent_withdrawn / 100) * length_of_travel
        upper_can = z0 + 108.58 + (percent_withdrawn / 100) * length_of_travel
        top_of_rod = z0 + 111.12 + (percent_withdrawn / 100) * length_of_travel

        if cr_type == 'afcr':
            '''

            If the control rod is to be air-followed, the geometry is as follows:
            
            From the bottom up,

                (1)     The bottoming steel fitting (end cap)
                (2)     A 15" air-filled can (air-follower)
                (3)     A magneformed spacer
                (4)     A 15"-section of boron carbide.
                (5)     Another expansion gap
                (6)     The top steel fitting (end cap)

                The steel cladding begins at the top of the bottom fitting and terminates at the bottom of the top 
                fitting.

            The entire control rod is surrounded by water (the boundary)

            '''

            # Surfaces for each component
            s_bound = RCC(center_base=(x0, y0, bottom_of_rod), height=magneform2 - bottom_of_rod, radius=1.7145,
                          axis='z', name=f"AFCR Bound in Position {position}")
            s_end_cap = RCC(center_base=(x0, y0, bottom_of_rod), height=end_cap - bottom_of_rod, radius=1.7145, axis='z',
                            name=f"Bottom End Cap of AFCR in Position {position}")
            s_air_follower = RCC(center_base=(x0, y0, end_cap), height=fuel_void - end_cap, radius=1.6637, axis='z',
                                 name=f"Air in AFCR in Position {position}")
            s_magneform1 = RCC(center_base=(x0, y0, fuel_void), height=magneform - fuel_void, radius=1.6637, axis='z',
                               name=f"Magneform Between Follower and Poison in AFCR in Position {position}")
            s_poison = RCC(center_base=(x0, y0, magneform), height=poison - magneform, radius=1.6637, axis='z',
                           name=f"Poison in FFCR in Position {position}")
            s_poison_void = RCC(center_base=(x0, y0, poison), height=poison_void - poison, radius=1.6637, axis='z',
                                name=f"Expansion Gap Above Poison in AFCR in Position {position}")
            s_top_end_cap = RCC(center_base=(x0, y0, poison_void), height=magneform2 - poison_void, radius=1.7145,
                                axis='z', name=f"Top End Cap of AFCR in Position {position}")
            s_cladding = RCC(center_base=(x0, y0, end_cap), height=poison_void - end_cap, radius=1.7145, axis='z',
                             name=f"Cladding of FFCR in Position {position}")

            # Cells for each component
            c_end_cap = openmc.Cell(fill=mats.steel, region=(-s_end_cap),
                                    name=f"Bottom End Cap of AFCR in Pos. {position}")
            c_air_follower = openmc.Cell(fill=mats.air, region=(-s_air_follower),
                                         name=f"Air in AFCR in Pos. {position}")
            c_magneform1 = openmc.Cell(fill=mats.steel, region=(-s_magneform1),
                                       name=f"Magneform Between Fuel and Poison in AFCR in Pos. {position}")
            c_poison = openmc.Cell(fill=mats.b4c, region=(-s_poison),
                                   name=f"Poison in AFCR in Pos. {position}")
            c_poison_void = openmc.Cell(fill=mats.air, region=(-s_poison_void),
                                        name=f"Void Above Poison in AFCR in Pos. {position}")
            c_top_end_cap = openmc.Cell(fill=mats.steel, region=(-s_top_end_cap),
                                        name=f"Magneform Above Poison in AFCR in Pos. {position}")
            c_cladding = openmc.Cell(fill=mats.steel, region=(-s_cladding & +s_air_follower & +s_magneform1
                                                              & +s_poison & +s_poison_void),
                                     name=f"Cladding of FFCR in Pos. {position}")
            c_bound = openmc.Cell(fill=mats.water, region=(+s_bound))

            # List of cells in FFCR
            afcr_cells = [c_end_cap,
                          c_air_follower,
                          c_magneform1,
                          c_poison,
                          c_poison_void,
                          c_top_end_cap,
                          c_cladding,
                          c_bound]

            # Assign temperatures to each cell in the rod
            for component in afcr_cells:
                component.temperature = core_temp

            # Generate the AFCR universe
            afcr = openmc.Universe(universe_id=position, name=f"AFCR in Position {position}", cells=afcr_cells)

            return afcr

        elif cr_type == 'ffcr':
            '''
            
            If the control rod is to be fuel-followed, the geometry is a bit different:
            
            From the bottom up, the geometry is:
                
                (1)     The bottoming steel fitting (end cap)
                (2)     The bottom "voided" displacement can
                (3)     A double-magneformed spacer
                (4a)    A 15"-section of TRIGA fuel. The fuel is slightly narrower in the control rods than a standard 
                        fuel-moderator element.
                (4b)    A 15" pure zirconium rod through the center of the TRIGA fuel.
                (5)     An "air-filled" expansion gap
                (6)     Another magneformed spacer
                (7)     A 15"-section of boron carbide.
                (8)     Another expansion gap
                (9)     Another magneformed spacer
                (10)    The top "voided" displacement can
                (11)    The top steel fitting (end cap)
                
                The steel cladding begins at the top of the bottom fitting and terminates at the bottom of the top 
                fitting.
            
            The entire control rod is surrounded by water (the boundary)
            
            '''
            # Surfaces for each component

            s_bound = RCC(center_base=(x0, y0,bottom_of_rod), height=top_of_rod - bottom_of_rod, radius=1.7145, axis='z',
                          name=f"FFCR Bound in Position {position}")
            s_end_cap = RCC(center_base=(x0, y0,bottom_of_rod), height=end_cap-bottom_of_rod, radius=1.7145, axis='z',
                            name=f"Bottom End Cap of FFCR in Position {position}")
            s_bottom_can = RCC(center_base=(x0, y0,end_cap), height=bottom_can-end_cap, radius=1.6637, axis='z',
                               name=f"Bottom Void Can of FFCR in Position {position}")
            s_magneform0 = RCC(center_base=(x0, y0,bottom_can), height=double_magneform-bottom_can, radius=1.6637,
                               axis='z', name=f"Magneform Below Fuel in FFCR in Position {position}")
            s_fuel_meat = RCC(center_base=(x0, y0,double_magneform), height=fuel-double_magneform, radius=1.6637,
                              axis='z', name=f"Fuel in FFCR in Position {position}")
            s_zirc_rod = RCC(center_base=(x0, y0,double_magneform), height=fuel-double_magneform, radius=0.3175,
                             axis='z', name=f"Zirconium Rod in Fuel in FFCR in Position {position}")
            s_fuel_void = RCC(center_base=(x0, y0,fuel), height=fuel_void-fuel, radius=1.6637, axis='z',
                              name=f"Expansion Gap Above Fuel in FFCR in Position {position}")
            s_magneform1 = RCC(center_base=(x0, y0,fuel_void), height=magneform-fuel_void, radius=1.6637, axis='z',
                               name=f"Magneform Between Fuel and Poison in FFCR in Position {position}")
            s_poison = RCC(center_base=(x0, y0,magneform), height=poison-magneform, radius=1.6637, axis='z',
                           name=f"Poison in FFCR in Position {position}")
            s_poison_void = RCC(center_base=(x0, y0,poison), height=poison_void-poison, radius=1.6637, axis='z',
                                name=f"Expansion Gap Above Poison in FFCR in Position {position}")
            s_magneform2 = RCC(center_base=(x0, y0,poison_void), height=magneform2-poison_void, radius=1.6637, axis='z',
                               name=f"Magneform Above Poison in FFCR in Position {position}")
            s_upper_can = RCC(center_base=(x0, y0,magneform2), height=upper_can-magneform2, radius=1.6637, axis='z',
                              name=f"Upper Void Can of FFCR in Position {position}")
            s_top_end_cap = RCC(center_base=(x0, y0,upper_can), height=top_of_rod-upper_can, radius=1.7145, axis='z',
                                name=f"Top End Cap of FFCR in Position {position}")
            s_cladding = RCC(center_base=(x0, y0,end_cap), height=upper_can-end_cap, radius=1.7145, axis='z',
                             name=f"Cladding of FFCR in Position {position}")

            # Cells for each component
            c_end_cap = openmc.Cell(fill=mats.steel, region=(-s_end_cap),
                                    name=f"Bottom End Cap of FFCR in Pos. {position}")
            c_bottom_can = openmc.Cell(fill=mats.air, region=(-s_bottom_can),
                                       name=f"Bottom Void Can of FFCR in Pos. {position}")
            c_magneform0 = openmc.Cell(fill=mats.steel, region=(-s_magneform0),
                                       name=f"Magneform Below Fuel in FFCR in Pos. {position}")
            c_fuel_meat = openmc.Cell(fill=mats.fuel, region=(-s_fuel_meat & +s_zirc_rod),
                                      name=f"Fuel in FFCR in Pos. {position}")
            c_zirc_rod = openmc.Cell(fill=mats.zirconium, region=(-s_zirc_rod),
                                     name=f"Zirconium Rod in FFCR in Pos. {position}")
            c_fuel_void = openmc.Cell(fill=mats.air, region=(-s_fuel_void),
                                      name=f"Expansion Gap Above Fuel of FFCR in Pos. {position}")
            c_magneform1 = openmc.Cell(fill=mats.steel, region=(-s_magneform1),
                                       name=f"Magneform Between Fuel and Poison in FFCR in Pos. {position}")
            c_poison = openmc.Cell(fill=mats.b4c, region=(-s_poison),
                                   name=f"Poison in FFCR in Pos. {position}")
            c_poison_void = openmc.Cell(fill=mats.air, region=(-s_poison_void),
                                        name=f"Void Above Poison in FFCR in Pos. {position}")
            c_magneform2 = openmc.Cell(fill=mats.steel, region=(-s_magneform2), cell_id=int(f"{position}" + "09"),
                                       name=f"Magneform Above Poison in FFCR in Pos. {position}")
            c_upper_can = openmc.Cell(fill=mats.air, region=(-s_upper_can),
                                      name=f"Upper Void Can of FFCR in Pos. {position}")
            c_top_end_cap = openmc.Cell(fill=mats.steel, region=(-s_top_end_cap),
                                        name=f"Top End Cap of FFCR in Pos. {position}")
            c_cladding = openmc.Cell(fill=mats.steel,
                                     region=(-s_cladding & +s_bottom_can & +s_magneform0
                                             & +s_fuel_meat & +s_fuel_void & +s_magneform1
                                             & +s_poison & +s_poison_void & +s_magneform2
                                             & +s_upper_can),
                                     name=f"Cladding of FFCR in Pos. {position}")
            c_bound = openmc.Cell(fill=mats.water, region=(+s_bound))

            # List of cells in FFCR
            ffcr_cells = [c_end_cap,
                          c_bottom_can,
                          c_magneform0,
                          c_fuel_meat,
                          c_zirc_rod,
                          c_fuel_void,
                          c_magneform1,
                          c_poison,
                          c_poison_void,
                          c_magneform2,
                          c_upper_can,
                          c_top_end_cap,
                          c_cladding,
                          c_bound]

            # Assign temperatures to each cell in the rod
            for component in ffcr_cells:
                component.temperature = core_temp

            # The fuel temperature is much higher than the rest of the rod
            ffcr_cells[3].temperature = fuel_temp
            ffcr_cells[4].temperature = fuel_temp

            # Generate the FFCR universe
            ffcr = openmc.Universe(universe_id=position, name=f"FFCR in Position {position}", cells=ffcr_cells)

            return ffcr

        else:
            raise Exception(f"Control Rod Type '{cr_type}' Not Supported")


    def gen_beam_port_can(center_base, axis, radius):
        """
        Generates a cylinder for the void cans in the reflector that are axially aligned
        with the beam ports.

        OpenMC currently cannot define cylinders other than those parallel to the 'x', 'y', or 'z' axes.
        As such, any cylinder in a different orientation must be adjusted after the standard cylinder is generated.

        :param center_base: (tuple) center of the bottom face of the cylinder (also the pivot point)
        :param axis:        (tuple) a vector describing the direction and length of the cylinder
        :param radius:      (float) radius of the cylinder

        :return:            (class) openmc.Surface() instance of the rotated RCC
        """

        vx = center_base[0]
        vy = center_base[1]
        vz = center_base[2]

        hx = axis[0]
        hy = axis[1]
        hz = axis[2]

        def rotation_matrix(v1, v2):
            """
            This function is a method pulled from the OpenMC-MCNP adapter to rotate cylinders.

            Compute rotation matrix that would rotate v1 into v2.

            Parameters
            ----------
            v1 : numpy.ndarray
                Unrotated vector
            v2 : numpy.ndarray
                Rotated vector

            Returns
            -------
            3x3 rotation matrix

            """

            # Normalize vectors
            u1 = v1 / np.linalg.norm(v1)
            u2 = v2 / np.linalg.norm(v2)

            # Calculate axis of rotation
            axis = np.cross(u1, u2)
            axis /= np.linalg.norm(axis)

            I = np.identity(3)

            # Handle special case where vectors are parallel or anti-parallel
            if abs(np.dot(u2, axis) - 1.0) < 1e-8:
                return I
            elif abs(np.dot(u2, axis) + 1.0) < 1e-8:
                return -I
            else:
                # Calculate rotation angle
                cos_angle = np.dot(u1, u2)
                sin_angle = np.sqrt(1 - cos_angle * cos_angle)

                # Create cross-product matrix K
                kx, ky, kz = axis
                K = np.array([
                    [0.0, -kz, ky],
                    [kz, 0.0, -kx],
                    [-ky, kx, 0.0]
                ])

                # Create rotation matrix using Rodrigues' rotation formula
                return I + K * sin_angle + (K @ K) * (1 - cos_angle)

        # Create vectors for Z-axis and cylinder orientation
        u = np.array([0., 0., 1.])
        h = np.array([hx, hy, hz])

        # Determine rotation matrix to transform u -> h
        rotation = rotation_matrix(u, h)

        # Create RCC aligned with Z-axis
        height = float(np.linalg.norm(h))
        surf = RCC((vx, vy, vz), height, radius, axis='z')

        # Rotate the RCC
        surf = surf.rotate(rotation, pivot=(vx, vy, vz))
        return surf


    ####################################################################################################################
    #                                                                                                                  #
    #  -------------------------------------------------- SURFACES --------------------------------------------------  #
    #                                                                                                                  #
    ####################################################################################################################

    print("Generating Surfaces...")

    origin = (0, 0, 0)

    # Problem Boundary
    s_bound = openmc.Sphere(x0=0, y0=0, z0=0, r=200, boundary_type='vacuum')

    #################################################
    #  ========== MAJOR CORE COMPONENTS ==========  #
    #################################################

    # Grid Plates and Core Bound
    s_lower_grid_plate = RCC(center_base=origin, height=1.27, radius=24.844375, axis='z')
    s_upper_grid_plate = RCC(center_base=(0, 0, 65.78), height=1.27, radius=26.749375, axis='z')
    s_core = RCC(center_base=origin, height=67.05, radius=26.749375, axis='z')

    # Reflector
    s_refl_od = RCC(center_base=(0, 0, 5.81), height=59.97, radius=56.75, axis='z')
    s_refl_lead_id = RCC(center_base=(0, 0, 5.81), height=59.97, radius=51.67, axis='z')

    # Lazy Susan
    s_lazy_susan = RCC(center_base=(0, 0, 39.11), height=28.67, radius=34.369375, axis='z')
    s_ls_air_od = RCC(center_base=(0, 0, 39.61), height=27.67, radius=33.869375, axis='z')
    s_ls_air_id = openmc.ZCylinder(x0=0, y0=0, r=27.249375)

    # Beam Port Cans In Reflector
    s_bp1 = gen_beam_port_can(center_base=(22.2154, 11.4662, 20.6375), axis=(35.5447, 18.346, 0), radius=7.62)
    s_bp2 = gen_beam_port_can(center_base=(22.2154,  -11.4662, 20.6375), axis=(35.5447,  -18.346, 0), radius=7.62)
    s_bp3 = gen_beam_port_can(center_base=(-33.1391, 16.3782, 20.6375), axis=(-30.1067,  -30.2663, 0), radius=7.62)
    s_bp4 = gen_beam_port_can(center_base=(-22.2154, 11.4662, 20.6375), axis=(-35.5447, 18.346, 0), radius=8.89)


    ##########################################
    # =========== GRID LOCATIONS =========== #
    ##########################################

    # These surfaces were generated by an Excel spreadsheet using the "CONCAT" function

    # A-Ring
    s_A1 = RCC(center_base=(0, 0, 1.27), height=85.78, radius=1.87325)                      # A1 CENTRAL THIMBLE

    # B-Ring
    s_B1 = RCC(center_base=(1.0414, -3.9116, 1.27), height=85.78, radius=1.87325, axis='z')           # B1 CLICIT
    s_B2 = RCC(center_base=(-2.8702, -2.8702, 1.27), height=65.78, radius=1.87325, axis='z')
    s_B3 = RCC(center_base=(-3.9116, 1.0414, 1.27), height=65.78, radius=1.87325, axis='z')
    s_B4 = RCC(center_base=(-1.0414, 3.9116, 1.27), height=65.78, radius=1.87325, axis='z')
    s_B5 = RCC(center_base=(2.8702, 2.8702, 1.27), height=65.78, radius=1.87325, axis='z')
    s_B6 = RCC(center_base=(3.9116, -1.0414, 1.27), height=65.78, radius=1.87325, axis='z')

    # C-Ring
    s_C1 = RCC(center_base=(0, -7.9756, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C2 = RCC(center_base=(-3.9878, -6.9088, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C3 = RCC(center_base=(-6.9088, -3.9878, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C4 = RCC(center_base=(-7.9756, 0, -44.07), height=150, radius=1.87325, axis='z')      # TRANSIENT ROD
    s_C5 = RCC(center_base=(-6.9088, 3.9878, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C6 = RCC(center_base=(-3.9878, 6.9088, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C7 = RCC(center_base=(0, 7.9756, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C8 = RCC(center_base=(3.9878, 6.9088, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C9 = RCC(center_base=(6.9088, 3.9878, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C10 = RCC(center_base=(0, -11.938, -44.07), height=150, radius=1.87325, axis='z')     # SHIM ROD
    s_C11 = RCC(center_base=(6.9088, -3.9878, 1.27), height=65.78, radius=1.87325, axis='z')
    s_C12 = RCC(center_base=(3.9878, -6.9088, 1.27), height=65.78, radius=1.87325, axis='z')

    # D-Ring
    s_D1 = RCC(center_base=(0, 11.938, -44.07), height=150, radius=1.87325, axis='z')       # SAFETY ROD
    s_D2 = RCC(center_base=(-4.0894, -11.2268, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D3 = RCC(center_base=(-7.6708, -9.1694, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D4 = RCC(center_base=(-10.3378, -5.969, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D5 = RCC(center_base=(-11.7602, -2.2352, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D6 = RCC(center_base=(-11.303, 2.0574, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D7 = RCC(center_base=(-10.3378, 5.969, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D8 = RCC(center_base=(-7.6708, 9.144, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D9 = RCC(center_base=(-4.0894, 11.2268, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D10 = RCC(center_base=(7.9756, 0, -44.07), height=150, radius=1.87325, axis='z')      # REGULATING ROD
    s_D11 = RCC(center_base=(4.0894, 11.2268, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D12 = RCC(center_base=(7.6708, 9.144, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D13 = RCC(center_base=(10.3378, 5.969, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D14 = RCC(center_base=(11.7602, 2.2098, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D15 = RCC(center_base=(11.303, -2.0828, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D16 = RCC(center_base=(10.3378, -5.969, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D17 = RCC(center_base=(7.6708, -9.1694, 1.27), height=65.78, radius=1.87325, axis='z')
    s_D18 = RCC(center_base=(4.0894, -11.2268, 1.27), height=65.78, radius=1.87325, axis='z')

    # E-Ring
    s_E1 = RCC(center_base=(0, -15.9258, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E2 = RCC(center_base=(-4.1148, -15.367, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E3 = RCC(center_base=(-7.9502, -13.7922, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E4 = RCC(center_base=(-11.2522, -11.2522, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E5 = RCC(center_base=(-13.7668, -7.9502, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E6 = RCC(center_base=(-15.1892, -4.2164, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E7 = RCC(center_base=(-15.1892, -0.254, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E8 = RCC(center_base=(-15.367, 4.1148, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E9 = RCC(center_base=(-13.7668, 7.9502, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E10 = RCC(center_base=(-11.2522, 11.2522, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E11 = RCC(center_base=(-7.9502, 13.7668, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E12 = RCC(center_base=(-4.1148, 15.367, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E13 = RCC(center_base=(0, 15.9004, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E14 = RCC(center_base=(4.1148, 15.367, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E15 = RCC(center_base=(7.9502, 13.7668, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E16 = RCC(center_base=(11.2522, 11.2522, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E17 = RCC(center_base=(13.7922, 7.9502, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E18 = RCC(center_base=(15.2146, 4.191, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E19 = RCC(center_base=(15.2146, 0.2286, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E20 = RCC(center_base=(15.367, -4.1148, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E21 = RCC(center_base=(13.7922, -7.9502, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E22 = RCC(center_base=(11.2522, -11.2522, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E23 = RCC(center_base=(7.9502, -13.7922, 1.27), height=65.78, radius=1.87325, axis='z')
    s_E24 = RCC(center_base=(4.1148, -15.367, 1.27), height=65.78, radius=1.87325, axis='z')

    # F-Ring
    s_F1 = RCC(center_base=(0, -19.8882, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F2 = RCC(center_base=(-4.1402, -19.4564, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F3 = RCC(center_base=(-8.0772, -18.161, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F4 = RCC(center_base=(-11.684, -16.1036, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F5 = RCC(center_base=(-14.7828, -13.3096, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F6 = RCC(center_base=(-17.1196, -10.1346, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F7 = RCC(center_base=(-18.923, -6.1468, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F8 = RCC(center_base=(-19.7866, -2.0828, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F9 = RCC(center_base=(-19.7866, 2.0828, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F10 = RCC(center_base=(-18.923, 6.1468, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F11 = RCC(center_base=(-17.2212, 9.9314, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F12 = RCC(center_base=(-14.7828, 13.3096, 1.27), height=85.78, radius=1.87325, axis='z')        # F12 ICIT
    s_F13 = RCC(center_base=(-11.684, 16.0782, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F14 = RCC(center_base=(-8.0772, 18.161, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F15 = RCC(center_base=(-4.1402, 19.4564, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F16 = RCC(center_base=(0, 19.8882, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F17 = RCC(center_base=(4.1402, 19.4564, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F18 = RCC(center_base=(8.1026, 18.161, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F19 = RCC(center_base=(11.684, 16.0782, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F20 = RCC(center_base=(14.7828, 13.3096, 1.27), height=85.78, radius=1.87325, axis='z')         # F20 CLOCIT
    s_F21 = RCC(center_base=(17.145, 10.1092, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F22 = RCC(center_base=(18.923, 6.1468, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F23 = RCC(center_base=(19.7866, 2.0828, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F24 = RCC(center_base=(19.7866, -2.0828, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F25 = RCC(center_base=(18.923, -6.1468, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F26 = RCC(center_base=(17.2212, -9.9568, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F27 = RCC(center_base=(14.7828, -13.3096, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F28 = RCC(center_base=(11.684, -16.1036, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F29 = RCC(center_base=(8.1026, -18.161, 1.27), height=65.78, radius=1.87325, axis='z')
    s_F30 = RCC(center_base=(4.1402, -19.4564, 1.27), height=65.78, radius=1.87325, axis='z')

    # G-Ring
    s_G1 = RCC(center_base=(0, -23.876, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G2 = RCC(center_base=(-4.1402, -23.495, 1.27), height=85.78, radius=1.87325, axis='z')          # G2 RABBIT
    s_G3 = RCC(center_base=(-8.1534, -22.4282, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G4 = RCC(center_base=(-11.938, -20.6756, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G5 = RCC(center_base=(-15.3416, -18.288, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G6 = RCC(center_base=(-18.2626, -15.3416, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G7 = RCC(center_base=(-20.6502, -11.938, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G8 = RCC(center_base=(-22.606, -8.0772, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G9 = RCC(center_base=(-23.495, -4.1402, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G10 = RCC(center_base=(-23.8506, 0, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G11 = RCC(center_base=(-23.495, 4.1402, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G12 = RCC(center_base=(-22.4282, 8.1534, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G13 = RCC(center_base=(-20.6502, 11.938, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G14 = RCC(center_base=(-18.2626, 15.3416, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G15 = RCC(center_base=(-15.3416, 18.2626, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G16 = RCC(center_base=(-11.938, 20.6502, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G17 = RCC(center_base=(-8.1534, 22.4282, 1.27), height=65.78, radius=1.87325, axis='z')         # AmBe SOURCE
    s_G18 = RCC(center_base=(-4.1402, 23.495, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G19 = RCC(center_base=(0, 23.8506, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G20 = RCC(center_base=(4.1402, 23.495, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G21 = RCC(center_base=(8.1534, 22.4282, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G22 = RCC(center_base=(11.938, 20.6502, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G23 = RCC(center_base=(15.3416, 18.2626, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G24 = RCC(center_base=(18.288, 15.3416, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G25 = RCC(center_base=(20.6756, 11.938, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G26 = RCC(center_base=(22.606, 8.0772, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G27 = RCC(center_base=(23.495, 4.1402, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G28 = RCC(center_base=(23.876, 0, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G29 = RCC(center_base=(23.495, -4.1402, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G30 = RCC(center_base=(22.4282, -8.1534, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G31 = RCC(center_base=(20.6756, -11.938, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G32 = RCC(center_base=(18.288, -15.3416, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G33 = RCC(center_base=(15.3416, -18.288, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G34 = RCC(center_base=(11.938, -20.6756, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G35 = RCC(center_base=(8.1534, -22.4282, 1.27), height=65.78, radius=1.87325, axis='z')
    s_G36 = RCC(center_base=(4.1402, -23.495, 1.27), height=65.78, radius=1.87325, axis='z')

    # The upper grid plate region is within the bounding cylinder and outside all grid positions.
    upper_grid_plate_region = (-s_upper_grid_plate & +s_A1 &
                               +s_B1 & +s_B2 & +s_B3 & +s_B4 & +s_B5 & +s_B6 &
                               +s_C1 & +s_C2 & +s_C3 & +s_C4 & +s_C5 & +s_C6 &
                               +s_C7 & +s_C8 & +s_C9 & +s_C10 & +s_C11 & +s_C12 &
                               +s_D1 & +s_D2 & +s_D3 & +s_D4 & +s_D5 & +s_D6 &
                               +s_D7 & +s_D8 & +s_D9 & +s_D10 & +s_D11 & +s_D12 &
                               +s_D13 & +s_D14 & +s_D15 & +s_D16 & +s_D17 & +s_D18 &
                               +s_E1 & +s_E2 & +s_E3 & +s_E4 & +s_E5 & +s_E6 &
                               +s_E7 & +s_E8 & +s_E9 & +s_E10 & +s_E11 & +s_E12 &
                               +s_E13 & +s_E14 & +s_E15 & +s_E16 & +s_E17 & +s_E18 &
                               +s_E19 & +s_E20 & +s_E21 & +s_E22 & +s_E23 & +s_E24 &
                               +s_F1 & +s_F2 & +s_F3 & +s_F4 & +s_F5 & +s_F6 &
                               +s_F7 & +s_F8 & +s_F9 & +s_F10 & +s_F11 & +s_F12 &
                               +s_F13 & +s_F14 & +s_F15 & +s_F16 & +s_F17 & +s_F18 &
                               +s_F19 & +s_F20 & +s_F21 & +s_F22 & +s_F23 & +s_F24 &
                               +s_F25 & +s_F26 & +s_F27 & +s_F28 & +s_F29 & +s_F30 &
                               +s_G1 & +s_G2 & +s_G3 & +s_G4 & +s_G5 & +s_G6 &
                               +s_G7 & +s_G8 & +s_G9 & +s_G10 & +s_G11 & +s_G12 &
                               +s_G13 & +s_G14 & +s_G15 & +s_G16 & +s_G17 & +s_G18 &
                               +s_G19 & +s_G20 & +s_G21 & +s_G22 & +s_G23 & +s_G24 &
                               +s_G25 & +s_G26 & +s_G27 & +s_G28 & +s_G29 & +s_G30 &
                               +s_G31 & +s_G32 & +s_G33 & +s_G34 & +s_G35 & +s_G36)

    # The water region is within the bounding cylinder and outside all grid positions and the grid plates.
    core_water_region =   (-s_core & +s_lower_grid_plate & +s_upper_grid_plate &
                           +s_A1 &
                           +s_B1 & +s_B2 & +s_B3 & +s_B4 & +s_B5 & +s_B6 &
                           +s_C1 & +s_C2 & +s_C3 & +s_C4 & +s_C5 & +s_C6 &
                           +s_C7 & +s_C8 & +s_C9 & +s_C10 & +s_C11 & +s_C12 &
                           +s_D1 & +s_D2 & +s_D3 & +s_D4 & +s_D5 & +s_D6 &
                           +s_D7 & +s_D8 & +s_D9 & +s_D10 & +s_D11 & +s_D12 &
                           +s_D13 & +s_D14 & +s_D15 & +s_D16 & +s_D17 & +s_D18 &
                           +s_E1 & +s_E2 & +s_E3 & +s_E4 & +s_E5 & +s_E6 &
                           +s_E7 & +s_E8 & +s_E9 & +s_E10 & +s_E11 & +s_E12 &
                           +s_E13 & +s_E14 & +s_E15 & +s_E16 & +s_E17 & +s_E18 &
                           +s_E19 & +s_E20 & +s_E21 & +s_E22 & +s_E23 & +s_E24 &
                           +s_F1 & +s_F2 & +s_F3 & +s_F4 & +s_F5 & +s_F6 &
                           +s_F7 & +s_F8 & +s_F9 & +s_F10 & +s_F11 & +s_F12 &
                           +s_F13 & +s_F14 & +s_F15 & +s_F16 & +s_F17 & +s_F18 &
                           +s_F19 & +s_F20 & +s_F21 & +s_F22 & +s_F23 & +s_F24 &
                           +s_F25 & +s_F26 & +s_F27 & +s_F28 & +s_F29 & +s_F30 &
                           +s_G1 & +s_G2 & +s_G3 & +s_G4 & +s_G5 & +s_G6 &
                           +s_G7 & +s_G8 & +s_G9 & +s_G10 & +s_G11 & +s_G12 &
                           +s_G13 & +s_G14 & +s_G15 & +s_G16 & +s_G17 & +s_G18 &
                           +s_G19 & +s_G20 & +s_G21 & +s_G22 & +s_G23 & +s_G24 &
                           +s_G25 & +s_G26 & +s_G27 & +s_G28 & +s_G29 & +s_G30 &
                           +s_G31 & +s_G32 & +s_G33 & +s_G34 & +s_G35 & +s_G36)

    # The water surrounding is outside of the reflector, core, lazy susan, control rods, and experiment facilities.
    water_outside_region = (-s_bound & +s_refl_od & +s_core &
                            +s_C4 & +s_C10 & +s_D1 & +s_D10 &
                            +s_A1 & +s_B1 & +s_F12 & +s_F20 &
                            +s_G2 & (+s_lazy_susan | -s_core.cyl))

    print("Surface Generation Complete. \n")

    ####################################################################################################################
    #                                                                                                                  #
    #  ---------------------------------------------------- CELLS ---------------------------------------------------  #
    #                                                                                                                  #
    ####################################################################################################################

    print("Generating Cells...")

    ##################################################
    # ========== MAIN ASSEMBLY COMPONENTS ========== #
    ##################################################

    # Grid Plates
    c_lower_grid_plate = openmc.Cell(fill=mats.aluminum, region=(-s_lower_grid_plate & +s_C4 & +s_C10 & +s_D1 & +s_D10 ),
                                     name="Lower Grid Plate")
    c_upper_grid_plate = openmc.Cell(fill=mats.aluminum, region=upper_grid_plate_region, name="Upper Grid Plate")

    # Water In Core
    c_core_water = openmc.Cell(fill=mats.water, region=core_water_region, name="Water In Core")

    # Reflector
    c_refl = openmc.Cell(fill=mats.graphite,
                         region=(-s_refl_lead_id & +s_core.cyl & +s_bp1 & +s_bp2 & +s_bp3 & +s_bp4 & +s_lazy_susan),
                         name="Reflector Graphite")
    c_refl_lead = openmc.Cell(fill=mats.lead,
                              region=(-s_refl_od & +s_refl_lead_id.cyl & +s_bp1 & +s_bp2 & +s_bp3 & +s_bp4),
                              name="Reflector Lead")

    # Lazy Susan
    c_lazy_susan_air = openmc.Cell(fill=mats.air, region=(-s_ls_air_od & +s_ls_air_id), name="Air in Lazy Susan")
    c_lazy_susan = openmc.Cell(fill=mats.aluminum, region=(-s_lazy_susan & ~c_lazy_susan_air.region & +s_core.cyl))


    # Beam Port Cans In Reflector
    c_bp1 = openmc.Cell(fill=mats.air, region=(-s_bp1 & +s_core & -s_refl_od), name="Beam Port #1 Void Can")
    c_bp2 = openmc.Cell(fill=mats.air, region=(-s_bp2 & +s_core & -s_refl_od), name="Beam Port #2 Void Can")
    c_bp4 = openmc.Cell(fill=mats.air, region=(-s_bp4 & +s_core & -s_refl_od), name="Beam Port #4 Void Can")
    c_bp3 = openmc.Cell(fill=mats.air, region=(-s_bp3 & +s_bp4 & +s_core & -s_refl_od),
                        name="Beam Port #3 Void Can")


    assembly_cells = [c_lower_grid_plate, c_upper_grid_plate, c_core_water,
                      c_refl, c_refl_lead, c_lazy_susan, c_lazy_susan_air,
                      c_bp1, c_bp2, c_bp4, c_bp3, c_bp4]

    for cell in assembly_cells:
        cell.temperature = core_water_temperature

    ##############################################
    # =========== CORE CONFIGURATION =========== #
    ##############################################

    '''
    This is the core configuration in use since 2017, called the two-CLICIT core. This configuration
    moved the G-Ring ICIT (GRICIT) from G14 to the F-Ring in F12. Fuel was moved to the G-Ring to level
    the flux shape around the tube, and provide more flux to BP4 and BP3. This configuration also added a CLICIT
    tube in F20, called the CLOCIT (Cadmium-Lined OUTER Core Irradiation Tube). From this point, only a CLICIT tube is
    used in B1 (before it could be swapped for an ICIT or a fuel element).
    '''

    # Control Rods

    # TRANSIENT
    c_C4 = openmc.Cell(fill=make_control_rod(304, s_C4, core_water_temperature, fuel_temperature,
                                             trans_pos, 'afcr'),
                       region=-s_C4, name='Position C4')
    # SAFETY
    c_D1 = openmc.Cell(fill=make_control_rod(401, s_D1, core_water_temperature, fuel_temperature,
                                             safety_pos, 'ffcr'),
                       region=-s_D1, name='Position D1')
    # SHIM
    c_D10 = openmc.Cell(fill=make_control_rod(410, s_D10, core_water_temperature, fuel_temperature,
                                              shim_pos, 'ffcr'),
                        region=-s_D10, name='Position D10')
    # REGULATING
    c_C10 = openmc.Cell(fill=make_control_rod(310, s_C10, core_water_temperature, fuel_temperature,
                                              reg_pos, 'ffcr'),
                        region=-s_C10, name='Position C10')

    # Everything Else

    # A-Ring
    c_A1 = openmc.Cell(fill=make_element(101, s_A1, core_water_temperature, 'ct'),
                       region=(-s_A1), name='Position A1')

    # B-Ring
    c_B1 = openmc.Cell(fill=make_element(201, s_B1, core_water_temperature, 'clicit'),
                       region=(-s_B1), name='Position B1')
    c_B2 = openmc.Cell(fill=make_element(202, s_B2, fuel_temperature, 'fuel'),
                       region=(-s_B2), name='Position B2')
    c_B3 = openmc.Cell(fill=make_element(203, s_B3, fuel_temperature, 'fuel'),
                       region=(-s_B3), name='Position B3')
    c_B4 = openmc.Cell(fill=make_element(204, s_B4, fuel_temperature, 'fuel'),
                       region=(-s_B4), name='Position B4')
    c_B5 = openmc.Cell(fill=make_element(205, s_B5, fuel_temperature, 'fuel'),
                       region=(-s_B5), name='Position B5')
    c_B6 = openmc.Cell(fill=make_element(206, s_B6, fuel_temperature, 'fuel'),
                       region=(-s_B6), name='Position B6')

    # C-Ring
    c_C1 = openmc.Cell(fill=make_element(301, s_C1, fuel_temperature, 'fuel'),
                       region=(-s_C1), name='Position C1')
    c_C2 = openmc.Cell(fill=make_element(302, s_C2, fuel_temperature, 'fuel'),
                       region=(-s_C2), name='Position C2')
    c_C3 = openmc.Cell(fill=make_element(303, s_C3, fuel_temperature, 'fuel'),
                       region=(-s_C3), name='Position C3')
    c_C5 = openmc.Cell(fill=make_element(305, s_C5, fuel_temperature, 'fuel'),
                       region=(-s_C5), name='Position C5')
    c_C6 = openmc.Cell(fill=make_element(306, s_C6, fuel_temperature, 'fuel'),
                       region=(-s_C6), name='Position C6')
    c_C7 = openmc.Cell(fill=make_element(307, s_C7, fuel_temperature, 'fuel'),
                       region=(-s_C7), name='Position C7')
    c_C8 = openmc.Cell(fill=make_element(308, s_C8, fuel_temperature, 'fuel'),
                       region=(-s_C8), name='Position C8')
    c_C9 = openmc.Cell(fill=make_element(309, s_C9, fuel_temperature, 'fuel'),
                       region=(-s_C9), name='Position C9')
    c_C11 = openmc.Cell(fill=make_element(311, s_C11, fuel_temperature, 'fuel'),
                        region=(-s_C11), name='Position C11')
    c_C12 = openmc.Cell(fill=make_element(312, s_C12, fuel_temperature, 'fuel'),
                        region=(-s_C12), name='Position C12')

    # D-Ring
    c_D2 = openmc.Cell(fill=make_element(402, s_D2, fuel_temperature, 'fuel'),
                       region=(-s_D2), name='Position D2')
    c_D3 = openmc.Cell(fill=make_element(403, s_D3, fuel_temperature, 'fuel'),
                       region=(-s_D3), name='Position D3')
    c_D4 = openmc.Cell(fill=make_element(404, s_D4, fuel_temperature, 'fuel'),
                       region=(-s_D4), name='Position D4')
    c_D5 = openmc.Cell(fill=make_element(405, s_D5, fuel_temperature, 'fuel'),
                       region=(-s_D5), name='Position D5')
    c_D6 = openmc.Cell(fill=make_element(406, s_D6, fuel_temperature, 'fuel'),
                       region=(-s_D6), name='Position D6')
    c_D7 = openmc.Cell(fill=make_element(407, s_D7, fuel_temperature, 'fuel'),
                       region=(-s_D7), name='Position D7')
    c_D8 = openmc.Cell(fill=make_element(408, s_D8, fuel_temperature, 'fuel'),
                       region=(-s_D8), name='Position D8')
    c_D9 = openmc.Cell(fill=make_element(409, s_D9, fuel_temperature, 'fuel'),
                       region=(-s_D9), name='Position D9')
    c_D11 = openmc.Cell(fill=make_element(411, s_D11, fuel_temperature, 'fuel'),
                        region=(-s_D11), name='Position D11')
    c_D12 = openmc.Cell(fill=make_element(412, s_D12, fuel_temperature, 'fuel'),
                        region=(-s_D12), name='Position D12')
    c_D13 = openmc.Cell(fill=make_element(413, s_D13, fuel_temperature, 'fuel'),
                        region=(-s_D13), name='Position D13')
    c_D14 = openmc.Cell(fill=make_element(414, s_D14, fuel_temperature, 'fuel'),
                        region=(-s_D14), name='Position D14')
    c_D15 = openmc.Cell(fill=make_element(415, s_D15, fuel_temperature, 'fuel'),
                        region=(-s_D15), name='Position D15')
    c_D16 = openmc.Cell(fill=make_element(416, s_D16, fuel_temperature, 'fuel'),
                        region=(-s_D16), name='Position D16')
    c_D17 = openmc.Cell(fill=make_element(417, s_D17, fuel_temperature, 'fuel'),
                        region=(-s_D17), name='Position D17')
    c_D18 = openmc.Cell(fill=make_element(418, s_D18, fuel_temperature, 'fuel'),
                        region=(-s_D18), name='Position D18')

    # E-Ring
    c_E1 = openmc.Cell(fill=make_element(501, s_E1, fuel_temperature, 'fuel'),
                       region=(-s_E1), name='Position E1')
    c_E2 = openmc.Cell(fill=make_element(502, s_E2, fuel_temperature, 'fuel'),
                       region=(-s_E2), name='Position E2')
    c_E3 = openmc.Cell(fill=make_element(503, s_E3, fuel_temperature, 'fuel'),
                       region=(-s_E3), name='Position E3')
    c_E4 = openmc.Cell(fill=make_element(504, s_E4, fuel_temperature, 'fuel'),
                       region=(-s_E4), name='Position E4')
    c_E5 = openmc.Cell(fill=make_element(505, s_E5, fuel_temperature, 'fuel'),
                       region=(-s_E5), name='Position E5')
    c_E6 = openmc.Cell(fill=make_element(506, s_E6, fuel_temperature, 'fuel'),
                       region=(-s_E6), name='Position E6')
    c_E7 = openmc.Cell(fill=make_element(507, s_E7, fuel_temperature, 'fuel'),
                       region=(-s_E7), name='Position E7')
    c_E8 = openmc.Cell(fill=make_element(508, s_E8, fuel_temperature, 'fuel'),
                       region=(-s_E8), name='Position E8')
    c_E9 = openmc.Cell(fill=make_element(509, s_E9, fuel_temperature, 'fuel'),
                       region=(-s_E9), name='Position E9')
    c_E10 = openmc.Cell(fill=make_element(510, s_E10, fuel_temperature, 'fuel'),
                        region=(-s_E10), name='Position E10')
    c_E11 = openmc.Cell(fill=make_element(511, s_E11, fuel_temperature, 'fuel'),
                        region=(-s_E11), name='Position E11')
    c_E12 = openmc.Cell(fill=make_element(512, s_E12, fuel_temperature, 'fuel'),
                        region=(-s_E12), name='Position E12')
    c_E13 = openmc.Cell(fill=make_element(513, s_E13, fuel_temperature, 'fuel'),
                        region=(-s_E13), name='Position E13')
    c_E14 = openmc.Cell(fill=make_element(514, s_E14, fuel_temperature, 'fuel'),
                        region=(-s_E14), name='Position E14')
    c_E15 = openmc.Cell(fill=make_element(515, s_E15, fuel_temperature, 'fuel'),
                        region=(-s_E15), name='Position E15')
    c_E16 = openmc.Cell(fill=make_element(516, s_E16, fuel_temperature, 'fuel'),
                        region=(-s_E16), name='Position E16')
    c_E17 = openmc.Cell(fill=make_element(517, s_E17, fuel_temperature, 'fuel'),
                        region=(-s_E17), name='Position E17')
    c_E18 = openmc.Cell(fill=make_element(518, s_E18, fuel_temperature, 'fuel'),
                        region=(-s_E18), name='Position E18')
    c_E19 = openmc.Cell(fill=make_element(519, s_E19, fuel_temperature, 'fuel'),
                        region=(-s_E19), name='Position E19')
    c_E20 = openmc.Cell(fill=make_element(520, s_E20, fuel_temperature, 'fuel'),
                        region=(-s_E20), name='Position E20')
    c_E21 = openmc.Cell(fill=make_element(521, s_E21, fuel_temperature, 'fuel'),
                        region=(-s_E21), name='Position E21')
    c_E22 = openmc.Cell(fill=make_element(522, s_E22, fuel_temperature, 'fuel'),
                        region=(-s_E22), name='Position E22')
    c_E23 = openmc.Cell(fill=make_element(523, s_E23, fuel_temperature, 'fuel'),
                        region=(-s_E23), name='Position E23')
    c_E24 = openmc.Cell(fill=make_element(524, s_E24, fuel_temperature, 'fuel'),
                        region=(-s_E24), name='Position E24')

    # F-Ring
    c_F1 = openmc.Cell(fill=mats.water, region=(-s_F1),
                       name='Position F1')
    c_F2 = openmc.Cell(fill=mats.water, region=(-s_F2),
                       name='Position F2')
    c_F3 = openmc.Cell(fill=make_element(603, s_F3, fuel_temperature, 'fuel'),
                       region=(-s_F3), name='Position F3')
    c_F4 = openmc.Cell(fill=make_element(604, s_F4, fuel_temperature, 'fuel'),
                       region=(-s_F4), name='Position F4')
    c_F5 = openmc.Cell(fill=make_element(605, s_F5, fuel_temperature, 'fuel'),
                       region=(-s_F5), name='Position F5')
    c_F6 = openmc.Cell(fill=make_element(606, s_F6, fuel_temperature, 'fuel'),
                       region=(-s_F6), name='Position F6')
    c_F7 = openmc.Cell(fill=make_element(607, s_F7, fuel_temperature, 'fuel'),
                       region=(-s_F7), name='Position F7')
    c_F8 = openmc.Cell(fill=make_element(608, s_F8, fuel_temperature, 'fuel'),
                       region=(-s_F8), name='Position F8')
    c_F9 = openmc.Cell(fill=make_element(609, s_F9, fuel_temperature, 'fuel'),
                       region=(-s_F9), name='Position F9')
    c_F10 = openmc.Cell(fill=make_element(610, s_F10, fuel_temperature, 'fuel'),
                        region=(-s_F10), name='Position F10')
    c_F11 = openmc.Cell(fill=make_element(611, s_F11, fuel_temperature, 'fuel'),
                        region=(-s_F11), name='Position F11')
    c_F12 = openmc.Cell(fill=make_element(612, s_F12, core_water_temperature, 'icit'),
                        region=(-s_F12), name='Position F12')
    c_F13 = openmc.Cell(fill=make_element(613, s_F13, fuel_temperature, 'fuel'),
                        region=(-s_F13), name='Position F13')
    c_F14 = openmc.Cell(fill=make_element(614, s_F14, fuel_temperature, 'fuel'),
                        region=(-s_F14), name='Position F14')
    c_F15 = openmc.Cell(fill=make_element(615, s_F15, fuel_temperature, 'fuel'),
                        region=(-s_F15), name='Position F15')
    #no
    c_F16 = openmc.Cell(fill=make_element(616, s_F16, core_water_temperature, 'refl'),
                        region=(-s_F16), name='Position F16')
    #no
    c_F17 = openmc.Cell(fill=make_element(617, s_F17, core_water_temperature, 'refl'),
                        region=(-s_F17), name='Position F17')
    c_F18 = openmc.Cell(fill=make_element(618, s_F18, fuel_temperature, 'fuel'),
                        region=(-s_F18), name='Position F18')
    c_F19 = openmc.Cell(fill=make_element(619, s_F19, fuel_temperature, 'fuel'),
                        region=(-s_F19), name='Position F19')
    c_F20 = openmc.Cell(fill=make_element(620, s_F20, core_water_temperature, 'clicit'),
                        region=(-s_F20), name='Position F20')
    c_F21 = openmc.Cell(fill=make_element(621, s_F21, fuel_temperature, 'fuel'),
                        region=(-s_F21), name='Position F21')
    c_F22 = openmc.Cell(fill=make_element(622, s_F22, fuel_temperature, 'fuel'),
                        region=(-s_F22), name='Position F22')
    c_F23 = openmc.Cell(fill=make_element(623, s_F23, fuel_temperature, 'fuel'),
                        region=(-s_F23), name='Position F23')
    c_F24 = openmc.Cell(fill=make_element(624, s_F24, fuel_temperature, 'fuel'),
                        region=(-s_F24), name='Position F24')
    c_F25 = openmc.Cell(fill=make_element(625, s_F25, fuel_temperature, 'fuel'),
                        region=(-s_F25), name='Position F25')
    c_F26 = openmc.Cell(fill=make_element(626, s_F26, fuel_temperature, 'fuel'),
                        region=(-s_F26), name='Position F26')
    c_F27 = openmc.Cell(fill=make_element(627, s_F27, fuel_temperature, 'fuel'),
                        region=(-s_F27), name='Position F27')
    c_F28 = openmc.Cell(fill=make_element(628, s_F28, fuel_temperature, 'fuel'),
                        region=(-s_F28), name='Position F28')
    c_F29 = openmc.Cell(fill=make_element(629, s_F29, fuel_temperature, 'fuel'),
                        region=(-s_F29), name='Position F29')
    c_F30 = openmc.Cell(fill=make_element(630, s_F30, fuel_temperature, 'fuel'),
                        region=(-s_F30), name='Position F30')

    # G-Ring
    c_G1 = openmc.Cell(fill=mats.water, region=(-s_G1),
                       name='Position G1')
    c_G2 = openmc.Cell(fill=make_element(702, s_G2, core_water_temperature, 'rabbit'),
                       region=(-s_G2), name='Position G2')
    c_G3 = openmc.Cell(fill=mats.water, region=(-s_G3),
                       name='Position G3')
    c_G4 = openmc.Cell(fill=make_element(704, s_G4, core_water_temperature, 'target'),
                       region=(-s_G4), name='Position G4')
    c_G5 = openmc.Cell(fill=make_element(705, s_G5, core_water_temperature, 'target'),
                       region=(-s_G5), name='Position G5')
    c_G6 = openmc.Cell(fill=make_element(706, s_G6, core_water_temperature, 'target'),
                       region=(-s_G6), name='Position G6')
    c_G7 = openmc.Cell(fill=make_element(707, s_G7, core_water_temperature, 'target'),
                       region=(-s_G7), name='Position G7')
    c_G8 = openmc.Cell(fill=make_element(708, s_G8, core_water_temperature, 'target'),
                       region=(-s_G8), name='Position G8')
    c_G9 = openmc.Cell(fill=make_element(709, s_G9, core_water_temperature, 'target'),
                       region=(-s_G9), name='Position G9')
    c_G10 = openmc.Cell(fill=make_element(710, s_G10, core_water_temperature, 'target'),
                        region=(-s_G10), name='Position G10')
    c_G11 = openmc.Cell(fill=make_element(711, s_G11, core_water_temperature, 'target'),
                        region=(-s_G11), name='Position G11')
    c_G12 = openmc.Cell(fill=make_element(712, s_G12, fuel_temperature, 'fuel'),
                        region=(-s_G12), name='Position G12')
    c_G13 = openmc.Cell(fill=make_element(713, s_G13, fuel_temperature, 'fuel'),
                        region=(-s_G13), name='Position G13')
    c_G14 = openmc.Cell(fill=make_element(714, s_G14, fuel_temperature, 'fuel'),
                        region=(-s_G14), name='Position G14')
    c_G15 = openmc.Cell(fill=make_element(715, s_G15, fuel_temperature, 'fuel'),
                        region=(-s_G15), name='Position G15')
    #no
    c_G16 = openmc.Cell(fill=make_element(716, s_G16, core_water_temperature, 'target'),
                        region=(-s_G16), name='Position G16')
    c_G17 = openmc.Cell(fill=make_element(717, s_G17, core_water_temperature, 'source'),
                        region=(-s_G17), name='Position G17')
    #no
    c_G18 = openmc.Cell(fill=make_element(718, s_G18, core_water_temperature, 'target'),
                        region=(-s_G18), name='Position G18')
    #no
    c_G19 = openmc.Cell(fill=make_element(719, s_G19, core_water_temperature, 'refl'),
                        region=(-s_G19), name='Position G19')
    #no
    c_G20 = openmc.Cell(fill=make_element(720, s_G20, core_water_temperature, 'refl'),
                        region=(-s_G20), name='Position G20')
    #no
    c_G21 = openmc.Cell(fill=make_element(721, s_G21, core_water_temperature, 'refl'),
                        region=(-s_G21), name='Position G21')
    #no
    c_G22 = openmc.Cell(fill=make_element(722, s_G22, core_water_temperature, 'refl'),
                        region=(-s_G22), name='Position G22')
    c_G23 = openmc.Cell(fill=make_element(723, s_G23, fuel_temperature, 'fuel'),
                        region=(-s_G23), name='Position G23')
    c_G24 = openmc.Cell(fill=make_element(724, s_G24, fuel_temperature, 'fuel'),
                        region=(-s_G24), name='Position G24')
    c_G25 = openmc.Cell(fill=make_element(725, s_G25, fuel_temperature, 'fuel'),
                        region=(-s_G25), name='Position G25')
    c_G26 = openmc.Cell(fill=make_element(726, s_G26, fuel_temperature, 'fuel'),
                        region=(-s_G26), name='Position G26')
    c_G27 = openmc.Cell(fill=make_element(727, s_G27, core_water_temperature, 'refl'),
                        region=(-s_G27), name='Position G27')
    c_G28 = openmc.Cell(fill=make_element(728, s_G28, core_water_temperature, 'refl'),
                        region=(-s_G28), name='Position G28')
    c_G29 = openmc.Cell(fill=make_element(729, s_G29, core_water_temperature, 'refl'),
                        region=(-s_G29), name='Position G29')
    c_G30 = openmc.Cell(fill=make_element(730, s_G30, core_water_temperature, 'refl'),
                        region=(-s_G30), name='Position G30')
    c_G31 = openmc.Cell(fill=make_element(731, s_G31, core_water_temperature, 'refl'),
                        region=(-s_G31), name='Position G31')
    c_G32 = openmc.Cell(fill=make_element(732, s_G32, core_water_temperature, 'refl'),
                        region=(-s_G32), name='Position G32')
    c_G33 = openmc.Cell(fill=make_element(733, s_G33, core_water_temperature, 'refl'),
                        region=(-s_G33), name='Position G33')
    c_G34 = openmc.Cell(fill=make_element(734, s_G34, core_water_temperature, 'refl'),
                        region=(-s_G34), name='Position G34')
    c_G35 = openmc.Cell(fill=make_element(735, s_G35, core_water_temperature, 'refl'),
                        region=(-s_G35), name='Position G35')
    c_G36 = openmc.Cell(fill=make_element(736, s_G36, fuel_temperature, 'fuel'),
                        region=(-s_G36), name='Position G36')

    # List of the Grid Position Cells
    grid_cells = [c_A1,
                  c_B1,  c_B2,  c_B3,  c_B4,  c_B5,  c_B6,
                  c_C1,  c_C2,  c_C3,  c_C4,  c_C5,  c_C6,
                  c_C7,  c_C8,  c_C9,  c_C10, c_C11, c_C12,
                  c_D1,  c_D2,  c_D3,  c_D4,  c_D5,  c_D6,
                  c_D7,  c_D8,  c_D9,  c_D10, c_D11, c_D12,
                  c_D13, c_D14, c_D15, c_D16, c_D17, c_D18,
                  c_E1,  c_E2,  c_E3,  c_E4,  c_E5,  c_E6,
                  c_E7,  c_E8,  c_E9,  c_E10, c_E11, c_E12,
                  c_E13, c_E14, c_E15, c_E16, c_E17, c_E18,
                  c_E19, c_E20, c_E21, c_E22, c_E23, c_E24,
                  c_F1,  c_F2,  c_F3,  c_F4,  c_F5,  c_F6,
                  c_F7,  c_F8,  c_F9,  c_F10, c_F11, c_F12,
                  c_F13, c_F14, c_F15, c_F16, c_F17, c_F18,
                  c_F19, c_F20, c_F21, c_F22, c_F23, c_F24,
                  c_F25, c_F26, c_F27, c_F28, c_F29, c_F30,
                  c_G1,  c_G2,  c_G3,  c_G4,  c_G5,  c_G6,
                  c_G7,  c_G8,  c_G9,  c_G10, c_G11, c_G12,
                  c_G13, c_G14, c_G15, c_G16, c_G17, c_G18,
                  c_G19, c_G20, c_G21, c_G22, c_G23, c_G24,
                  c_G25, c_G26, c_G27, c_G28, c_G29, c_G30,
                  c_G31, c_G32, c_G33, c_G34, c_G35, c_G36]

    # Water Around The Assembly
    c_water_around = openmc.Cell(fill=mats.water, region=water_outside_region, name="Water Surrounding Assembly")
    c_water_around.temperature = bulk_water_temperature

    print("Cell Generation Complete.\n")

    # Create a root universe for the problem geometry. Every OpenMC model needs a root universe.
    root_univ = openmc.Universe()

    # Add all cells to the root universe
    root_univ.add_cells(assembly_cells)
    root_univ.add_cells(grid_cells)
    root_univ.add_cells([c_water_around])

    # Create the Geometry class
    geometry = openmc.Geometry(surface_precision=6)

    # Select the root universe
    geometry.root_universe = root_univ

    # Apply the geometry to the model (for the "geometry.xml" input file)
    model.geometry = geometry
    print("Geometry Generated.\n")

    ##########################################################################################################
    #                                                                                                        #
    # ---------------------------------------------- SETTINGS ---------------------------------------------- #
    #                                                                                                        #
    ##########################################################################################################

    # Create the Settings class
    settings = openmc.Settings()

    # Simulation settings (similar to MCNP KCODE)
    settings.particles = n_particles
    settings.inactive = n_inactive
    settings.batches = n_inactive + n_active
    settings.max_lost_particles = 10
    settings.rel_max_lost_particles = 1e-6

    # Simulation temperature specifications
    settings.temperature['default'] = core_water_temperature
    settings.temperature['method'] = 'interpolation'

    # Source Definition
    lower_left = (-27, -27, -45)
    upper_right = (27, 27, 150)
    source_region = openmc.stats.Box(lower_left, upper_right)
    source = openmc.IndependentSource(space=source_region, particle='neutron', constraints={'fissionable': True})
    settings.source = source

    # Apply the settings to the model (for the "settings.xml" input file)
    model.settings = settings

    print("Settings Generated.\n")

    ##########################################################################################################
    #                                                                                                        #
    # ------------------------------------------------ PLOTS ----------------------------------------------- #
    #                                                                                                        #
    ##########################################################################################################

    # Get color map for materials from mats
    cmap = mats.cmap

    # For all slices, ensure the plot.width is the same aspect ratio as plot.pixels, or the image will be stretched

    # A slice through the XZ plane, colored by cell
    plot1 = openmc.Plot()
    plot1.filename = 'Center_XZ_Slice_Cells'
    plot1.width = (200, 150)
    plot1.basis = 'xz'
    plot1.origin = (0, 0, 34.305)
    plot1.pixels = (4000, 3000)
    plot1.color_by = 'cell'
    plot1.id = 1000

    # A slice through the YZ plane, colored by cell
    plot2 = openmc.Plot()
    plot2.filename = 'Center_YZ_Slice_Cells'
    plot2.width = (200, 150)
    plot2.basis = 'yz'
    plot2.origin = (0, 0, 34.305)
    plot2.pixels = (4000, 3000)
    plot2.color_by = 'cell'
    plot2.id = 1001

    # A slice through a plane parallel to the XY plane, axially aligned with the beam port cans
    plot3 = openmc.Plot()
    plot3.filename = 'XY_Slice_BPs_Cells'
    plot3.width = (120, 120)
    plot3.basis = 'xy'
    plot3.origin = (0, 0, 20.6375)
    plot3.pixels = (4000, 4000)
    plot3.color_by = 'cell'
    plot3.id = 1003

    # A slice through the XZ plane, colored by cell
    plot11 = openmc.Plot()
    plot11.filename = 'Center_XZ_Slice_Mats'
    plot11.width = (200, 150)
    plot11.basis = 'xz'
    plot11.origin = (0, 0, 34.305)
    plot11.pixels = (4000, 3000)
    plot11.colors = cmap
    plot11.color_by = 'material'
    plot11.id = 1011

    # A slice through the YZ plane, colored by cell
    plot12 = openmc.Plot()
    plot12.filename = 'Center_YZ_Slice_Mats'
    plot12.width = (200, 150)
    plot12.basis = 'yz'
    plot12.origin = (0, 0, 34.305)
    plot12.pixels = (4000, 3000)
    plot12.colors = cmap
    plot12.color_by = 'material'
    plot12.id = 1012

    # A slice through a plane parallel to the XY plane, axially aligned with the beam port cans
    plot13 = openmc.Plot()
    plot13.filename = 'XY_Slice_BPs_Mats'
    plot13.width = (120, 120)
    plot13.basis = 'xy'
    plot13.origin = (0, 0, 20.6375)
    plot13.pixels = (4000, 4000)
    plot13.colors = cmap
    plot13.color_by = 'material'
    plot13.id = 1013

    plot_list = [plot1, plot2, plot3, plot11, plot12, plot13]

    plots = openmc.Plots(plot_list)

    # Specify the plots for the "plots.xml" input file
    model.plots = plots

    print("Plot Definitions Generated.\n")

    ##########################################################################################################
    #                                                                                                        #
    # ---------------------------------------------- TALLIES ----------------------------------------------- #
    #                                                                                                        #
    ##########################################################################################################

    # Any tallies would go here and be added to an openmc.Tallies() class and then applied to the model

    tallies = openmc.Tallies()
    
    
    np_tally = openmc.Tally(name='Np237_Analysis')
    np_tally.filters = [openmc.MaterialFilter(mats.npo2)] 
    
    # 2. Specify nuclides to differentiate reactions
    # This gives you the Pu238 production rate (via Np237 n,gamma)
    np_tally.nuclides = ['Np237'] 
    np_tally.scores = ['(n,gamma)', 'fission', 'absorption'] # Fission not revelant but included just cause
    
    tallies.append(np_tally)

    return model


def main():
    """

    The following define some arguments for use when running in the command line:

    -n              = (int)         Number of particle histories to track per cycle
    -i              = (int)         Number of inactive cycles to run before tracking statistics
    -a              = (int)         Number of active cycles for which statistics are tracked
    -m              = (bool)        Generate the model with temperatures at full power (1 MW)
    --run           = (bool)        Run the OpenMC model after it is generated
    --plot          = (bool)        Plot any specified plots from the model after it is generated
    --rods          = (float) x4    Four floats for the height of the control rods, i.e. <TR> <SA> <SH> <REG>
    --num_threads   = (int)         Number of OpenMP threads to run OpenMC with
    --dir           = (str)         Path for which to generate AND run OpenMC (if specified)

    An example command:

    python OSTR_OpenMC_Model.py -n 25000 -i 30 -a 470 -m --run --plot --rods 30 40 50 60 --dir ./folder --num_threads 16

    This will generate a model of the OSTR with the Transient Rod 30% withdrawn, Safety Rod 40% withdrawn, Shim Rod 50%
    withdrawn, and the Regulating Rod 60% withdrawn with temperatures set to those at 1 MW (full power). This model is
    set to track 25000 particles per cycle, with 30 inactive cycles and 470 active cycles (total of 500). The input
    files will be generated in "./folder". OpenMC will then run in that same directory using 16 threads, and plots will
    be generated in that same directory.

    """

    ap = ArgumentParser()

    ap.add_argument('-n', dest='n_particles', type=int, default=10000,
                    help='Number of particle histories per generation')
    ap.add_argument('-i', dest='n_inactive', type=int, default=50,
                    help='Number of inactive cycles')
    ap.add_argument('-a', dest='n_active', type=int, default=250,
                    help='Number of active cycles')
    ap.add_argument('-m', dest='megawatt', action='store_true',
                    help='Generate model with temperatures at full power (1 MW)')
    ap.add_argument('--run', action='store_true',
                    help='Run OpenMC')
    ap.add_argument('--plot', action='store_true',
                    help='Plot Geometry')
    ap.add_argument('--rods', dest='rods', nargs=4, type=float, default=[0, 0, 0, 0],
                    help='Rod Heights')
    ap.add_argument('--num_threads', dest='n_threads', type=int, default=os.cpu_count(),
                    help='Number of OpenMP threads for openmc.run()')
    ap.add_argument('--dir', dest='directory', type=str, default='.',
                    help='Directory for generating XML files and running OpenMC')
    ap.add_argument('--deplete', action='store_true',
                    help='Run depletion calculation')
    ap.add_argument('--days', type=float, default=30.0,
                    help='Depletion time in days')
    ap.add_argument('--power', type=float, default=1.0e6,
                    help='Reactor power in Watts')

    args = ap.parse_args()

    rod_heights = tuple(args.rods)
    for rod_height in rod_heights:
        if rod_height < 0 or rod_height > 100:
            raise ValueError('Rod height must be between 0 and 100')

    # Generate the Model
    model = ostr_model(args.n_particles, args.n_inactive, args.n_active, rod_heights, args.megawatt)

    # Generate the Input Files from the Model
    model.export_to_xml(directory=args.directory)
    print("*** Input Files Generated. ***\n")
    

    if args.run and not args.deplete:
        print("Beginning OpenMC Run\n")
        openmc.run(threads=args.n_threads, cwd=args.directory)
        print("\n\n-- RUN COMPLETE --\n\n")


    if args.deplete:


        chain_file = openmc.config['chain_file']

        openmc.lib.init()

        print("\n--- MATERIALS IN OPENMC LIB ---")

        for mat_id in openmc.lib.materials:
            print(mat_id)

        openmc.lib.finalize()
        
        # Change your operator to use fission-q
        op = openmc.deplete.CoupledOperator(
            model, 
            chain_file, 
            normalization_mode='fission-q', # ADD THIS LINE
            #directory=args.directory
        )


        # 30 daily steps
        time_steps = [1.0] * int(args.days)

        integrator = openmc.deplete.CECMIntegrator(
            operator=op,
            timesteps=time_steps,
            power=args.power,
            timestep_units='d'
        )

        integrator.integrate()

        print("\n\n-- DEPLETION COMPLETE --\n\n")


    if args.plot and not args.run:
        print("Plotting...\n")
        openmc.plot_geometry(cwd=args.directory)
        print("\n\n-- PLOTTING COMPLETE --\n\n")

    if args.plot and args.run:
        print("Plotting...\n")
        openmc.plot_geometry(cwd=args.directory)
        print("\n\n-- PLOTTING COMPLETE --\n\n")

        print("Beginning OpenMC Run\n")
        openmc.run(threads=args.n_threads, cwd=args.directory)
        print("\n\n-- RUN COMPLETE --\n\n")


if __name__ == '__main__':
    # This section allows the methods from this python script to be imported to other scripts without running main()
    main()
    print("/////      END OF SCRIPT      /////\n\n")
