###############################################################################
# This file is part of HumanMobility which is fork of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durhama.ac.uk)
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
# Copyright (c) 2024 Dmitry Nikolaenko (dmitry.nikolaenko@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

import numpy as np
import h5py as h5
import write_gadget as wg


class Humans(object):
    """
    Holder class for human properties. Also contains some methods to e.g.
    set their ... velocities based on their positions. These properties
    are set using the 'generationmethod' functions below.
    """

    def __init__(self, meta):
        self.gravitymass = meta["gravitymass"]
        self.nhumans = meta["nhumans"] ** 2
        self.humanmass = meta["humanmass"]
        self.boxsize = meta["boxsize"]

        self.smoothing = np.zeros(self.nhumans) + meta["smoothing"]
        self.internalenergy = np.zeros(self.nhumans) + meta["internalenergy"]

        self.positions = np.array([])
        self.radii = np.array([])
        self.theta = np.array([])
        self.phi = np.array([])
        self.velocities = np.array([])
        self.ids = np.array([])
        self.densities = np.array([])
        self.masses = np.array([])

        return

    def calculate_velocities(self, angle=0):
        """
        Calculates velocities of humans.
        """
        v_x = (np.random.rand(self.nhumans)-0.5)*1000 #np.zeros(self.nhumans) #
        v_y = (np.random.rand(self.nhumans)-0.5)*1000 #np.zeros(self.nhumans) #
        v_z = np.zeros(self.nhumans)

        self.velocities = np.array([v_x, v_y, v_z]).T

        return self.velocities

    def calculate_masses(self):
        """
        Calculate the individual masses for the humans,
        currently, just 1 mass units.
        """
        mass_factor = self.humanmass
        self.masses = self.densities * mass_factor

        return self.masses

    def generate_ids(self):
        """
        Generate consecutive IDs from 0 based on the number of humans
        currently in the object.
        """
        self.ids = np.arange(self.nhumans)

        return self.ids

    def save_to_gadget(self, filename, boxsize=10000, type="particles"):
        """
        Save the human data to a GADGET .hdf5 file.

        Uses the internal options, but you must specify a filename.
    
        Parameters:
            filename (str): The output HDF5 filename
            boxsize (float): Size of the simulation box
            type (str): Type of particles to write - either "gas" or "particles"
        """
        # Validate particle type
        if type not in ["gas", "particles"]:
            raise ValueError('type must be either "gas" or "particles"')
            
        # Map particle type to index
        type_index = 0 if type == "gas" else 1
        
        # Create arrays with the correct particle type index
        np_total = np.zeros(6, dtype=int)
        np_total[type_index] = self.nhumans
        
        mass_table = np.zeros(6)
        mass_table[type_index] = self.humanmass

        with h5.File(filename, "w") as handle:
            wg.write_header(
                handle,
                boxsize=[boxsize, boxsize],
                flag_entropy=0,
                np_total=np_total,
                np_total_hw=np.array([0, 0, 0, 0, 0, 0]),
                other={
                    "MassTable": mass_table,
                    "Time": 0,
                    "Dimension": 2,
                    "Flag_Entropy_ICs": 0,
                },
            )

            wg.write_runtime_pars(handle, periodic_boundary=1)

            wg.write_units(
                handle, current=1.0, length=1.0, mass=1, temperature=1.0, time=1.0
            )

            wg.write_block(
                handle,
                type_index,  # gas, dark matter particles
                self.positions,
                self.velocities,
                self.ids,
                mass=self.masses,
                int_energy=self.internalenergy,
                smoothing=self.smoothing,
#                other={"Density": self.densities},
            )

        return


def gen_humans_grid(meta):
    """
    Generates humans on a grid and returns a filled Humans object.
    """
    humans = Humans(meta)
    positions = (0, meta["boxsize"]) # -0.5*meta["boxsize"]
    # centre_of_ring = [meta["boxsize"] * 0.0] * 3

    # Because we are using a uniform grid we actually use the same x and y
    # positions for the initial human setup.
    width = positions[1] - positions[0]
    step = width / meta["nhumans"]

    x_values = np.arange(0, width, step, dtype=float) # -0.5*width

    # These are 2d arrays which isn't actually that helpful.
    x, y = np.meshgrid(x_values, x_values)
    x = x.flatten() + np.random.rand(humans.nhumans) * 5
    y = y.flatten() + np.random.rand(humans.nhumans) * 5
    z = np.zeros(humans.nhumans)
    # z = np.zeros_like(x) + meta["boxsize"] / 2

    humans.positions = np.array([x, y, z]).T
    humans.densities = np.ones(humans.nhumans)
    humans.calculate_velocities()
    humans.calculate_masses()

    humans.generate_ids()

    return humans


if __name__ == "__main__":
    import argparse as ap

    PARSER = ap.ArgumentParser(
        description="""
                    Initial conditions generator for the Human Mobility
                    example. It has <..> defaults for <..>, but if you
                    wish to run the example with <..> you sould use
                    --generationmethod <..>.
                    """
    )

    PARSER.add_argument(
        "-p",
        "--gravitymass",
        help="""
             GM for the central point mass. Default: 1.
             """,
        required=False,
        default=1.0,
    )

    PARSER.add_argument(
        "-f",
        "--filename",
        help="""
             Filename for your initial conditions.
             Default: initial_conditions.hdf5.
             """,
        required=False,
        default="initial_conditions.hdf5",
    )

    PARSER.add_argument(
        "-n",
        "--nhumans",
        help="""
             Square-root of the number of humans, i.e. the default
             nhumans=100 leads to a square with 100^2 humans in it.
             """,
        required=False,
        default=100,
    )

    PARSER.add_argument(
        "-m",
        "--humanmass",
        help="""
             Mass of the humans. Default: 1.
             """,
        required=False,
        default=1e5,
    )

    PARSER.add_argument(
        "-s",
        "--smoothing",
        help="""
             Initial smoothing length for all of the humans.
             Default: Boxsize/N
             """,
        required=False,
        default=-1,
    )

    PARSER.add_argument(
        "-i",
        "--internalenergy",
        help="""
             Initial internal energy for all of the humans. Not used currently
             for the human population.
             Default: 1.
             """,
        required=False,
        default=1.0,
    )

    PARSER.add_argument(
        "-g",
        "--generationmethod",
        help="""
             Generation method for the humans. Currently, only grid is implemented
             where the humans are generated
             in a way that minimises the energy in SPH. For more details on
             this method see Cartwright, Stamatellos & Whitworth (2009).
             Default: grid.
             """,
        required=False,
        default="grid",
    )

    PARSER.add_argument(
        "-b",
        "--boxsize",
        help="""
             The box size.
             Default: 10000 m (10 km)
             """,
        required=False,
        default=10000,
    )

    PARSER.add_argument(
        "-t",
        "--type",
        help="""
            Type of particles to use - either 'gas' or 'particles'.
            'gas' will save humans as SPH particles (PartType0),
            'particles' will save them as dark matter particles (PartType1).
            Default: gas
            """,
        required=False,
        choices=['gas', 'particles'],
        default='gas',
    )

    ### --- ### --- Argument Parsing --- ### --- ###

    ARGS = vars(PARSER.parse_args())

    if ARGS["generationmethod"] == "grid":
        gen_humans = gen_humans_grid
    else:
        print(
            "ERROR: {} is an invalid generation method. Exiting.".format(
                ARGS["generationmethod"]
            )
        )
        exit(1)

    if ARGS["smoothing"] == -1:
        smoothing = float(ARGS["boxsize"]) / int(ARGS["nhumans"])
    else:
        smoothing = float(ARGS["smoothing"])

    META = {
        "gravitymass": float(ARGS["gravitymass"]),
        "nhumans": int(ARGS["nhumans"]),
        "humanmass": float(ARGS["humanmass"]),
        "smoothing": smoothing,
        "internalenergy": float(ARGS["internalenergy"]),
        "boxsize": float(ARGS["boxsize"]),
    }

    HUMANS = gen_humans(META)

    # For SPH gas particles (PartType0) or dark matter particles (PartType1)
    HUMANS.save_to_gadget(
        filename=ARGS["filename"], 
        boxsize=ARGS["boxsize"], 
        type=ARGS["type"]
    )
    print("Initial condition generated")

