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
        v_x = np.random.rand(self.nhumans, 1) * 1
        v_y = np.random.rand(self.nhumans, 1) * 1
        # v_z = 0

        self.velocities = np.array([v_x, v_y]).T #, v_z

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

    def save_to_gadget(self, filename, boxsize=10000):
        """
        Save the human data to a GADGET .hdf5 file.

        Uses the internal options, but you must specify a filename.
        """
        with h5.File(filename, "w") as handle:
            wg.write_header(
                handle,
                boxsize=boxsize,
                flag_entropy=0,
                np_total=np.array([self.nparts, 0, 0, 0, 0, 0]),
                np_total_hw=np.array([0, 0, 0, 0, 0, 0]),
                other={
                    "MassTable": np.array([self.particlemass, 0, 0, 0, 0, 0]),
                    "Time": 0,
                },
            )

            wg.write_runtime_pars(handle, periodic_boundary=1)

            wg.write_units(
                handle, current=1.0, length=1.0, mass=1, temperature=1.0, time=1.0
            )

            wg.write_block(
                handle,
                0,  # gas
                self.positions,
                self.velocities,
                self.ids,
                mass=self.masses,
                int_energy=self.internalenergy,
                smoothing=self.smoothing,
                other={"Density": self.densities},
            )

        return


def gen_humans_grid(meta):
    """
    Generates humans on a grid and returns a filled Humans object.
    """
    humans = Humans(meta)
    range = (0, meta["boxsize"])
    centre_of_ring = [meta["boxsize"] / 2.0] * 3

    # Because we are using a uniform grid we actually use the same x and y
    # range for the initial human setup.
    step = (range[1] - range[0]) / meta["nhumans"]

    x_values = np.arange(0, range[1] - range[0], step)

    # These are 2d arrays which isn't actually that helpful.
    x, y = np.meshgrid(x_values, x_values)
    x = x.flatten() + centre_of_ring[0] - (range[1] - range[0]) / 2
    y = y.flatten() + centre_of_ring[1] - (range[1] - range[0]) / 2
    # z = np.zeros_like(x) + meta["boxsize"] / 2

    humans.positions = np.array([x, y]).T #, z
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
        "-m",
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
        "-h",
        "--humanmass",
        help="""
             Mass of the humans. Default: 1.
             """,
        required=False,
        default=1.0,
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
        smoothing = float(ARGS["boxsize"]) / int(ARGS["nparts"])
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

    HUMANS.save_to_gadget(filename=ARGS["filename"], boxsize=ARGS["boxsize"])



#################################################

# Generates a SWIFT IC file with ...

# Parameters
periodic = 0  # 1 For periodic box
boxSize = 10  # 1 km
rho = 200  # Population density in code units ?
T = 1  # Initial intensity of human motion (how many people are moving in a box or percentage?)
gamma = 5.0 / 3.0  # Gas adiabatic index
fileName = "mobilityBox.hdf5"
# ---------------------------------------------------

# defines some constants
# need to be changed in plotTemperature.py too
h_frac = 0.76
mu = 4.0 / (1.0 + 3.0 * h_frac)

m_h_cgs = 1.67e-24
k_b_cgs = 1.38e-16

# defines units
unit_length = 1  # 1m
unit_mass = 1  # a unit mass of 1 particle-human
unit_time = 1  # 1 s ?

# Read id, position and h from glass
glass = h5.File("humans.hdf5", "r")
ids = glass["/PartType0/ParticleIDs"][:]
pos = glass["/PartType0/Coordinates"][:, :] * boxSize
h = glass["/PartType0/SmoothingLength"][:] * boxSize

# Compute basic properties

# need to define `pos`

numHum = np.size(pos) // 2
mass = boxSize ** 2 * rho # number of humans in a unit of territory
internalEnergy = k_b_cgs * T * mu / ((gamma - 1.0) * m_h_cgs)
internalEnergy *= (unit_time / unit_length) ** 2

# File
f = h5.File(fileName, "w")

# Header
grp = f.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumHum_Total"] = [numHum, 0, 0, 0, 0, 0]
grp.attrs["NumHum_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumHum_ThisFile"] = [numHum, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# grp.attrs["Flag_Entropy_ICs"] = 0

# Runtime parameters
grp = f.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

# Units
grp = f.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_length
grp.attrs["Unit mass in cgs (U_M)"] = unit_mass
grp.attrs["Unit time in cgs (U_t)"] = unit_time
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = f.create_group("/PartType0") # humans

v = np.zeros((numHum, 2))
ds = grp.create_dataset("Velocities", (numHum, 2), "f")
ds[()] = v

m = np.full((numHum, 1), mass)
ds = grp.create_dataset("Masses", (numHum, 1), "f")
ds[()] = m

h = np.reshape(h, (numHum, 1))
ds = grp.create_dataset("SmoothingLength", (numHum, 1), "f")
ds[()] = h

u = np.full((numHum, 1), internalEnergy)
ds = grp.create_dataset("InternalEnergy", (numHum, 1), "f")
ds[()] = u

ids = np.reshape(ids, (numHum, 1))
ds = grp.create_dataset("ParticleIDs", (numHum, 1), "L")
ds[()] = ids

ds = grp.create_dataset("Coordinates", (numHum, 2), "d")
ds[()] = pos

f.close()

print("Initial condition generated")
