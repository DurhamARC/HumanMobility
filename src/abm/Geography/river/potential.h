/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_POTENTIAL_RIVER_H
#define SWIFT_POTENTIAL_RIVER_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"
#include "abm/HumanMobility/abm_utils.h"
#include "common_io.h"

/**
 * @brief External Potential Properties - River
 *        (a lake could be another option)
 */
struct external_potential {

  /*! Position of the river (horizontal on the map) */
  double y[2]; // the northern and southern bank of the river

  /*! Mass */
  double mass;

  /*! Cached acceleration field data */
  float *ax;
  float *ay;
  double box_size[2];
  int grid_size[2];
};


/**
 * @brief Computes the time-step due to the acceleration from river
 *
 * @param time The current time.
 * @param potential The properties of the externa potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  return FLT_MAX;
}

/**
 * @brief Computes the gravitational acceleration of a human near river
 *
 * We change acceleration of humans if they go inside the river.
 *
 * @param time The current time.
 * @param potential The properties of the external potential (representing a river).
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data (representing a human).
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {
  
    double x, y;
    double box_size[2];
    int accel_size[2];
    float *accel_x, *accel_y;
    int xi, yi;
    float rx, ry;
    float ax_ll, ay_ll, ax_lr, ay_lr, ax_ul, ay_ul, ax_ur, ay_ur;
    float ax, ay;

    // get the location of a human
    x = g->x[0];
    y = g->x[1];
    box_size[0] = potential->box_size[0];
    box_size[1] = potential->box_size[1];
    accel_size[0] = potential->grid_size[0];
    accel_size[1] = potential->grid_size[1];
    accel_x = potential->ax;
    accel_y = potential->ay;

    // find the nearest indices of the location of a human
    find_nearest_indices_2D(x, y, box_size, accel_size, &xi, &yi, &rx, &ry); //+0.5*box_size[0]

    // correct the indices if they are out of the range
    if (xi < 0) xi = 0;
    if (xi >= accel_size[0]-1) xi = accel_size[0]-2;
    if (yi < 0) yi = 0;
    if (yi >= accel_size[1]-1) yi = accel_size[1]-2;

    // get x- and y-acceleration at the 4 anchor points
    ax_ll = accel_x[xi * accel_size[0] + yi];
    ay_ll = accel_y[xi * accel_size[0] + yi];
    ax_lr = accel_x[(xi+1) * accel_size[0] + yi];
    ay_lr = accel_y[(xi+1) * accel_size[0] + yi];
    ax_ul = accel_x[xi * accel_size[0] + (yi+1)];
    ay_ul = accel_y[xi * accel_size[0] + (yi+1)];
    ax_ur = accel_x[(xi+1) * accel_size[0] + (yi+1)];
    ay_ur = accel_y[(xi+1) * accel_size[0] + (yi+1)];

    // interpolate acceleration for a human at the location
    // (be careful of periodic boundary condition if applicable)
    bilinear_interpolation(ax_ll, ay_ll,
                           ax_lr, ay_lr,
                           ax_ul, ay_ul,
                           ax_ur, ay_ur,
                           rx,    ry,
                           &ax,   &ay);

  // if (ax > 1e+10 || ay > 1e+10)
  //   error("Acceleration is too high (human is inside river): %f %f", ax, ay);

  g->a_grav[0] += ax;
  g->a_grav[1] += ay;
  // g->a_grav[2] = az;

  // gravity_add_comoving_potential(g, value); // value ?
}

/**
 * @brief Computes the gravitational potential energy of a human near river.
 *
 * We return 0.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  return 0.f;
}

/**
 * 
*/
static INLINE void geography_read_acceleration_field(
    struct swift_params* parameter_file, struct external_potential* potential) {
#if defined(HAVE_HDF5)

  /*! Acceleration field filename */
  char filename[DESCRIPTION_BUFFER_SIZE];
  
  /* Read acceleration field file path */
  parser_get_param_string(parameter_file, "RiverPotential:parameter_file", filename);

  /* Load acceleration field data */

  /* Open file */
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("Unable to open file %s", filename);

  /* Read grid size and box size */

  /* Open group */
  hid_t group_id = H5Gopen(file_id, "Header", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group Header.\n");

  /* Read the arrays */
  io_read_array_attribute(group_id, "BoxSize", DOUBLE, potential->box_size, 2);
  io_read_array_attribute(group_id, "GridSize", INT, potential->grid_size, 2);

  /* Close group */
  hid_t status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  /* Allocate and read acceleration fields */
  const int size = potential->grid_size[0] * potential->grid_size[1];
  potential->ax = (float*)malloc(size * sizeof(float));
  potential->ay = (float*)malloc(size * sizeof(float));
  printf("size: %d\n", size);

  /* Open group */
  group_id = H5Gopen(file_id, "AccelerationField", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group AccelerationField.\n");

  /* Read the datasets */
  io_read_array_dataset(group_id, "ax", FLOAT, potential->ax, size);
  io_read_array_dataset(group_id, "ay", FLOAT, potential->ay, size);

  /* Close group */
  status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing file.");

#else
  message("Cannot read the acceleration field without HDF5");
#endif
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct space* s,
    struct external_potential* potential) {

  /* Read river banks position */
  parser_get_param_double_array(parameter_file,
                                "RiverPotential:position", // northern and southern bank coordinatates
                                2, potential->y);

  /* Read mass parameter */
  potential->mass =
      parser_get_param_double(parameter_file, "RiverPotential:mass");

  /* Load acceleration field data */
  geography_read_acceleration_field(parameter_file, potential);
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message("External potential is 'River'.");
}

/**
 * @brief Cleans up the external potential by freeing allocated memory.
 *
 * @param potential The external potential to clean up.
 */
static INLINE void potential_cleanup_backend(
  struct external_potential* potential) {

#if defined(HAVE_HDF5)
/* Free acceleration field arrays if they were allocated */
if (potential->ax != NULL) {
  free(potential->ax);
  potential->ax = NULL;
}

if (potential->ay != NULL) {
  free(potential->ay);
  potential->ay = NULL;
}
#endif

}
#endif /* SWIFT_POTENTIAL_RIVER_H */
