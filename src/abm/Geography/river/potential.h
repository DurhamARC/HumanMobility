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

  /*! Mass */
  double mass;

  /*! Cached acceleration field data */
  float *ax;
  float *ay;
  double box_size[2];
  int grid_size[2];

  /*! River geometry data for multiple rivers */
  int num_rivers;             // Number of rivers
  float **left_bank_x;        // [num_rivers][bank_points[i]]
  float **left_bank_y;
  float **right_bank_x;
  float **right_bank_y;
  int *bank_points;           // Number of points for each river
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

  g->a_grav[0] += ax;
  g->a_grav[1] += ay;
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

  char filename[DESCRIPTION_BUFFER_SIZE];
  parser_get_param_string(parameter_file, "RiverPotential:parameter_file", filename);

  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("Unable to open file %s", filename);

  // Read grid size, box size, and mass
  hid_t group_id = H5Gopen(file_id, "Header", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group Header.\n");
  io_read_array_attribute(group_id, "BoxSize", DOUBLE, potential->box_size, 2);
  io_read_array_attribute(group_id, "GridSize", INT, potential->grid_size, 2);
  io_read_attribute(group_id, "Mass", DOUBLE, &potential->mass);
  hid_t status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  // Allocate and read acceleration fields
  const long long size = potential->grid_size[0] * potential->grid_size[1];
  potential->ax = (float*)malloc(size * sizeof(float));
  potential->ay = (float*)malloc(size * sizeof(float));
  printf("acceleration grid size: %lld\n", size);

  group_id = H5Gopen(file_id, "AccelerationField", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group AccelerationField.\n");
  io_read_array_dataset(group_id, "ax", FLOAT, potential->ax, size);
  io_read_array_dataset(group_id, "ay", FLOAT, potential->ay, size);
  status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  // Open RiverGeometry group and read number of rivers
  group_id = H5Gopen(file_id, "RiverGeometry", H5P_DEFAULT);
  if (group_id < 0) error("unable to open group RiverGeometry.\n");

  // Read num_rivers attribute
  int num_rivers = 1;
  if (H5Aexists(group_id, "num_rivers") > 0) {
    hid_t attr = H5Aopen(group_id, "num_rivers", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &num_rivers);
    H5Aclose(attr);
  }
  potential->num_rivers = num_rivers;

  // Allocate arrays for each river
  potential->left_bank_x  = (float**)malloc(num_rivers * sizeof(float*));
  potential->left_bank_y  = (float**)malloc(num_rivers * sizeof(float*));
  potential->right_bank_x = (float**)malloc(num_rivers * sizeof(float*));
  potential->right_bank_y = (float**)malloc(num_rivers * sizeof(float*));
  potential->bank_points  = (int*)malloc(num_rivers * sizeof(int));

  // Read each river's geometry
  for (int i = 0; i < num_rivers; i++) {
    char river_name[32];
    snprintf(river_name, sizeof(river_name), "river_%d", i);
    hid_t river_group = H5Gopen(group_id, river_name, H5P_DEFAULT);
    if (river_group < 0) error("unable to open group %s.\n", river_name);

    // Get size of bank arrays
    hid_t dataset = H5Dopen(river_group, "left_bank_x", H5P_DEFAULT);
    if (dataset < 0) error("unable to open dataset left_bank_x.\n");
    hid_t space = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space, dims, NULL);
    int npoints = dims[0];
    H5Sclose(space);
    H5Dclose(dataset);

    potential->bank_points[i] = npoints;
    potential->left_bank_x[i]  = (float*)malloc(npoints * sizeof(float));
    potential->left_bank_y[i]  = (float*)malloc(npoints * sizeof(float));
    potential->right_bank_x[i] = (float*)malloc(npoints * sizeof(float));
    potential->right_bank_y[i] = (float*)malloc(npoints * sizeof(float));

    io_read_array_dataset(river_group, "left_bank_x", FLOAT,  potential->left_bank_x[i],  npoints);
    io_read_array_dataset(river_group, "left_bank_y", FLOAT,  potential->left_bank_y[i],  npoints);
    io_read_array_dataset(river_group, "right_bank_x", FLOAT, potential->right_bank_x[i], npoints);
    io_read_array_dataset(river_group, "right_bank_y", FLOAT, potential->right_bank_y[i], npoints);

    status = H5Gclose(river_group);
    if (status < 0) error("error closing river group.");
  }

  status = H5Gclose(group_id);
  if (status < 0) error("error closing RiverGeometry group.");

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

  /* Load acceleration field data and mass parameter */
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

/* Free river geometry arrays if they were allocated */
if (potential->left_bank_x != NULL) {
  for (int i = 0; i < potential->num_rivers; i++) {
    if (potential->left_bank_x[i] != NULL) free(potential->left_bank_x[i]);
    if (potential->left_bank_y[i] != NULL) free(potential->left_bank_y[i]);
    if (potential->right_bank_x[i] != NULL) free(potential->right_bank_x[i]);
    if (potential->right_bank_y[i] != NULL) free(potential->right_bank_y[i]);
  }
  free(potential->left_bank_x);
  free(potential->left_bank_y);
  free(potential->right_bank_x);
  free(potential->right_bank_y);
  potential->left_bank_x = NULL;
  potential->left_bank_y = NULL;
  potential->right_bank_x = NULL;
  potential->right_bank_y = NULL;
}

if (potential->bank_points != NULL) {
  free(potential->bank_points);
  potential->bank_points = NULL;
}
#endif

}
#endif /* SWIFT_POTENTIAL_RIVER_H */
