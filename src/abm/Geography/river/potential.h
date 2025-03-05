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

/**
 * @brief External Potential Properties - River
 *        (a lake could be another option)
 */
struct external_potential {

  /*! Position of the river (horizontal on the map) */
  double y[2]; // the northern and southern bank of the river

  /*! Mass */
  double mass;

  /*! Width */
  // double width;

  /*! Time-step condition pre-factor */
  // float timestep_mult;
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
 * @param potential The properties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data (representing a human).
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {
  
    float x, y;
    float dx, dy;
    int box_size[2], accel_size[2];
    float* accel_x, accel_y;
    int index_low[2], index_high[2];
    float *rx, *ry;
    float ax_ll, ay_ll, ax_lr, ay_lr, ax_ul, ay_ul, ax_ur, ay_ur;
    float *ax, *ay;

    const double box_size[2] = {e->s->dim[0], e->s->dim[1]};

    read_acceleration_field(filename, accel_x, accel_y, box_size); // only once!

    // within the 4 anchor points, where are we:
    // x_rel = ((x / box_size_x) * accel_size[0]) % 1;
    // y_rel = ((y / box_size_y) * accel_size[1]) % 1;
  
    // get the location of a human
    x = g->x[0];
    y = g->x[1];

    // find the nearest indices of the location of a human
    find_nearest_indices_2D(x, y, box_size, accel_size, index_low, index_high, rx, ry);

    // get x- and y-acceleration at the 4 anchor points
    ax_ll = accel_x[index_low[0] + index_low[1] * accel_size[0]];
    ay_ll = accel_y[index_low[0] + index_low[1] * accel_size[0]];
    ax_lr = accel_x[index_high[0] + index_low[1] * accel_size[0]];
    ay_lr = accel_y[index_high[0] + index_low[1] * accel_size[0]];
    ax_ul = accel_x[index_low[0] + index_high[1] * accel_size[0]];
    ay_ul = accel_y[index_low[0] + index_high[1] * accel_size[0]];
    ax_ur = accel_x[index_high[0] + index_high[1] * accel_size[0]];
    ay_ur = accel_y[index_high[0] + index_high[1] * accel_size[0];

    // interpolate acceleration for a human at the location
    // (be careful of periodic boundary condition if applicable)
    bilinear_interpolation(ax_ll, ay_ll,
                           ax_lr, ay_lr,
                           ax_ul, ay_ul,
                           ax_ur, ay_ur,
                           rx,    ry,
                           ax,    ay);

    // particle_i.velocity.x += particle_i_accel_x * dt
    // particle_i.velocity.y += particle_i_accel_y * dt

    dy = get_distance_to_river(g->x[0], g->x[1], potential->y[0], potential->y[1]);

/*     
    if(g->x[1] < potential->y[0])
      dy = g->x[1] - potential->y[0];
    else if(g->x[1] > potential->y[1])
      dy = g->x[1] - potential->y[1];
    else {
      dy = 0;
      message("Human is inside the river!!!");
    }

  const float rinv = 1.f / sqrtf(dy * dy);
  const float rinv3 = rinv * rinv * rinv;

  // The acceleration must change when the human is in the river
  // g->a_grav[0] = ax + potential->mass * dx * rinv3;
  g->a_grav[1] = ay + potential->mass * dy * rinv3;
  // g->a_grav[2] = az + potential->mass * dz * rinv3;
 */
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

  /* Read in the position of the centre of potential */
  parser_get_param_double_array(parameter_file, "RiverPotential:position", // northern and southern bank coordinatates
                                2, potential->y);

  /* Read the other parameters of the model */
  potential->mass =
      parser_get_param_double(parameter_file, "RiverPotential:mass");

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

#endif /* SWIFT_POTENTIAL_RIVER_H */
