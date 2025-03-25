/*******************************************************************************
 * This file is part of SWIFT_ABM.
 * Copyright (c) 2025 Dmitry Nikolaenko (dmitry.nikolaenko@durham.ac.uk)
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
#ifndef SWIFT_HUMANMOBILITY_ABM_UTILS_H
#define SWIFT_HUMANMOBILITY_ABM_UTILS_H

/* Standard headers */
#include <math.h>

/**
 * @brief A utility to set some initialisation parametres.
 *
 * @param h The human whose parameters are to be set.
 */
// __attribute__((always_inline)) INLINE static void _abm_rand_params(
//     struct part *restrict h, const float max_rand, const float div_rand) {

//   h->max_rand = max_rand;
//   h->div_rand = div_rand;
//   h->sub_rand = max_rand / div_rand / 2;
// }


/**
 * @brief Search for indices of the nearest grid points for a given position (x, y)
 *      in the 2D box of `accel_size`
 *      describing the area of size `box_size` in which the human is located
 *
 * @param x The x-coordinate of the human location (x, y).
 * @param y The y-coordinate of the human location (x, y).
 * @param box_size The size of the box in which the human is located.
 * @param accel_size The size of the acceleration grid for the box.
 * @param xi The lower x-index of the human location (x, y) in the grid.
 * @param yi The lower y-index of the human location (x, y) in the grid.
 * @param rx The relative anchor `x` of the human location (x, y) in the grid.
 * @param ry The relative anchor `y` of the human location (x, y) in the grid.
 */
__attribute__((always_inline))
INLINE static void find_nearest_indices_2D(double x, double y,
                                          const double* box_size,
                                          const int* accel_size,
                                          int* xi, int* yi,
                                          float *rx, float *ry) {
  *xi = (int)floor(x / box_size[0] * (accel_size[0]-1));
  *yi = (int)floor(y / box_size[1] * (accel_size[1]-1));
  *rx = fmod((x / box_size[0]) * (accel_size[0]-1), 1.0f);
  *ry = fmod((y / box_size[1]) * (accel_size[1]-1), 1.0f);
}



/**
 * @brief Bilinear interpolation of the acceleration grid for a given position (x, y)
 *
 * @param h The human whose parameters are to be set.
 */
__attribute__((always_inline)) INLINE static void bilinear_interpolation(
    const float ax_ll, const float ay_ll,   // Lower-left accelerations
    const float ax_lr, const float ay_lr,   // Lower-right accelerations  
    const float ax_ul, const float ay_ul,   // Upper-left accelerations
    const float ax_ur, const float ay_ur,   // Upper-right accelerations
    const float rx, const float ry,         // Relative position (0-1)
    float* ax_out, float* ay_out) {         // Output accelerations

    /* Weights for the four corners */
    const float w1 = (1.0f - rx) * (1.0f - ry);  // Lower-left
    const float w2 = rx * (1.0f - ry);           // Lower-right
    const float w3 = (1.0f - rx) * ry;           // Upper-left  
    const float w4 = rx * ry;                    // Upper-right

    /* Interpolate x-acceleration */
    *ax_out = ax_ll * w1 + ax_lr * w2 + ax_ul * w3 + ax_ur * w4;

    /* Interpolate y-acceleration */
    *ay_out = ay_ll * w1 + ay_lr * w2 + ay_ul * w3 + ay_ur * w4;
}


/**
 * @brief Check if a human is in the river (between the banks for a given position (x, y)).
 *
 * @param 
 */
__attribute__((always_inline)) INLINE static void check_human_in_river( ) {

}

/**
 * @brief A utility to advance random walk of a human.
 *
 * @param h The human.
 */
__attribute__((always_inline)) INLINE static void _kick_random_walk(
    struct part *restrict h, const float dx[3], const float max_rand, const float div_rand) {

  const float sub_rand = max_rand / div_rand / 2;

  /* Change a_hydro (only 2D for humans) */
  h->a_hydro[0] = (fmod((float)rand(), max_rand) / div_rand - sub_rand) * dx[0];
  h->a_hydro[1] = (fmod((float)rand(), max_rand) / div_rand - sub_rand) * dx[1];
  h->a_hydro[2] = 0.0; //((rand()%max_rand)/div_rand - sub_rand) * dx[2];
}


#endif /* SWIFT_HUMANMOBILITY_ABM_UTILS_H */
