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
#ifndef SWIFT_ABM_H
#define SWIFT_ABM_H

/* Local headers. */
// #include "abm_properties.h"

/* Import the right functions */
#if defined(NONE_ABM)
#include "./None/abm.h"
#include "./None/abm_iact.h"
#define ABM_IMPLEMENTATION "No ABM scheme"
#elif defined(HUMANMOBILITY_ABM)
#include "./HumanMobility/abm.h"
#include "./HumanMobility/abm_iact.h"
#define ABM_IMPLEMENTATION "HUMANMOBILITY ABM"
#else
#error "Invalid choice of ABM variant"
#endif

#endif /* SWIFT_ABM_H */
