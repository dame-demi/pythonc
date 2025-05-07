/*******************************************************************************
* CONVERGENT SCIENCE CONFIDENTIAL                                              *
* All rights reserved.                                                         *
* All information contained herein is the property of Convergent Science.      *
* The intellectual and technical concepts contained herein are                 *
* proprietary to Convergent Science.                                           *
* Dissemination of this information or reproduction of this material           *
* is strictly forbidden unless prior written permission is obtained from       *
* Convergent Science.                                                          *
*******************************************************************************/
#include "lagrangian/env.h"

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**********************************************************************/
/*                                                                    */
/* Name: user_parcel_prop                                             */
/*                                                                    */
/* Description: user_parcel_prop sets values for the     custom       */
/* parcel properties when new parcels are created. This is            */
/* accomplished through the four functions: user_parcel_inject,       */
/* user_parcel_child, user_parcel_film, user_parcel_splash            */
/*                                                                    */
/*                                                                    */
/**********************************************************************/


// initialize values for the custom parcel properties when new parcels
// are first injected into the domain
CONVERGE_UDF(parcel_inject,
             IN(VALUE(CONVERGE_index_t, passed_parcel_idx), VALUE(CONVERGE_cloud_t, passed_spray_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&parcel_cloud, passed_spray_cloud);

   parcel_semi_mass_old = parcel_cloud.density[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
                          parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx];

   parcel_cloud.density[passed_parcel_idx]     = parcel_cloud.density[passed_parcel_idx] * 1.0001;
   parcel_cloud.radius[passed_parcel_idx]      = parcel_cloud.radius[passed_parcel_idx];
   parcel_cloud.density_tm1[passed_parcel_idx] = parcel_cloud.density_tm1[passed_parcel_idx];

   parcel_semi_mass_new = parcel_cloud.density[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
                          parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx];

   parcel_cloud.num_drop[passed_parcel_idx] =
      parcel_cloud.num_drop[passed_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;
}

// initialize values for the custom parcel properties when new child parcels
// are created from the Kelvin-Helmholtz stripping breakup mechanism

// NOTE: use the following syntax to set the child parcel value equal to its parent's value

CONVERGE_UDF(parcel_child,
             IN(VALUE(CONVERGE_index_t, passed_child_parcel_idx),
                VALUE(CONVERGE_index_t, passed_parent_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&parcel_cloud, passed_spray_cloud);

   parcel_semi_mass_old = parcel_cloud.density[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx] *
                          parcel_cloud.radius[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx];

   parcel_cloud.density[passed_child_parcel_idx]     = parcel_cloud.density[passed_parent_parcel_idx];
   parcel_cloud.radius[passed_child_parcel_idx]      = parcel_cloud.radius[passed_parent_parcel_idx];
   parcel_cloud.density_tm1[passed_child_parcel_idx] = parcel_cloud.density_tm1[passed_parent_parcel_idx];

   parcel_semi_mass_new = parcel_cloud.density[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx] *
                          parcel_cloud.radius[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx];

   parcel_cloud.num_drop[passed_child_parcel_idx] =
      parcel_cloud.num_drop[passed_child_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;
}
// set values for the custom parcel properties when film parcels separate
// from a surface and are converted to spray parcels

// NOTE: set the properties equal to themselves if they are to retain their
//       pre-separation values


CONVERGE_UDF(parcel_splash,
             IN(VALUE(CONVERGE_index_t, passed_spray_parcel_idx),
                VALUE(CONVERGE_index_t, passed_film_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud),
                VALUE(CONVERGE_cloud_t, passed_film_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud spray_parcel_cloud;
   struct ParcelCloud film_parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&spray_parcel_cloud, passed_spray_cloud);
   load_user_cloud(&film_parcel_cloud, passed_film_cloud);

   parcel_semi_mass_old =
      spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
      spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.density[passed_spray_parcel_idx]     = film_parcel_cloud.density[passed_film_parcel_idx];
   spray_parcel_cloud.radius[passed_spray_parcel_idx]      = film_parcel_cloud.radius[passed_film_parcel_idx];
   spray_parcel_cloud.density_tm1[passed_spray_parcel_idx] = film_parcel_cloud.density_tm1[passed_film_parcel_idx];

   parcel_semi_mass_new =
      spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
      spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.num_drop[passed_spray_parcel_idx] =
      spray_parcel_cloud.num_drop[passed_spray_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;
}
// initialize values for the custom parcel properties when new parcels are created
// from film stripping

// NOTE: use the following syntax to set the splashed parcel value equal to the impinged parcel's value

CONVERGE_UDF(parcel_strip,
             IN(VALUE(CONVERGE_index_t, passed_spray_parcel_idx),
                VALUE(CONVERGE_index_t, passed_film_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud),
                VALUE(CONVERGE_cloud_t, passed_film_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud spray_parcel_cloud;
   struct ParcelCloud film_parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&spray_parcel_cloud, passed_spray_cloud);
   load_user_cloud(&film_parcel_cloud, passed_film_cloud);

   parcel_semi_mass_old =
      spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
      spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.density[passed_spray_parcel_idx]     = film_parcel_cloud.density[passed_film_parcel_idx];
   spray_parcel_cloud.radius[passed_spray_parcel_idx]      = film_parcel_cloud.radius[passed_film_parcel_idx];
   spray_parcel_cloud.density_tm1[passed_spray_parcel_idx] = film_parcel_cloud.density_tm1[passed_film_parcel_idx];

   parcel_semi_mass_new =
      spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
      spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.num_drop[passed_spray_parcel_idx] =
      spray_parcel_cloud.num_drop[passed_spray_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;
}
