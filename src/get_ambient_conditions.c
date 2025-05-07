#include <CONVERGE/udf.h>
// #include "user_header.h"

CONVERGE_precision_t ambient_pres, ambient_temp;
CONVERGE_precision_t MW_Fuel;
CONVERGE_precision_t fuel_den, fuel_hvap, fuel_pvap, fuel_temp, fuel_cond, fuel_visc, noz_dia, noz_area, noz_len, vap_den, mdot, sigma, sound_speed;

// Transport
/************************************************************************************/
/* name: CONVERGE_BEFORE_TRANSPORT(<name>, IN(...), OUT(...))                       */
/*                                                                                  */
/* description: Routine that is called before transport is called.                  */
/*                                                                                  */
/* inputs: Field variables                                                          */
/*                                                                                  */
/* outputs: Field variables                                                         */
/*                                                                                  */
/************************************************************************************/
#define VCELL_VARS_IN \
	VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_precision_t *, volume), FIELD(CONVERGE_precision_t *, pressure), FIELD(CONVERGE_precision_t *, temperature)

#define VCELL_VARS_OUT \
       CONVERGE_VOID

CONVERGE_BEFORE_TRANSPORT(get_ambient_conditions,
		IN(VCELL_VARS_IN),
		OUT(VCELL_VARS_OUT))
{

	// Skip solid region
	if(CONVERGE_mesh_is_rigid(mesh)!=0) return;

	int rank;
	CONVERGE_mpi_comm_rank(&rank);

	CONVERGE_id_t reg_id = 0;
	CONVERGE_region_t chamber_region = CONVERGE_mesh_region_from_id(mesh, reg_id);


	/*Get name of chamber region */
	CONVERGE_string_literal_t name_chamber_region = CONVERGE_region_name(chamber_region);

	int count = 0;
	CONVERGE_precision_t vol_1, vol_tot, pres_1, pres_tot, temp_1, temp_tot;
	vol_1    = 0.0;
	pres_1   = 0.0;
	temp_1   = 0.0;
	vol_tot  = 0.0;
	pres_tot = 0.0;
	temp_tot = 0.0;

	int num_regions = CONVERGE_mesh_num_regions(mesh);

	/* Create an iterator over the cells of a region */
	CONVERGE_iterator_t rit;
	CONVERGE_region_iterator_create(chamber_region, &rit);
	for(CONVERGE_index_t kk = CONVERGE_iterator_first(rit); kk != -1; kk = CONVERGE_iterator_next(rit))
	{
		vol_1 += volume[kk];
		pres_1 += pressure[kk]*volume[kk];
		temp_1 += temperature[kk]*volume[kk];
	}

	CONVERGE_mpi_reduce(&vol_1, &vol_tot, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);
	CONVERGE_mpi_reduce(&pres_1, &pres_tot, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);
	CONVERGE_mpi_reduce(&temp_1, &temp_tot, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);

	if (rank==0)
	{
		ambient_pres = pres_tot/vol_tot;
		ambient_temp = temp_tot/vol_tot;
	}
	CONVERGE_mpi_bcast(&ambient_pres,1,CONVERGE_DOUBLE,0);
	CONVERGE_mpi_bcast(&ambient_temp,1,CONVERGE_DOUBLE,0);

	printf (" rank = %d ambient_pres = %.12e ambient_temp = %.12e \n", rank, ambient_pres, ambient_temp);
//
        // SPECIES
        const CONVERGE_species_t species = CONVERGE_mesh_species(mesh);
	const int num_total_species = CONVERGE_species_num_tot(species);
        const int num_gas_species   = CONVERGE_species_num_gas(species);

	printf (" number of total species = %d \n", num_total_species);

	if (rank ==0)
	{

		for(int isp = 0; isp < num_total_species; isp++)
		{
			if(isp == CONVERGE_species_index(species, "IC8H18"))
			{
				MW_Fuel = CONVERGE_species_mw(species, isp);
			}
		}
	}

	CONVERGE_mpi_bcast(&MW_Fuel,1,CONVERGE_DOUBLE,0);

	printf (" MW_Fuel = %.12e \n", MW_Fuel);
	printf("This is the end of get ambient conditions \n");



	return;
}

