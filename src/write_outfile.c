#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>

#define CONVERGE_OUTPUT_VARS \
   VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_precision_t*, volume), FIELD(CONVERGE_precision_t*, pressure)

/**********************************************************************************************************/
/* macro: CONVERGE_OUTPUT                                                                                 */
/* name: user defined                                                                                     */
/* notes: Allows for use of setup, init, final, and close stages                                          */
/*                                                                                                        */
/* inputs:                                                                                                */
/*          Mesh & Field variables                                                                        */
/*                                                                                                        */
/*          void* passed_file (optional, only if passed_file is not an output)                            */
/*             User buffer that holds the pointer to the file object.                                     */
/*                                                                                                        */
/* outputs:                                                                                               */
/*          Mesh & Field variables                                                                        */
/*                                                                                                        */
/*          void* passed_file (optional, only if passed_file is not an input)                             */
/*             User buffer that holds the pointer to the file object.                                     */
/*                                                                                                        */
/**********************************************************************************************************/

CONVERGE_OUTPUT(outfile_method_setup,
                IN(VALUE(CONVERGE_mesh_t, mesh)),
                OUT(REF(CONVERGE_outfile_handle_t*, passed_file)))
{
   CONVERGE_logger_debug("Calling user_output setup");
   CONVERGE_index_t user_num_regions = CONVERGE_mesh_num_regions(mesh);

   int rank;
   CONVERGE_mpi_comm_rank(&rank);


   if(rank == 0)
   {
      // allocate udf data for output outfile handles of each of the regions in the mesh and the whole mesh
      CONVERGE_outfile_handle_t* outputs_file = malloc((user_num_regions + 1) * sizeof(CONVERGE_outfile_handle_t));

      *passed_file = outputs_file;
      char filename[31];
      CONVERGE_outfile_header_t header;
      CONVERGE_outfile_row_template_t row;

      CONVERGE_outfile_header_create(2, &header);
      CONVERGE_outfile_row_template_create(&row);

      // One universal row header for each user_outputs file
      CONVERGE_outfile_header_add_column(header,
                                         0,
                                         "Pressure",
                                         "Max_pres",
                                         "Min_pres",
                                         CONVERGE_OUTFILE_HEADER_ARGS_END);
      CONVERGE_outfile_header_add_column(header, 1, "(MPa)", "(MPa)", "(MPa)", CONVERGE_OUTFILE_HEADER_ARGS_END);
      CONVERGE_outfile_row_template_add_column(row,
                                               CONVERGE_outfile_DBL,
                                               CONVERGE_outfile_DBL,
                                               CONVERGE_outfile_DBL,
                                               CONVERGE_outfile_TEMPLATE_ARGS_END);

      // Iterating over each region file
      CONVERGE_iterator_t it;
      CONVERGE_mesh_region_index_iterator_create(mesh, &it);
      CONVERGE_index_t ii = 0;
      for(CONVERGE_index_t ireg = CONVERGE_iterator_first(it); ireg != -1; ireg = CONVERGE_iterator_next(it))
      {
         CONVERGE_outfile_handle_t* file_handle = &(outputs_file[ii]);

         sprintf(filename, "user_outfile_method_region%d", CONVERGE_mesh_region_id(mesh, ireg));

         CONVERGE_outfile_header_t dup_header;
         CONVERGE_outfile_row_template_t dup_row;

         CONVERGE_outfile_header_dup(header, &dup_header);
         CONVERGE_outfile_row_template_dup(row, &dup_row);

         CONVERGE_outfile_create(filename,
                                 CONVERGE_OUTFILE_NORMAL_COL_WIDTH,
                                 dup_header,
                                 dup_row,
                                 CONVERGE_OUTFILE_FEATURES_TIMESTEP_COL,
                                 file_handle);

         CONVERGE_outfile_row_template_destroy(&dup_row);
         CONVERGE_outfile_header_destroy(dup_header);
         ii++;
      }
      CONVERGE_iterator_destroy(&it);

      // Writing header in the file for the whole mesh
      sprintf(filename, "user_outfile_method_%s", CONVERGE_mesh_name(mesh));
      CONVERGE_outfile_create(filename,
                              CONVERGE_OUTFILE_NORMAL_COL_WIDTH,
                              header,
                              row,
                              CONVERGE_OUTFILE_FEATURES_TIMESTEP_COL,
                              &(outputs_file[user_num_regions]));
      CONVERGE_outfile_row_template_destroy(&row);
      CONVERGE_outfile_header_destroy(header);
   }

   return;
}

CONVERGE_OUTPUT(outfile_method,
                IN(CONVERGE_OUTPUT_VARS, VALUE(CONVERGE_outfile_handle_t*, passed_file)),
                OUT(CONVERGE_VOID))
{
   CONVERGE_logger_debug("Calling user_output");
   // Variables available through CONVERGE_OUTPUT_VARS
   const CONVERGE_precision_t* global_volume = volume;
   const CONVERGE_precision_t* global_pres   = pressure;

   // local vars
   CONVERGE_index_t user_num_regions = CONVERGE_mesh_num_regions(mesh);
   CONVERGE_precision_t vol_tot, vol_tot_1;
   CONVERGE_precision_t pres_avg, pres_avg_1;
   CONVERGE_precision_t pres_peak, pres_peak_1;
   CONVERGE_precision_t pres_min, pres_min_1;
   CONVERGE_precision_t vol_tot_reg, vol_tot_reg_1;
   CONVERGE_precision_t pres_avg_reg, pres_avg_reg_1;
   CONVERGE_precision_t pres_peak_reg, pres_peak_reg_1;
   CONVERGE_precision_t pres_min_reg, pres_min_reg_1;
   const CONVERGE_outfile_handle_t* outputs_file = passed_file;

   int rank = CONVERGE_mesh_get_mpi_rank(mesh);
   /*int rank;*/
   /*MPI_Comm_rank(MPI_COMM_WORLD, &rank);*/

   // Intializing to 0 for complete domain calculation
   vol_tot_1   = 0.0;
   pres_avg_1  = 0.0;
   pres_peak_1 = 0.0;
   pres_min_1  = 1e20;

   // Iterating over each region
   CONVERGE_index_t ii = 0;
   CONVERGE_iterator_t it;
   CONVERGE_mesh_region_index_iterator_create(mesh, &it);
   for(CONVERGE_index_t ireg = CONVERGE_iterator_first(it); ireg != -1; ireg = CONVERGE_iterator_next(it))
   {
      // Initializing to 0 for region index ireg
      vol_tot_reg_1   = 0.0;
      pres_avg_reg_1  = 0.0;
      pres_peak_reg_1 = 0.0;
      pres_min_reg_1  = 1e20;
      // Iterating over cells in region index ireg
      CONVERGE_iterator_t rit;
      CONVERGE_mesh_region_iterator_create(mesh, ireg, &rit);
      for(CONVERGE_index_t kk = CONVERGE_iterator_first(rit); kk != -1; kk = CONVERGE_iterator_next(rit))
      {
         vol_tot_1 += global_volume[kk];
         vol_tot_reg_1 += global_volume[kk];
         pres_avg_1 += global_pres[kk] * global_volume[kk];
         pres_avg_reg_1 += global_pres[kk] * global_volume[kk];

         if(global_pres[kk] > pres_peak_1)
            pres_peak_1 = global_pres[kk];
         if(global_pres[kk] > pres_peak_reg_1)
            pres_peak_reg_1 = global_pres[kk];
         if(global_pres[kk] < pres_min_1)
            pres_min_1 = global_pres[kk];
         if(global_pres[kk] < pres_min_reg_1)
            pres_min_reg_1 = global_pres[kk];
      }
      CONVERGE_mesh_region_iterator_destroy(&rit);

      // MPI operations for the region
      CONVERGE_mpi_reduce(&vol_tot_reg_1, &vol_tot_reg, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);
      CONVERGE_mpi_reduce(&pres_avg_reg_1, &pres_avg_reg, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);
      if(vol_tot_reg > 0.0)
         pres_avg_reg = pres_avg_reg / vol_tot_reg;
      CONVERGE_mpi_reduce(&pres_peak_reg_1, &pres_peak_reg, 1, CONVERGE_DOUBLE, CONVERGE_MPI_MAX, 0);
      CONVERGE_mpi_reduce(&pres_min_reg_1, &pres_min_reg, 1, CONVERGE_DOUBLE, CONVERGE_MPI_MIN, 0);

      // Writing region based output file
      if(rank == 0)
      {
         CONVERGE_outfile_data_row_t row;
         CONVERGE_outfile_data_row_create(outputs_file[ii], &row);
         CONVERGE_outfile_data_row_add_data(row, 3, pres_avg_reg, pres_peak_reg, pres_min_reg);
         CONVERGE_outfile_add_row(outputs_file[ii], row);
         CONVERGE_outfile_data_row_destroy(row);
      }
      ii++;
   }
   CONVERGE_iterator_destroy(&it);

   // MPI operations for the whole mesh
   CONVERGE_mpi_reduce(&vol_tot_1, &vol_tot, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);
   CONVERGE_mpi_reduce(&pres_avg_1, &pres_avg, 1, CONVERGE_DOUBLE, CONVERGE_MPI_SUM, 0);
   if(vol_tot > 0.0 && rank == 0)
      pres_avg = pres_avg / vol_tot_reg;
   CONVERGE_mpi_reduce(&pres_peak_1, &pres_peak, 1, CONVERGE_DOUBLE, CONVERGE_MPI_MAX, 0);
   CONVERGE_mpi_reduce(&pres_min_1, &pres_min, 1, CONVERGE_DOUBLE, CONVERGE_MPI_MIN, 0);

   // Writing output file for the whole mesh
   if(rank == 0)
   {
      CONVERGE_outfile_data_row_t row;
      CONVERGE_outfile_data_row_create(outputs_file[user_num_regions], &row);
      CONVERGE_outfile_data_row_add_data(row, 3, pres_avg, pres_peak, pres_min);
      CONVERGE_outfile_add_row(outputs_file[user_num_regions], row);
      CONVERGE_outfile_data_row_destroy(row);
   }

   return;
}

CONVERGE_OUTPUT(outfile_method_close,
                IN(VALUE(CONVERGE_mesh_t, mesh)),
                OUT(VALUE(CONVERGE_outfile_handle_t*, passed_file)))
{
   CONVERGE_logger_debug("Calling user_output close");
   CONVERGE_index_t user_num_regions       = CONVERGE_mesh_num_regions(mesh);
   CONVERGE_outfile_handle_t* outputs_file = passed_file;
   int rank;
   CONVERGE_mpi_comm_rank(&rank);

   if(rank == 0)
   {
      for(CONVERGE_index_t ii = 0; ii < user_num_regions + 1; ii++)
      {
         CONVERGE_outfile_destroy(outputs_file[ii]);
      }
   }
   free(outputs_file);
   return;
}
