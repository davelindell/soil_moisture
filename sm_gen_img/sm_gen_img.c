/*
 * sm_gen_img.c
 *
 *  Created on: Jun 2, 2015
 *      Author: lindell
 */

#include <stdlib.h>
#include <argp.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <pthread.h>

#include <netcdf.h>

#define NUM_DAYS 365
#define YEAR_START 2009

#define NUM_THREADS 24
#define NDIMS 2
#define YEAR_START 2009
#define YEAR_END 2014
#define NUM_YEARS 6

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

// some global mutexes
pthread_mutex_t fopen_lock;
pthread_mutex_t fclose_lock;
pthread_mutex_t netcdfop_lock;

/* Program documentation. */
static char doc[] =
  "sm_gen_img.c-- Program compile images from the SWI/ms files\n\
 Region must be a defined type 'NAm','SAm', etc.";

/* A description of the arguments we accept. */
static char args_doc[] = "Region";

/* The options we understand. */
static struct argp_option options[] = {
  {"verbose",  'v', 0,      0,  "Produce verbose output" },
  {"grd",  'g', 0,      0,  "Generate grd images" },
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *region;                /* Region */
  int verbose;
  int grd;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'v':
      arguments->verbose = 1;
      break;
    case 'g':
      arguments->grd = 1;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 1) {
        /* Too many arguments. */
        argp_usage (state);
      } else if (arguments->region == NULL) {
          arguments->region = arg;
      }


      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1) {
        /* Not enough arguments. */
        argp_usage (state);
      } else if (arguments->region == NULL) {
          argp_failure(state, 1, 0, "ERROR, region not defined!");
      } else if (
      strcmp(arguments->region,"Ama") != 0 &&
      strcmp(arguments->region,"Aus") != 0 &&
      strcmp(arguments->region,"Ber") != 0 &&
      strcmp(arguments->region,"CAm") != 0 &&
      strcmp(arguments->region,"ChJ") != 0 &&
      strcmp(arguments->region,"Eur") != 0 &&
      strcmp(arguments->region,"Ind") != 0 &&
      strcmp(arguments->region,"NAf") != 0 &&
      strcmp(arguments->region,"NAm") != 0 &&
      strcmp(arguments->region,"SAf") != 0 &&
      strcmp(arguments->region,"SAm") != 0 &&
      strcmp(arguments->region,"SAs") != 0) {
          argp_failure(state, 1, 0, "ERROR, inputted region not defined!");
      }
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

void free_array(float **array, int num_rows, int num_columns) {
    int i;
    for (i = 0; i < num_rows; i++) {
        free(array[i]);
    }
    free(array);
    return;
}

/* Check to see if a directory exists */
void checkdir(char* dirname) {
// http://stackoverflow.com/questions/9314586/c-faster-way-to-check-if-a-directory-exists
    struct stat s;
    int err = stat(dirname, &s);
    if(-1 == err) {
        if(ENOENT == errno) {
            /* does not exist-- creating dir */
            mkdir(dirname, 0700);
        } else {
            perror("stat error");
            exit(1);
        }
    } else {
        if(S_ISDIR(s.st_mode)) {
            /* it's a dir */
        } else {
            /* exists but is no dir */
            perror("Error creating directory (already exists as file?)");
        }
    }
    return;
}

typedef struct {
    float ****swi_ts;
    float ****ms_ts;
    float ****dry_ts;
    float ****tmp;
    int start_i;
    int stop_i;
    int num_columns;
    char *region;
    pthread_mutex_t *fopen_lock;
    pthread_mutex_t *fclose_lock;
    pthread_mutex_t *netcdfop_lock;
} thread_args;

void *mthreadGetImgData(void *arg) {
    thread_args *t_args = (thread_args*)arg;
    float ****swi_ts = t_args->swi_ts;
    float ****ms_ts = t_args->ms_ts;
    float ****dry_ts = t_args->dry_ts;
    int start_row = t_args->start_i;
    int stop_row = t_args->stop_i;
    char *region = t_args->region;
    pthread_mutex_t *fopen_lock = t_args->fopen_lock;
    char fname[100];
    int i;
    int retval, ncid1, swi_varid, ms_varid, dry_varid;

    // for all pixel files, grab the value and store in the image arrays
    for (i = start_row; i <= stop_row; i++) {
        setvbuf (stdout, NULL, _IONBF, 0);
        printf("Processing Row: %04d\n",i+1);

        /* Open netcdf file */
        pthread_mutex_lock(fopen_lock);
        sprintf(fname,"/auto/temp/lindell/soilmoisture/swi/swi_%s_%04d.nc",region,i+1);

        if ((retval = nc_open(fname, NC_NOWRITE, &ncid1)))
            ERR(retval);

        /* Get the varid of the data variable, based on its name. */
        if ((retval = nc_inq_varid(ncid1, "swi", &swi_varid)))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid1, "dry", &dry_varid)))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid1, "ms", &ms_varid)))
            ERR(retval);

        /* read values from netCDF variable */
        if ((retval = nc_get_var_float(ncid1, swi_varid, &swi_ts[i][0][0][0])))
            ERR(retval);

        /* read values from netCDF variable */
        if ((retval = nc_get_var_float(ncid1, ms_varid, &ms_ts[i][0][0][0])))
            ERR(retval);

        /* read values from netCDF variable */
        if ((retval = nc_get_var_float(ncid1, dry_varid, &dry_ts[i][0][0][0])))
            ERR(retval);

        if ((retval = nc_close(ncid1)))
            ERR(retval);

        pthread_mutex_unlock(fopen_lock);
    }

    return NULL;
}

void *mthreadSwiRearrArray(void *arg) {
    thread_args *t_args = (thread_args*)arg;
    float ****swi_ts = t_args->swi_ts;
    float ****tmp = t_args->tmp;
    int start_row = t_args->start_i;
    int stop_row = t_args->stop_i;
    int num_columns = t_args->num_columns;
    int i,j,k,m;

    for (i = start_row; i <= stop_row; i++) {
        for (j = 0; j < num_columns; j++) {
            for (k = 0; k < NUM_YEARS; k++) {
                for (m = 0; m < NUM_DAYS; m++) {
                    tmp[k][m][i][j] = swi_ts[i][j][k][m];
                }
            }
        }
    }
    return NULL;
}

void *mthreadMsRearrArray(void *arg) {
    thread_args *t_args = (thread_args*)arg;
    float ****ms_ts = t_args->ms_ts;
    float ****tmp = t_args->tmp;
    int start_row = t_args->start_i;
    int stop_row = t_args->stop_i;
    int num_columns = t_args->num_columns;
    int i,j,k,m;

    for (i = start_row; i <= stop_row; i++) {
        for (j = 0; j < num_columns; j++) {
            for (k = 0; k < NUM_YEARS; k++) {
                for (m = 0; m < NUM_DAYS; m++) {
                    tmp[k][m][i][j] = ms_ts[i][j][k][m];
                }
            }
        }
    }
    return NULL;
}

void *mthreadRearrDryArray(void *arg) {
    thread_args *t_args = (thread_args*)arg;
    float ****dry_ts = t_args->dry_ts;
    float ****tmp = t_args->tmp;
    int start_row = t_args->start_i;
    int stop_row = t_args->stop_i;
    int num_columns = t_args->num_columns;
    int i,j,k,m;

    for (i = start_row; i <= stop_row; i++) {
        for (j = 0; j < num_columns; j++) {
            for (k = 0; k < NUM_YEARS; k++) {
                for (m = 0; m < NUM_DAYS; m++) {
                    tmp[k][m][i][j] = dry_ts[i][j][k][m];
                }
            }
        }
    }
    return NULL;
}

int main (int argc, char **argv)
{
    struct arguments arguments;

    /* Default values. */
    arguments.verbose = 0;
    arguments.region = NULL;
    arguments.grd = 0;

    /* Parse our arguments; every option seen by parse_opt will
     be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    /* Set up variables */
    int num_columns;
    int num_rows;
    int grd = arguments.grd;
    char* region = arguments.region;
    char swi_fname[100];
    int i,j,k;

    // Initialize NETCDF Variables
    int ncid, row_dimid, col_dimid;
    int swi_varid, ms_varid, dry_varid;
    int retval;
    int dimids[NDIMS];

    // multithread args
    thread_args t_args[NUM_THREADS];
    pthread_t thread_id[NUM_THREADS];
    int ind_per_thread;
    int start_index;
    int stop_index;

    printf ("GEN_C0\n---------------\nBeginning processing with options:\n");

    printf ("Region = %s\nVERBOSE = %s\nGRD = %s\n---------------\n",
      arguments.region,
      arguments.verbose ? "yes" : "no",
      arguments.grd ? "yes" : "no");

    // define image areas based on region
    if (!grd) {
        if (strcmp(region,"Ama") == 0) {
            num_columns = 1128;
            num_rows = 744;
        } else if (strcmp(region,"Ber") == 0) {
            num_columns = 1350;
            num_rows = 750;
        } else if (strcmp(region,"CAm") == 0) {
              num_columns = 1440;
              num_rows = 700;
        } else if (strcmp(region,"ChJ") == 0) {
              num_columns = 1980;
              num_rows = 950;
        } else if (strcmp(region,"Eur") == 0) {
              num_columns = 1530;
              num_rows = 1040;
        } else if (strcmp(region,"Ind") == 0) {
              num_columns = 1800;
              num_rows = 680;
        } else if (strcmp(region,"NAf") == 0) {
              num_columns = 2120;
              num_rows = 1130;
        } else if (strcmp(region,"NAm") == 0) {
              num_columns = 1890;
              num_rows = 1150;
        } else if (strcmp(region,"SAf") == 0) {
              num_columns = 1220;
              num_rows = 1260;
        } else if (strcmp(region,"SAm") == 0) {
              num_columns = 1310;
              num_rows = 1850;
        }  else if (strcmp(region,"SAs") == 0) {
              num_columns = 1760;
              num_rows = 720;
        } else {
            printf("ERROR SETTING REGION SIZES!");
            exit(-1);
        }
    } else {
        if (strcmp(region,"NAm") == 0) {
            num_columns = 672;
            num_rows = 410;
        } else {
            printf("ERROR SETTING REGION SIZES!");
            exit(-1);
        }
    }

    // Allocate memory for 4D image timeseries array
    setvbuf (stdout, NULL, _IONBF, 0);
    printf("Allocating Memory...");

    float ****row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    float ***column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    float **year_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_YEARS);
    float *day_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****swi_ts = row_ptr;

    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        swi_ts[i] = column_ptr;
        for (j = 0; j < num_columns; j++, year_ptr += NUM_YEARS) {
            swi_ts[i][j] = year_ptr;
            for (k = 0; k < NUM_YEARS; k++, day_ptr += NUM_DAYS) {
                swi_ts[i][j][k] = day_ptr;
            }
        }
    }
    memset(swi_ts[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));

    row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    year_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_YEARS);
    day_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****ms_ts = row_ptr;

    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        ms_ts[i] = column_ptr;
        for (j = 0; j < num_columns; j++, year_ptr += NUM_YEARS) {
            ms_ts[i][j] = year_ptr;
            for (k = 0; k < NUM_YEARS; k++, day_ptr += NUM_DAYS) {
                ms_ts[i][j][k] = day_ptr;
            }
        }
    }
    memset(ms_ts[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));

    row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    year_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_YEARS);
    day_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****dry_ts = row_ptr;

    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        dry_ts[i] = column_ptr;
        for (j = 0; j < num_columns; j++, year_ptr += NUM_YEARS) {
            dry_ts[i][j] = year_ptr;
            for (k = 0; k < NUM_YEARS; k++, day_ptr += NUM_DAYS) {
                dry_ts[i][j][k] = day_ptr;
            }
        }
    }
    memset(dry_ts[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));

    float ****year_ptr_r = (float****)malloc(sizeof(float ***)*NUM_YEARS);
    float ***day_ptr_r = (float***)malloc(sizeof(float **)*NUM_YEARS*NUM_DAYS);
    float **row_ptr_r = (float**)malloc(sizeof(float *)*NUM_DAYS * NUM_YEARS*num_rows);
    float *column_ptr_r = (float*)malloc(sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);
    float ****tmp = year_ptr_r;

    for (i = 0; i < NUM_YEARS; i++, day_ptr_r += NUM_DAYS) {
        tmp[i] = day_ptr_r;
        for (j = 0; j < NUM_DAYS; j++, row_ptr_r += num_rows) {
            tmp[i][j] = row_ptr_r;
            for (k = 0; k < num_rows; k++, column_ptr_r += num_columns) {
                tmp[i][j][k] = column_ptr_r;
            }
        }
    }
    memset(tmp[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));



    printf("Done.\n");

    // get threads ready
    // split up column processing based on number of threads/rows
    ind_per_thread = num_rows / NUM_THREADS -1;
    start_index = 0;
    stop_index = 0;
    for (i = 0; i < NUM_THREADS; i++) {
        if (i == NUM_THREADS - 1) {
            stop_index = num_rows-1;
        } else {
            stop_index = start_index + ind_per_thread;
        }
        t_args[i].swi_ts = swi_ts;
        t_args[i].ms_ts = ms_ts;
        t_args[i].dry_ts = dry_ts;
        t_args[i].tmp = tmp;
        t_args[i].start_i = start_index;
        t_args[i].stop_i = stop_index;
        t_args[i].num_columns = num_columns;
        t_args[i].region = region;
        t_args[i].fopen_lock = &fopen_lock;
        t_args[i].fclose_lock = &fclose_lock;
        t_args[i].netcdfop_lock = &netcdfop_lock;
        start_index = stop_index + 1;
    }

    // call multithreaded function
    // submit threads
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadGetImgData, &t_args[i]);
    }

    // join threads
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }


    printf("Rearranging SWI...");
    // call multithreaded function
    // submit threads

    // ***** REARRANGE SWI *****
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadSwiRearrArray, &t_args[i]);
    }

    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }

    // Free old array
    free(swi_ts[0][0][0]);
    free(swi_ts[0][0]);
    free(swi_ts[0]);
    free(swi_ts);

    // Allocate new array
    year_ptr_r = (float****)malloc(sizeof(float ***)*NUM_YEARS);
    day_ptr_r = (float***)malloc(sizeof(float **)*NUM_YEARS*NUM_DAYS);
    row_ptr_r = (float**)malloc(sizeof(float *)*NUM_DAYS * NUM_YEARS*num_rows);
    column_ptr_r = (float*)malloc(sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);
    float ****swi_ts_rearr = year_ptr_r;

    for (i = 0; i < NUM_YEARS; i++, day_ptr_r += NUM_DAYS) {
        swi_ts_rearr[i] = day_ptr_r;
        for (j = 0; j < NUM_DAYS; j++, row_ptr_r += num_rows) {
            swi_ts_rearr[i][j] = row_ptr_r;
            for (k = 0; k < num_rows; k++, column_ptr_r += num_columns) {
                swi_ts_rearr[i][j][k] = column_ptr_r;
            }
        }
    }

    memcpy((void *)&swi_ts_rearr[0][0][0][0], (void *)&tmp[0][0][0][0], sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);

    printf("Done.\n");

    // ***** REARRANGE MS *****
    printf("Rearranging ms...");
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadMsRearrArray, &t_args[i]);
    }

    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }

    // Free old array
    free(ms_ts[0][0][0]);
    free(ms_ts[0][0]);
    free(ms_ts[0]);
    free(ms_ts);

    // Allocate new array
    year_ptr_r = (float****)malloc(sizeof(float ***)*NUM_YEARS);
    day_ptr_r = (float***)malloc(sizeof(float **)*NUM_YEARS*NUM_DAYS);
    row_ptr_r = (float**)malloc(sizeof(float *)*NUM_DAYS * NUM_YEARS*num_rows);
    column_ptr_r = (float*)malloc(sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);
    float ****ms_ts_rearr = year_ptr_r;

    for (i = 0; i < NUM_YEARS; i++, day_ptr_r += NUM_DAYS) {
       ms_ts_rearr[i] = day_ptr_r;
       for (j = 0; j < NUM_DAYS; j++, row_ptr_r += num_rows) {
           ms_ts_rearr[i][j] = row_ptr_r;
           for (k = 0; k < num_rows; k++, column_ptr_r += num_columns) {
               ms_ts_rearr[i][j][k] = column_ptr_r;
           }
       }
    }

    memcpy((void *)&ms_ts_rearr[0][0][0][0], (void *)&tmp[0][0][0][0], sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);
    printf("Done.\n");

    // ***** REARRANGE DRY *****
    printf("Rearranging sigma0-dry...");
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadRearrDryArray, &t_args[i]);
    }

    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }

    // Free old array
    free(dry_ts[0][0][0]);
    free(dry_ts[0][0]);
    free(dry_ts[0]);
    free(dry_ts);

    // Allocate new array
    year_ptr_r = (float****)malloc(sizeof(float ***)*NUM_YEARS);
    day_ptr_r = (float***)malloc(sizeof(float **)*NUM_YEARS*NUM_DAYS);
    row_ptr_r = (float**)malloc(sizeof(float *)*NUM_DAYS * NUM_YEARS*num_rows);
    column_ptr_r = (float*)malloc(sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);
    float ****dry_ts_rearr = year_ptr_r;

    for (i = 0; i < NUM_YEARS; i++, day_ptr_r += NUM_DAYS) {
        dry_ts_rearr[i] = day_ptr_r;
        for (j = 0; j < NUM_DAYS; j++, row_ptr_r += num_rows) {
            dry_ts_rearr[i][j] = row_ptr_r;
            for (k = 0; k < num_rows; k++, column_ptr_r += num_columns) {
                dry_ts_rearr[i][j][k] = column_ptr_r;
            }
        }
    }

    memcpy((void *)&dry_ts_rearr[0][0][0][0], (void *)&tmp[0][0][0][0], sizeof(float)*NUM_DAYS*NUM_YEARS*num_rows*num_columns);

    // Free temporary holding array
    free(tmp[0][0][0]);
    free(tmp[0][0]);
    free(tmp[0]);
    free(tmp);

    printf("Done.\n");

    // write netcdf files
    printf("Writing NetCDF Files!\n");
    for (i = 0; i < NUM_YEARS; i++) {
        for (j = 0; j < NUM_DAYS; j+=2) {
            printf("    Day: %03d Year: %04d\n",j+1,i+YEAR_START);
            sprintf(swi_fname,"/auto/temp/lindell/soilmoisture/swi/combined/swi_%s_%04d_%03d.nc",region,i+YEAR_START,j+1);
            if ((retval = nc_create(swi_fname, NC_NETCDF4, &ncid)))
                ERR(retval);

            /* Define the dimensions. */
            if ((retval = nc_def_dim(ncid, "row", num_rows, &row_dimid)))
                ERR(retval);
            if ((retval = nc_def_dim(ncid, "column", num_columns, &col_dimid)))
                ERR(retval);

            /* Define the netCDF variables. The dimids array is used to pass
                the dimids of the dimensions of the variables.*/
            dimids[0] = row_dimid;
            dimids[1] = col_dimid;

            /* define the variable */
            if ((retval = nc_def_var(ncid, "swi", NC_FLOAT, NDIMS, dimids, &swi_varid)))
                ERR(retval);

            if ((retval = nc_def_var(ncid, "ms", NC_FLOAT, NDIMS, dimids, &ms_varid)))
                ERR(retval);

            if ((retval = nc_def_var(ncid, "dry", NC_FLOAT, NDIMS, dimids, &dry_varid)))
                ERR(retval);

            /* End define mode. */
            if ((retval = nc_enddef(ncid)))
                ERR(retval);

            /* Write the data. */
            if ((retval = nc_put_var_float(ncid, swi_varid, &swi_ts_rearr[i][j][0][0])))
                ERR(retval);

            if ((retval = nc_put_var_float(ncid, ms_varid, &ms_ts_rearr[i][j][0][0])))
                ERR(retval);

            if ((retval = nc_put_var_float(ncid, dry_varid, &dry_ts_rearr[i][j][0][0])))
                ERR(retval);

            /* Close the file. */
            if ((retval = nc_close(ncid)))
                ERR(retval);
        }
    }

    printf("Freeing memory.\n");

    free(swi_ts_rearr[0][0][0]);
    free(swi_ts_rearr[0][0]);
    free(swi_ts_rearr[0]);
    free(swi_ts_rearr);

    free(ms_ts_rearr[0][0][0]);
    free(ms_ts_rearr[0][0]);
    free(ms_ts_rearr[0]);
    free(ms_ts_rearr);

    free(dry_ts_rearr[0][0][0]);
    free(dry_ts_rearr[0][0]);
    free(dry_ts_rearr[0]);
    free(dry_ts_rearr);

    printf("Finished.\n");
    exit (0);
}
