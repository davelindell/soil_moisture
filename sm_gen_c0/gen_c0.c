/*
 * sm_gen_c0.c
 *
 *  Created on: May 26, 2015
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

#define NUM_TS_DAYS 2500
#define NUM_THREADS 24
#define NUM_YEARS 6
#define YEAR_START 2009
#define YEAR_END 2014

#define NDIMS 2

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

/* some global mutexes */
pthread_mutex_t fopen_lock;

/* Program documentation. */
static char doc[] =
  "sm_gen_c0.c-- Program to estimate the c0 values from the files\
 generated by the sm_gen_time_series program.\n\
 Region must be a defined type 'NAm','SAm', etc.";

/* A description of the arguments we accept. */
static char args_doc[] = "Region Image_type";

/* The options we understand. */
static struct argp_option options[] = {
  {"verbose",  'v', 0,      0,  "Produce verbose output" },
  {"grd",  'g', 0,      0,  "Generate c0 over grd files" },
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *region;                /* Region */
  char *type;
  int grd;
  int verbose;
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
      if (state->arg_num >= 2) {
        /* Too many arguments. */
        argp_usage (state);
      }
      else if (arguments->region == NULL) {
          arguments->region = arg;
      }
      else if (arguments->type == NULL) {
          arguments->type = arg;
      }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2) {
        /* Not enough arguments. */
        argp_usage (state);
      }
      else if (arguments->region == NULL) {
          argp_failure(state, 1, 0, "ERROR, region not defined!");
      }
      else if (
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
      else if (strcmp(arguments->type,"a") != 0 && strcmp(arguments->type,"b") != 0) {
          argp_failure(state, 1, 0, "ERROR, inputted image type not defined!");
      }
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/* Check to see if a directory exists */
void checkdir(char* dirname) {
/* http://stackoverflow.com/questions/9314586/c-faster-way-to-check-if-a-directory-exists */
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
    float a;
    float b;
} sigma0_value;

int wetcmpfunc (const void * a, const void * b)
{
   float result = ( *(float*)a - *(float*)b );
   if (result > 0)
       return 1;
   else if (result == 0)
       return 0;
   else
       return -1;
}

int drycmpfunc (const void * a, const void * b)
{
   float result = ( (*(sigma0_value*)a).a - (*(sigma0_value*)b).a );
   if (result > 0)
       return 1;
   else if (result == 0)
       return 0;
   else
       return -1;
}

float mean(float *arr, int size) {
    int i;
    float avg = 0;
    for (i = 0; i < size; i++) {
        avg += arr[i];
    }
    avg = avg / size;
    return avg;
}

float dry_mean(sigma0_value *arr, int size, int opt) {
    int i;
    float avg = 0;
    if (opt == 0) {
        for (i = 0; i < size; i++) {
            avg += arr[i].a;
        }
    }
    else {
        for (i = 0; i < size; i++) {
            avg += arr[i].b;
        }
    }
    avg = avg / size;
    return avg;
}

void find_min_max(float ****tseries, float ****tseriesb, int row, int col , float *min, float *max, float *slope, char *type) {
    int i, year, day;
    int cur_ind = 0;
    int num_days = 365*(YEAR_END-YEAR_START+1);
    float cur_data, cur_slope;
    float dry_iqr, wet_iqr, dry_min, dry_max, wet_min, wet_max;
    int q1_loc,q3_loc, dry_start, dry_stop, wet_start, wet_stop;
    int found_dry_start, found_wet_start;

    /* Allocate memory for the filtered, ordered array */
    float *filt_tseries = (float*)malloc(num_days*sizeof(float));
    if (!filt_tseries) {
            fprintf(stderr, "Memory Error!\n");
            exit(-1);
    }
    memset(filt_tseries,0,num_days*sizeof(float));

    /* Allocate memory for the array adjusted to 25 deg. inc angle for theta dry */
    sigma0_value *tseries_dry = (sigma0_value*)malloc(num_days*sizeof(sigma0_value));

    if (!tseries_dry) {
            fprintf(stderr, "Memory Error!\n");
            exit(-1);
    }
    memset(tseries_dry,0,num_days*sizeof(sigma0_value));

    year = -1;
    day = 0;
    for (i = 0; i < num_days; i++) {
        day = i%366;
        if (day == 0) {
            year++;
        }

        cur_data = tseries[row][col][year][day];
        cur_slope = tseriesb[row][col][year][day];

        if (!strcmp(type,"a")) {
            if (cur_data != 0 &&
                    abs(abs(cur_data) - 33) > .01) {
                filt_tseries[cur_ind] = cur_data;

                /* Fill in the 25 deg reference values */
                /* sigma0_25 = A + B * (25 - 40)  where 25 is theta_dry, 40 is current inc angle*/
                tseries_dry[cur_ind].a = cur_data + cur_slope * (-15);
                tseries_dry[cur_ind].b = cur_slope;
                cur_ind++;
            }
        } else {
            if (cur_slope != 0 &&
                    abs(abs(cur_slope) - 3) > .01) {
                filt_tseries[cur_ind] = cur_slope;
                tseries_dry[cur_ind].b = cur_slope;
                cur_ind++;
            }
        }
    }

    if (cur_ind == 0) {
        *min = 0;
        *max = 0;
        return;
    }
    qsort(filt_tseries,cur_ind,sizeof(float),wetcmpfunc);
    qsort(tseries_dry,cur_ind,sizeof(sigma0_value),drycmpfunc);

    /* use filt_tseries to get c0wet, tseries_dry for c0dry */
    /* c0wet from max values, c0dry from min values*/
    /* First find interquartile range -- q1_loc = (N+1)/4 */
    q1_loc = (cur_ind + 1) / 4;
    q3_loc = 3*(cur_ind + 1) / 4;

    wet_iqr = filt_tseries[q3_loc] - filt_tseries[q1_loc];
    dry_iqr = tseries_dry[q3_loc].a - tseries_dry[q1_loc].a;

    /* remove values greater than 3*IQR away from mean */
    /* values must be greater than dry cap and less than wet cap */
    wet_min = mean(filt_tseries, cur_ind) - (3 * wet_iqr);
    wet_max = mean(filt_tseries, cur_ind) + (3 * wet_iqr);

    dry_min = dry_mean(tseries_dry, cur_ind, 0) - (3 * dry_iqr);
    dry_max = dry_mean(tseries_dry, cur_ind, 0) + (3 * dry_iqr);

    dry_start = 0;
    dry_stop = 0;
    wet_start = 0;
    wet_stop = 0;
    found_dry_start = 0;
    found_wet_start = 0;

    for (i = 0; i < cur_ind; i++) {
        if (!found_dry_start && tseries_dry[i].a > dry_min) {
            dry_start = i;
            found_dry_start = 1;
        }
        if (!found_wet_start && filt_tseries[i] > wet_min) {
            wet_start = i;
            found_wet_start = 1;
        }
        if (found_dry_start && tseries_dry[i].a < dry_max) {
            dry_stop = i;
        }
        if (found_wet_start && filt_tseries[i] < wet_max) {
            wet_stop = i;
        }
    }

    // now calculate the mean again and remove any outliers 1.5 IQR away from mean
    // before I used an alternate more ad hoc method, I think this is the actual method
    // described in naeimi2009, but their description is somewhat ambiguous.
    q1_loc = (wet_stop - wet_start + 1) / 4 + wet_start;
    q3_loc = 3*(wet_stop - wet_start + 1) / 4 + wet_start;
    wet_iqr = filt_tseries[q3_loc] - filt_tseries[q1_loc];

    q1_loc = (dry_stop - dry_start + 1) / 4 + dry_start;
    q3_loc = 3*(dry_stop - dry_start + 1) / 4 + dry_start;
    dry_iqr = tseries_dry[q3_loc].a - tseries_dry[q1_loc].a;

    /* remove values greater than 1.5*IQR away from mean */
    /* values must be greater than dry cap and less than wet cap */
    wet_min = mean(filt_tseries+wet_start, wet_stop - wet_start + 1) - (1.5 * wet_iqr);
    wet_max = mean(filt_tseries+wet_start, wet_stop - wet_start + 1) + (1.5 * wet_iqr);

    dry_min = dry_mean(tseries_dry+dry_start, dry_stop - dry_start + 1, 0) - (1.5 * dry_iqr);
    dry_max = dry_mean(tseries_dry+dry_start, dry_stop - dry_start + 1, 0) + (1.5 * dry_iqr);

    dry_start = 0;
    dry_stop = 0;
    wet_start = 0;
    wet_stop = 0;
    found_dry_start = 0;
    found_wet_start = 0;

    for (i = 0; i < cur_ind; i++) {
        if (!found_dry_start && tseries_dry[i].a > dry_min) {
            dry_start = i;
            found_dry_start = 1;
        }
        if (!found_wet_start && filt_tseries[i] > wet_min) {
            wet_start = i;
            found_wet_start = 1;
        }
        if (found_dry_start && tseries_dry[i].a < dry_max) {
            dry_stop = i;
        }
        if (found_wet_start && filt_tseries[i] < wet_max) {
            wet_stop = i;
        }
    }

    /* Now average top/bottom 10% values minus what we skimmed off */
    int num_dry_avg = (cur_ind-dry_start)*.05;
    int num_wet_avg = wet_stop*.05;
    *min = dry_mean(tseries_dry + dry_start, num_dry_avg, 0);
    // There are wet_stop total measurements used, so avg top 10% of those
    *max = mean(filt_tseries + wet_stop - (num_wet_avg-1), num_wet_avg);
    *slope = dry_mean(tseries_dry + dry_start, num_dry_avg,1);


//    This is the old method where I used the separate groups of high/low values
//    dry_size = 50;
//    wet_size = 20;
//
//    q1_loc = (dry_size + 1) / 4;
//    q3_loc = 3*(dry_size + 1) / 4;
//    dry_iqr = tseries_dry[q3_loc + dry_start].a - tseries_dry[q1_loc + dry_start].a;
//
//    q1_loc = (wet_size + 1) / 4;
//    q3_loc = 3*(wet_size + 1) / 4;
//    wet_iqr = filt_tseries[q3_loc + wet_stop - wet_size] - filt_tseries[q1_loc + wet_stop - wet_size];
//
//    /* now take the top and bottom values, get rid of 1.5 irq away from mean, take avg */
//    /* take 50 for dry and 20 for wet */
//    /* now take the average and remove values greater than 1.5 IQR away from mean */
//    dry_min = dry_mean(tseries_dry + dry_start, dry_size, 0) - (1.5 * dry_iqr);
//    wet_max = mean(filt_tseries + wet_stop - (wet_size - 1), wet_size) + (1.5 * wet_iqr);
//
//    dry_stop = dry_start + dry_size - 1;
//    dry_start = 0;
//    wet_start = wet_stop - wet_size + 1;
//    wet_stop = 0;
//    found_dry_start = 0;
//    found_wet_start = 0;
//
//    for (i = 0; i < cur_ind; i++) {
//        if (!found_dry_start && tseries_dry[i].a > dry_min) {
//            dry_start = i;
//            found_dry_start = 1;
//        }
//        if (filt_tseries[i] < wet_max) {
//            wet_stop = i;
//        }
//    }
//
//    /* Now average the 50 and 20 values minus what we skimmed off */
//    *min = dry_mean(tseries_dry + dry_start, dry_stop-dry_start+1, 0);
//    *max = mean(filt_tseries + wet_start, wet_stop-wet_start+1);
//    *slope = dry_mean(tseries_dry, cur_ind,1);

    free(filt_tseries);
    free(tseries_dry);
    return;
}

typedef struct {
    float ****tseries;
    float ****tseriesb;
    float **c0_wet;
    float **c0_dry;
    float **dry_slope;
    int start_i;
    int stop_i;
    int num_columns;
    char *region;
    char *type;
} thread_args;

void *mthreadGenC0(void *arg) {
    thread_args *t_args = (thread_args*)arg;
    float ****tseries = t_args->tseries;
    float ****tseriesb = t_args->tseriesb;
    float **c0_wet = t_args->c0_wet;
    float **c0_dry = t_args->c0_dry;
    float **dry_slope = t_args->dry_slope;
    int start_row = t_args->start_i;
    int stop_row = t_args->stop_i;
    int num_columns = t_args->num_columns;
    char *type = t_args->type;
    int i,j;
    float min,max,slope;

    /* for all pixel files find min/max, store */
    for (i = start_row; i <= stop_row; i++) {
        setvbuf (stdout, NULL, _IONBF, 0);
        printf("Processing Row: %04d\n",i+1);

        for (j = 0; j < num_columns-1; j++) {

            /* find min and max */
            find_min_max(tseries, tseriesb, i,j, &min, &max, &slope, type);

            /* store in 2d array */
            c0_dry[i][j] = min;
            c0_wet[i][j] = max;
            dry_slope[i][j] = slope;

        }
    }
    return NULL;
}


int main (int argc, char **argv)
{
    struct arguments arguments;

    /* Default values. */
    arguments.verbose = 0;
    arguments.grd = 0;
    arguments.region = NULL;
    arguments.type = NULL;
    static const int NUM_DAYS = 365;

    /* Parse our arguments; every option seen by parse_opt will
     be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    /* Set up variables */
    int num_columns;
    int num_rows;
    char* region = arguments.region;
    char* type = arguments.type;
    int grd = arguments.grd;
    int i,j,k;

    /* Initialize NETCDF Variables */
    int ncid, row_dimid, col_dimid;
    int varid, wet_varid, dry_varid, slope_varid;
    int retval;
    char FILE_NAME[100];
    int dimids[NDIMS];

    /* multithread args */
    thread_args t_args[NUM_THREADS];
    pthread_t thread_id[NUM_THREADS];
    int ind_per_thread;
    int start_index;
    int stop_index;

    printf ("GEN_C0\n---------------\nBeginning processing with options:\n");

    printf ("Region = %s\nVERBOSE = %s\nTYPE = %s\n---------------\n",
      arguments.region,
      arguments.verbose ? "yes" : "no",
      arguments.type);

    /* define image areas based on region */
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

    /* allocate memory for 2d arrays */
    printf("Allocating Memory...");
    float **c0_dry = (float**)malloc(sizeof(float *)*num_rows);
    float **c0_wet = (float**)malloc(sizeof(float *)*num_rows);
    float **dry_slope = (float**)malloc(sizeof(float *)*num_rows);
    c0_dry[0] = (float*)malloc(sizeof(float)*num_rows*num_columns);
    c0_wet[0] = (float*)malloc(sizeof(float)*num_rows*num_columns);
    dry_slope[0] = (float*)malloc(sizeof(float)*num_rows*num_columns);
    for (i = 1; i < num_rows; i++) {
        c0_dry[i] = c0_dry[0] + i * num_columns;
        c0_wet[i] = c0_wet[0] + i * num_columns;
        dry_slope[i] = dry_slope[0] + i * num_columns;
    }
    memset(c0_dry[0],0,num_rows*num_columns*sizeof(float));
    memset(c0_wet[0],0,num_rows*num_columns*sizeof(float));
    memset(dry_slope[0],0,num_rows*num_columns*sizeof(float));

    /* allocate memory for NetCDF File */
    float ****row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    float ***column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    float **year_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_YEARS);
    float *day_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****tseries = row_ptr;

    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        tseries[i] = column_ptr;
        for (j = 0; j < num_columns; j++, year_ptr += NUM_YEARS) {
            tseries[i][j] = year_ptr;
            for (k = 0; k < NUM_YEARS; k++, day_ptr += NUM_DAYS) {
                tseries[i][j][k] = day_ptr;
            }
        }
    }

    row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    year_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_YEARS);
    day_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****tseriesb = row_ptr;

    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        tseriesb[i] = column_ptr;
        for (j = 0; j < num_columns; j++, year_ptr += NUM_YEARS) {
            tseriesb[i][j] = year_ptr;
            for (k = 0; k < NUM_YEARS; k++, day_ptr += NUM_DAYS) {
                tseriesb[i][j][k] = day_ptr;
            }
        }
    }

    printf("done\n");

    setvbuf (stdout, NULL, _IONBF, 0);
    printf("Reading NetCDF Files...");
    /* Open the netCDF time series a file*/
    sprintf(FILE_NAME,"/auto/temp/lindell/soilmoisture/ts/ts_%s_%s.nc",region,"a");
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "data", &varid)))
        ERR(retval);

    /* read values from netCDF variable */
    if ((retval = nc_get_var_float(ncid, varid, &tseries[0][0][0][0])))
       ERR(retval);

    /* Open the netCDF time series b file*/
    sprintf(FILE_NAME,"/auto/temp/lindell/soilmoisture/ts/ts_%s_%s.nc",region,"b");
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "data", &varid)))
        ERR(retval);

    /* read values from netCDF variable */
    if ((retval = nc_get_var_float(ncid, varid, &tseriesb[0][0][0][0])))
       ERR(retval);

    printf("done\n");

    /* get threads ready */
    /* split up column processing based on number of threads/rows */
    ind_per_thread = num_rows / NUM_THREADS - 1;
    start_index = 0;
    stop_index = 0;
    for (i = 0; i < NUM_THREADS; i++) {
        if (i == NUM_THREADS - 1) {
            stop_index = num_rows-1;
        } else {
            stop_index = start_index + ind_per_thread;
        }
        t_args[i].tseries = tseries;
        t_args[i].tseriesb = tseriesb;
        t_args[i].c0_wet = c0_wet;
        t_args[i].c0_dry = c0_dry;
        t_args[i].dry_slope = dry_slope;
        t_args[i].start_i = start_index;
        t_args[i].stop_i = stop_index;
        t_args[i].num_columns = num_columns;
        t_args[i].region = region;
        t_args[i].type = type;
        start_index = stop_index + 1;
    }

    /* call multithreaded function */
    /* submit threads */
    printf("Starting Processing\n");
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadGenC0, &t_args[i]);
    }

    /* join threads */
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }

    /* save min/max 2d arrays to netcdf file */
    /* Create the file. */
    sprintf(FILE_NAME,"/auto/temp/lindell/soilmoisture/c0/c0_%s.nc",region);
    if ((retval = nc_create(FILE_NAME, NC_NETCDF4, &ncid)))
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
    if ((retval = nc_def_var(ncid, "dry", NC_FLOAT, NDIMS, dimids, &dry_varid)))
        ERR(retval);

    if ((retval = nc_def_var(ncid, "wet", NC_FLOAT, NDIMS, dimids, &wet_varid)))
        ERR(retval);

    if ((retval = nc_def_var(ncid, "dry_slope", NC_FLOAT, NDIMS, dimids, &slope_varid)))
        ERR(retval);

    /* End define mode. */
    if ((retval = nc_enddef(ncid)))
        ERR(retval);

    /* Write the data. */
    if ((retval = nc_put_var_float(ncid, dry_varid, &c0_dry[0][0])))
        ERR(retval);
    if ((retval = nc_put_var_float(ncid, wet_varid, &c0_wet[0][0])))
        ERR(retval);
    if ((retval = nc_put_var_float(ncid, slope_varid, &dry_slope[0][0])))
        ERR(retval);

    /* Close the file. */
    if ((retval = nc_close(ncid)))
        ERR(retval);

    /* Free memory for 3D image timeseries array */
    printf("Finishing up...");
    free(tseries[0][0][0]);
    free(tseries[0][0]);
    free(tseries[0]);
    free(tseries);

    free(c0_wet[0]);
    free(c0_dry[0]);
    free(dry_slope[0]);
    free(c0_wet);
    free(c0_dry);
    free(dry_slope);

    printf("done\n");

    exit (0);
}
