/*
 * gen_time_series.c
 *
 *  Created on: May 22, 2015
 *      Author: lindell
 */

#include <stdlib.h>
#include <stdio.h>
#include <argp.h>
#include <string.h>
#include <sir_ez.h>
#include <sir3.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <fcntl.h>
#include <errno.h>

#include <pthread.h>

#include <netcdf.h>

/* This is the name of the data file we will read. */
#define NUM_THREADS 24
#define NDIMS 2
#define YEAR_START 2007
#define YEAR_END 2014
#define NUM_YEARS 8
#define NUM_DAYS 183

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e));}


pthread_mutex_t fopen_lock;

/* Program documentation. */
static char doc[] =
  "gen_warp_images.c-- Program to generate images from the warp files.";

/* A description of the arguments we accept. */
static char args_doc[] = "Region";

/* The options we understand. */
static struct argp_option options[] = {
  {"verbose",  'v', 0,      0,  "Produce verbose output" },
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *region;                /* Region */
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
    case ARGP_KEY_ARG:
      if (state->arg_num >= 1) {
        /* Too many arguments. */
        argp_usage (state);
      }
      else if (arguments->region == NULL) {
    	  arguments->region = arg;
      }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1) {
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
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/* Thread arg struct */
typedef struct {
    float ****storage;
    float ****count;
	int start_i;
	int stop_i;
	int num_rows;
	int num_columns;
	char *region;
	char **warp_list;
	sir_head *head;
	pthread_mutex_t *fopen_lock;
} thread_args;

void remove_newline(char *text) {
    int last_char = strlen(text) - 1;
    if (text[last_char] == '\n')
        text[last_char] = '\0';
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

void *mthreadLoadImg(void *arg) {
    thread_args *t_args = (thread_args*)arg;
    float ****storage = t_args->storage;
    float ****count = t_args->count;
    int start_i = t_args->start_i;
    int stop_i = t_args->stop_i;
    int num_rows = t_args->num_rows;
    int num_columns = t_args->num_columns;
    pthread_mutex_t *fopen_lock = t_args->fopen_lock;
    char **warp_list = t_args->warp_list;
    int retval;
    int ncid;
    size_t i,j;
    int file_i;
    sir_head *head = t_args->head;
    // netcdf vars
    int lon_varid, lat_varid, time_varid, sm_varid, rsize_varid;
    int loc_dimid, obs_dimid;
    size_t loc_len,obs_len;
    // netcdf storage
    double *lon, *lat, *time;
    long long *rsize;
    int8_t *sm;
    // time storage
    struct tm *tm;



    for (file_i = start_i; file_i <= stop_i; file_i++) {
        printf("Processing file %d of %d\n",file_i,935);
        pthread_mutex_lock(fopen_lock);

        if ((retval = nc_open(warp_list[file_i], NC_NOWRITE, &ncid))) {
            ERR(retval);
            continue;
        }
        // get the dimids and the lengths of the dimensions
        if ((retval = nc_inq_dimid (ncid, "locations", &loc_dimid))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_dimid (ncid, "obs", &obs_dimid))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_dimlen (ncid, loc_dimid, &loc_len))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_dimlen (ncid, obs_dimid, &obs_len))){
            ERR(retval);
            continue;
        }
        // allocate space for all of the variables based on dim length
        lon = malloc(sizeof(double)*loc_len);
        lat = malloc(sizeof(double)*loc_len);
        time = malloc(sizeof(double)*obs_len);
        rsize = malloc(sizeof(int64_t)*loc_len);
        sm = malloc(sizeof(int8_t)*obs_len);

        // where to store the time conversions
        tm = malloc(sizeof(struct tm)*obs_len);

        memset(lon,0,sizeof(double)*loc_len);
        memset(lat,0,sizeof(double)*loc_len);
        memset(time,0,sizeof(double)*obs_len);
        memset(rsize,0,sizeof(int64_t)*loc_len);
        memset(sm,0,sizeof(int8_t)*obs_len);
        memset(tm,0,sizeof(struct tm)*obs_len);

        /* Get the varid of the data variable, based on its name. */
        if ((retval = nc_inq_varid(ncid, "lon", &lon_varid))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_varid(ncid, "lat", &lat_varid))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_varid(ncid, "time", &time_varid))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_varid(ncid, "row_size", &rsize_varid))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_inq_varid(ncid, "sm", &sm_varid))){
            ERR(retval);
            continue;
        }
        /* read values from netCDF variable */
        if ((retval = nc_get_var_double(ncid, lon_varid, lon))){
            ERR(retval);
            continue;
        }
        /* read values from netCDF variable */
        if ((retval = nc_get_var_double(ncid, lat_varid, lat))){
            ERR(retval);
            continue;
        }
        /* read values from netCDF variable */
        if ((retval = nc_get_var_double(ncid, time_varid, time))){
            ERR(retval);
            continue;
        }
        /* read values from netCDF variable */
        if ((retval = nc_get_var_longlong(ncid, rsize_varid, rsize))){
            ERR(retval);
            continue;
        }
        /* read values from netCDF variable */
        if ((retval = nc_get_var(ncid, sm_varid, sm))){
            ERR(retval);
            continue;
        }
        if ((retval = nc_close(ncid))){
            ERR(retval);
            continue;
        }

        pthread_mutex_unlock(fopen_lock);

        // now load the image into the storage array if necessary
        int use_file = 0;
        float x,y;
        int x_int, y_int;
        for (i = 0; i < loc_len; i++) {
            sir_latlon2pix(lon[i],lat[i],&x,&y,head);
            x_int = (int)x;
            y_int = (int)y;
            if ( x_int < num_columns && x_int >= 0 && y_int < num_rows && y_int >= 0) {
                use_file = 1;
                break;
            }
        }

        if (!use_file) {
            // free data
            free(lon);
            free(lat);
            free(time);
            free(rsize);
            free(sm);
            free(tm);
            continue;
        }

        // if I'm here, then there are measurements that need to be stored
        // first convert the time to doy and year
        // time is stored as a double -- number of days since 1900-1-1 00:00:00
        double cur_time;
        for (i = 0; i < obs_len; i++){
            cur_time = time[i];
            // first convert to seconds
            // (1 day * 24 hours/day * 60 min/hr * 60 sec/min) = 84600
            cur_time = cur_time * 86400;
            // convert from seconds since 1900 to unix time
            // note that accd to time protocol in RFC 868, this is 2208988800
            cur_time = cur_time - 2208988800;

            // here, cast from double to time_t and use gmtime_r to convert to year/doy
            time_t cur_time_conv = (time_t) cur_time;
            gmtime_r(&cur_time_conv, &tm[i]);
            // note that tm_year is years since 1900, so adjust to true year
            tm[i].tm_year = tm[i].tm_year + 1900;
            assert(tm[i].tm_year-(YEAR_START)>=0);
        }
        // now all the times should be correct in the tm struct
        // parse through the measurements

        size_t ind = 0;
        int affected_days[3] = {-1,-1,-1};
        int k;
        for (i = 0; i < loc_len; i++){
            sir_latlon2pix(lon[i],lat[i],&x,&y,head);
            x_int = (int)x;
            y_int = (int)y;
            // if lat/lon are in im bounds
            if ( x_int < num_columns && x_int >= 0 && y_int < num_rows && y_int >= 0) {
                for (j = 0; j < rsize[i]; j++){

                    if (tm[ind].tm_yday%2 == 0){
                        affected_days[0] = tm[ind].tm_yday-3;
                        affected_days[1] = tm[ind].tm_yday-1;
                        affected_days[2] = -1;
                    }
                    else {
                        affected_days[0] = tm[ind].tm_yday-4;
                        affected_days[1] = tm[ind].tm_yday-2;
                        affected_days[2] = tm[ind].tm_yday;
                    }
                    // convert affected days to a 1:2:365 index

                    for (k = 0; k < 3; k++) {
                        if (affected_days[k] < 0)
                            continue;
                        affected_days[k] = affected_days[k]/2;
                        // shorthand
                        float *st = &storage[y_int][x_int][affected_days[k]][tm[ind].tm_year-YEAR_START];
                        float *ct = &count[y_int][x_int][affected_days[k]][tm[ind].tm_year-YEAR_START];
                        // do running average
                        *st = *st + (((float)sm[ind]) - *st) / *ct;
                        (*ct)++;
                    }
                    ind++;
                }
            }
            else {
                ind += rsize[i];
            }
        }
        // free data
        free(lon);
        free(lat);
        free(time);
        free(rsize);
        free(sm);
        free(tm);
    }
    return NULL;
}

int main (int argc, char **argv)
{
    struct arguments arguments;

    // Set up variables
    int num_columns;
    int num_rows;

    // multithreading params
    thread_args t_args[NUM_THREADS];
    pthread_t thread_id[NUM_THREADS];
    int ind_per_thread;
    int start_index;
    int stop_index;
    size_t i,j,k,m;

    // Initialize NETCDF Variables
    int ncid, row_dimid, col_dimid;
    int varid;
    int retval;
    char FILE_NAME[100];
    int dimids[NDIMS];

    /* Default values. */
    arguments.verbose = 0;
    arguments.region = NULL;

    /* Parse our arguments; every option seen by parse_opt will
     be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    char* region = arguments.region;

    printf ("GEN_WARP_IMAGES\n---------------\nBeginning processing with options:\n");

    printf ("Region = %s\nVERBOSE = %s\n---------------\n",
          arguments.region,
          arguments.verbose ? "yes" : "no");

    // define image areas based on region
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

    // open example sir file
    sir_head head;
    char sir_fname[150];
    sprintf(sir_fname,"/home/lindell/workspace/soil_moisture/sm_gen_warp_images/sir/%s.sir",region);
    FILE *sir_fid = fopen(sir_fname,"r");
    if (sir_fid == NULL) {
        fprintf(stderr,"*** could not open list file %s\n",sir_fname);
        exit(-1);
    }
    get_sir_head_file(sir_fid,&head);
    head.ascale = head.ascale/2.809;
    head.bscale = head.bscale/2.809;

//    printf("NUM_ROWS: %d\nNUM_COL: %d\nNUM_DAYS: %d\nNUM_YEARS: %d\nREG:%s\n",num_rows,num_columns,NUM_DAYS,NUM_YEARS,region);
//    exit(-1);
    // Allocate memory for 4D image timeseries array
    setvbuf (stdout, NULL, _IONBF, 0);
    printf("Allocating Memory...");
    float ****row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    float ***column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    float **day_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_DAYS);
    float *year_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****storage = row_ptr;

    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        storage[i] = column_ptr;
        for (j = 0; j < num_columns; j++, day_ptr += NUM_DAYS) {
            storage[i][j] = day_ptr;
            for (k = 0; k < NUM_DAYS; k++, year_ptr += NUM_YEARS) {
                storage[i][j][k] = year_ptr;
            }
        }
    }
    memset(storage[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));

    row_ptr = (float****)malloc(sizeof(float ***)*num_rows);
    column_ptr = (float***)malloc(sizeof(float **)*num_rows * num_columns);
    day_ptr = (float**)malloc(sizeof(float *)*num_rows * num_columns*NUM_DAYS);
    year_ptr = (float*)malloc(sizeof(float)*num_rows*num_columns*NUM_YEARS*NUM_DAYS);
    float ****count = row_ptr;
    for (i = 0; i < num_rows; i++, column_ptr += num_columns) {
        count[i] = column_ptr;
        for (j = 0; j < num_columns; j++, day_ptr += NUM_DAYS) {
            count[i][j] = day_ptr;
            for (k = 0; k < NUM_DAYS; k++, year_ptr += NUM_YEARS) {
                count[i][j][k] = year_ptr;
            }
        }
    }
//    memset(count[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));
    for (i = 0; i < num_rows; i++){
        for (j = 0; j < num_columns; j++) {
            for (k = 0; k < NUM_DAYS; k++) {
                for (m = 0; m < NUM_YEARS; m++) {
                    count[i][j][k][m] = 1;
                }
            }
        }
    }
    printf("Done\n");

    // get list of files to be opened
    char warp_fname[] = "/home/lindell/workspace/soil_moisture/sm_gen_warp_images/warp.list";
    FILE* file_id = fopen(warp_fname,"r");
    if (file_id == NULL) {
        fprintf(stderr,"*** could not open list file %s\n",warp_fname);
        exit(-1);
    }

    // count lines in file
    char ch;
    size_t list_len = 0;
    while(!feof(file_id)){
        ch = fgetc(file_id);
        if (ch == '\n')
            ++list_len;
    }

    // reset file pointer and put filenames into array buffer
    rewind(file_id);
    int fname_i = 0;
    char fname[150];

    char **warp_list;
    warp_list = malloc(list_len * sizeof(char*));
    for (i = 0; i < list_len; i++)
        warp_list[i] = malloc((150) * sizeof(char));

    while (fgets(fname,sizeof(fname),file_id)!=NULL) {
        remove_newline(fname);
        strcpy(warp_list[fname_i],fname);
        ++fname_i;
    }

    printf("Preparing for processing...\n");

    // sort through pixels/images
    // split up row processing based on number of threads/files
    ind_per_thread = (list_len-1) / NUM_THREADS;

    if ((ind_per_thread % 2) == 0) {
        ind_per_thread--;
    }

    start_index = 0;
    stop_index = 0;
    for (i = 0; i < NUM_THREADS; i++) {
        if (i == NUM_THREADS - 1) {
            stop_index = list_len-1;
        } else {
            stop_index = start_index + ind_per_thread;
        }
        t_args[i].storage = storage;
        t_args[i].count = count;
        t_args[i].start_i = start_index;
        t_args[i].stop_i = stop_index;
        t_args[i].num_rows = num_rows;
        t_args[i].num_columns = num_columns;
        t_args[i].region = region;
        t_args[i].fopen_lock = &fopen_lock;
        t_args[i].warp_list = warp_list;
        t_args[i].head = &head;
        start_index = stop_index + 1;
    }

    // submit threads
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadLoadImg, &t_args[i]);
    }

    // join threads
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }

    // save storage array to netcdf file
    printf("Done processing, preparing to save files\n");

    // malloc a temp array;
    float **tmp_arr = malloc(sizeof(float**)*num_rows);
    tmp_arr[0] = malloc(sizeof(float)*num_rows*num_columns);
    for (i = 0; i < num_rows; i++) {
        tmp_arr[i] = tmp_arr[0] + i*num_columns;
    }
    memset(tmp_arr[0],0,num_rows*num_columns*sizeof(float));

    int year_i;
    int day_i;
    for (year_i = 0; year_i < NUM_YEARS; year_i++) {
        for (day_i = 0; day_i < NUM_DAYS; day_i++){
            printf("Saving data for %03d %d...\n",day_i*2+1,year_i+2007);

            // copy data from storage into temp array
            for (i = 0; i < num_rows; i++) {
                for (j = 0; j < num_columns; j++) {
                    tmp_arr[i][j] = storage[i][j][day_i][year_i];
                    if (count[i][j][day_i][year_i] == 1) {
                        tmp_arr[i][j] = -1; // set nodata flag
                    }
                }
            }

            sprintf(FILE_NAME,"/auto/temp/lindell/soilmoisture/warp/%s_%04d_%03d.nc",region,year_i+2007,day_i*2+1);

            if ((retval = nc_create(FILE_NAME, NC_NETCDF4, &ncid))){
                ERR(retval);
                continue;
            }
            /* Define the dimensions. */
            if ((retval = nc_def_dim(ncid, "row", num_rows, &row_dimid))) {
                ERR(retval);
                continue;
            }
            if ((retval = nc_def_dim(ncid, "column", num_columns, &col_dimid))){
                ERR(retval);
                continue;
            }

            /* Define the netCDF variables. The dimids array is used to pass
               the dimids of the dimensions of the variables.*/
            dimids[0] = row_dimid;
            dimids[1] = col_dimid;

            /* define the variable */
            if ((retval = nc_def_var(ncid, "sm", NC_FLOAT, NDIMS, dimids, &varid)))
                ERR(retval);

            /* End define mode. */
            if ((retval = nc_enddef(ncid)))
                ERR(retval);

            /* Write the data. */
            if ((retval = nc_put_var_float(ncid, varid, &tmp_arr[0][0])))
                ERR(retval);

            /* Close the file. */
            if ((retval = nc_close(ncid)))
                ERR(retval);

            // clear the tmp array
            memset(tmp_arr[0],0,num_rows*num_columns*sizeof(float));
        }
    }

    free(tmp_arr[0]);
    free(tmp_arr);

    // Free memory for 3D image timeseries array
    printf("Finishing up...");
    free(storage[0][0][0]);
    free(storage[0][0]);
    free(storage[0]);
    free(storage);

    free(count[0][0][0]);
    free(count[0][0]);
    free(count[0]);
    free(count);

    for (i = 0; i < list_len; i++)
        free((void*)warp_list[i]);
    free((void*)warp_list);

    printf("done\n");

    exit (0);
}
