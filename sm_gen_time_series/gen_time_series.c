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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <fcntl.h>
#include <errno.h>

#include <pthread.h>

#include <netcdf.h>

/* This is the name of the data file we will read. */
#define NUM_THREADS 24
#define NDIMS 4
#define YEAR_START 2009
#define YEAR_END 2014
#define NUM_YEARS 6

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}


pthread_mutex_t fopen_lock;

/* Program documentation. */
static char doc[] =
  "gen_time_series.c-- Program to parse sir file pixels into time\
 series, stored in binary files in a predetermined directory tree.\n\
 Region must be a defined type 'NAm','SAm', etc.\n\
 Image_Type must be 'a' or 'b'";

/* A description of the arguments we accept. */
static char args_doc[] = "Region Image_Type";

/* The options we understand. */
static struct argp_option options[] = {
  {"verbose",  'v', 0,      0,  "Produce verbose output" },
  {"grd",  'g', 0,      0,  "Parse grd files" },
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
      else if (strcmp(arguments->type,"a") != 0 && strcmp(arguments->type,"b") != 0)
    	  argp_failure(state, 1, 0, "ERROR, inputted image type not defined!");
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
    float ****tseries;
	int start_i;
	int stop_i;
	int num_rows;
	int num_columns;
	char *region;
	int grd;
	char *type;
	int NUM_DAYS;
	pthread_mutex_t *fopen_lock;
	pthread_mutex_t *store_lock;
	float ****thread_buffer;
} thread_args;

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

// Copy file utility
// http://stackoverflow.com/questions/2180079/how-can-i-copy-a-file-on-unix-using-c
int cp(const char *to, const char *from)
{
    int fd_to, fd_from;
    char buf[4096];
    ssize_t nread;
    int saved_errno;

    fd_from = open(from, O_RDONLY);
    if (fd_from < 0)
        return -1;

    fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd_to < 0)
        goto out_error;

    while (nread = read(fd_from, buf, sizeof buf), nread > 0)
    {
        char *out_ptr = buf;
        ssize_t nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if (nwritten >= 0)
            {
                nread -= nwritten;
                out_ptr += nwritten;
            }
            else if (errno != EINTR)
            {
                goto out_error;
            }
        } while (nread > 0);
    }

    if (nread == 0)
    {
        if (close(fd_to) < 0)
        {
            fd_to = -1;
            goto out_error;
        }
        close(fd_from);

        /* Success! */
        return 0;
    }

  out_error:
    saved_errno = errno;

    close(fd_from);
    if (fd_to >= 0)
        close(fd_to);

    errno = saved_errno;
    return -1;
}

void *mthreadParseImg(void *arg) {
	thread_args *t_args = (thread_args*)arg;
	float ****tseries = t_args->tseries;
	int start_day = t_args->start_i;
	int stop_day = t_args->stop_i;
	int num_rows = t_args->num_rows;
	int num_columns = t_args->num_columns;
	int grd = t_args->grd;
	char *region = t_args->region;
	char *type = t_args->type;
	pthread_mutex_t *fopen_lock = t_args->fopen_lock;
	int row, column;
	int year;
	int day;
	char ascat_path[100];
	char ascat_internet_path[100];
	char tmpdir[100];
	char cmd[100];
	float pix_val;
	sir_head head;

      for (year = YEAR_START; year <= YEAR_END; ++year) {
          for (day = start_day; day <= stop_day; day+=2) {
              setvbuf (stdout, NULL, _IONBF, 0);
              printf("    Day: %03d of %04d\n", day, year);

              if (strcmp("a",type) == 0) {
                  if (!grd) { // here we do sir files
                      sprintf(ascat_path,
                      "/auto/temp/lindell/soilmoisture/msfa-%s-%s%02d-%03d-%03d.sir",
                      type,region,year-2000,day,day+4);

                      // check if sir file exists in temp dir
                      if( access( ascat_path, F_OK ) == -1 ) {
                          // not unzipped in temp dir, copy over unzipped version and use that
                          // check if exists on internet
                          sprintf(ascat_internet_path,
                          "/auto/internet/ftp/data/ascat/%d/sir/msfa/%s/%03d/%s/msfa-%s-%s%02d-%03d-%03d.sir.gz",
                          year,region,day,type,type,region,year-2000,day,day+4);
                          if( access( ascat_internet_path, F_OK ) != -1 ) {
                              sprintf(tmpdir,"%s.gz",ascat_path);
                              cp(tmpdir,ascat_internet_path);
                              sprintf(cmd, "gunzip %s.gz",ascat_path);
                              system(cmd);
                          } else { //then we don't have the file
                              printf("Error reading file for day %03d, year %d\n",day,year);
                              continue;
                          }
                      }
                  }
                  else { // grab the grd files
                    sprintf(ascat_path,
                        "/auto/temp/lindell/soilmoisture/grd/%04d/%03d-%03d-%04d/msfa-%s-%s%02d-%03d-%03d.grd",
                        year,day,day+4,year,type,region,year-2000,day,day+4);
                    if( access( ascat_path, F_OK ) == -1 ) { // we don't have custom grd files anywhere else
                        printf("Error reading file for day %03d, year %d:\n%s",day,year,ascat_path);
                        continue;
                    }
                  }
              } else { // type is 'b'
                  int day1 = day-14;
                  int day2 = day+15;
                  int year_tmp = year;
                  if (day2 > 365) {
                      day2 = (day2%366)+1;
                  }
                  // year boundary corrections
                  if (day1 < 1) {
                      day1 = day1+366;
                      day2 = day+16;
                      year_tmp = year-1;
                  }
                  if (!grd) {
                      sprintf(ascat_path,
                         "/auto/temp/lindell/soilmoisture/ave/%04d/%03d-%03d-%04d/msfa-%s-%s%02d-%03d-%03d.ave",
                         year_tmp,day1,day2,year_tmp,type,region,year_tmp-2000,day1,day2);
                  }
                  else { // grab grd files
                      sprintf(ascat_path,
                           "/auto/temp/lindell/soilmoisture/grd/%04d/%03d-%03d-%04d/msfa-%s-%s%02d-%03d-%03d.grd",
                            year_tmp,day1,day2,year_tmp,type,region,year_tmp-2000,day1,day2);
                  }
                  if( access( ascat_path, F_OK ) == -1 ) {
                      printf("Error reading file %s for day %03d, year %d\n",ascat_path,day,year);
                      continue;
                  }
              }
              pthread_mutex_lock(fopen_lock);
              FILE *sir = fopen(ascat_path,"r");
              pthread_mutex_unlock(fopen_lock);

              if (!sir) {
                  printf("SIR FILE ERROR day %03d, year %d\n",day,year);
                  continue;
              }

              sir_init_head(&head);
              get_sir_head_file(sir, &head);

              for (row = 1; row <= num_rows; row++) {
                    for (column = 1; column <= num_columns; column++) {
                      // retrieve pixel
                      if (get_sir_data_block(sir, &pix_val, &head, column, row, column, row) < 0)
                          printf("ERROR READING SIR DATA BLOCK!\n");

                      // store pixel in buffer
                      tseries[row-1][column-1][year-YEAR_START][day-1] = pix_val;
                  }
              }

              // close file
              fclose(sir);
          }
	}
	return NULL;
}

int main (int argc, char **argv)
{
    struct arguments arguments;

    // Set up variables
    static const int NUM_DAYS = 365;
    int num_columns;
    int num_rows;

    // multithreading params
    thread_args t_args[NUM_THREADS];
    pthread_t thread_id[NUM_THREADS];
    int ind_per_thread;
    int start_index;
    int stop_index;
    int i,j,k;

    // Initialize NETCDF Variables
    int ncid, row_dimid, col_dimid, year_dimid, day_dimid;
    int varid;
    int retval;
    char FILE_NAME[100];
    int dimids[NDIMS];

    /* Default values. */
    arguments.verbose = 0;
    arguments.grd = 0;
    arguments.region = NULL;
    arguments.type = NULL;

    /* Parse our arguments; every option seen by parse_opt will
     be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    char* region = arguments.region;
    char* type = arguments.type;
    int grd = arguments.grd;

    printf ("GEN_TIME_SERIES\n---------------\nBeginning processing with options:\n");

    printf ("Region = %s\nVERBOSE = %s\nIMAGE_TYPE = %s\n---------------\n",
          arguments.region,
          arguments.verbose ? "yes" : "no",
          arguments.type);

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
    memset(tseries[0][0][0],0,num_rows*num_columns*NUM_YEARS*NUM_DAYS*sizeof(float));

    printf("Done\n");

    // sort through pixels/images
    // split up row processing based on number of threads/rows
    ind_per_thread = 364 / NUM_THREADS;

    if ((ind_per_thread % 2) == 0) {
        ind_per_thread--;
    }

    start_index = 1;
    stop_index = 0;
    for (i = 0; i < NUM_THREADS; i++) {
        if (i == NUM_THREADS - 1) {
            stop_index = 364;
        } else {
            stop_index = start_index + ind_per_thread;
        }
        t_args[i].tseries = tseries;
        t_args[i].start_i = start_index;
        t_args[i].stop_i = stop_index;
        t_args[i].num_rows = num_rows;
        t_args[i].num_columns = num_columns;
        t_args[i].region = region;
        t_args[i].grd = grd;
        t_args[i].type = type;
        t_args[i].NUM_DAYS = NUM_DAYS;
        t_args[i].fopen_lock = &fopen_lock;
        start_index = stop_index + 1;
    }

    // submit threads
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_create(&thread_id[i], NULL, mthreadParseImg, &t_args[i]);
    }

    // join threads
    for (i = 0; i < NUM_THREADS; i++) {
        pthread_join(thread_id[i], NULL);
    }

    printf("Saving NetCDF File...");
    /* Create the file. */
    sprintf(FILE_NAME,"/auto/temp/lindell/soilmoisture/ts/ts_%s_%s.nc",region,type);
    if ((retval = nc_create(FILE_NAME, NC_NETCDF4, &ncid)))
        ERR(retval);

    /* Define the dimensions. */
    if ((retval = nc_def_dim(ncid, "row", num_rows, &row_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "column", num_columns, &col_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "year", NUM_YEARS, &year_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "day", NUM_DAYS, &day_dimid)))
          ERR(retval);

    /* Define the netCDF variables. The dimids array is used to pass
        the dimids of the dimensions of the variables.*/
    dimids[0] = row_dimid;
    dimids[1] = col_dimid;
    dimids[2] = year_dimid;
    dimids[3] = day_dimid;

    /* define the variable */
    if ((retval = nc_def_var(ncid, "data", NC_FLOAT, NDIMS, dimids, &varid)))
        ERR(retval);

    /* End define mode. */
    if ((retval = nc_enddef(ncid)))
        ERR(retval);

    /* Write the data. */
    if ((retval = nc_put_var_float(ncid, varid, &tseries[0][0][0][0])))
        ERR(retval);

    /* Close the file. */
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("done\n");

    // Free memory for 3D image timeseries array
    printf("Finishing up...");
    free(tseries[0][0][0]);
    free(tseries[0][0]);
    free(tseries[0]);
    free(tseries);

    printf("done\n");

    exit (0);
}
