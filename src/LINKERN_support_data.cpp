#include <stdio.h>
#include <stdlib.h>

#include "LINKERN_support_data.h"


static int
    matrix_edgelen (int i, int j, compass_data *data);

typedef struct compass_data compass_data;

static int edgelen_nonorm (int i, int j, compass_data *data)
{
    fprintf(stderr, "erroor %d\n", data->adj[i][j]);
    fprintf(stderr, "OPutil_dat_edgelen has been called with no norm set\n");
    fprintf(stderr, "This is a FATAL ERROR\n");
    if (i != 0 || j != 0 || data != (compass_data *) NULL) {
        /* so the compiler won't complain about unused variables */
        fprintf(stderr, "This is a FATAL ERROR\n");
        exit (1);
    }
    exit (1);
    return -1;
}

static void init_data ( compass_data *data)
{ data->n = 0;
  data->x = (double *) NULL;
  data->y = (double *) NULL;
  data->z = (double *) NULL;
  data->adj = (int **) NULL;
  data->adjspace = (int *) NULL;
  data->len      = (int **) NULL;
  data->lenspace = (int *) NULL;
  data->degree   = (int *) NULL;
  data->norm = 0;
  data->dsjrand_param = 1;
  data->dsjrand_factor = 1.0;
  data->default_len = 100000;
  data->sparse_ecount = 0;
  data->edgelen = edgelen_nonorm;

  data->ndepot = 0;
  data->orig_ncount = 0;
  data->depotcost = (int *) NULL;
  data->orig_names = (int *) NULL;
  return;
}

void compass_init_data( compass_data *data)
{ init_data(data);
  return;
}

int compass_get_edge_len (int i, int j, compass_data *data)
{
    return (data->edgelen)(i, j, data);
}

static void delete_data(compass_data *data);

void compass_erase_data(compass_data *data)
{
      delete_data(data);
      init_data(data);
      return;
}

static void delete_data(compass_data *data)
{
  if (data->adj != (int **) NULL) free (data->adj);
  if (data->adjspace != (int *) NULL) free (data->adjspace);
  return;
}

void compass_delete_data(compass_data *data)
{
  if (data != (compass_data*) NULL){
    delete_data(data);
    free(data);
  }
  return;
}

int compass_data_set_norm (compass_data *data, int norm)
{
    switch (norm) {
    case CC_MATRIXNORM:
        data->edgelen = matrix_edgelen;
        break;
	default:
        fprintf(stderr, "ERROR:  Unknown NORM %d.\n", norm);
        return 1;
    }
    data->norm = norm;

    return 0;
}

static int matrix_edgelen (int i, int j, compass_data *data)
{
    if (i > j)
        return (data->adj[i])[j];
    else
        return (data->adj[j])[i];
}

