#include <stdio.h>
#include <stdlib.h>
#include "NSSTools.h"

void freeSpBlMat(SparseBlockStructuredMatrix *blmat)
{

  /*    int i; */

  free(blmat->RowIndex);
  free(blmat->ColumnIndex);
  free(blmat->blocksize);
  /*    for (i = 0 ; i < blmat->nbblocks ; i++) free(blmat->block[i]);  */
  free(blmat->block);

}

