/* generate_random_data_file.c

   James S. Plank
   plank@cs.utk.edu
   http://web.eecs.utk.edu/~plank

   Professor
   EECS Department
   University of Tennessee
   Knoxville, TN 37996

   January, 2013.

   Part of the SD coding library.  Please see Technical Report UT-CS-13-704
   http://web.eecs.utk.edu/~plank/plank/papers/FAST-2013-SD.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

void usage(char *s)
{
  fprintf(stderr, "usage: generate_random_data_file n m s r blocksize seed(-1 for time(0))\n");
  if (s != NULL) fprintf(stderr, "%s\n", s);
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, n, m, s, r, size;
  long seed;
  int *erased;

  if (argc != 7) usage("Wrong number of arguments");
  if (sscanf(argv[1], "%d", &n) == 0 || n <= 0) usage("Bad n");  // d: decimal, the variable follow must be a pointer to integer &n.  sscanf returns the number of input items successfully matched.
      // different from which simply reads the input in C. 
  if (sscanf(argv[2], "%d", &m) == 0 || m <= 0 || m >= n) usage("Bad m");
  if (sscanf(argv[3], "%d", &s) == 0 || s <= 0) usage("Bad s");
  if (sscanf(argv[4], "%d", &r) == 0 || r <= 0) usage("Bad r");
  if (m*r+s >= n*r) usage("m*r+s is greater than or equal to n*r");
  if (sscanf(argv[5], "%d", &size) == 0 || size <= 0) usage("Bad size");
  if (sscanf(argv[6], "%ld", &seed) == 0) usage("Bad seed");
  if (size % 8 != 0) usage("Size has to be a multiple of 8\n");

  if (seed == -1) seed = time(0);
  srand48(seed);

  erased = (int *) malloc(sizeof(int)*n*r);
  for (i = 0; i < n*r; i++) erased[i] = 0;
  for (i = 0; i < m; i++) for (j = 0; j < r; j++) {
    erased[n-i-1+j*n] = 1;
  }
  j = n*r-1;
  for (i = 0; i < s; i++) {
    while (erased[j]) j--;
    erased[j] = 1;
  }
  for (i = 0; i < n*r; i++) {
    if (erased[i]) {
      for (j = 0; j < size; j++) printf("00 ");
    } else {
      for (j = 0; j < size; j++) printf("%02lx ", lrand48()%256);
    }
    printf("\n");
  }
  return 0;
}
