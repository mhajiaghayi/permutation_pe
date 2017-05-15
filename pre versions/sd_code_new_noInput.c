/* sd_codec.c
 
 Revision 2.0
 
 James S. Plank
 plank@cs.utk.edu
 http://web.eecs.utk.edu/~plank
 
 Professor
 EECS Department
 University of Tennessee
 Knoxville, TN 37996
 
 May, 2013.
 
 Part of the SD coding library.  Please see XXXX
 */

#include <strings.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h> 

#include "gf_complete.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

void usage(char *s)
{
    fprintf(stderr, "usage: sd_code n m s r w size IO(0|1) PCM Input Output d_0 .. d_m-1 s_0 .. s_x-1\n\n");
    //stderr stands for standard error device. In console programming it is the console -- the screen. It is essentially the same as stdout. it shows this message in the console instead of in a file. 
    fprintf(stderr, "The Parity Check Matrix should be in the file PCM.\n");
    fprintf(stderr, "If IO = 0, Input and Output are ignored.\n");
    fprintf(stderr, "Otherwise, Input contains the codeword -- size bytes per block,\n");
    fprintf(stderr, "      The failed blocks will be ignored, but you must include them in the file.\n");
    fprintf(stderr, "Output contains the decoded codeword -- size bytes per block.\n");
    fprintf(stderr, "You must have m disk failures and s sector failures.\n");
    fprintf(stderr, "\n");
    if (s != NULL) fprintf(stderr, "%s\n", s);
    exit(1);
}

void timer_start (double *t)
{
    struct timeval  tv;
    
    gettimeofday (&tv, NULL);
    *t = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

double
timer_split (const double *t)
{
    struct timeval  tv;
    double  cur_t;
    
    gettimeofday (&tv, NULL);
    cur_t = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
    return (cur_t - *t);
}

void print_data(int n, int r, int size, uint8_t **array, char *file)
{
    int i, j;
    FILE *fout;
    
    fout = fopen(file, "w");
    if (fout == NULL) { perror(file); exit(1); }
    
    for (i = 0; i < n*r; i++) {
        for (j = 0; j < size; j++) {
            if (j != 0) fprintf(fout, " ");
            fprintf(fout, "%02x", array[i][j]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}

int invert_matrix(int *mat, int *inv, int rows, gf_t *gf)
{
    int cols, i, j, k, x, rs2;
    int row_start, tmp, inverse;
    
    cols = rows;
    
    k = 0;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            inv[k] = (i == j) ? 1 : 0;
            k++;
        }
    }
    
    /* First -- convert into upper triangular  */
    for (i = 0; i < cols; i++) {
        row_start = cols*i;
        
        /* Swap rows if we ave a zero i,i element.  If we can't swap, then the
         matrix was not invertible  */
        
        if (mat[row_start+i] == 0) {
            for (j = i+1; j < rows && mat[cols*j+i] == 0; j++) ;
            if (j == rows) return -1;
            rs2 = j*cols;
            for (k = 0; k < cols; k++) {
                tmp = mat[row_start+k];
                mat[row_start+k] = mat[rs2+k];
                mat[rs2+k] = tmp;
                tmp = inv[row_start+k];
                inv[row_start+k] = inv[rs2+k];
                inv[rs2+k] = tmp;
            }
        }
        
        /* Multiply the row by 1/element i,i  */
        tmp = mat[row_start+i];
        if (tmp != 1) {
            inverse = gf->divide.w32(gf, 1, tmp);
            for (j = 0; j < cols; j++) {
                mat[row_start+j] = gf->multiply.w32(gf, mat[row_start+j], inverse);
                inv[row_start+j] = gf->multiply.w32(gf, inv[row_start+j], inverse);
            }
        }
        
        /* Now for each j>i, add A_ji*Ai to Aj  */
        k = row_start+i;
        for (j = i+1; j != cols; j++) {
            k += cols;
            if (mat[k] != 0) {
                if (mat[k] == 1) {
                    rs2 = cols*j;
                    for (x = 0; x < cols; x++) {
                        mat[rs2+x] ^= mat[row_start+x];
                        inv[rs2+x] ^= inv[row_start+x];
                    }
                } else {
                    tmp = mat[k];
                    rs2 = cols*j;
                    for (x = 0; x < cols; x++) {
                        mat[rs2+x] ^= gf->multiply.w32(gf, tmp, mat[row_start+x]);
                        inv[rs2+x] ^= gf->multiply.w32(gf, tmp, inv[row_start+x]);
                    }
                }
            }
        }
    }
    
    /* Now the matrix is upper triangular.  Start at the top and multiply down  */
    
    for (i = rows-1; i >= 0; i--) {
        row_start = i*cols;
        for (j = 0; j < i; j++) {
            rs2 = j*cols;
            if (mat[rs2+i] != 0) {
                tmp = mat[rs2+i];
                mat[rs2+i] = 0;
                for (k = 0; k < cols; k++) {
                    inv[rs2+k] ^= gf->multiply.w32(gf, tmp, inv[row_start+k]);
                }
            }
        }
    }
    return 0;
}


main(int argc, char **argv)
{
/*    gf_t gfs, gfm;
    int i, j, n, m, r, s, w, c, size, nf;
    int *matrix;
    int *encoder, *inverse;
    int *erased, *eindex;
    int k, symbol, l;
    uint32_t coef;
    uint8_t **array;
    uint8_t **syndromes;
    double timer, split;
    int IO;
    FILE *fin, *fout, *fpcm;
    // test params 
    n = 6 ;    m = 2 ; s = 2 ; r = 4 ; w = 8 ; size = 8 ; IO =0;
    int L= pow(3,5);
    int d_x[2] = { 4 , 5} ;  // erased disks
    int s_x[] = {20,21};
    double LL = pow(n,m); */ 

    gf_t gfs, gfm;
    int i, j, n, k, m, r, s, w, c, size, nf;
    double LL;
    uint8_t **array;
    int *matrix;
    int *lam;
    
    
    int *encoder, *inverse;
    int *erased, *eindex;
    //    int symbol, l;
    uint32_t coef;
    
    uint8_t **syndromes;
    double timer, split;
    int IO;
    FILE *fin, *fpcm;  // *fout
    // Initialization. test params
    n = 5 ;  k=2;  m = (n-k) ;  s = 0 ; r = 4 ; w = 8 ; size = 1 ; IO =0;
    LL = pow (m,k);
    //   printf("L %02d", l) ;
    int d_x[2] = { 4 , 5} ;  // erased disks
    int s_x[] = {20,21};
    fprintf(stderr, "size L %02f/n", LL);
    lam = talloc(int, k);
    for (i=0; i<k; i++)
        lam[i] = i+1 ;
    
    
    
/*
    if (argc <= 8) usage(NULL);
    if (sscanf(argv[1], "%d", &n) == 0 || n <= 0) usage("Bad n");
    if (sscanf(argv[2], "%d", &m) == 0 || m <= 0) usage("Bad m");
    if (sscanf(argv[3], "%d", &s) == 0 || s <  0) usage("Bad s");
    if (sscanf(argv[4], "%d", &r) == 0 || r <= 0) usage("Bad r");
    if (sscanf(argv[5], "%d", &w) == 0 || w <= 0) usage("Bad w");
    if (sscanf(argv[6], "%d", &size) == 0 || size <= 0) usage("Bad size");
    if (sscanf(argv[7], "%d", &IO) == 0 || IO < 0 || IO > 1) usage("Bad IO");
    if (w != 8 && w != 16 && w != 32) usage("W has to be 8, 16 or 32\n");
    if (argc != 11 + (m+s)) usage("Wrong number of arguments");
    
    if (w == 16 && size % 2 != 0) usage("When w=16, size must be a multiple of 2");
    if (w == 32 && size % 4 != 0) usage("When w=32, size must be a multiple of 4");
    */
 //   fpcm = fopen(argv[8], "r");  // read module 
 //   if (fpcm == NULL) { perror(argv[8]); exit(1); }
    fpcm = fopen("All-Coefficients.txt", "r");
    if (fpcm == NULL) {
        usage("bad coef matrix");
    }
    if (w == 16 || w == 32) {
        if (!gf_init_hard(&gfm, w, GF_MULT_SPLIT_TABLE, GF_REGION_ALTMAP | GF_REGION_SSE, GF_DIVIDE_DEFAULT, 0, 4, w, NULL, NULL)) {
            printf("Bad gf spec\n");
            exit(1);
        }
    } else if (w == 8) {
        if (!gf_init_easy(&gfm, w)) {
            printf("Bad gf spec\n");
            exit(1);
        }
    } else {
        printf("Not supporting w = %d\n", w);
    }
    
    if (!gf_init_easy(&gfs, w)) {
        printf("Bad gf spec\n");
        exit(1);
    }
    
    
    fpcm = fopen("All-Coefficients.txt", "r");
    if (fpcm == NULL) {
        usage("bad coef matrix");
    }
    
    
    if (!IO) {
        srand48(time(0));
    } else {
        fin = fopen(argv[9], "r");
        if (fin == NULL) { perror(argv[9]); exit(1); }
    }
    
    array = talloc(uint8_t *, n*LL);   // 2-d matrix with size ( n*r x size) holds the data.
    
    for (i = 0; i < n*LL; i++) array[i] = talloc(uint8_t, size);
    
    // data loading to all n*r blocks!
    for (i = 0; i < n*LL; i++) {
        for (j = 0; j < size; j++) {
            if (!IO) {
                array[i][j] = lrand48()%256;    // generate it randomly or
            } else {
                if (fscanf(fin, "%x", &c) != 1) {  // read it from a file
                    fprintf(stderr, "Bad input at block %d byte %d\n", i, j);
                } else array[i][j] = c;
            }
        }
    }
    if (IO) fclose(fin);
    
    /* Setting the parity to zero */
    
    nf = LL*m + s;
    
    erased = talloc(int, n*r);
    for (i = 0; i < n*LL; i++) erased[i] = 0;
    eindex = talloc(int, nf);
    
    // the number of erased disks should be m
    for (i = k; i < n; i++) {
        //    if (sscanf(argv[11+i], "%d", &j) != 1 || j < 0 || j >= n) usage("Bad d_x");
        j = i ;
        if (erased[j]) usage("Duplicate failed disk");
        while (j < n*LL) {
            erased[j] = 1;      // erased vector with size (n*r) indicates which sector is erased.
            // if erased[i] = 1 => sector i is erased.
            j += n;
        }
    }
    for (i = 0; i < n*LL; i++) {
        if (erased[i]) {
            bzero(array[i], size);  // turn them into zero
        }
    }
     print_data(n, LL, size, array, "test.txt");
    
    
    matrix = talloc(int, n*r*(m*r+s));
    encoder = talloc(int, (m*r+s)*(m*r+s));
    inverse = talloc(int, (m*r+s)*(m*r+s));
    
    /* Encoding */
    
    
    
    
    
    
    /*
     for (i = 0; i < s; i++) {
     //   if (sscanf(argv[11+m+i], "%d", &j) != 1 || j < 0 || j >= n*r) usage("Bad s_y");
     j = s_x[i];
     if (erased[j]) usage("Duplicate failed sector");
     erased[j] = 1;
     }
     
     j = 0;
     for (i = 0; i < n*r; i++) {
     if (erased[i]) {
     eindex[j] = i;  // save the index of erased sectors in eindex vector.
     j++;
     }
     }
     */
    
    
    
    
    
    
    
    /* parity loading */
    syndromes = talloc(uint8_t *, (m*r+s));
    for (i = 0; i < m*r+s; i++) {
        syndromes[i] = talloc(uint8_t, size);
        bzero(syndromes[i], size);
    }
    
    /* Read the parity check matrix. */
    
    for (i = 0; i < n*r*(m*r+s); i++) {
        
        //       printf( "Hello %d \n",fscanf(fpcm,"%u", &coef)) ;
        //     fscanf(fpcm,"%u", &coef);
        //    fprintf(stderr, "coef %d", coef);
        if (fscanf(fpcm, "%u", &coef) != 1) {
            fprintf(stderr, "Too few elements in the parity check matrix.\n");
            exit(1);
        }
        if (w != 32 && coef >= (1 << w)) {
            fprintf(stderr, "Bad Parity check element row %d col %d: %u\n", i/n/r, i %(n*r), coef);
            exit(1);
        }
        matrix[i] = coef;
    }
    fclose(fpcm);
    
    /* Erase the bad blocks (turn them into zeros, and create a decoding matrix
     from the columns of the parity check matrix that correspond to the failed
     blocks. */
    
    
    
    timer_start(&timer);
    if (invert_matrix(encoder, inverse, r*m+s, &gfs) == -1) {   // r*m+s unknowns for parity sectors. we made an encoder matrix A and calculate the inverse.
        printf("Can't invert\n");
        exit(0);
    }
    
    for (i = 0; i < n*r; i++) {
        if (!erased[i]) {
            for (j = 0; j < m*r+s; j++) {
                coef = matrix[j*n*r+i];
                if (coef != 0) {
                    gfm.multiply_region.w32(&gfm, array[i], syndromes[j], coef, size, 1);
                }
            }
        }
    }
    /* ultimately save the syndroms[j] into the data */
    for (i = 0; i < m*r+s; i++) {
        for (j = 0; j < m*r+s; j++) {
            coef = inverse[j*(m*r+s)+i];
            if (coef != 0) {
                gfm.multiply_region.w32(&gfm, syndromes[i], array[eindex[j]], coef, size, 1);
            }
        }
    }
    
    split = timer_split(&timer);
    
    fprintf(stderr, "Timer: %.6lf\n", split);
    if (IO) print_data(n, r, size, array, argv[10]);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    return 0 ;
    
}