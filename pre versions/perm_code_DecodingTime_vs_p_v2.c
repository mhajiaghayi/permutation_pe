/* Permutation code implementation with modified recoverying scheme. 
 
 
 coder: Mahdi Hajiaghayi
 Program summary: it implements the encoding and failure recovery of the permutation code based on the Jafar's paper on permutational code and our modified scheme. 
 
 plot: It plots the recovery time vs the number of participating nodes. The total 
 size of data is fixed while changin
 
 
 May, 2013.
 
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
// printing a matrix containing the data. 
void print_data(int n, long int r, int size, int **array, char *file)
{
    int i, j;
    FILE *fout;
    
    fout = fopen(file, "w");
    if (fout == NULL) { perror(file); exit(1); }
    
    for (i = 0; i < n*r; i++) {
        for (j = 0; j < size; j++) {
            if (j != 0) fprintf(fout, " ");
            fprintf(fout, "%04x", array[i][j]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}



/* --------- My functions ------------------ */
void dec2basep(int n, int p, int *v){
    int i = 0;
    for (int k = n  ; k> 0 ; k/=p) {
        v[i] = k % p ;
        i++ ; }
    
}
int basep2dec( int *v, int p, int size) {
    int n = 0 ;
    for (int i =0 ; i < size; i++)
        n = n + pow(p,i)*v[i] ;
    return n ; 
}

int XoRbasep(int a,int b, int p) {
    return ((a +b) % p) ;
}

void vec_initialize(int *vec, int*vec0, int size){
    for (int i=0 ; i<size; i++)
        vec[i] = vec0[i];
}
void Init_zero(int * vec, int size ){
    for (int i=0; i<size; i++)
        vec[i] = 0 ;
}
// printing array

void printArray(double anArray[], int left, int right, char *ArrayName)
{	printf("Printing %s [ %d ... %d ] = ", ArrayName, left, right); // ,char MatrixName
	for (int i = left; i<=right; i++ )
    {
        printf("%f, ",anArray[i]) ;
    }
    printf("\n") ;
    
}

main(int argc, char **argv)
{


    gf_t gfm;  // galois field objects
    void *scratch ; 
    int i, j, n, k,  m,   w, c, size, symb_size, p ,q , pnum, itr;
    int itr_max; 
    
    long double M ; long int L;

    int **array;
    int **parity ;
    int ** new_node ;
    int temp ; 
    
    
    int *erased; 
    int *helper_vec;    //  binary list of participating parity nodes 



    int indx, indx_hat, indx_q, indx_f ;
    int x,x_q, x_hat ;  // x: vector number.
    int *x_vec, *x_vec_p, *x_vec_q, *x_vec_f ;  // row numbers in a vector form base (n-k)

    double *RecoveryTime ;
    

    int *lam;
    uint32_t coef;
    
    double timer, split;
    int IO_read, IO_write;
    FILE *fin ;  // *fout 
    
    /* ----------- Initialization. test params ---- */
    n = 10 ;  k =3;  m = (n-k) ; w = 16; IO_read =0; IO_write =0  ;
    itr_max = 1;
    L = (int) round(pow (m,k));
    int max_data = ((int) round(pow(2,w))); 
    /* ------- file size M to block size, size ---- */  
    
    M = 32 ;   // in Meg.
    size = (int) round(M * pow(1024,2)/k/L/(w/8));
    
    /* ------ block size, size, to file size M ------ */
    
   // size = 1 ;
   // M = L *L *size *(w/8)/pow(1024,2);  //k*L ;  // in Meg */
   // fprintf(stderr, "M%Lf",M) ;
    
    symb_size = size * ((int) round(w/8)); 
    lam = talloc(int, k);
    for (i=0; i<k; i++)
        lam[i] = 1 ;  // i+1 initalize the coefficients

    //   printf("L %02d", l) ;
   // printf("size L %02d\n", L);
    

        
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
    
 /*   fpcm = fopen("All-Coefficients.txt", "r");
    if (fpcm == NULL) {
        usage("bad coef matrix");
    }  */
    
    scratch =  NULL ;  // (void *) malloc(symb_size);
    if (w == 16 || w == 32) {
        if (!gf_init_hard(&gfm, w, GF_MULT_SPLIT_TABLE, GF_REGION_ALTMAP | GF_REGION_SSE, GF_DIVIDE_DEFAULT, 0, 4, w, NULL, scratch)) {
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
    
    
    
    RecoveryTime = talloc(double, m);
    for (i=0; i<m; i++) RecoveryTime[i] = 0 ;
    
 
    
for (itr = 0; itr< itr_max; itr++) {   // monte Carlo simulation
    
    fprintf(stderr, "RecoveryTime[0] = %f", RecoveryTime[0]);
    
    helper_vec = talloc(int,m) ;
    for (i =0 ; i<m; i++)
        helper_vec[i] = 0 ;
    
/* --------- Data loading to all n*L blocks!  --- */

    srand48(time(0));
 // srand48(1) ; // generates the same value.
    
    array = talloc(int *, n*L);   // 2-d matrix with size ( n*r x size) holds the data.
    for (i = 0; i < n*L; i++) array[i] = talloc(int, size);
    

    for (i = 0; i < n*L; i++) {
        for (j = 0; j < size; j++) {
            if (!IO_read) {
                array[i][j] = lrand48()% max_data;    // generate it randomly or
 //               fprintf(stderr,"array %d, temp %d", array[i][j],temp) ;
            } else {
                if (fscanf(fin, "%x", &c) != 1) {  // read it from a file
                    fprintf(stderr, "Bad input at block %d byte %d\n", i, j);
                } else array[i][j] = c;
            }
        }
    }
    if (IO_read) fclose(fin);
    
    /* ------------- Setting the parity to zero----------*/
    
    
    erased = talloc(int, n*L);
    for (i = 0; i < n*L; i++) erased[i] = 0;  // set erased vector to zero.

    
    // the number of erased disks should be m
    for (i = k; i < n; i++) {
        //    if (sscanf(argv[11+i], "%d", &j) != 1 || j < 0 || j >= n) usage("Bad d_x");
        j = i ;
        if (erased[j]) usage("Duplicate failed disk");
        while (j < n*L) {
            erased[j] = 1;      // erased vector with size (n*r) indicates which sector is erased.
            j += n;             //going over all parity rows for column j
        }
    }
    for (i = 0; i < n*L; i++) {
        if (erased[i]) {            // if erased[i] = 1 => sector i is erased.
            bzero(array[i], symb_size);  // set the first size * ((int) round(w/8)) bytes to zero
        }
    }
    /*  print the loaded data in a file */
    if (IO_write)   print_data(n, L, size, array, "systematic.txt");
 
    
    
    /* ----------------- Encoding using permutation code ---------- */

    parity= talloc(int *, (m*L));
    for (i = 0; i < m*L; i++) {
        parity[i] = talloc(int, size);
        bzero(parity[i], symb_size);
    }


    x_vec = talloc(int, k) ;
    x_vec_p = talloc(int, k) ;
    x_vec_f = talloc(int, k) ;
    x_vec_q = talloc(int, k) ;
    // main loop 
    for (p =0; p < m;  p++) {  // over parity nodes 1:m 
        for (x =0 ; x< L; x++) {  // for each row 1:L => 0 .. L-1
            Init_zero(x_vec, k) ; // 
         //   printf("x_vec before %d ",x_vec[k-1]);
            dec2basep(x, m, x_vec) ; // convert it to base m(n-k). note that LSB is x_vec[0] 


            for (i = 0; i < k; i++) {
                vec_initialize(x_vec_p, x_vec, k);  // initalize x_vec_p with x_vec. next applies
                x_vec_p[k-i-1] = XoRbasep(x_vec_p[k-i-1],p,m); // x_vec_p[k-i-1] +p mod m
                x_hat = basep2dec(x_vec_p, m, k); 
              //  printArray(x_vec_p, 0, k-1, "x_vec_p");
                coef = round( pow(lam[i],i)) ; // to make it integer 
                indx_hat = x_hat * n + i ;
                indx = x*m + p;  // colum p, row x; 
                if (coef != 0) 
                    gfm.multiply_region.w32(&gfm, array[indx_hat], parity[indx], coef, size, 1);
                // it is summed over
           //     fprintf(stderr, "parity: %d, array: %d",parity[indx][0], array[indx_hat][0]);
           }
            
        }
    }
   if (IO_write)      print_data(m, L, size, parity, "parity.txt");
    
    
    /* ------------  Single failure recovery ---------- */ 

//    int helper_dec = basep2dec(helper_vec, 2, m) ;

    new_node= talloc(int *, L);
    
    
  for (pnum =1; pnum<= m; pnum++) {  // main loop for different number of participating nodes 
     helper_vec[pnum-1] = 1 ;

    
    for (i = 0; i < L; i++) {
        new_node[i] = talloc(int, size);
        bzero(new_node[i], symb_size);
    }
    
    
    timer_start(&timer);
    
    for ( p=0; p < m; p++ )  // over each parity nodes
    {
        if (helper_vec[p] ==1 ) {
            for (x=0; x< round(L/m); x++)  // focus on each row between 1 and L/m
            {   Init_zero(x_vec, k) ;
                dec2basep(x, m, x_vec) ;
                indx_f = p* round(L/m) + x ;  // failed node index ;
                indx = x*m + p;  // colum p, row x of parity nodes ;
                
                for (i =1; i<k; i++ ) { // cancell the effect of all systematic nodes on each row of parity nodes.
                    vec_initialize(x_vec_p, x_vec, k);  // initalize x_vec_p with x_vec. next applies
                    //   printf("x_vec %d ",x_vec[k-i-1]);
                    x_vec_p[k-i-1] = XoRbasep(x_vec_p[k-1-i],p,m);
                    x_hat = basep2dec(x_vec_p, m, k);
                    //       printArray(x_vec_p, 0, k-1, "x_vec_p");
                    coef = round( pow(lam[i],i)) ; // to make it integer
                    indx_hat = x_hat * n + i ;
                    new_node[indx_f] = parity[indx]; // load each sector of new_node with its corresponding parity sector, and subtract (add) the interference parts.
                    if (coef != 0)
                        gfm.multiply_region.w32(&gfm, array[indx_hat], new_node[indx_f], coef, size, 1);
                    //
                    
                }
                
            }
            
        }
        else {
            for ( q=0; q < m; q++ ) {
                if (helper_vec[q] ==1 ) {  // find the first available parity node
                    for (x =0; x< round(L/m); x++) {
                        indx_f = p*round(L/m) + x ;  // index of sectors missing from the first phase of recovery 
                        Init_zero(x_vec_f, k);
                        dec2basep(indx_f, m, x_vec_f);  // 
                        x_vec_f[k-1] = XoRbasep(x_vec_f[k-1],m-q,m); // location of that sector in parity node q in a binary vector form.
                        x_q = basep2dec(x_vec_f, m, k);
                        indx_q = x_q*m + q ; // in parity matrix 
                        for (i=1; i<k ; i++) {
                            vec_initialize(x_vec_q,x_vec_f,k);
                            x_vec_q[k-i-1] = XoRbasep(x_vec_q[k-1-i],q,m);
                            x_hat = basep2dec(x_vec_q, m, k) ;
                            coef = round(pow(lam[i],i));
                            indx_hat = x_hat *n + i ;
                            new_node[indx_f] = parity[indx_q];
                            if (coef != 0)
                                gfm.multiply_region.w32(&gfm, array[indx_hat], new_node[indx_f], coef, size, 1);

                        }
                        
                    }
                    break;
                }
                else {
                    if (q==m-1) {
                        printf("recovery is impossible with the selected parity nodes \n");
                        break;
                    }
                }
                    
            }
        }
    
    }
    RecoveryTime[pnum-1] = 1/itr_max * timer_split(&timer) + RecoveryTime[pnum-1];
    split = RecoveryTime[pnum-1];
    fprintf(stderr, "RecoveryTime[0] = %f", RecoveryTime[0]);
  } // end of main loop
   fprintf(stderr, "RecoveryTime[0] = %f", RecoveryTime[0]);
}
    printArray(RecoveryTime, 0, m-1, "Recovery Time");
    fprintf(stderr, "Timer: %.6lf\n", split);
    fprintf(stderr, "Total filze size: %.10Lf Meg \n", M);
    /* ---- Print the new node ----  */
    print_data(1, L, size, new_node, "newnode.txt");
    
    return 0 ;
    
}