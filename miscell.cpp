// this cpp file implements all the header functions introduced in miscell. h 





#include <iostream> 
#include<string> 
//#include<conio.h>
#include<math.h> 
#include<cstdlib>
#include<ctime> 
#include "miscell.h"
#include<vector> 
using namespace std ; 
void sayHello()//Implementation
{
    cout << "Hello!";
}

int partition(int *a, int *index, int left, int right)
{	int pivot = a[left];
	while (true)
	{
		while (a[left] < pivot) left++;
		while (a[right] >= pivot) right--;
        //		cout << "left=" << left << endl ;
        //		cout << "right=" << right << endl ;
        //		system("pause>nul");   // waits until a key is pressed
		if (left < right)
		{
			swap(a[left], a[right]);
            //		printArray(a, 0, 6, "a") ;
            
			swap(index[left],index[right]);
		}
		else
		{
			return left;
		}
	}
    
}

void quicksort(int *a, int *index, int left, int right)
{
	if (left < right)
	{
        int pivot = partition(a, index, left, right);
		quicksort(a, index, left, pivot-1);
		quicksort(a, index, pivot+1, right);
	}
}

// pringitng matrix

void printMatrix(int** anMatrix, int ro, int co, string MatrixName)
{
	cout<<"Printing "<<MatrixName <<"["<<ro<<","<<co<<"]"<<endl; // ,char MatrixName
	for (int i = 0; i<ro; i++ ) {
		for (int j = 0; j<co; j++ ){
			cout<<anMatrix[i][j]<<"\t";
		}
		cout<<endl;
	}
}
// printing array
void printArray(int anArray[], int left, int right, string ArrayName)
{	cout<<"Printing "<<ArrayName <<"["<<left<<"..."<<right<<"]"<<endl; // ,char MatrixName
	for (int i = left; i<=right; i++ )
    {
        cout<<anArray[i]<<"\t";
    }
}

// computing factorial

// -----------factorial ----------
int factorial( int n) {
    int result ;
    if (n > 1)
        result = n * factorial(n-1) ;
    else
        result = 1 ;
    return  result ;
    
}


// ----------------- Reverse order -------------
/*  ------------------- Reverse the sentence word by word ---------------- */
bool reverseWords( char str[]){
    // char * buffer;
    int tokenReadPos, wordReadPos, wordEnd, writePos = 0;
    
    /* Position of the last character is length - 1 */
    tokenReadPos =  strlen(str) - 1;
    char *buffer ;
    buffer = new char [ tokenReadPos+2] ;
    // buffer = (char *) malloc(tokenReadPos + 2);   //allocate a block (tokenreadspos +2) of memory
    if( !buffer )
        return false; /* reverseWords failed */
    
    while( tokenReadPos >= 0 ){
        
        if( str[tokenReadPos] == ' ' ){ /* Non-word characters */
            /* Write character */
            buffer[writePos++] = str[tokenReadPos--];
            
        } else {  /* Word characters */
            /* Store position of end of word */
            wordEnd = tokenReadPos;
            
            /* Scan to next non-word character */
            while( tokenReadPos >= 0 && str[tokenReadPos] != ' ' )
                tokenReadPos--;
            
            /* tokenReadPos went past the start of the word */
            wordReadPos = tokenReadPos + 1;
            
            /* Copy the characters of the word */
            while( wordReadPos <= wordEnd ){
                buffer[writePos++] = str[wordReadPos++];
            }
        }
    }
    /* null terminate buffer and copy over str */
    buffer[writePos] = '\0';
    strcpy(str, buffer);  // copy two strings to each other.
    //str = buffer ;
    free(buffer);
    
    return true; /* ReverseWords successful */
}



/* ----------obtain position of words */

void charindexarray( char str [], vector<int>  D[] ){
    int positionindex = 0 ;
    int len =  strlen(str) ;
    for (positionindex = 0; positionindex < len ; positionindex++){
        int chartoint = str[positionindex] - 'a' ;
        D[chartoint].push_back(positionindex);
        
    }
}


// -----------  shortest sequence containing a sequence----- 

/*void StartEndSequence ( vector<int> D[], char L[], int start, end )
{
    int len_L = strlen(L) -1 ;
    int endofL = L(len_L);
    int sizeofDofEndofL = D[endofL].size() - 1 ;
    temp = 1000; 
    for (int i=sizeofDofEndofL ; i>-1 ; i--)
    {
        end_temp = D[endofL][i] ;
        char buffer[] ;
        strncpy(buffer, L, len_L);  // L-1 elements were copied to buffer
        start_temp = startseq(buffer, end)
        if (end-)
    }
}*/


















