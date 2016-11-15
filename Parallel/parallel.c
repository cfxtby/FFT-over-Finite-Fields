#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include "DataStruct.h"
#include "tables.h"
#include "time.hpp"
#include <mpi.h>
#include <limits.h>

void taylorExpansion(polynomial *p,long n, long t, uint64_t *g, long rowIndex);
uint64_t* fourierTransform(polynomial* p, long n, long s, int field);
uint64_t WofI(polynomial * p,long i, int m);
void initialTaylor(polynomial *p, long n, long t, long MPI_Index, uint64_t *mine);


void naiveTranspose(void *matrix, void *d, long nRows, long nCols, size_t size,long aR,long aC) {
	long i,j;
	for(i=0;i<nRows;i++) {
		for(j=0;j<nCols;j++) {
			memcpy( (d + (size*( ( j*aR ) + i ))) , (matrix + (size*((i*aC)+j))), size);
		}
	}

}

void transpose(void *matrix, void *d, long nRows, long nCols, long lowCurrRow, long highCurrRow, long lowCurrCol, long highCurrCol, size_t size, int level,int rec) {


	long rowDiff = highCurrRow - lowCurrRow;
	long colDiff = highCurrCol - lowCurrCol;

	if(level==rec) {
		naiveTranspose( (matrix + (size*((lowCurrRow*nCols) + lowCurrCol))), (d + (size*((lowCurrCol*nRows) + lowCurrRow))), rowDiff+1,colDiff+1,size,nRows,nCols);
	}

	else {

		if(rowDiff > colDiff) {

			transpose(matrix,d,nRows,nCols,lowCurrRow,(highCurrRow + lowCurrRow)/2,lowCurrCol,highCurrCol,size,level+1,rec);
			transpose(matrix,d,nRows,nCols,((highCurrRow + lowCurrRow)/2)+1,highCurrRow,lowCurrCol,highCurrCol,size,level+1,rec);
		}
		else if(colDiff > rowDiff) {

			transpose(matrix,d,nRows,nCols,lowCurrRow,highCurrRow,lowCurrCol,(highCurrCol + lowCurrCol)/2,size,level+1,rec);
			transpose(matrix,d,nRows,nCols,lowCurrRow,highCurrRow,((highCurrCol + lowCurrCol)/2)+1,highCurrCol,size,level+1,rec);
		}
		else if((highCurrCol - lowCurrCol) == 0 && (highCurrRow - lowCurrRow) == 0) { //base case
			memcpy( (d + (size*((lowCurrCol * nRows) + lowCurrRow))), (matrix + (size*((lowCurrRow * nCols) + lowCurrCol))), size);
		}
		else {
			transpose(matrix,d,nRows,nCols,lowCurrRow,(highCurrRow + lowCurrRow)/2,lowCurrCol,highCurrCol,size,level+1,rec);
			transpose(matrix,d,nRows,nCols,((highCurrRow + lowCurrRow)/2)+1,highCurrRow,lowCurrCol,highCurrCol,size,level+1,rec);
		}
	
	}

}

int numProcesses,rank,flag;
uint64_t *taylorRow;
long share;
//Assumption: All polynomials are represented in descending order of degrees
int main(int argc, char **argv) {

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	buildTable();
	buildPols();

	long T;
	int gf;
	long N;

	flag = atoi(argv[1]);

	if(rank==0) {

		polynomial *p = (polynomial *) malloc(sizeof(polynomial));
		initialize(p);
		int c = -5;
	
		int  counter = 1;

		uint64_t theCoff;

		FILE *fp = fopen("input.txt","r");

		fscanf(fp,"%d",&gf);
		fscanf(fp,"%ld",&N);
		while(true) {
		
			counter = 1;
			theCoff = 0;

			fscanf(fp,"%d",&c);
			if(c==-1) break;
			else {
				while(true) {
	
					fscanf(fp,"%d",&c);
					if(c==-1) break;
					else {
						setBit(&theCoff,c);
					}
			
				}

				fscanf(fp,"%d",&c);
				insertElement(p,theCoff,c);
			}
		}
		fclose(fp);	


		Time timer;
		timer.start();

		long k=1;
		int m= 0;
	// n= 2^m, get m
		while(k != N) {
			k=k<<1;
			m++;
		}

		T = 1;
		m=m>>1;

		T=T<<m;
	 	share = T/numProcesses;
		uint64_t *rowMatrix;
		if(flag==1 || flag==2) rowMatrix = (uint64_t *) malloc(sizeof(uint64_t) *T*share);
		else rowMatrix = (uint64_t *) malloc(sizeof(uint64_t) *T*T);
		taylorRow = (uint64_t *) malloc(sizeof(uint64_t)*T);

		initialTaylor(p,N,T,0,rowMatrix);

		destroyPolynomial(p);

		free(taylorRow);
		long i,j;

		uint64_t *colMatrix;
		if(flag==1 || flag==2) {
			colMatrix= (uint64_t *) malloc(sizeof(uint64_t)*T*share);
			if(flag==1) transpose(rowMatrix,colMatrix,T,share,0,T-1,0,share-1,sizeof(uint64_t),0,16);
			else naiveTranspose(rowMatrix,colMatrix,T,share,sizeof(uint64_t),T,share);	
			free(rowMatrix);
		}

		uint64_t *transforms;
		polynomial *p1;
		uint64_t *swap;		
		if(flag==1 || flag==2) {
			for(i=0;i<share;i++) {
				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);
				for(j=0;j<T;j++) {
					if(colMatrix[i*T+j]!=0) {
						insertElement(p1,colMatrix[i*T+j],T-1-j);
					}
				}
				transforms = fourierTransform(p1,T,0,gf);

				memcpy(&colMatrix[i*T],transforms,sizeof(uint64_t)*T);
			
				free(transforms);
				//free the polynomial
				destroyPolynomial(p1);
			}
		
			rowMatrix = (uint64_t *) malloc(sizeof(uint64_t) *N);
			if(flag==1) transpose(colMatrix,rowMatrix,share,T,0,share-1,0,T-1,sizeof(uint64_t),0,16);
			else naiveTranspose(colMatrix,rowMatrix,share,T,sizeof(uint64_t),share,T);
		}
	
		else {
			for(i=0;i<share;i++) {
				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);
				for(j=0;j<T;j++) {
					if(rowMatrix[j*share+i]!=0) {
						insertElement(p1,rowMatrix[j*share+i],T-1-j);
					}
				}
				transforms = fourierTransform(p1,T,0,gf);
			
				for(j=0;j<T;j++) {
					rowMatrix[j*share+i] = transforms[j];
				}
			
				free(transforms);
				//free the polynomial
				destroyPolynomial(p1);
			}
			colMatrix= (uint64_t *) malloc(sizeof(uint64_t)*T*share);

		}
		if(T==65536 && numProcesses<16) MPI_Buffer_attach(colMatrix,INT_MAX);
		else MPI_Buffer_attach(colMatrix,T*share);

		for(i=0;i<numProcesses;i++) {
			if(i!=rank) {
				for(j=0;j<share;j++) {
					MPI_Bsend(&rowMatrix[((i*share)+j)*share],share, MPI_UINT64_T, i,(j*T) + (rank*share),MPI_COMM_WORLD);
				}
			}
		}
				
		swap = (uint64_t *) malloc(sizeof(uint64_t)*share*share);

		memcpy(swap,&rowMatrix[rank*share*share],sizeof(uint64_t)*share*share);

		for(i=0;i<share;i++) {
			memcpy(&rowMatrix[i*T + (rank*share)], &swap[(i*share)],sizeof(uint64_t)*share);
		}
		
		free(swap);
		
		for(i=0;i<numProcesses;i++) {
			if(i!=rank) {
				for(j=0;j<share;j++) {
				MPI_Recv(&rowMatrix[j*T + (i*share)],share, MPI_UINT64_T , MPI_ANY_SOURCE,(j*T)+(i*share),MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
            }
		}
  
		int h;
		MPI_Buffer_detach(colMatrix,&h);
		free(colMatrix);   
		
		for(i=0;i<share;i++) {
			p1 = (polynomial *) malloc(sizeof(polynomial));
			initialize(p1);
			for(j=0;j<T;j++) {
				if(rowMatrix[i*T+j]!=0) {
					insertElement(p1,rowMatrix[(i*T)+j],T-1-j);
				}
			}
			//call the fft on p1
			transforms = fourierTransform(p1,T,i*T,gf);

			memcpy(&rowMatrix[i*T],transforms,sizeof(uint64_t)*T);
			free(transforms);
			destroyPolynomial(p1);
		}

		MPI_Barrier(MPI_COMM_WORLD);	
		MPI_Gather(MPI_IN_PLACE,share*T,MPI_UINT64_T,rowMatrix,T*share,MPI_UINT64_T,0,MPI_COMM_WORLD);

		timer.stop();
		timer.printTime();

		FILE *out = fopen("output.txt","w");
		for(i=0;i<N;i++) {
			fprintf(out,"Point %ld: ",i);
			filePrintBeta(out,&rowMatrix[i]);
			fprintf(out,"\n\n");
		}

		free(rowMatrix);
	}
	else {
		polynomial *p = (polynomial *) malloc(sizeof(polynomial));
		initialize(p);
		int c = -5;
	
		int  counter = 1;

		uint64_t theCoff;

		FILE *fp = fopen("input.txt","r");

		fscanf(fp,"%d",&gf);
		fscanf(fp,"%ld",&N);
		while(true) {
		
			counter = 1;
			theCoff = 0;

			fscanf(fp,"%d",&c);
			if(c==-1) break;
			else {
				while(true) {
	
					fscanf(fp,"%d",&c);
					if(c==-1) break;
					else {
						setBit(&theCoff,c);
					}
			
				}

				fscanf(fp,"%d",&c);
				insertElement(p,theCoff,c);
			}
		}
		fclose(fp);


		long k=1;
		int m= 0;
	// n= 2^m, get m
		while(k != N) {
			k=k<<1;
			m++;
		}

		T = 1;
		m=m>>1;

		T=T<<m;

		share = T/numProcesses;

		uint64_t *rowMatrix = (uint64_t *) malloc(sizeof(uint64_t) *T*share);

		taylorRow = (uint64_t *) malloc(sizeof(uint64_t)*T);

		initialTaylor(p,N,T,0,rowMatrix);
		free(taylorRow);
		destroyPolynomial(p);

		long i,j;
		//now build the polynomials and call the fft on them

		uint64_t *transforms;
		polynomial *p1;

		uint64_t *colMatrix;
		if(flag==1 || flag==2) {
 			colMatrix = (uint64_t *) malloc(sizeof(uint64_t)*T*share);
			if(flag==1) transpose(rowMatrix,colMatrix,T,share,0,T-1,0,share-1,sizeof(uint64_t),0,16);
			else naiveTranspose(rowMatrix,colMatrix,T,share,sizeof(uint64_t),T,share);
			free(rowMatrix);
		}
		uint64_t *swap;

		if(flag==1 || flag==2) {
			for(i=0;i<share;i++) {
				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);
				for(j=0;j<T;j++) {
					if(colMatrix[i*T+j]!=0) {
						insertElement(p1,colMatrix[i*T+j],T-1-j);
					}
				}

				transforms = fourierTransform(p1,T,0,gf);			
				//then update the coefficients
				memcpy(&colMatrix[i*T],transforms,sizeof(uint64_t)*T);
				free(transforms);
				destroyPolynomial(p1);
			}

			rowMatrix = (uint64_t *) malloc(sizeof(uint64_t) *T*share);
			if(flag==1) transpose(colMatrix,rowMatrix,share,T,0,share-1,0,T-1,sizeof(uint64_t),0,16);
			else naiveTranspose(colMatrix,rowMatrix,share,T,sizeof(uint64_t),share,T);

		}	
		else {
			for(i=0;i<share;i++) {
				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);
				for(j=0;j<T;j++) {
					if(rowMatrix[j*share+i]!=0) {
						insertElement(p1,rowMatrix[j*share+i],T-1-j);
					}
				}
				transforms = fourierTransform(p1,T,0,gf);
			
				for(j=0;j<T;j++) {
					rowMatrix[j*share+i] = transforms[j];
				}
			
				free(transforms);
				//free the polynomial
				destroyPolynomial(p1);
			}
			colMatrix= (uint64_t *) malloc(sizeof(uint64_t)*T*share);
		}
		
		if(T==65536 && numProcesses<16) MPI_Buffer_attach(colMatrix,INT_MAX);
		else MPI_Buffer_attach(colMatrix,T*share);

                for(i=0;i<numProcesses;i++) {
                        if(i!=rank) {
                                for(j=0;j<share;j++) {
                                        MPI_Bsend(&rowMatrix[i*share*share+(j*share)],share, MPI_UINT64_T, i,(j*T) + (rank*share),MPI_COMM_WORLD);
                                }
                        }
                }

		swap = (uint64_t *) malloc(sizeof(uint64_t)*share*share);

		memcpy(swap,&rowMatrix[rank*share*share],sizeof(uint64_t)*share*share);

                for(i=0;i<share;i++) {
                        memcpy(&rowMatrix[i*T + (rank*share)], &swap[(i*share)],sizeof(uint64_t)*share);
                }

	
		free(swap);
	
                for(i=0;i<numProcesses;i++) {
                        if(i!=rank) {
                                for(j=0;j<share;j++) {
                                MPI_Recv(&rowMatrix[j*T + (i*share)],share, MPI_UINT64_T , MPI_ANY_SOURCE,(j*T)+(i*share), MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                }
                        }
                }

		int h;
		MPI_Buffer_detach(colMatrix,&h);
		free(colMatrix);
		
		for(i=0;i<share;i++) {
			p1 = (polynomial *) malloc(sizeof(polynomial));	
			initialize(p1);
			for(j=0;j<T;j++) {
				if(rowMatrix[(i*T)+j]) {
					insertElement(p1,rowMatrix[(i*T)+j],T-1-j);
				}
			}

			long st = ((rank*share)+i)*T;

			transforms = fourierTransform(p1,T,st,gf);

			memcpy(&rowMatrix[i*T],transforms,sizeof(uint64_t)*T);
			free(transforms);
			destroyPolynomial(p1);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(rowMatrix,T*share,MPI_UINT64_T,NULL,T*share,MPI_UINT64_T,0,MPI_COMM_WORLD);
		free(rowMatrix);

	}

	freeTable();
	freePols();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;

}


void initialTaylor(polynomial *p, long n, long t, long MPI_Index, uint64_t *mine) {
	if(n<=t) {
		long j,curr;
		curr = 0;
		for(j=0;j<t;j++) {
			taylorRow[j] = 0;
			if(curr<=p->curr) {
				if(p->terms[curr].degree == t-1-j) {
					taylorRow[j] = p->terms[curr].coefficient;	
					curr++;
				}
			}
	
		}
		memcpy(mine + (MPI_Index*share), &taylorRow[rank*share],sizeof(uint64_t)*share);
	}
	else {
		long k = 0;			//Step 1
		while(n > t * (1 <<(k+1)) ) {		
			k++;
		}

		long factor = (t << k);	//Step 2
	
		polynomial* f0 = (polynomial *) malloc(sizeof(polynomial));
		polynomial* f1 = (polynomial *) malloc(sizeof(polynomial));
		polynomial* f2 = (polynomial *) malloc(sizeof(polynomial));

		initialize(f0);
		initialize(f1);
		initialize(f2);

		long deg;
		int w = 0;

	//Can make this a binary search problem
		long i;
		for(i=0;i<=p->curr;i++) {
			if(p->terms[i].degree < factor) break;
		}

	//insert all elements starting from i to the right in f0

		long j;
		int count;
	
		for(j=i;j<=p->curr;j++) {
			insertElement(f0,p->terms[j].coefficient,p->terms[j].degree);
		}


		for(j=0;j<i;j++) {
			if(p->terms[j].degree < (factor << 1) - (1<<k) ) break;
		}
	//j is marker2

		long u;
		for(u=j;u<i;u++) {
			insertElement(f1,p->terms[u].coefficient,p->terms[u].degree - factor);
		}

		for(u=0;u<j;u++) {
			insertElement(f2,p->terms[u].coefficient,p->terms[u].degree - ( factor + ((t-1)<<k) ) );
		}

		polynomial* h1;
		polynomial* h2 = (polynomial *) malloc(sizeof(polynomial));
		polynomial* g0;
		polynomial* g1;

		initialize(h2);

		h1 = addPolynomials(f1,f2);

		for(i=0;i<=h1->curr;i++) {
			insertElement(h2,h1->terms[i].coefficient,h1->terms[i].degree + (1<<k) );
		}

		g0 = addPolynomials(f0,h2);

		for(i=0;i<=f2->curr;i++) {
			f2->terms[i].degree += (t-1) * (1<<(k));
		}

		g1 = addPolynomials(h1,f2);

		destroyPolynomial(f0);
		destroyPolynomial(f1);
		destroyPolynomial(f2);
		destroyPolynomial(h1);
		destroyPolynomial(h2);

		initialTaylor(g0,factor,t,2*MPI_Index,mine);
		initialTaylor(g1,n-factor,t,(2*MPI_Index) +1,mine);

		destroyPolynomial(g0);
		destroyPolynomial(g1);
	}
}


void taylorExpansion(polynomial *p, long n, long t,uint64_t *g, long rowIndex) {
	if(n<=t) {
		long i,curr;
		curr = 0;
		for(i=0;i<t;i++) {
			g[rowIndex*t + i] = 0;
			if(curr<=p->curr) {
				if(p->terms[curr].degree == t-1-i) {
					g[rowIndex*t + i] = p->terms[curr].coefficient;
					curr++;
				}
			}
		}
	}
	else {
		long k = 0;			//Step 1
		while(n > t * (1 <<(k+1)) ) {		
			k++;
		}

		long factor = (t << k);	//Step 2
	
		polynomial* f0 = (polynomial *) malloc(sizeof(polynomial));
		polynomial* f1 = (polynomial *) malloc(sizeof(polynomial));
		polynomial* f2 = (polynomial *) malloc(sizeof(polynomial));

		initialize(f0);
		initialize(f1);
		initialize(f2);

		long deg;
		int w = 0;

	//Can make this a binary search problem
		long i;
		for(i=0;i<=p->curr;i++) {
			if(p->terms[i].degree < factor) break;
		}

	//insert all elements starting from i to the right in f0

		long j;
		int count;
	
		for(j=i;j<=p->curr;j++) {
			insertElement(f0,p->terms[j].coefficient,p->terms[j].degree);
		}


		for(j=0;j<i;j++) {
			if(p->terms[j].degree < (factor << 1) - (1<<k) ) break;
		}
	//j is marker2

		long u;
		for(u=j;u<i;u++) {
			insertElement(f1,p->terms[u].coefficient,p->terms[u].degree - factor);
		}

		for(u=0;u<j;u++) {
			insertElement(f2,p->terms[u].coefficient,p->terms[u].degree - ( factor + ((t-1)<<k) ) );
		}

		polynomial* h1;
		polynomial* h2 = (polynomial *) malloc(sizeof(polynomial));
		polynomial* g0;
		polynomial* g1;

		initialize(h2);

		h1 = addPolynomials(f1,f2);

		for(i=0;i<=h1->curr;i++) {
			insertElement(h2,h1->terms[i].coefficient,h1->terms[i].degree + (1<<k) );
		}

		g0 = addPolynomials(f0,h2);

		for(i=0;i<=f2->curr;i++) {
			f2->terms[i].degree += (t-1) * (1<<(k));
		}

		g1 = addPolynomials(h1,f2);

		destroyPolynomial(f0);
		destroyPolynomial(f1);
		destroyPolynomial(f2);
		destroyPolynomial(h1);
		destroyPolynomial(h2);

		taylorExpansion(g0,factor,t,g,2*rowIndex);
		taylorExpansion(g1,n-factor,t,g,2*rowIndex +1);

		destroyPolynomial(g0);
		destroyPolynomial(g1);
	}
}

uint64_t* fourierTransform(polynomial* p, long n, long s, int field) {

	if(n==2) {
		uint64_t* ret = (uint64_t *) malloc(sizeof(uint64_t) *2);

		uint64_t first = WofI(p,2*s,field);
		uint64_t second = WofI(p,2*s+1,field);

		ret[0] = first;
		ret[1] = second;

		return ret;

	}

	long k=1;
	int m= 0;
// n= 2^m, get m
	while(k != n) {
		k=k<<1;
		m++;
	}
	long t = 1;
	m=m>>1;

	t=t<<m;
	uint64_t *g = (uint64_t *) malloc(sizeof(uint64_t)*n);

	taylorExpansion(p,n,t,g,0);
	long i,j;

	polynomial *p1;
	uint64_t * fft;

	if(flag==1 || flag==2) {
	
		uint64_t *g1 = (uint64_t *) malloc(sizeof(uint64_t)*n);

		if(flag==1) {
			if(t>=65536) transpose(g,g1,t,t,0,t-1,0,t-1,sizeof(uint64_t),0,16);
			else transpose(g,g1,t,t,0,t-1,0,t-1,sizeof(uint64_t),0,0);	
		}

		else naiveTranspose(g,g1,t,t,sizeof(uint64_t),t,t);

		for(i=0;i<t;i++) {

			p1 = (polynomial *) malloc(sizeof(polynomial));
			initialize(p1);
			for(j=0;j<t;j++) {
				if(g1[i*t+j] !=0) {
					insertElement(p1,g1[i*t+j],t-1-j);
				}
			}

			fft = fourierTransform(p1,t,s*t,field);

			memcpy(&g1[i*t],fft,sizeof(uint64_t)*t);

			free(fft);
			destroyPolynomial(p1);
		}

		if(flag==1) {
			if(t>=65536) transpose(g1,g,t,t,0,t-1,0,t-1,sizeof(uint64_t),0,16);
			else transpose(g1,g,t,t,0,t-1,0,t-1,sizeof(uint64_t),0,0);
		}
		else naiveTranspose(g1,g,t,t,sizeof(uint64_t),t,t);
		free(g1);
	}

	else {
		for(i=0;i<t;i++) {

			p1 = (polynomial *) malloc(sizeof(polynomial));
			initialize(p1);
			for(j=0;j<t;j++) {
				if(g[j*t+i] !=0) {
					insertElement(p1,g[j*t+i],t-1-j);
				}
			}

			fft = fourierTransform(p1,t,s*t,field);		
			for(j=0;j<t;j++) {
				g[j*t+i]= fft[j];
			}

			free(fft);
			destroyPolynomial(p1);
		}

	}

	for(i=0;i<t;i++) {

		p1 = (polynomial *) malloc(sizeof(polynomial));
		initialize(p1);

		for(j=0;j<t;j++) {
			if(g[i*t+j] !=0) {
				insertElement(p1,g[i*t+j],t-1-j);
			}
		}
		//now compute fft of p1;
		fft = fourierTransform(p1,t,(s*n) + (i*t),field);
		memcpy(&g[i*t],fft,sizeof(uint64_t)*t);
		free(fft);
		destroyPolynomial(p1);
	}

	return g;
}

///

uint64_t WofI(polynomial* p,long i, int m) {
	long x = i;
	uint64_t* table = getTable(m);

	uint64_t eval = 0;
	long j,k;

	for(j=0;j<m;j++) {
		if(x%2 ==1) {	
			eval = eval ^ table[j];
		}
		x=(x>>1);
	}

	uint64_t ans = 0;

	uint64_t prim = getPrimPol(m);

	uint64_t prod,y;

	for(j = 0;j<=p->curr;j++) {
	
		prod= 0;	

		if(p->terms[j].degree == 0) y = 1;
		else y = eval;

		if(y==1) ans = ans ^ p->terms[j].coefficient;

		else {
			int carry;
			uint64_t temp = p->terms[j].coefficient;
			for(k=(m-1);k>=0;k--) {
				if(getBit(&y,0)==1) prod = prod ^ temp;	//Step1
				y = (y>>1);												//Step2
				carry = getBit(&temp,(m-1));			//Step3
				clearBit(&temp,(m-1));				//Step4
				temp = (temp << 1);	//Step4
				if(carry==1) temp = temp ^ prim;
			}
			ans = ans ^ prod;
		}
	}	
	
//////
	return ans;

}


