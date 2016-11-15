#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include "DataStruct.h"
#include "tables.h"
#include "time.hpp"


void taylorExpansion(polynomial *p,long n, long t, uint64_t *g, long rowIndex);
uint64_t* fourierTransform(polynomial* p, long n, long s, int field);
uint64_t WofI(polynomial * p,long i, int m);


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
int limit,flag;
//Assumption: All polynomials are represented in descending order of degrees
int main(int argc, char **argv) {

	buildTable();
	buildPols();

	limit = 23100;
	flag = atoi(argv[1]);

	polynomial *p = (polynomial *) malloc(sizeof(polynomial));
	initialize(p);
	int c = -5;
	long deg;	

	int counter = 1;

	uint64_t theCoff;

	FILE *fp = fopen("input.txt","r");
	long n;
	int gf;

	gf = 32;
	n = 4294967296;


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

			fscanf(fp,"%ld",&deg);
		
			insertElement(p,theCoff,deg);
		}
	}
	fclose(fp);	



	Time timer;
	timer.start();
	uint64_t *l = fourierTransform(p,n,0,gf);

	timer.stop();
	timer.printTime();

	destroyPolynomial(p);

	free(l);

	freeTable();
	freePols();

	return 0;

}

void taylorExpansion(polynomial *p, long n, long t,uint64_t *g, long rowIndex) {
	if(n<=t) {
		if(t==65536) {
		//	printf("Tay %ld\n",rowIndex);
			long i,curr;
			curr = 0;

			if(rowIndex<limit) {
				for(i=0;i<limit;i++) {
					g[rowIndex*limit + i] = 0;
					if(curr<=p->curr) {
						if(p->terms[curr].degree == t-1-i) {
							g[rowIndex*limit + i] = p->terms[curr].coefficient;
							curr++;
						}
					}
				}
			}
		}
		else {
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
	uint64_t *g;
	if(t==65536) {
		g = (uint64_t *) malloc(sizeof(uint64_t)*(limit*limit));
	}
	else {
		g = (uint64_t *) malloc(sizeof(uint64_t)*n);
	}
	taylorExpansion(p,n,t,g,0);

	long i,j,f;


	uint64_t *g1;
	polynomial *p1;
	if(flag==1 || flag==2) {
		if(t==65536) { 
			g1 = (uint64_t *) malloc(sizeof(uint64_t)*limit*limit);
			if(flag==1) transpose(g,g1,limit,limit,0,limit-1,0,limit-1,sizeof(uint64_t),0,16);
			else naiveTranspose(&g,g1,limit,limit,sizeof(uint64_t),limit,limit);
			free(g);
		}

		else {
			g1 = (uint64_t *) malloc(sizeof(uint64_t)*n);
			naiveTranspose(g,g1,t,t,sizeof(uint64_t),t,t);
		}


		if(t==65536) {

			for(i=0;i<limit;i++) {
				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);

				for(j=0;j<limit;j++) {
					if(g1[i*limit+j] !=0) {
						insertElement(p1,g1[i*limit+j],t-1-j);
					}
				}

				uint64_t *fft = fourierTransform(p1,t,s*t,field);

				memcpy(&g1[i*limit],fft,sizeof(uint64_t)*limit);

				free(fft);
				destroyPolynomial(p1);
			}


		}
		else {

			for(i=0;i<t;i++) {

				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);
				for(j=0;j<t;j++) {
					if(g1[i*t+j] !=0) {
						insertElement(p1,g1[i*t+j],t-1-j);
					}
				}

				uint64_t *fft = fourierTransform(p1,t,s*t,field);
				memcpy(&g1[i*t],fft,sizeof(uint64_t)*t);

				free(fft);
				destroyPolynomial(p1);
			}

		}

		if(t==65536) { 
			g = (uint64_t *) malloc(sizeof(uint64_t)*(limit*limit));
			if(flag==1) transpose(g1,g,limit,limit,0,limit-1,0,limit-1,sizeof(uint64_t),0,16);
			else naiveTranspose(g1,g,limit,limit,sizeof(uint64_t),limit,limit);
		}
		else naiveTranspose(g1,g,t,t,sizeof(uint64_t),t,t);
		free(g1);

	}

	else {
		if(t==65536) {
			for(i=0;i<limit;i++) {

				p1 = (polynomial *) malloc(sizeof(polynomial));
				initialize(p1);

				for(j=0;j<limit;j++) {
					if(g[j*limit+i] !=0) {
						insertElement(p1,g[j*limit+i],t-1-j);
					}
				}

				uint64_t *fft = fourierTransform(p1,t,s*t,field);

				for(j=0;j<limit;j++) {
					g[j*limit+i] = fft[j];
				}

				free(fft);
				destroyPolynomial(p1);
			}

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

				uint64_t *fft = fourierTransform(p1,t,s*t,field);
				for(j=0;j<t;j++) {
					g[j*t+i] = fft[j];
				}

				free(fft);
				destroyPolynomial(p1);
			}

		}

	}

	if(t==65536) {

		for(i=0;i<limit;i++) {
			p1 = (polynomial *) malloc(sizeof(polynomial));
			initialize(p1);

			for(j=0;j<limit;j++) {
				if(g[i*limit+j] !=0) {
					insertElement(p1,g[i*limit+j],t-1-j);
				}
			}

			//now compute fft of p1;
			uint64_t* fft = fourierTransform(p1,t,(s*n) + (i*t),field);

			memcpy(&g[i*limit],fft,sizeof(uint64_t)*limit);
			free(fft);
			destroyPolynomial(p1);
		}

	}

	else {
		for(i=0;i<t;i++) {

			p1 = (polynomial *) malloc(sizeof(polynomial));
			initialize(p1);

			for(j=0;j<t;j++) {
				if(g[i*t+j] !=0) {
					insertElement(p1,g[i*t+j],t-1-j);
				}
			}

			//now compute fft of p1;
			uint64_t* fft = fourierTransform(p1,t,(s*n) + (i*t),field);
			memcpy(&g[i*t],fft,sizeof(uint64_t)*t);
			free(fft);
			destroyPolynomial(p1);
		}
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

