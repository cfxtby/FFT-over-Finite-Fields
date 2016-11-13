#include <stdint.h>
//Assuming a GF of maximum size 64


typedef struct primitivePolynomial {
	int* terms;
	int size;
} primitivePolynomial;

primitivePolynomial* getPrimitive() {
	primitivePolynomial *p = (primitivePolynomial *) malloc(sizeof(primitivePolynomial));
	p->size = 0;
	return p;
}

void insertPrimitiveElement(primitivePolynomial *p, int d) {
	if(p->size==0) {
		p->terms = (int *)malloc(sizeof(int));
		p->size =1;
		p->terms[0] = d;
	}
	else {
		p->terms = (int *) realloc(p->terms,sizeof(int) * (++p->size));
		p->terms[p->size-1] = d;
	}

}

typedef struct pol {
	int *terms;
	long curr;
	long size;
} pol;

void init(pol *p) {
	p->curr = -1;
	p->size = 1;
	p->terms = (int *) malloc(sizeof(int));
}

void destroyPol(pol *p) {
	free(p->terms);
	free(p);
}
//num is how many elements in a
void insertTerm(pol *p, int num, int *a) {
	while(p->curr + num >= p->size) {
		p->size = (p->size +32);
		int *tmp =  (int *)realloc(p->terms,sizeof(int) * p->size);
		if(tmp) {
			p->terms = tmp;
		}
		else printf("realloc error\n");
	}


	long i;
	for(i=p->curr+1;i<=p->curr + num;i++) {
		p->terms[i] = a[i-(p->curr+1)];
	}
	p->curr = p->curr + num;
}


void filePrintPol(pol *p,FILE *fp) {

	long i;
	long j;

	int n = 0;
	long pDeg;

	for(i=0;i<=p->curr;i+=n) {
		n = p->terms[i]+1;
		pDeg = p->terms[i]+i;
		fprintf(fp,"(");
		for(j=i+1;j<pDeg;j++) {
			fprintf(fp,"B^%d",p->terms[j]);
			if(j!=p->terms[i]+i -1) fprintf(fp," + ");
		}
		fprintf(fp,")x^%d",p->terms[pDeg]);
		if(i+n <=p->curr) fprintf(fp," + ");
	}
}

void printPol(pol *p) {

	long i;
	long j;

	int n = 0;
	long pDeg;

	for(i=0;i<=p->curr;i+=n) {
		n = p->terms[i]+1;
		pDeg = p->terms[i]+i;
		printf("(");
		for(j=i+1;j<pDeg;j++) {
			printf("B^%d",p->terms[j]);
			if(j!=p->terms[i]+i -1) printf (" + ");
		}
		printf(")x^%d",p->terms[pDeg]);
		if(i+n <=p->curr) printf(" + ");
	}
}


pol* addPols(pol *p, pol *q) {

	pol *ans = (pol *) malloc(sizeof(pol));
	init(ans);


	long i=0;
	long j=0;
	int t;
	int n;
	int *nums = (int *) malloc(sizeof(int));
	long pDeg;
	long qDeg;

	while(i<=p->curr && j<=q->curr) {
//make i and j jump to the start of the next term
		pDeg = p->terms[i] +i;
		qDeg = q->terms[j] +j;
		if(p->terms[pDeg] > q->terms[qDeg]) {
			n = p->terms[i]+1;
			nums = (int *)realloc(nums,sizeof(int)*n);
			for(t=0;t<n;t++) {
				nums[t] = p->terms[i+t];
			}
			i=i+n;
			insertTerm(ans,n,nums);
		}
	
		else if(p->terms[pDeg] < q->terms[qDeg]) {
			n = q->terms[j]+1;
			nums = (int *)realloc(nums,sizeof(int)*n);
			for(t=0;t<n;t++) {
				nums[t] = q->terms[j+t];
			}
			j=j+n;
			insertTerm(ans,n,nums);
		}

		else {
//equal degrees, must merge the coefficients
			int k = p->terms[i]+1;
			int l = q->terms[j]+1;
			int counter = 1;
			int size = 1;

			long t1 = j+1;
			long t0 = i+1;

			while(t0<pDeg && t1<qDeg) {
				if(counter==size) {
					size++;
					nums = (int *)realloc(nums,(sizeof(int)*size));
				}
				if(p->terms[t0] > q->terms[t1]) {
					nums[counter] = p->terms[t0];
					t0++;
					counter++;
				}
				else if(p->terms[t0] < q->terms[t1]) {
					nums[counter] = q->terms[t1];
					t1++;
					counter++;	
				}
				else {
					t0++;
					t1++;
				}
			}
			while(t0<pDeg) {
				if(counter==size) {
					size++;
					nums = (int *)realloc(nums,(sizeof(int)*size));
				}
				nums[counter] = p->terms[t0];
				counter++;
				t0++;
			}
			while(t1<qDeg) {
				if(counter==size) {
					size++;
					nums = (int *)realloc(nums,(sizeof(int)*size));
				}
				nums[counter] = q->terms[t1];
				counter++;
				t1++;
			}
			if(counter==size) {
				size++;
				nums = (int *)realloc(nums,(sizeof(int)*size));
			}
		
			if(counter!=1) {
				nums[counter] = p->terms[pDeg];
				n = counter+1;
				nums = (int *)realloc(nums,sizeof(int)*n);
				nums[0] = n-1;

				insertTerm(ans,n,nums);
			}	
			i = i+k;
			j = j+l;

		}
	}

	while(i<=p->curr) {
		n = p->terms[i]+1;
		nums = (int *)realloc(nums,sizeof(int)*n);
		for(t=0;t<n;t++) {
			nums[t] = p->terms[i+t];
		}

		insertTerm(ans,n,nums);
		i=i+n;
	}
	while(j<=q->curr) {
		n = q->terms[j]+1;
		nums = (int *)realloc(nums,sizeof(int)*n);
		for(t=0;t<n;t++) {
			nums[t] = q->terms[j+t];
		}

		insertTerm(ans,n,nums);
		j=j+n;
	}
	free(nums);
	
	return ans;

}

typedef struct GenPolTerm {
	int coefficient;
	int power;
} GenPolTerm;

typedef struct GenPolList {
	long size;
	GenPolTerm *terms;
} GenPolList;

void printGen(GenPolList *g) {
	long i;
	for(i=0;i<g->size-1;i++) {
		printf("%d B^%d + ",g->terms[i].coefficient,g->terms[i].power);
	}
	printf("%d B^%d\n",g->terms[i].coefficient,g->terms[i].power);

}

void add(GenPolList *g, GenPolTerm *t) {
	if(g->size==0) {
		g->terms = (GenPolTerm *) malloc(sizeof(GenPolTerm));
		g->size++;
	}
	else {
		g->size++;
		g->terms = (GenPolTerm *)realloc(g->terms,sizeof(GenPolTerm)*g->size);
	}
	g->terms[g->size-1].coefficient = t->coefficient;
	g->terms[g->size-1].power = t->power;
}

typedef struct ArrayList {
	long size;
	int *array;
}ArrayList;

void initializeList(ArrayList *a) {
	a->array = (int *) malloc(sizeof(int));
	a->size = 0;
}

void addNum(ArrayList *a,int n) {
	a->size++;
	a->array = (int *)realloc(a->array,sizeof(int)*a->size);
	
	a->array[a->size-1] = n;

}

int* getArray(ArrayList *a) {

	int *ret = (int *)malloc(sizeof(int)*a->size);
	long i;
	for(i=0;i<a->size;i++) {
		ret[i] = a->array[i];
	}
	return ret;

}

typedef struct PolyList2 {
	pol *p;
	long n;
} PolyList2;

void destroyPolList(PolyList2 *list) {
	long i;
	for(i=0;i<list->n;i++) {
		free(list->p[i].terms);
	}
	free(list->p);
	free(list);

}

void push(PolyList2 *list, long i, pol *p) {

	long j;
	int t;
	int n = 0;
	int *nums = (int *) malloc(sizeof(int));
	for(j=0;j<=p->curr;j+=n) {
		n = p->terms[j]+1;
		nums = (int *)realloc(nums,sizeof(int) * n);
		for(t=0;t<n;t++) {
			nums[t] = p->terms[j+t];
		}
		insertTerm(&list->p[i],n,nums);
	}
	free(nums);
	
}

PolyList2* concatenate2(PolyList2 *p1, PolyList2 *p2) {

	PolyList2 *q = (PolyList2 *) malloc(sizeof(PolyList2));
	q->n = p1->n + p2->n;
	q->p = (pol *) malloc(sizeof(pol) * q->n);

	long i;
	long j;
	int t;
	int n = 0;

	for(i=0;i<p1->n;i++) {
		init(&q->p[i]);
		push(q,i,&p1->p[i]);
	} 
	for(i=p1->n;i<q->n;i++) {
		init(&q->p[i]);
		push(q,i,&p2->p[i-p1->n]);
	} 

	return q;

}


unsigned int getBit(uint64_t *l, int i) {
//if i>63, error
	unsigned int bit = (*l >> i) & 1;
	return bit;
}

void setBit(uint64_t *l, int i) {
//if i>63, error
	*l = *l | (((uint64_t) 1) << i);
}

void clearBit(uint64_t *l, int i) {
	*l = *l & ~(((uint64_t)1) << i);
}

void toggleBit(uint64_t *l, int i) {
	*l = *l ^ (((uint64_t) 1) << i);
}

void filePrintBeta(FILE *o, uint64_t *b) {

	int i;
	int k = 65;
	for(i=63;i>0;i--) {
		if(getBit (b, i)==1 ) {
			if(k==65) {
				fprintf(o,"B^%d",i);
				k = 3;
			}
			else {
				fprintf(o," + B^%d",i);
			}
		}
	}
	if (getBit (b, i)==1){
		if(k==3) fprintf(o," + B^0");
		else {
			fprintf(o,"B^0");
			k=3;
		}
	}
	if(k==65) fprintf(o,"0");

}

void printBeta(uint64_t *b) {

	int i;
	int k = 65;
	for(i=63;i>0;i--) {
		if(getBit (b, i)==1 ) {
			if(k==65) {
				printf("B^%d ",i);
				k = 3;
			}
			else {
				printf ("+ B^%d ",i);
			}
		}
	}
	if (getBit (b, i)==1) printf (" + 1");

}


typedef struct element {
	uint64_t coefficient;
	long degree;
} element;

//when you initialize, size = 1, curr = -1
//Maybe implement an initializer function?
typedef struct polynomial {
	long size;
	long curr;
	element *terms;
} polynomial;

void initialize(polynomial *p) {
	p->size = 1;
	p->curr = -1;
	p->terms = (element *) malloc(sizeof(element));
}

void destroyPolynomial(polynomial *p) {
	free(p->terms);
	free(p);
}

void printPolynomial(polynomial *p) {

	if(p->curr < 0) {
		printf("Polynomial is empty\n");
		return;
	}
	
	long i;
	for(i=0;i<=p->curr;i++) {
		
		printf("(");
		printBeta(&(p->terms[i].coefficient));
		printf(")");
		printf("x^%ld",p->terms[i].degree);
		
		if(i!=p->curr) printf(" + ");
	}

}


void filePrintPolynomials(FILE *o, polynomial *p) {

	if(p->curr < 0) {
		fprintf(o,"Polynomial is empty\n");
		return;
	}
	
	long i;
	for(i=0;i<=p->curr;i++) {
		
		fprintf(o,"(");
		filePrintBeta(o,&(p->terms[i].coefficient));
		fprintf(o,")");
		fprintf(o,"x^%ld",p->terms[i].degree);
		
		if(i!=p->curr) fprintf(o," + ");
	}
}


void insertElement(polynomial *p, uint64_t c, long d) {
	if(p->curr == p->size-1) {
//use realloc to double the size of array;
		p->size = (p->size + 32);		//p->size << 1
		p->terms = (element *) realloc(p->terms,sizeof(element) * p->size);

	};

	p->curr = p->curr + 1;

	p->terms[p->curr].coefficient = c;
	p->terms[p->curr].degree = d;

}

uint64_t addBetas1(uint64_t c1, uint64_t c2) {
	return c1^c2;
}

polynomial* addPolynomials(polynomial *p, polynomial *q) {

	polynomial *n = (polynomial *) malloc(sizeof(polynomial));
	
	initialize(n);
	
	uint64_t c;
	long d;


	long i=0;
	long j=0;
	while(i<=p->curr && j<=q->curr) {
		if(p->terms[i].degree > q->terms[j].degree) {
			c = p->terms[i].coefficient;
			d = p->terms[i].degree;
			i++;
		}
		else if(q->terms[j].degree > p->terms[i].degree) {
			c = q->terms[j].coefficient;
			d = q->terms[j].degree;
			j++;
		}
		else {
			c = p->terms[i].coefficient ^ q->terms[j].coefficient;
			d = p->terms[i].degree;
			i++;
			j++;
		}
		if(c!=0) insertElement(n,c,d);
	}

	while(i<=p->curr) {
		insertElement(n,p->terms[i].coefficient,p->terms[i].degree);
		i++;
	}
	while(j<=q->curr) {
		insertElement(n,q->terms[j].coefficient,q->terms[j].degree);
		j++;
	}
	return n;

}
typedef struct PolyList {
	polynomial *p;
	long n;
} PolyList;

void destroyList(PolyList *list) {
	long i;
	for(i=0;i<list->n;i++) {
		free(list->p[i].terms);
	}
	free(list->p);
	free(list);
}

PolyList* concatenate(PolyList *p1, PolyList *p2) {

	PolyList *q = (PolyList *) malloc(sizeof(PolyList));
	q->n = p1->n + p2->n;
	q->p = (polynomial *) malloc(sizeof(polynomial) * q->n);

	long i;
	long j;
//fill q with all elements of p1
	for(i=0;i<p1->n;i++) {
		initialize(&(q->p[i]));
		for(j=0;j<=p1->p[i].curr;j++) {
			insertElement(&q->p[i],p1->p[i].terms[j].coefficient,p1->p[i].terms[j].degree);
		}
	}


	for(i=p1->n;i< q->n;i++) {
		initialize(&(q->p[i]));

		for(j=0;j<=p2->p[i-p1->n].curr;j++) {
			insertElement(&q->p[i],p2->p[i-p1->n].terms[j].coefficient,p2->p[i-p1->n].terms[j].degree);
		}
	}
	
	return q;
}

uint64_t* listConcatenate(uint64_t* a, uint64_t* b) {

	uint64_t *g = (uint64_t *) malloc(sizeof(uint64_t) * (a[0]+b[0]+1));

	g[0] = a[0]+b[0];

	long i;
//fill q with all elements of p1
	for(i=1;i<=a[0];i++) {
		g[i] = a[i];
	}


	for(i=a[0]+1;i<= g[0];i++) {
		g[i] = b[i-a[0]];
	}

	return g;


}
