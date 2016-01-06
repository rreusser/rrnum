#include "sparse.h"
int assert(int i, int n);


rr_sparse_mat* rr_sparse_mat_alloc(int n) {
    //fprintf(stderr,"Allocating sparse_mat\n");
    int i;
    rr_sparse_mat* A;
    A = malloc(sizeof(rr_sparse_mat));
    A->ind = (int*)NULL;
    A->permutation = (int*)NULL;
    A->length = (int*)NULL;
    A->dat = (double*)NULL;
    A->diagonal = (double*)NULL;
    A->scale = (double*)NULL;
    A->row = (node**)NULL;
    A->cind = (int*)NULL;
    A->cdat = (int*)NULL;
    A->col_indexed = 0;
    A->dim = n;
    A->row = malloc(n*sizeof(node*));
    A->length = malloc(n*sizeof(int));
    A->diagonal = malloc(n*sizeof(double));
    A->scale = malloc(n*sizeof(double));
    A->ind=NULL;
    A->dat=NULL;
    for(i=0; i<n; i++) {
	A->row[assert(i,A->dim)] = NULL;
	A->length[assert(i,A->dim)] = 0;
	A->diagonal[assert(i,A->dim)] = 0.0;
    }
    return A;
}
void rr_sparse_mat_dealloc(rr_sparse_mat* A) {
    rr_sparse_mat_unlink( A );
    if(A->dat!=NULL) {free(A->dat);}
    //for(i=0; i<A->dim; i++) {
	//printf("A->ind[%i] = %i\n",i,A->ind[i]);
    //}
    if(A->permutation!=NULL) {free(A->permutation);}
    if(A->ind!=NULL) {free(A->ind);}
    if(A->length!=NULL) {free(A->length);}
    if(A->cind!=NULL) {free(A->cind);}
    if(A->cdat!=NULL) {free(A->cdat);}
    if(A->diagonal!=NULL) {free(A->diagonal);}
    if(A->scale!=NULL) {free(A->scale);}
    if(A->row!=NULL) {free(A->row);}
    free(A);
}
void rr_sparse_mat_build(rr_sparse_mat* A) {
    int c, i;
    //get the total size of the array
    A->size=A->dim+1;
    for(i=0; i<A->dim; i++) {
	A->size+=A->length[i];
    }
    //allocate the arrays
    A->ind = malloc(sizeof(int)*A->size);
    A->dat = malloc(sizeof(double)*A->size);
    //transfer the diagonal elements
    for(i=0; i<A->dim; i++) {
	A->dat[i]=A->diagonal[i];
    }
    //transfer the rest, row by row
    c=A->dim+1;
    node* current;
    for(i=0; i<A->dim; i++) {
	A->ind[i] = c;
	current=A->row[i];
	while(current!=NULL) {
	    A->ind[c]=current->index;
	    A->dat[c]=current->data;
	    current=current->next;
	    c++;
	}
    }
    A->ind[A->dim]=c;

    rr_sparse_mat_unlink(A);
}
void rr_sparse_mat_insert(rr_sparse_mat* A, int i, int j, double data) {
    node* last=NULL;
    if(i>A->dim-1||j>A->dim-1) {printf("sparse_insert: error: out of bounds: %i,%i vs %i\n",i,j,A->dim);}
    if(i==j) {
	A->diagonal[assert(i,A->dim)] = data;
    } else {
	node* current = A->row[assert(i,A->dim)];
	if(current==NULL) {
	    current=A->row[assert(i,A->dim)]=malloc(sizeof(node));
	    current->next=NULL;
	    current->data=data;
	    current->index=j;
	    A->length[assert(i,A->dim)]++;
	} else if(j<current->index) {
	    A->row[assert(i,A->dim)] = malloc(sizeof(node));
	    A->row[assert(i,A->dim)]->next = current;
	    A->row[assert(i,A->dim)]->data = data;
	    A->row[assert(i,A->dim)]->index = j;
	    A->length[assert(i,A->dim)]++;
	} else {
	    int repeat=0;
	    while(current != NULL) {
		if(j==current->index) {
		    repeat = 1;
		    break;
		}
		if(j<current->index) break;
		last=current;
		current=current->next;
	    }
	    if(repeat) {
		current->data = data;
	    } else {
		last->next = malloc(sizeof(node));
		last->next->next = current;
		last->next->data = data;
		last->next->index = j;
		A->length[assert(i,A->dim)]++;
	    }
	}
    }
}
void rr_sparse_mat_calculate_scale(rr_sparse_mat* A) {
    int i;
    double *scale = A->scale;
    node* current;
    if(A->row==NULL) {
	fprintf(stderr,"Error: must pivot after inserting elements and before building matrix.\n");
    } else {
	for(i=0;i<A->dim;i++) {
	    scale[i] = ABS(A->diagonal[i]);
	    current = A->row[i];
	    while(current!=NULL) {
		scale[i] = MAX(scale[i],ABS(current->data));
		current = current->next;
	    }
	}
    }
}

// THIS IS BASICALLY O(n^2) WHICH MAKES IT USELESS NEXT
// TO AN O(n*log n) ALGORITHM... HMM...
void rr_sparse_mat_partial_pivot( rr_sparse_mat* A ) {
    int i,j,k, imax;
    double dmax, dcomp, sdcomp=0.0;
    node* current, *last;
    int rows_switched=0;
    if(A->row==NULL) {
	fprintf(stderr,"Error: must pivot after inserting elements and before building matrix.\n");
    } else {
	if(A->permutation==NULL) {
	    A->permutation = malloc(sizeof(int)*A->dim);
	    for(i=0;i<A->dim;i++) {
		A->permutation[i] = i;
	    }
	}
	for(i=0;i<A->dim;i++) {
	    // Look for a pivot element for this element:
	    dmax=A->diagonal[i]/A->scale[i];
	    imax=i;
	    for(j=i+1;j<A->dim;j++) {
		current = A->row[j];
		while(current!=NULL) {
		    if(current->index==i) {
			dcomp=ABS(current->data/A->scale[j]);
			if(dcomp > dmax) {
			    imax=j;
			    dmax=dcomp;
			    sdcomp=current->data;
			}
		    }
		    current = current->next;
		}
	    }
	    if(imax!=i) {
		// Make the switch:
		//fprintf(stderr,"We should switch rows %i and %i\n",i,imax);
		rows_switched++;
		j=A->permutation[i];
		A->permutation[i]=A->permutation[imax];
		A->permutation[imax]=j;
		current=A->row[i];
		A->row[i]=A->row[imax];
		A->row[imax]=current;
		// Now sort out the diagonal:
		// insert the old diagonal elements into the rows
		if(A->diagonal[i] != 0.0) {
		    rr_sparse_mat_insert( A, imax, i, A->diagonal[i] );
		}
		if(A->diagonal[imax] != 0.0) {
		    rr_sparse_mat_insert( A, i, imax, A->diagonal[imax] );
		}
		// take the old off-diagonal elements and place them on the diag:
		A->diagonal[i]=sdcomp;
		A->diagonal[imax]=0.0;
		current=A->row[imax];
		while(current!=NULL) {
		    if(current->index==imax) {
			A->diagonal[imax]=current->data;
		    }
		    current=current->next;
		}
		// newly-diagonal elements from the rows:
		current=A->row[i];
		last=NULL;
		while(current!=NULL) {
		    if(current->index==i) {
			if(last==NULL) {
			    A->row[i]=current->next;
			} else {
			    last->next=current->next;
			}
			free(current);
			break;
		    }
		    last=current;
		    current=current->next;
		}
		current=A->row[imax];
		last=NULL;
		while(current!=NULL) {
		    if(current->index==imax) {
			if(last==NULL) {
			    A->row[i]=current->next;
			} else {
			    last->next=current->next;
			}
			free(current);
			break;
		    }
		    last=current;
		    current=current->next;
		}
	    }
	}
	fprintf(stderr,"Pivot: swapped %i rows\n",rows_switched);
    }
}

void rr_sparse_mat_unlink(rr_sparse_mat* A) {
    int i;
    node* next;
    if(A->row != NULL) {
	for(i=0; i<A->dim; i++) {
	    while(A->row[assert(i,A->dim)] != NULL) {
		next = A->row[assert(i,A->dim)]->next;
		free(A->row[assert(i,A->dim)]);
		A->row[assert(i,A->dim)] = next;
	    }
	}
	if(A->row!=NULL) free(A->row);
	if(A->diagonal!=NULL) free(A->diagonal);
	if(A->length!=NULL) free(A->length);
	A->row=NULL;
	A->diagonal=NULL;
	A->length=NULL;
    }
}

int rr_sparse_mat_is_col_indexed(rr_sparse_mat* A) {
    return A->col_indexed;
}
void rr_sparse_mat_copy(rr_sparse_mat* B, rr_sparse_mat* A) { //A=source, B=destination
    int i;
    B->size = A->size;
    B->dim = A->dim;
    if(B->dat!=NULL) free(B->dat);
    if(B->ind!=NULL) free(B->ind);
    B->ind = malloc(sizeof(int)*A->size);
    B->dat = malloc(sizeof(double)*A->size);
    for(i=0; i<A->size; i++) {
	B->ind[assert(i,B->size)] = A->ind[assert(i,A->size)];
	B->dat[assert(i,B->size)] = A->dat[assert(i,A->size)];
    }
}
void rr_sparse_mat_copyT(rr_sparse_mat* B, rr_sparse_mat* A) {
    int i;
    B->size = A->size;
    B->dim = A->dim;
    if(!(rr_sparse_mat_is_col_indexed(A))) {
	rr_sparse_mat_col_index(A);
    } 
    for(i=0; i<A->size;i++) {
	B->ind[i] = A->cind[assert(i,A->size)];
	B->dat[i] = A->dat[assert(A->cdat[assert(i,A->size)],A->size)];
    }
}
int rr_sparse_mat_nonzeros(rr_sparse_mat* A) {
    return A->ind[A->dim]-1;
}

void rr_sparse_mat_output(rr_sparse_mat* A) {
    int i;
    //for(i=0; i<A->dim; i++) {
	//printf("A->diagonal[%i] = %f\n",i,A->diagonal[i]);
    //}
    printf("A->size = %i\n",A->size);
    for(i=0; i<A->size; i++) {
	printf("%i\t%i\t%f\n",i,A->ind[i],A->dat[i]);
    }
}
/*void vector_mult(rr_sparse_mat* A, double* x, double* b) {
    int i,k;
    for(i=0; i<A->dim; i++) {
	b[i]=A->dat[assert(i,A->size)]*x[i];
    }
    for(i=0; i<A->dim; i++) {
	for(k=A->ind[assert(i,A->size)]; k<A->ind[assert(i+1,A->size)]; k++) {
	    b[assert(i,A->dim)] += A->dat[k]*x[assert(A->ind[assert(k,A->size)],A->dim)];
	}
    }
}*/
void vector_mult(rr_sparse_mat* A, double* x, double* b) {
    int i,k;
    for(i=0; i<A->dim; i++) {
	b[i]=A->dat[i]*x[i];
    }
    for(i=0; i<A->dim; i++) {
	for(k=A->ind[i]; k<A->ind[i+1]; k++) {
	    b[i] += A->dat[k]*x[A->ind[k]];
	}
    }
}
void vector_multT(rr_sparse_mat* A, double* x, double* b) {
    int i,k;
    if(!rr_sparse_mat_is_col_indexed(A)) rr_sparse_mat_col_index(A);
    for(i=0;i<A->dim;i++) {
	b[i]=A->dat[i]*x[i];
	for(k=A->cind[i];k<A->cind[i+1];k++) {
	    b[i]+=A->dat[A->cdat[k]]*x[A->cind[k]];
	}
    }
}
void rr_sparse_mat_cleanup( rr_sparse_mat* A) {
    if(rr_sparse_mat_is_col_indexed(A)) {
	free(A->cind);
	A->cind = NULL;
	free(A->cdat);
	A->cdat = NULL;
	A->col_indexed=0;
    }
}
void rr_sparse_mat_col_output(rr_sparse_mat* A) {
    fprintf(stderr,"Column index output:\n");
    int i;
    if(rr_sparse_mat_is_col_indexed(A) == 1) {
	for(i=0; i<A->size; i++) {
	    printf("%i\t%i\t%i\n",i,A->cind[i],A->cdat[i]);
	}
    } else {
	printf("Error: not column indexed.\n");
    }
}
void rr_sparse_mat_ind_output(rr_sparse_mat* A) {
    int i,j;
    for(i=0;i<A->dim;i++) {
	printf("(%i,%i) = %f\n",i,i,A->dat[i]);
    }
    for(i=0; i<A->dim; i++) {
	for(j=A->ind[i]; j<A->ind[i+1];j++) {
	    printf("(%i,%i) = %f\n",i,A->ind[j],A->dat[j]);
	}
    }
}
void rr_sparse_mat_col_index(rr_sparse_mat* A) {
    int h, *ch, icol, i,k;
    A->col_indexed = 1;
    A->cind = malloc(A->size*sizeof(int));
    A->cdat = malloc(A->size*sizeof(int));
    for(i=0;i<A->dim+1;i++) {
	assert(i,A->size);
	A->cind[i]=A->dim+1;
	A->cdat[i]=i;
    }
    for(i=A->dim+1;i<A->size;i++) {
	A->cind[i]=A->cdat[i]=-1;
	//fprintf(stderr,"A->cind[%i] = A->cdat[%i] = -1\n",i,i);
    }
    ch = calloc(A->dim,sizeof(int));
    //fprintf(stderr,"i from %i to %i\n",A->ind[0],A->ind[A->dim]);
    for(i=A->ind[0];i<A->ind[A->dim];i++) {
	//fprintf(stderr,"ch[A->ind[%i]] = ch[%i] ++\n",i,A->ind[i]);
	ch[A->ind[i]]++;
    }
    for(i=1;i<A->dim+1;i++) {
	//fprintf(stderr,"A->cind[%i] = %i+%i\n",i,A->ind[i-1],ch[i-1]);
	A->cind[i]=A->cind[i-1]+ch[i-1];
	ch[i-1]=0;
    }
    ch[A->dim-1]=0;
    
    for(i=0;i<A->dim;i++) {
	A->cdat[i]=i;
	for(k=A->ind[i];k<A->ind[i+1];k++) {
	    icol=A->ind[k];
	    h=A->cind[icol]+ch[icol];
	    assert(h,A->size);
	    A->cind[h]=i;
	    A->cdat[h]=k;
	    ch[icol]++;
	}
    }
    free(ch);
}
void rr_sparse_mat_ppm_output(rr_sparse_mat* A, char fname[], float scale) {
    int i,j,k,l;
    int ni = (int)(A->dim*scale);
    unsigned char* mat = calloc(ni*ni,sizeof(unsigned char));
    if(mat==NULL) {
	fprintf(stderr,"Error allocating image array.  Too big?\n");
    } else {
	for(i=0;i<A->dim;i++) {
	    l = (int)((A->dim-i-1)*scale);
	    for(k=A->ind[i]; k<A->ind[i+1]; k++) {
		j = (int)(A->ind[k]*scale);
		mat[j*ni+l] = 1;
	    }
	    j = (int)(i*scale);
	    mat[j*ni+l] = 2;
	}
	FILE* fptr = fopen(fname,"w");
	if(fptr!=NULL) {
	    fprintf(fptr,"P3\n");
	    fprintf(fptr,"# rrnum output\n");
	    fprintf(fptr,"%i %i\n",ni,ni);
	    fprintf(fptr,"%i\n",255);
	    for(j=ni-1;j>=0;j--) {
		for(i=0;i<ni;i++) {
		    switch(mat[j+ni*i]) {
			case 0: fprintf(fptr,"255 255 255\n"); break;
			case 1: fprintf(fptr,"0 0 0\n"); break;
			case 2: fprintf(fptr,"255 0 0\n"); break;
		    }
		}
	    }
	    fclose(fptr);
	} else {
	    fprintf(stderr,"Error opening file `%s' for writing.\n",fname);
	}
	free(mat);
    }
}

void rr_sparse_mat_stats(rr_sparse_mat* A) {
    printf("number of nonzero elements = %i\n",rr_sparse_mat_nonzeros(A));
    printf("fill = %f %%\n",(double)(rr_sparse_mat_nonzeros(A))/(double)(A->dim)/(double)(A->dim)*100.0);
}
