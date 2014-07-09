#ifndef _MATRIX_MATH__H__
#define _MATRIX_MATH__H__


struct MMatrix
{
	int nrows, ncols;
	float data[];

	float& operator()(int i, int j) { return data[i]; }
};

static void MMatrixDump(MMatrix* m);


static MMatrix* MMatrixCreate(int nrows, int ncols, int bZero)
{
	MMatrix* mat = (MMatrix*)malloc(sizeof(MMatrix)+(nrows*ncols)*sizeof(float));
	mat->nrows = nrows;
	mat->ncols = ncols;
	if (bZero)
		memset(&mat->data, 0, sizeof(float)*(nrows*ncols));
	return mat;
}

static void MMatrixDelete(MMatrix* m)
{
	free(m);
}

static void MMatrixMul(MMatrix* matC, const MMatrix* matA, const MMatrix* matB)
{
	assert(matC!=matA);
	assert(matC!=matB);
	assert(matA->ncols == matB->nrows);
	assert(matC->nrows == matA->nrows);
	assert(matC->ncols == matB->ncols);

	const int bInc = matB->ncols;
	const int cInc = matC->ncols;
	const float *pa=(const float*)&matA->data;
	float *pc = (float*)&matC->data;

	for (int r=0; r<matC->nrows; r++, pa+=matA->ncols)
	{
		for (int c=0; c<matC->ncols; c++, pc++)
		{
			float sum = 0.f;
			const float *pb = (float*)(matB->data+c);
			for (int i=0; i<matA->ncols; i++, pb+=bInc) sum += pa[i] * (*pb);
			*pc = sum;
		}
	}
}

// Slow hack
static MMatrix* MMatrixCreateAugmented(const MMatrix* matA, const MMatrix* matB)
{
	assert(matA->nrows==matB->nrows);
	MMatrix* aug = MMatrixCreate(matA->nrows, matA->ncols + matB->ncols, 0);
	for (int r=0; r<matA->nrows; r++)
	{
		memcpy(&aug->data[r*aug->ncols], &matA->data[r*matA->ncols], sizeof(float)*matA->ncols);
		memcpy(&aug->data[r*aug->ncols+matA->ncols], &matB->data[r*matB->ncols], sizeof(float)*matB->ncols);
	}
	return aug;
}

static MMatrix* MMatrixCreateSquareId(int n)
{
	MMatrix* m = MMatrixCreate(n,n,1);
	float *p = (float*)&m->data;
	for (int i=0; i<n; i++)
	{
		*p = 1.f;
		p += n+1;
	}
	return m;
}

///////////////////////////////////////////////////////////////////////////////////
// Gauuss-Jordan Elimination : create zeros above and below each diaganol pivot
// 0) Add ID to create the augmented matrix
// 1) Can switch rows
// 2) Can scale a row by any value
// 3) Can add any other scaled row
//
// First augment to create the following
//
//    | A A A A 1 0 0 0 | 
//    | A A A A 0 1 0 0 | 
//    | A A A A 0 0 1 0 | 
//    | A A A A 0 0 0 1 | 
// 
// then use row switching to achive
//
//    | 1 0 0 0 B B B B | 
//    | 0 1 0 0 B B B B | 
//    | 0 0 1 0 B B B B | 
//    | 0 0 0 1 B B B B | 
//
//  B contains the inverse (hopefully)
///////////////////////////////////////////////////////////////////////////////////

static bool MMatrixGaussJordanInvert(MMatrix* solution, const MMatrix* m)
{
	// Square Matrices Only!
	assert(m->ncols == m->nrows);
	assert(solution->ncols == m->nrows);

	MMatrix* id = MMatrixCreateSquareId(m->nrows);
	MMatrix* aug = MMatrixCreateAugmented(m, id);

	int n = m->nrows;
	int nrows = aug->nrows;
	int ncols = aug->ncols;
	float *p = (float*)&aug->data;

	// Sort so that largest is at the bottom
	if(0)for(int i=(n-1); i>0; i--)
    {
		// Swap a row?
        if(p[(i-1)*ncols]>p[i*ncols])
		{
			for(int j=0; j<ncols; j++)
			{
				float tmp = p[i*ncols + j];
				p[i*ncols + j] = p[(i-1)*ncols + j];
				p[(i-1)*ncols + j] = tmp;
			}
		}
    }

#define P(i,j) (p[(i*ncols)+j])

	// Lower tri
	for (int r=0; r<n; r++)
	{
		float d = P(r,r);
		if (d!=0.f)
		{
			// First scale the row
			for(int j=r; j<ncols; j++) P(r,j)/=d;
			// Eliminate the underneath rows
			for (int r2=r+1; r2<n; r2++)
			{
				float d2 = P(r2,r);
				for (int j=r; j<ncols; j++)
					P(r2,j) = P(r2,j) - d2*P(r,j);
			}
		}
		else
		{
			return false;
		}
		//printf("Debug:\n");
		//MMatrixDump(aug);
	}

	// Upper tri
	for (int r=n-1; r>=0; r--)
	{
		// Eliminate the underneath rows
		for (int r2=r-1; r2>=0; r2--)
		{
			float d2 = P(r2,r);
			for (int j=r; j<ncols; j++)
				P(r2,j) = P(r2,j) - d2*P(r,j);
		}
		//printf("Debug:\n");
		//MMatrixDump(aug);
	}

	//printf("Debug:\n");
	//MMatrixDump(aug);

	// Copy out the solution
	float* q = (float*)&solution->data;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++, q++)
			*q = P(i,j+n);
	}

	MMatrixDelete(id);
	MMatrixDelete(aug);
	
#undef P
	return true;
}

// Insert a into b at (row,col)
static void MMatrixInsert(const MMatrix* a, MMatrix* b, int row0, int col0)
{
	// assert((a->ncols+col0) <= b->ncols);
	// assert((a->nrows+row0) <= b->nrows);
	// for (int i=0; i<a->nrows; i++)
	// {
	// 	for (int j=0; j<a->ncols; j++)
	// 	{
	// 		b(i+row0, 0);//[j+col0] = 0;// a[i][j];
	// 	}
	// }
}


static void MMatrixInit(MMatrix* m, float* data)
{
	for (int i=0; i<m->nrows*m->ncols; i++)
		m->data[i] = data[i];
}

static void MMatrixDump(MMatrix* m)
{
	float *p = (float*)&m->data;
	for (int r=0; r<m->nrows; r++)
	{
		for (int c=0; c<m->ncols; c++)
			printf("%f,", *p++);
		printf("\n");
	}
}


#endif
