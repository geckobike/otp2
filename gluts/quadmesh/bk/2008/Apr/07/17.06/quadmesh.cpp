#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine/engine.h"

#include <GL/gl.h>
#include <GL/glut.h>


#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#define rnd rand
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define assert(x) { if (!(x)) { printf("ABORT: from file %s, line %d: expr failed: %s\n", __FILE__, __LINE__, #x); fflush(stdout); exit(0); } }

#define veclist(v) (v)[0], (v)[1], (v)[2], (v)[3]

typedef unsigned short u16;
typedef unsigned char u8;

static int s_width = 800;
static int s_height = 600;

static bool s_paused = false;
static float s_averageDt = 0.f;

static float s_elev=0.f;
static float s_yaw=0.f;
static float s_cameraRadius=1.f;
static vec3 s_cameraPos={0.f, 1.f, 0.f};
static vec3 s_cameraDir;
static vec3 s_cameraForward = {0.f, 0.f, 1.f};	// normalised direction in XZ plane
static vec3 s_cameraLeft = {-1.f, 0.f, 0.f};	// normalised direction in XZ plane
static int s_mouseDown = false;
static int s_lastMouseX=0;
static int s_lastMouseY=0;

#include "cloths.h"
#include "cloths.cpp"

#define ARRAY_SIZEOF(a) ((int)( sizeof((a)) / sizeof((a)[0]) ))

static int s_numCloths=ARRAY_SIZEOF(g_cloths);
static int s_cloth = 0;

/*
=======================================================================================
    Utils
=======================================================================================
*/


// Utils
template<class T>
struct List
{
	int count;
	T* array;
	List(T* _array)
	{
		array=_array;
		count = 0;
	}

	T* Append(T& t)
	{
		array[count]=t;
		T* tmp = array + count;
		count++;
		return tmp;
	}

	T* Append()
	{
		T* t = array + count;
		count++;
		return t;
	}
	T* ptr(int idx) { return array + idx; }
	T& ref(int idx) { return array[idx]; }
	const T& operator()(int idx) { return array[idx]; }
	T& operator[](int idx) { return array[idx]; }
	operator T* () { return array; }
};

#define LIST(type, x, size) type x##Buffer[size]; List<type> x(x##Buffer);

template <class T>
static void sortQuickSort(T* array, int l, int r)	// use: sortQuickSort(array, 0, N-1), need to implement "operator < ()"
{
	int i,j, mi;
	T k, m;

	i = l; j = r; mi = (l + r)/2;
	m = array[mi];
	while (i <= j) {
		while((i<=j) && array[i] < m) i++;
		while((j>l) && m < array[j]) j--;
		if (i <= j) {
			k = array[i]; array[i] = array[j]; array[j] = k;
			i++; j--;
		}
	}
	if (l < j) sortQuickSort<T>(array, l, j);
	if (i < r) sortQuickSort<T>(array, i, r);
}
void dlDrawLine(const vec3* from, const vec3* to, const vec4* colour)
{
	glBegin(GL_LINES);
	glColor3f(colour->x, colour->y, colour->z);
	glVertex3f(from->x, from->y, from->z);
	glColor3f(colour->x, colour->y, colour->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}

void dlDrawLine(const vec3* from, const vec3* to)
{
	glBegin(GL_LINES);
	glVertex3f(from->x, from->y, from->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}



static double averageTime()
{
    static bool init = false;
    static struct timeval lastTime;

    if (!init) { gettimeofday(&lastTime,NULL); init = true; }
    
    struct timeval now;
 
    long long sec;
    long long microsec;
    
    gettimeofday(&now,NULL);
    
    sec = now.tv_sec - lastTime.tv_sec;
    microsec = now.tv_usec - lastTime.tv_usec;
    
    lastTime = now;
   
    // Filter the average time to smooth it out

    double newDt = (double)sec + ( (double)microsec / (1000.0*1000.0) );

    return s_averageDt += 0.5*(newDt - s_averageDt);
}

static void wait(long long delay) // in micro secs
{
    struct timeval last;
    struct timeval now;
    
    gettimeofday(&last,NULL);
    gettimeofday(&now,NULL);

    while((now.tv_usec - last.tv_usec) < delay)
    {
	gettimeofday(&now,NULL);
	if (now.tv_sec > last.tv_sec)
	{
	    now.tv_usec += 1000000;
	}
    }
}

static void drawCircle(float x, float y, float r, bool filled = true)
{
    u32 N = 20;
    float dAlpha = PIx2 / (float)N;
    float a = 0.0f;
    filled ? glBegin(GL_POLYGON) : glBegin(GL_LINE_STRIP);
	for (u32 i=0; i<(N+1); i++)
	{
	    glVertex3f(x+r*cosf(a), y+r*sinf(a), 0.f);
	    a += dAlpha;
	}
    glEnd();
}

static const float drawEps = 0.003f;

static void drawPoint(float x, float y, float size = drawEps)
{
    glBegin(GL_QUADS);
	glVertex3f(x-size, y-size, 0.f);
	glVertex3f(x-size, y+size, 0.f);
	glVertex3f(x+size, y+size, 0.f);
	glVertex3f(x+size, y-size, 0.f);
    glEnd();
}


void updateCamera(const vec3* zz, const vec3* yy)
{
    vec3 up={0.f, 1.f, 0.f};
    vec3 forward={0.f, 0.f, -1.f};

    mtx rot;
    matrixRotY(&rot, s_yaw);
    vec3mtx33mulvec3(&s_cameraForward, &rot, &forward);

    s_cameraForward.y = 0.f;
    s_cameraForward.z = -cosf(s_yaw);
    s_cameraForward.x = +sinf(s_yaw);

    s_cameraDir = s_cameraForward;

    veccross(&s_cameraLeft, &up, &s_cameraForward);
    matrixRot(&rot, &s_cameraLeft, s_elev);
    vec3mtx33mulvec3(&s_cameraDir, &rot, &s_cameraForward);

    vec3 target = s_cameraPos;
    vec3 pos;
    vecsubscale(&pos, &target, &s_cameraDir, s_cameraRadius);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(pos.x, pos.y, pos.z, target.x, target.y, target.z, up.x, up.y, up.z);
}


/*
========================================================================================================
	Mesh
========================================================================================================
*/


typedef Cloth Mesh;
//struct Mesh
//{
//    int numVerts;
//    int numTris;
//    vec3* verts;
//    u16* tris;
//};

static void meshGenerateUniform(Mesh* mesh, int numX, int numY)
{
	mesh->numVerts = numX*numY;
	mesh->numTris = (numX-1)*(numY-1)*2;
	mesh->tris = (u16*)malloc(sizeof(u16)*3*mesh->numTris);
	mesh->verts = (vec3*)malloc(sizeof(vec3)*mesh->numVerts);
	
	float x = 0.f;
	float y = 0.f;
	float dx = 1.f/(numX-1.f);
	float dy = 1.f/(numY-1.f);

	vec3* v = mesh->verts;

	for (int j=0; j<numY; j++)
	{
		x = 0;
		for (int i=0; i<numX; i++)
		{
			vecset(v, x, y, 0.f);
			v++;
			x += dy;
		}
		y += dx;
	}

#define V(i,j) (numX*(j)+(i))

	u16* t = mesh->tris;
	bool even;

	for (int j=0; j<numY-1; j++)
	{
		even = j%2==0;
		for (int i=0; i<numX-1; i++)
		{
			u16* tmp = t;
			bool e = even;
			if (i==j) e = !e;
			if (i*2==j) e = !e;
			if (e)
			{
				// First tri from quad
				*t = V(i,j); t++;
				*t = V(i+1,j); t++;
				*t = V(i,j+1); t++;

				//Second tri from quad
				*t = V(i,j+1); t++;
				*t = V(i+1,j); t++;
				*t = V(i+1,j+1); t++;
			}
			else
			{
				// First tri from quad
				*t = V(i,j); t++;
				*t = V(i+1,j+1); t++;
				*t = V(i,j+1); t++;

				//Second tri from quad
				*t = V(i,j); t++;
				*t = V(i+1,j); t++;
				*t = V(i+1,j+1); t++;
			}
			if (1)
			{
				#define DELETETRIS {tmp[0]=0;tmp[1]=0;tmp[2]=0;tmp[3]=0;tmp[4]=0;tmp[5]=0;}
				#define DELETETRI_A {tmp[0]=0;tmp[1]=0;tmp[2]=0;}
				#define DELETETRI_B {tmp[3]=0;tmp[4]=0;tmp[5]=0;}
				if (j==5 && i==4) DELETETRIS;
				if (j==5 && i==5) DELETETRIS;
				if (j==5 && i==6) DELETETRIS;//_A;
				if (j==6 && i==6) DELETETRIS;//_B;
				if (j==6 && i==5) DELETETRIS;
				if (j==7 && i==6) DELETETRIS;

				if (j==2 && i==4) DELETETRIS;
				if (j==3 && i==4) DELETETRIS;
				if (j==3 && i==5) DELETETRIS;
				if (j==4 && i==5) DELETETRIS;
				if (j==4 && i==6) DELETETRIS;
			}
			even = !even;
		}
	}

#undef V
}
	
static void meshRandomise(Mesh* mesh)
{
	const int N = mesh->numVerts;
	int vertMap[1024];
	for (int i=0; i<N; i++)
	{
		vertMap[i] = i;
	}

	for (int i=0; i<N; i++)
	{
		int idx0 = rnd() % N;
		int idx1 = rnd() % N;
		// Swap
		int tmp = vertMap[idx0];
		vertMap[idx0] = vertMap[idx1];
		vertMap[idx1] = tmp;
	}

	vec3 verts[1024];
	for (int i=0; i<mesh->numVerts; i++)
	{
		verts[vertMap[i]] = mesh->verts[i];
	}
	for (int i=0; i<mesh->numVerts; i++)
	{
		mesh->verts[i] = verts[i];
	}

	for (int i=0; i<mesh->numTris*3; i++)
	{
		mesh->tris[i] = vertMap[mesh->tris[i]];
	}
}

static void meshTranslate(Mesh* mesh, float dx, float dy)
{
	vec3 dv={dx, dy, 0.f};
	for (int i=0; i<(const int)mesh->numVerts; i++)
	{
		vecadd(&mesh->verts[i], &mesh->verts[i], &dv);
	}
}

static void meshScale(Mesh* mesh, float scale)
{
	for (int i=0; i<(const int)mesh->numVerts; i++)
	{
		vecscale(&mesh->verts[i], &mesh->verts[i], scale);
	}
}

static void meshDraw(const Mesh* mesh)
{
	u16* t = mesh->tris;
	for (int i=0; i<mesh->numTris; i++)
	{
		u16 i0 = *t; t++;
		u16 i1 = *t; t++;
		u16 i2 = *t; t++;
		const vec3* v0 = mesh->verts + i0;
		const vec3* v1 = mesh->verts + i1;
		const vec3* v2 = mesh->verts + i2;
		glBegin(GL_LINE_STRIP);
		glVertex3f(v0->x, v0->y, v0->z);
		glVertex3f(v1->x, v1->y, v1->z);
		glVertex3f(v2->x, v2->y, v2->z);
		glVertex3f(v0->x, v0->y, v0->z);
		glEnd();
	}
}

/*
========================================================================================
    Adjacency
========================================================================================
*/

struct Adjacency
{
    int num;
    int* idxs;
    
    int Find(int idx)
    {
	for (int i=0; i<num; i++)
	    if (idxs[i] == idx) return i;
	return -1;
    }
    void Remove(int idx)
    {
	for (int i=0; i<num; i++)
	{
	    if (idxs[i] == idx)
	    {
		idxs[i] = idxs[num-1];
		num--;
	    }
	}
    }
};

#define memAlloc(l,g) ::malloc(l)
#define memFree(p) ::free(p)

// The two verts are assumed to be the first two u16's in the structure
static Adjacency* adjacencyBuild(u32 numVerts, const void* links, u32 skipSize, u32 numLinks)
{
    // Allocate memory upfront. This avoids lots of little allocs
    Adjacency* ad = (Adjacency*) memAlloc(sizeof(Adjacency)*numVerts + sizeof(int)*numLinks*2, g_clothHeap);
    memset(ad, 0, sizeof(Adjacency)*numVerts);
    int* memory = (int*)(&ad[numVerts]); // Memory pool to hold the adjacency indices

    // First count the adjacencies for each vert
    for (u32 i=0; i<numLinks; i++)
    {
	const u16* verts = (u16*)((u8*)(links) + skipSize*i);
	u16 i0 = verts[0];
	u16 i1 = verts[1];
	assert(i0<numVerts);
	assert(i1<numVerts);
	ad[i0].num++;
	ad[i1].num++;
    }

    // Set up the pointers into the memory pool
    u32 count = 0;
    for (u32 i=0; i < numVerts; i++)
    {
	ad[i].idxs = &memory[count];
	count = count + ad[i].num;
	ad[i].num = 0;
    }
    assert(count == 2*numLinks);
    
    // Finally create the adjacency data
    for (u32 i=0; i<numLinks; i++)
    {
	const u16* verts = (u16*)((u8*)(links) + skipSize*i);
	u16 i0 = verts[0];
	u16 i1 = verts[1];
	
	ad[i0].idxs[ad[i0].num] = i1;
	ad[i1].idxs[ad[i1].num] = i0;
	
	ad[i0].num++;
	ad[i1].num++;
    }
    
    return ad;
}

static void adjacencyDump(Adjacency* adjacency, u32 numVerts)
{
    printf("adjacency:\n");
    for (u32 i=0; i<numVerts; i++)
    {
	printf("%d: ", i);
	u32 num = adjacency[i].num;
	for (u32 j=0; j<num; j++)
	{
	    printf(" %d", adjacency[i].idxs[j]);
	}
	printf("\n");
    }
    printf("--------------------------\n");
}

/*
========================================================================================================
	Links
========================================================================================================
*/

struct Edge
{
	union
	{
		u32 hash;
		struct
		{
			u16 v0;
			u16 v1;
		};
	};
	u16 third;
	u16 tri;
	bool deleted;
	bool deletable;
	bool isBorder;

	void create(int _v0, int _v1, int _third, int _tri)
	{
		if (_v0<_v1)
		{
			v0 = _v0;
			v1 = _v1;
		}
		else
		{
			v0 = _v1;
			v1 = _v0;
		}
		third = _third;
		tri = _tri;
	}

	bool operator==(const Edge& e)
	{
		return e.v0==v0 && e.v1==v1;
	}
	bool operator<(const Edge& e)
	{
		if (v0<e.v0) return true;
		if (v0>e.v0) return false;
		return v1<e.v1;
	}
};


static void edgesDraw(List<Edge>* edges, const vec3* verts)
{
	for (int i=0; i<edges->count; i++)
	{
		Edge* e = edges->ptr(i);
		if (1&&e->deletable)
		//if (e->isBorder)
		{
			glColor3f(0.f, 0.f, 1.f);
		}
		else
		{
			glColor3f(1.f, 1.f, 1.f);
		}
		dlDrawLine(&verts[e->v0], &verts[e->v1]);
	}
}

struct EdgeHandle : public Edge
{
    bool isBorder;
    bool deletable;
    bool deleted;
    int fourth; // fourth vert [only valid for an edge that joins two tris!]
    void assign(const Edge* e)
    {
	v0 = e->v0;
	v1 = e->v1;
    }
};

struct EdgeHandleTable
{
    enum {WIDTH=4096};
    EdgeHandle edges[WIDTH];
    int maxSearch;

    void init()
    {
	memset(edges, 0, sizeof(edges));
	maxSearch = 1;
    }

    int hash(const Edge* e)
    {
	return (((e->v0*31) ^ (e->v1*127)) ^ (e->v0*128) ^ (e->v1*128)) % WIDTH;
    }

    EdgeHandle* find(int v0, int v1)
    {
	Edge e;
	e.create(v0, v1,0,0);
	return find(&e);
    }

    EdgeHandle* find(const Edge* e)
    {
	int idx = hash(e);
	for (int i=0; i<maxSearch; i++)
	{
	    if (edges[idx].hash==e->hash) return edges + idx;
	    idx = (idx + 1)%WIDTH;
	}
	return NULL;
    }

    EdgeHandle* insert(const Edge* e)
    {
	int idx = hash(e);
	int c = 1;
	while (edges[idx].hash)
	{
	    idx = (idx + 1)%WIDTH;
	    c++;
	}
	maxSearch = MAX(c, maxSearch);
	edges[idx].assign(e);
	return edges + idx;
    }
};


static void edgeHandleDelete(EdgeHandleTable& table, int a, int b, int c, int d=-1)
{
    EdgeHandle* eh = table.find(a,b);
    eh->deleted = true;

    EdgeHandle* eh1 = table.find(a,c);
    EdgeHandle* eh2 = table.find(b,c);
    eh1->deletable = false;
    eh2->deletable = false;
    if (d>=0)
    {
	EdgeHandle* eh3 = table.find(a,d);
	EdgeHandle* eh4 = table.find(b,d);
	eh3->deletable = false;
	eh4->deletable = false;
    }
}

static bool edgeHandleCanDelete(EdgeHandleTable& table, int a, int b, int c, int d=-1)
{
	EdgeHandle* eh = table.find(a,b);
	bool canDelete = eh->deletable;
	
	EdgeHandle* eh1 = table.find(a,c);
	EdgeHandle* eh2 = table.find(b,c);
	canDelete = canDelete && !eh1->deleted && !eh2->deleted && !eh1->isBorder && !eh2->isBorder;
	canDelete = canDelete && !eh1->deletable && !eh2->deletable;

	if (canDelete && d>=0)
	{
		EdgeHandle* eh3 = table.find(a,d);
		EdgeHandle* eh4 = table.find(b,d);
		canDelete = canDelete && !eh3->deleted && !eh3->deleted;
		canDelete = canDelete && !eh3->isBorder && !eh4->isBorder;
		canDelete = canDelete && !eh3->deletable && !eh4->deletable;
	}
	return canDelete;
}

static bool edgeHandleShouldDelete(int pass, EdgeHandleTable& table, int a, int b, int c, int d, const vec3* verts, int countA, int countB, int countC, int countD, int borderA, int borderB, int borderC, int borderD)
{
	EdgeHandle* eh = table.find(a,b);
	bool canDelete = eh->deletable;
	
	if (canDelete)
	{
		EdgeHandle* eh1 = table.find(a,c);
		EdgeHandle* eh2 = table.find(b,c);
		EdgeHandle* eh3 = table.find(a,d);
		EdgeHandle* eh4 = table.find(b,d);
		int count = 0;
		count += !eh1->deletable ? 1 : 0;
		count += !eh2->deletable ? 1 : 0;
		count += !eh3->deletable ? 1 : 0;
		count += !eh4->deletable ? 1 : 0;
		
		if (!eh1->deletable && !eh2->deletable)
		{
		    return true;
		}

		if (!eh3->deletable && !eh4->deletable)
		{
		    return true;
		}

		if (count ==4)
		{
			return true;
		}
		if (pass>0 && (!eh1->isBorder && !eh2->isBorder && !eh3->isBorder && !eh4->isBorder))
		{
		    float sizeAB = vecdistsq(&verts[eh->v0], &verts[eh->v1]);

		    if (!eh1->deletable && !eh1->isBorder)
		    {
			float sizeBC = vecdistsq(&verts[eh2->v0], &verts[eh2->v1]);
			if (sizeAB > sizeBC) return true;
		    }

		    if (!eh2->deletable && !eh2->isBorder)
		    {
			float sizeAC = vecdistsq(&verts[eh1->v0], &verts[eh1->v1]);
			if (sizeAB > sizeAC) return true;
		    }

		    if (!eh3->deletable && !eh3->isBorder)
		    {
			float sizeBD = vecdistsq(&verts[eh4->v0], &verts[eh4->v1]);
			if (sizeAB > sizeBD) return true;
		    }

		    if (!eh4->deletable && !eh4->isBorder)
		    {
			float sizeAD = vecdistsq(&verts[eh3->v0], &verts[eh3->v1]);
			if (sizeAB > sizeAD) return true;
		    }
		}

		////if (eh1->isBorder || eh2->isBorder)// &&*/ !eh1->deletable && !eh2->deletable)
		//{
		//	float sizeA = vecdistsq(&verts[eh->v0], &verts[eh->v1]);
		//	float sizeB = vecdistsq(&verts[eh1->v0], &verts[eh1->v1]);
		//	float sizeC = vecdistsq(&verts[eh2->v0], &verts[eh2->v1]);
		//	if (sizeA>sizeB*1.1 && sizeA>sizeC*1.1) return true;
		//}
		////if (eh3->isBorder || eh4->isBorder)// &&*/ !eh3->deletable && !eh4->deletable)
		//{
		//	float sizeA = vecdistsq(&verts[eh->v0], &verts[eh->v1]);
		//	float sizeB = vecdistsq(&verts[eh3->v0], &verts[eh3->v1]);
		//	float sizeC = vecdistsq(&verts[eh4->v0], &verts[eh4->v1]);
		//	if (sizeA>sizeB*1.1 && sizeA>sizeC*1.1) return true;
		//}
		//if (!eh1->deletable && !eh4->deletable) return true;
		//if (!eh2->deletable && !eh3->deletable) return true;

		//if ((countA==3 || countA==3) && borderA && !borderB && !borderC && borderD) return true;
		//if ((countA==3 || countA==3) && borderA && !borderB && borderC && !borderD) return true;
		//if ((countB==3 || countB==3) && borderB && !borderA && !borderC && borderD) return true;
		//if ((countB==3 || countB==3) && borderB && !borderA && borderC && !borderD) return true;

	}
	return false;
}

static int gMaxToDelete=0;

static void createLinks(Mesh* mesh, List<Edge>* edgesOut)
{
    const vec3* verts = mesh->verts;
    enum {maxN=10240};
    LIST(Edge, edges, maxN);
    int count[maxN]={0};
    bool vertIsBorder[maxN]={0};
    int i;

    {
	// Insert each edge fom every triangle
	u16* t = mesh->tris;
	for (i=0; i<mesh->numTris; i++)
	{
	    u16 i0 = *t; t++;
	    u16 i1 = *t; t++;
	    u16 i2 = *t; t++;
	    Edge e;
	    e.create(i0,i1,i2,i); edges.Append(e);
	    e.create(i1,i2,i0,i); edges.Append(e);
	    e.create(i2,i0,i1,i); edges.Append(e);
	}
    }

    sortQuickSort(edges.array, 0, edges.count-1);

    EdgeHandleTable table;
    table.init();

    // Count the number of equal edges that are contiguous in the list
    // and setup the vert count and hash table
    for (i=0; i<edges.count; )
    {	
	int j=i;
	Edge e = edges[i];
	EdgeHandle* eh = table.insert(&e);
	while (j<edges.count && edges[j]==e)
	{
	    j++;
	    count[i]++;

	}
	if (count[i]==2)
	{
	    eh->fourth = edges[i+1].third;
	    eh->deletable = true;
	}
	else
	{
	    // borders and degenerate edges can't be deleted
	    eh->deletable= false;
	    eh->fourth = -1;
	}
	if (count[i]==1)
	{
	    vertIsBorder[e.v0]=true;
	    vertIsBorder[e.v1]=true;
	}
	eh->third = e.third;
	eh->isBorder = count[i]==1;
	eh->deleted = false;
	i = i + count[i];
    }
    
    //printf("maxSearch = %d\n", table.maxSearch);

    //======================
    // Create Adjacency Data
    //======================
    LIST(Edge, reducedEdges, maxN);
    for (i=0; i<edges.count; )
    {	
	Edge e = edges[i];
	reducedEdges.Append(e);
	i = i + count[i];
    }

    Adjacency* adjacency = adjacencyBuild(mesh->numVerts, reducedEdges.array, sizeof(reducedEdges.array[0]), reducedEdges.count);
    //adjacencyDump(adjacency, mesh->numVerts);

    {
	for (int i=0; i<mesh->numVerts; i++)
	{
	    Adjacency* ad = adjacency + i;
	    if (!vertIsBorder[i] && ad->num==4)
	    {
		// Then these struts can't be deleted
		for (int j=0; j<ad->num; j++)
		{
		    EdgeHandle* eh = table.find(i, ad->idxs[j]);
		    eh->deletable = false;
		}
	    }
	}

	// Reduce further
	for (int i=0; i<edges.count; )
	{
	    Edge* e = edges + i;
	    EdgeHandle* eh = table.find(e);
	    if (eh->deletable)
	    {
		int a = e->v0;
		int b = e->v1;
		int c = e->third;
		int countA = adjacency[a].num;
		int countB = adjacency[b].num;
		int countC = adjacency[c].num;
		if (count[i]==2)
		{
		    Edge* e1 = e + 1;
		    int d = e1->third;
		    int countD = adjacency[d].num;
		    if (edgeHandleShouldDelete(0, table, a, b, c, d, mesh->verts,countA,countB,countC,countD,vertIsBorder[a], vertIsBorder[b], vertIsBorder[c], vertIsBorder[d]))
		    {
			if (!e->deleted)
			{
			    //if (numDeleted>=gMaxToDelete) goto finished;
			    edgeHandleDelete(table, a, b, c, d);
			    e->deleted=true;
			    //numDeleted++;
			}
		    }
		}
	    }
	    i = i + count[i];
	}

	int numDeletedTotal=0;
	int numDeleted;
	if(0)do
	{
	    numDeleted = 0;
	    for (int i=0; i<(const int)mesh->numVerts; i++)
	    {
		float longest = 0.f;
		EdgeHandle* longestEdge = NULL;

		Adjacency* ad = adjacency + i;
		if (!vertIsBorder[i] && ad->num>4)
		{
		    // Find the longest edge and delete it
		    for (int j=0; j<(const int)ad->num; j++)
		    {
			const int other = ad->idxs[j];
			bool ok = (!vertIsBorder[other] && adjacency[other].num>4);
			ok = ok || (vertIsBorder[other] && adjacency[other].num>2);
			if(ok)
			{
			    EdgeHandle* eh = table.find(i,other);
			    if (eh->deletable)
			    {
				float l = vecdistsq(&verts[i], &verts[other]);
				if (l>longest)
				{
				    longest = l;
				    if (!eh->deleted) longestEdge = eh;
				}
			    }
			}
		    }
		}
		if (longestEdge)
		{
		    longestEdge->deleted = true;
		    adjacency[longestEdge->v0].Remove(longestEdge->v1);
		    adjacency[longestEdge->v1].Remove(longestEdge->v0);

		    if (longestEdge->fourth>=0)
		    {
			edgeHandleDelete(table, longestEdge->v0, longestEdge->v1, longestEdge->third, longestEdge->fourth);
		    }
		    numDeleted++;
		    numDeletedTotal++;
		    if (0&&numDeletedTotal>=gMaxToDelete)
		    {
			numDeleted=0;
			break;
		    }
		}
	    }

	} while(numDeleted>0);
    }

    for (int rep=0; rep<100; rep++)
    // Reduce further
    for (i=0; i<edges.count; )
    {
	Edge* e = edges + i;
	EdgeHandle* eh = table.find(e);
	if (eh->deletable)
	{
	    int a = e->v0;
	    int b = e->v1;
	    int c = e->third;
	    int countA = adjacency[a].num;
	    int countB = adjacency[b].num;
	    int countC = adjacency[c].num;
	    if (count[i]==2)
	    {
		Edge* e1 = e + 1;
		int d = e1->third;
		int countD = adjacency[d].num;
		if (edgeHandleShouldDelete(1, table, a, b, c, d, mesh->verts,countA,countB,countC,countD,vertIsBorder[a], vertIsBorder[b], vertIsBorder[c], vertIsBorder[d]))
		{
		    if (!e->deleted)
		    {
			//if (numDeleted>=gMaxToDelete) goto finished;
			edgeHandleDelete(table, a, b, c, d);
			e->deleted=true;
			//numDeleted++;
		    }
		}
	    }
	}
	i = i + count[i];
    }
    finished: ;

    // Copy out the edges

    edgesOut->count=0;
    for (i=0; i<edges.count; )
    {
      if (!edges[i].deleted)
      {
	  EdgeHandle* eh = table.find(&edges[i]);
	  if (!eh->deleted)
	  {
	      Edge* edge = edgesOut->Append(edges[i]);
	      edge->deletable = eh->deletable;
	      edge->isBorder = eh->isBorder;
	  }
      }
      i = i + count[i];
    }
}


#define gMesh (g_cloths[s_cloth])
LIST(Edge, gEdges, 10240);


static void animateScene ()
{
    createLinks(&gMesh, &gEdges);
}

void render()
{
    animateScene();
    
    updateCamera(&s_cameraPos, &s_cameraDir);
    
    u32 N = 10;
    float width=7.f;
    float dx = width/(float)(N-1);
    float dy = width/(float)(N-1);
    float y=-width*0.5f;
    float x=-width*0.5f;
    
    glColor3f(0.2f, 0.2f, 0.2f);
    for (u32 i=0; i<N; i++)
    {
	vec3 from = {x, 0.f, -y};
	vec3 to = {x, 0.f, +y};
	dlDrawLine(&from, &to);
	x += dx;
    }
    
    y=-width*0.5f;
    x=-width*0.5f;
    for (u32 i=0; i<N; i++)
    {
	vec3 from = {+x, 0.f, y};
	vec3 to = {-x, 0.f, y};
	dlDrawLine(&from, &to);
	y += dy;
    }


    glDisable(GL_DEPTH_TEST);
    glColor3f(0.f, 0.6f, 0.f);
    meshDraw(&gMesh);
    glColor3f(1.f, 1.f, 1.f);
    edgesDraw(&gEdges, gMesh.verts);
}


/*
gggg
*/


static void start()
{
    //meshGenerateUniform(&gMesh, 10, 10);
    //meshRandomise(&gMesh);
    //meshTranslate(&gMesh, -0.5f, -0.5f);
    //meshScale(&gMesh, 15.f);
    createLinks(&gMesh, &gEdges);

    GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    /*	light_position is NOT default value	*/
    GLfloat light_position0[] = { 1.0, 10.0, 1.0, 0.0 };
    GLfloat light_position1[] = { -1.0, -10.0, -1.0, 0.0 };
  
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glClearColor(0.3,0.3,0.3,0);

    //  glEnable(GL_CULL_FACE);
    //  glCullFace(GL_BACK);
}

static void keyboard(unsigned char key, int x, int y)
{
    if (key==27) exit(0);

    if (key==']')
    {
	s_cloth++;
	s_cloth = s_cloth % s_numCloths;
	printf("s_cloth = %d\n", s_cloth);
    }
    if (key=='[')
    {
	s_cloth--;
	s_cloth = s_cloth % s_numCloths;
	printf("s_cloth = %d\n", s_cloth);
    }

    vec3 up={0.f, 1.f, 0.f};

    if (key=='q') vecaddscale(&s_cameraPos, &s_cameraPos, &up, -0.1f);
    if (key=='e') vecaddscale(&s_cameraPos, &s_cameraPos, &up, +0.1f);

    if (key=='w')
    {
	s_cameraRadius -= 0.1f;
    }
    if (key=='s')
    {
	s_cameraRadius += 0.1f;
    }

    if (key=='a')
    {
	vecaddscale(&s_cameraPos, &s_cameraPos, &s_cameraLeft, +0.1f);
    }
    if (key=='d')
    {
	vecaddscale(&s_cameraPos, &s_cameraPos, &s_cameraLeft, -0.1f);
    }
    if (key=='p') gMaxToDelete++;
    if (key=='l') gMaxToDelete--;

    glutPostRedisplay();
}


static void display()
{
    // clear the window
    glClearColor (0.5,0.5,0.5,0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // go to GL_MODELVIEW matrix mode and set the camera
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    // leave openGL in a known state - flat shaded white, no textures
    glDisable (GL_TEXTURE_2D);
    glShadeModel (GL_FLAT);
    glDisable (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);
    glColor4f (1.f,1.f,1.f,1.f);

    render();

    glutSwapBuffers();
}

static void reshape(GLint width, GLint height)
{
    s_width = width;
    s_height = height;
    
    // setup viewport
    glViewport (0,0,s_width,s_height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    const float vnear = 0.01f;
    const float vfar = 100.0f;
    const float k = 0.8f;     // view scale, 1 = +/- 45 degrees
    if (s_width >= s_height) {
	float k2 = float(s_height)/float(s_width);
	glFrustum (-vnear*k,vnear*k,-vnear*k*k2,vnear*k*k2,vnear,vfar);
    }
    else {
	float k2 = float(s_width)/float(s_height);
	glFrustum (-vnear*k*k2,vnear*k*k2,-vnear*k,vnear*k,vnear,vfar);
    }
    glutPostRedisplay();
}

static void initGraphics()
{
}

static void mouseButton(int button, int state, int x, int y)
{
    // 0 == down, 1 == up
    int mouseDown = !state;

    if (!s_mouseDown && mouseDown)
    {
	s_lastMouseX = x;
	s_lastMouseY = y;
    }
    s_mouseDown = mouseDown;
}

static void mouseMotion(int x, int y)
{
    if (!s_mouseDown) return;

    float dx = (float)(x - s_lastMouseX)/(float)s_width;
    float dy = (float)(y - s_lastMouseY)/(float)s_height;

    s_lastMouseX = x;
    s_lastMouseY = y;
    
    s_yaw += dx*1.0f;
    s_elev += dy*1.0f;

    if (s_elev > 0.9*PI_D2) s_elev = 0.9*PI_D2;
    if (s_elev < -0.9*PI_D2) s_elev = -0.9*PI_D2;

    glutPostRedisplay();
}


int main (int argc, char **argv)
{
    start();

    // GLUT Window Initialization:
    glutInit (&argc, argv);
    glutInitWindowSize (s_width, s_height);
    glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow ("triangles to quad mesh");

    // Initialize OpenGL graphics state
    initGraphics();

    // Register callbacks:
    glutDisplayFunc (display);
    glutReshapeFunc (reshape);
    glutKeyboardFunc (keyboard);
    glutMouseFunc (mouseButton);
    glutMotionFunc (mouseMotion);
    glutIdleFunc (animateScene);

    // Turn the flow of control over to GLUT
    glutMainLoop ();

    //printf("average dt was = %f, and FPS = %f\n", s_averageDt, 1.f/s_averageDt);

    return 0;
}


