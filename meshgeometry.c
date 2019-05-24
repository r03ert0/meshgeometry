//char version[]="meshgeometry, version 1, roberto toro, 12 August 2009";
//char version[]="meshgeometry, version 2, roberto toro, 28 April 2010";
//char version[]="meshgeometry, version 3, roberto toro, 1 May 2010";     // added -laplace and -taubinLM
//char version[]="meshgeometry, version 4, roberto toro, 21 May 2010";    // commands are processed as a chain
//char version[]="meshgeometry, version 5, roberto toro, 28 May 2012";    // added several commands: foldLength, volume, absgi, texture threshold, countClusters, and includes meshconvert v8
//char version[]="meshgeometry, version 6, roberto toro, 10 November 2012"; // added randomverts, help, centre, normalise, normal, verbose, off mesh format (load and save), added to github
//char version[]="meshgeometry, version 7, roberto toro, 17 Decembre 2014"; // vtk support
//char version[]="meshgeometry, version 8, roberto toro, 26 Decembre 2015";
//char version[]="meshgeometry, version 9, roberto toro, 10 June 2017"; // add gii reader
char version[]="meshgeometry, version 10, roberto toro, 7 November 2017"; // add asc reader/writer

/*
    To use:
    
    ./meshgeometry_mac -i /Applications/_Neuro/freesurfer510/subjects/bert/surf/lh.inflated -i /Applications/_Neuro/freesurfer510/subjects/bert/surf/lh.curv -drawSurface bert.tif lat

    To compile:
    
    On Mac OS X:
    gcc -Wall meshgeometry.c -o meshgeometry_mac -framework Carbon -framework OpenGL -framework GLUT

    On Unix:
    gcc -Wall meshgeometry.c -o meshgeometry_unix -lGL -lGLU -lglut -lm -lz

    On Windows:
    gcc -Wall meshgeometry.c -o meshgeometry_win.exe -lopengl32 -lglut32
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

// OpenGL libraries
#ifdef __APPLE__
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
  #include <GLUT/glut.h>
#else
  #include <GL/gl.h>
  #include <GL/glut.h>
#endif

#define pi 3.14159265358979323846264338327950288419716939927510
#define EPSILON  1e-8 // small enough to avoid division overflow
#define min(a,b) (((a)<(b))?(a):(b))

#define kMAXNETRIS          100
#define kFreeSurferMesh     1
#define kFreeSurferData     2
#define kFreeSurferAnnot    3    
#define kBrainVisaMesh      4
#define kFloatData          5
#define kRawFloatData       6
#define kText               7
#define kTextData           8
#define kVRMLMesh           9
#define kObjMesh            10
#define kPlyMesh            11
#define kSTLMesh            12
#define kSmeshMesh          13
#define kBinMesh            14
#define kOffMesh            15
#define kMGHData            16
#define kVTKMesh            17
#define kDPVData            18
#define kCivetObjMesh       19
#define kGiiMesh            20
#define kGiiData            21
#define kAscMesh            22

typedef struct
{
    float   x,y,z;
}float3D;
typedef struct
{
    int    a,b,c;
}int3D;
typedef struct
{
    int    a,b;
}int2D;
#define SIZESTACK    64
typedef struct
{
    int     n;
    int     t[SIZESTACK];
}NTriRec;
typedef struct
{
    int     np;         // number of vertices
    int     nt;         // number of triangles
    int     ddim;       // data dimensions (default: 1)
    float3D *p;         // vertices
    int3D   *t;         // triangles
    float   *data;      // data
    char    *selection; // selection
    NTriRec *NT;        // neighbouring triangles
}Mesh;

float area(Mesh *m);
float volume(Mesh *m);
int smooth(Mesh *m);
int taubin(float lambda, float mu, int N, Mesh *m);
float minData(Mesh *m);
float maxData(Mesh *m);
int nonmanifold_tris(Mesh *mesh);

int g_gluInitFlag=0;

Mesh    mesh;
float   R;
int     verbose=0;

int smoothData(Mesh *m,float l,int niter);

#pragma mark -
float dot3D(float3D a, float3D b)
{
    return (float){a.x*b.x+a.y*b.y+a.z*b.z};
}
float3D cross3D(float3D a, float3D b)
{
    return (float3D){a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
float3D add3D(float3D a, float3D b)
{
    return (float3D){a.x+b.x,a.y+b.y,a.z+b.z};
}
float3D sub3D(float3D a, float3D b)
{
    return (float3D){a.x-b.x,a.y-b.y,a.z-b.z};
}
float3D sca3D(float3D a, float t)
{
    return (float3D){a.x*t,a.y*t,a.z*t};
}
float norm3D(float3D a)
{
    return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}
float norm3Dsqr(float3D a)
{
    return a.x*a.x+a.y*a.y+a.z*a.z;
}
float3D normal3D(int i, Mesh *m)
{
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D N;
    
    N=cross3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
    return sca3D(N,1/norm3D(N));
}

float determinant(float3D a, float3D b, float3D c)
{
    float   D= a.x*(b.y*c.z-c.y*b.z)+
                a.y*(b.z*c.x-c.z*b.x)+
                a.z*(b.x*c.y-c.x*b.y);
    return D;
}

int multMatVec(float *m, float3D v, float3D *result)
{
    /*
    Multiplies the vector v by the 4x4 matrix m, and puts the result into result
    */
    float3D r;
    r.x= m[0]*v.x +m[1]*v.y +m[2]*v.z+m[3];
    r.y= m[4]*v.x +m[5]*v.y +m[6]*v.z+m[7];
    r.z= m[8]*v.x +m[9]*v.y +m[10]*v.z+m[11];
    result->x=r.x;
    result->y=r.y;
    result->z=r.z;
    
    return 0;
}


#pragma mark -
#pragma mark [ Utilities ]
int     endianness;
#define kMOTOROLA   1
#define kINTEL      2
void checkEndianness(void)
{
    char    b[]={1,0,0,0};
    int     num=*(int*)b;
    
    if(num==16777216)
        endianness=kMOTOROLA;
    else
        endianness=kINTEL;
}
void swapint(int *n)
{
    char    *by=(char*)n;
    char    sw[4]={by[3],by[2],by[1],by[0]};
    
    *n=*(int*)sw;
}
void swapfloat(float *n)
{
    char    *by=(char*)n;
    char    sw[4]={by[3],by[2],by[1],by[0]};
    
    *n=*(float*)sw;
}
void swaptriangles(Mesh *m)
{
    int     nt=m->nt;
    int3D   *t=m->t;
    int     i;
    
    for(i=0;i<nt;i++)
    {
        swapint(&t[i].a);
        swapint(&t[i].b);
        swapint(&t[i].c);
    }
}
void swapvertices(Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p;
    int     i;
    
    for(i=0;i<np;i++)
    {
        swapfloat(&p[i].x);
        swapfloat(&p[i].y);
        swapfloat(&p[i].z);
    }
}
// triangle area using Heron's formula
float triangle_area(float3D p0, float3D p1, float3D p2)
{
    float   a,b,c;    // side lengths
    float   s;        // semiperimeter
    float   area;
    
    a=norm3D(sub3D(p0,p1));
    b=norm3D(sub3D(p1,p2));
    c=norm3D(sub3D(p2,p0));
    s=(a+b+c)/2.0;
    
    if(s*(s-a)*(s-b)*(s-c)<0)
        area=0;
    else
        area=sqrt(s*(s-a)*(s-b)*(s-c));
    
    return area;
}
// Adapted from intersect_RayTriangle()
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
//    Input:  vector "x", triangle index "it"
//    Output: *c0, *c1 = the triangle-based coordinates of the intersection (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 = disjoint (no intersect)
//             1 = intersect in unique point I1
//             2 = are in the same plane
// code from:http://geometryalgorithms.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
int intersect_VectorTriangle(float3D x, int i, float *c0, float *c1, Mesh *m)
{
    float3D *p=m->p;
    int3D   *t=m->t;
    int3D   T=t[i];
    double xx[3];
    double u[3], v[3], n[3];   // triangle vectors
    double dir[3],w0[3], w[3]; // ray vectors
    double  r, a, b;           // params to calc ray-plane intersect
    double  uu, uv, vv, wu, wv, D;
    double  ss,tt;

    u[0]=p[T.b].x-p[T.a].x;
    u[1]=p[T.b].y-p[T.a].y;
    u[2]=p[T.b].z-p[T.a].z;
    
    v[0]=p[T.c].x-p[T.a].x;
    v[1]=p[T.c].y-p[T.a].y;
    v[2]=p[T.c].z-p[T.a].z;
    
    n[0]=u[1]*v[2]-u[2]*v[1];
    n[1]=u[2]*v[0]-u[0]*v[2];
    n[2]=u[0]*v[1]-u[1]*v[0];

    if(sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])<1e-10)
    {
        //printf("%lf\n", sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]));        // triangle is degenerate, do not deal with this case
        return -1;
    }

    dir[0]=x.x;
    dir[1]=x.y;
    dir[2]=x.z;
    
    w0[0] = -p[T.a].x;
    w0[1] = -p[T.a].y;
    w0[2] = -p[T.a].z;
    
    a = n[0]*w0[0]+n[1]*w0[1]+n[2]*w0[2];    //a = dot3D(n,w0);
    b = n[0]*dir[0]+n[1]*dir[1]+n[2]*dir[2]; //b = dot3D(n,dir);
    
    if (b>-EPSILON && b<EPSILON) { // ray is parallel to triangle plane
        if (a == 0.0)              // ray lies in triangle plane
            return 2;
        else
            return 0;              // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = -a/b;
    if (r < 0.0)                    // ray goes away from triangle
        return 0;                   // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    xx[0]=dir[0]*r;
    xx[1]=dir[1]*r;
    xx[2]=dir[2]*r; // intersect point of ray and plane

    // is I inside T?
    uu=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
    uv=u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
    vv=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    w[0]=xx[0]-p[T.a].x;
    w[1]=xx[1]-p[T.a].y;
    w[2]=xx[2]-p[T.a].z;
    wu=w[0]*u[0]+w[1]*u[1]+w[2]*u[2];
    wv=w[0]*v[0]+w[1]*v[1]+w[2]*v[2];
    D = uv * uv - uu * vv;

    // get and test parametric coords
    ss = (uv * wv - vv * wu) / D;
    if(ss>-EPSILON && ss<EPSILON) ss=0;
    if((1-ss)>-EPSILON && (1-ss)<EPSILON) ss=1;

    tt = (uv * wu - uu * wv) / D;
    if(tt>-EPSILON && tt<EPSILON)   tt=0;
    if((1-tt)>-EPSILON && (1-tt)<EPSILON) tt=1;

    *c0=(float)ss;
    *c1=(float)tt;

    if (ss < 0.0 || tt < 0.0 || (ss + tt) > 1.0)  // I is outside T
        return 0;

    return 1;                      // I is in T
}
void neighbours(Mesh *m)
{
    // find incident triangles for every vertex
    int     np=m->np;
    int     nt=m->nt;
    int3D   *t=m->t;
    NTriRec **NT=&(m->NT);
    int     i;
        
    if(*NT)
        free(*NT);
    *NT=(NTriRec*)calloc(np,sizeof(NTriRec));
    if(*NT==NULL)
    {
        printf("ERROR: Cannot create NT structure in neighbours() function\n");
        return;
    }       
    for(i=0;i<nt;i++)
    {
        ((*NT)[t[i].a]).t[((*NT)[t[i].a]).n++] = i;
        ((*NT)[t[i].b]).t[((*NT)[t[i].b]).n++] = i;
        ((*NT)[t[i].c]).t[((*NT)[t[i].c]).n++] = i;
    }
}
#define TINY 1.0e-10    // A small number.
#define NMAX 500000     //Maximum allowed number of function evaluations.
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
float amotry(float *p, float y[], float psum[], int ndim,float (*funk)(float []), int ihi, float fac)
{
    int     j;
    float   fac1,fac2,ytry,*ptry;

    ptry=(float*)calloc(ndim,sizeof(float));
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=0;j<ndim;j++)
        ptry[j]=psum[j]*fac1-p[ndim*ihi+j]*fac2;
    ytry=(*funk)(ptry); // Evaluate the function at the trial point.
    if (ytry < y[ihi])
    {
        // If it's better than the highest, then replace the highest.
        y[ihi]=ytry;
        for (j=0;j<ndim;j++)
        {
            psum[j] += ptry[j]-p[ndim*ihi+j];
            p[ndim*ihi+j]=ptry[j];
        }
    }
    free(ptry);
    return ytry;
}
#define GET_PSUM  for(j=0;j<ndim;j++){for(sum=0.0,i=0;i<mpts;i++)sum+=p[ndim*i+j];psum[j]=sum;}
void amoeba(float *p, float y[], int ndim, float ftol,float (*funk)(float []),int *nfunk)
// Multidimensional minimization of the function funk(x) where x[0..ndim-1] is a vector in ndim
// dimensions, by the downhill simplex method of Nelder and Mead. From Numerical Recipes in C
{
    int     i,ihi,ilo,inhi,j,mpts=ndim+1;
    float   rtol,sum,swap,ysave,ytry,*psum;

    psum=(float*)calloc(ndim,sizeof(float));
    *nfunk=0;
    GET_PSUM
    for (;;)
    {
        ilo=0;
        ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (i=0;i<mpts;i++)
        {
            if (y[i] <= y[ilo])
                ilo=i;
            if (y[i] > y[ihi])
            {
                inhi=ihi;
                ihi=i;
            }
            else if (y[i] > y[inhi] && i != ihi)
                inhi=i;
        }
        rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        if (rtol < ftol)
        {
            SWAP(y[0],y[ilo])
            for (i=0;i<ndim;i++) SWAP(p[ndim*(0)+i],p[ndim*ilo+i])
            break;
        }
        if (*nfunk >= NMAX){break;};
        *nfunk += 2;
        ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
        if (ytry <= y[ilo])
            ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi])
        {
            ysave=y[ihi];
            ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave)
            {
                for (i=0;i<mpts;i++)
                if (i != ilo)
                {
                    for (j=0;j<ndim;j++)
                        p[ndim*i+j]=psum[j]=0.5*(p[ndim*i+j]+p[ndim*(ilo-1)+j]);
                    y[i]=(*funk)(psum);
                }
                *nfunk += ndim; // Keep track of function evaluations.
                GET_PSUM // Recompute psum.
            }
        }
        else
            --(*nfunk); // Correct the evaluation count.
    } // Go back for the test of doneness and the next
    free(psum);
}
#pragma mark -
#pragma mark [ Format conversion ]
int getformatindex(char *path)
{
    char    *formats[]={"orig",     "pial",  "white",    "mesh",    "sratio",
                        "float",    "curv",  "txt",      "inflated","sphere",
                        "sulc",     "reg",   "txt1",     "wrl",     "obj",
                        "ply",      "stl",   "smesh",    "off",     "bin",
                        "mgh",      "annot", "raw",      "vtk",     "dpv",
                        "civet_obj","gii",   "data_gii", "asc"};
    int     i,n=sizeof(formats)/sizeof(long); // number of recognised formats
    int     found,index;
    char    *extension;
    
    for(i=strlen(path);i>=0;i--)
        if(path[i]=='.')
            break;
    if(i==0)
    {
        printf("ERROR: Unable to find the format extension\n");
        return 0;
    }
    extension=path+i+1;
    
    for(i=0;i<n;i++)
    {
        found=(strcmp(formats[i],extension)==0);
        if(found)
            break;
    }
    
    index=-1;
    if(i==0 || i==1 || i==2 || i==8 || i==9 ||i==11)
    {
        index=kFreeSurferMesh;
        if(verbose)
            printf("Format: FreeSurfer mesh\n");
    }
    else
    if(i==21)
    {
        index=kFreeSurferAnnot;
        if(verbose)
            printf("Format: FreeSurfer annot\n");
    }
    else
    if(i==3)
    {
        index=kBrainVisaMesh;
        if(verbose)
            printf("Format: BrainVisa mesh\n");
    }
    else
    if(i==4 || i==6 || i==10)
    {
        index=kFreeSurferData;
        if(verbose)
            printf("Format: FreeSurfer Data\n");
    }
    else
    if(i==5)
    {
        index=kFloatData;
        if(verbose)
            printf("Format: Float Data\n");
    }
    else
    if(i==22)
    {
        index=kRawFloatData;
        if(verbose)
            printf("Format: Raw Float Data\n");
    }
    else
    if(i==12)
    {
        index=kTextData;
        if(verbose)
            printf("Format: Text Data\n");
    }
    else
    if(i==7)
    {
        index=kText;
        if(verbose)
            printf("Format: Text mesh\n");
    }
    else
    if(i==13)
    {
        index=kVRMLMesh;
        if(verbose)
            printf("Format: VRML mesh\n");
    }
    else
    if(i==14)
    {
        index=kObjMesh;
        if(verbose)
            printf("Format: Obj mesh\n");
    }
    else
    if(i==15)
    {
        index=kPlyMesh;
        if(verbose)
            printf("Format: Ply mesh\n");
    }
    else
    if(i==16)
    {
        index=kSTLMesh;
        if(verbose)
            printf("Format: STL mesh\n");
    }
    else
    if(i==17)
    {
        index=kSmeshMesh;
        if(verbose)
            printf("Format: Smesh mesh\n");
    }
    else
    if(i==18)
    {
        index=kOffMesh;
        if(verbose)
            printf("Format: Off mesh\n");
    }
    else
    if(i==19)
    {
        index=kBinMesh;
        if(verbose)
            printf("Format: Bin mesh\n");
    }
    else
    if(i==20)
    {
        index=kMGHData;
        if(verbose)
            printf("Format: MGH data\n");
    }
    else
    if(i==23)
    {
        index=kVTKMesh;
        if(verbose)
            printf("Format: VTK Mesh\n");
    }
    else
    if(i==24)
    {
        index=kDPVData;
        if(verbose)
            printf("Format: DPV Data\n");
    }
    else
    if(i==25)
    {
        index=kCivetObjMesh;
        if(verbose)
            printf("Format: Civet Obj Mesh\n");
    }
    else
    if(i==26)
    {
        index=kGiiMesh;
        if(verbose)
            printf("Format: Gii Mesh\n");
    }
    else
    if(i==27)
    {
        index=kGiiData;
        if(verbose)
            printf("Format: Gii Data\n");
    }
    else
    if(i==28)
    {
        index=kAscMesh;
        if(verbose)
            printf("Format: Asc Mesh\n");
    }
        
    return index;
}
#pragma mark -
int FreeSurfer_load_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     id,a,b,c;
    char    date[256],info[256];


    f=fopen(path,"r");
    
    if(f==NULL)
        return 1;

    // read triangle/quad identifier: 3 bytes
    a=((int)(u_int8_t)fgetc(f))<<16;
    b=((int)(u_int8_t)fgetc(f))<<8;
    c=(u_int8_t)fgetc(f);
    id=a+b+c;
    if(id==16777214)    // triangle mesh
    {
        fgets(date,256,f);
        fgets(info,256,f);
        fread(np,1,sizeof(int),f);    if(endianness==kINTEL) swapint(np);
        fread(nt,1,sizeof(int),f);    if(endianness==kINTEL) swapint(nt);
        // read vertices
        *p=(float3D*)calloc(*np,sizeof(float3D));
            if((*p)==NULL) printf("ERROR: Cannot allocate memory for points [FreeSurfer_load_mesh]\n");
            else
            {
        fread((char*)(*p),*np,3*sizeof(float),f);      if(endianness==kINTEL) swapvertices(m);
        // read triangles
        *t=(int3D*)calloc(*nt,sizeof(int3D));
            if((*t)==NULL) printf("ERROR: Cannot allocate memory for triangles [FreeSurfer_load_mesh]\n");
            else
            {
        fread((char*)(*t),*nt,3*sizeof(int),f);        if(endianness==kINTEL) swaptriangles(m);
            }
            }
    }
    fclose(f);
    
    return 0;
}
int FreeSurfer_load_data(char *path, Mesh *m)
{
    int     *np=&(m->np);
    float   **data=&(m->data);
    FILE    *f;
    int     i,j;
    int     id,a,b,c;
    char    byte4[4];

    if(verbose)
        printf("* FreeSurfer_load_data\n");

    f=fopen(path,"r");
    if(f==NULL)
        return 1;

    // read identifier: 3 bytes
    a=((int)(u_int8_t)fgetc(f))<<16;
    b=((int)(u_int8_t)fgetc(f))<<8;
    c=(u_int8_t)fgetc(f);
    id=a+b+c;
    if(id==16777215)    // triangle mesh
    {
        if(endianness==kINTEL)
            for(i=0;i<4;i++) byte4[3-i]=fgetc(f);
        else
            fread(byte4,4,sizeof(char),f);
        *np=*(int*)byte4;
        if(verbose)
            printf("FS #vertex_data %i\n",*np);
        
        *data=(float*)calloc(*np,sizeof(float));
        
        // disregard FaceCount and ValsPerVertex
        fgetc(f);fgetc(f);fgetc(f);fgetc(f);
        fgetc(f);fgetc(f);fgetc(f);fgetc(f);
        
        // read vertex data
        for(j=0;j<*np;j++)
        {
            if(endianness==kINTEL)
                for(i=0;i<4;i++) byte4[3-i]=fgetc(f);
            else
                fread(byte4,4,sizeof(char),f);
            (*data)[j]=*(float*)byte4;
        }
    }
    if(verbose)
        printf("FSData finished\n");
    
    fclose(f);
    
    return 0;
}
int FreeSurfer_load_annot(char *path, Mesh *m)
{
    if(verbose)
        printf("* FreeSurfer_load_annot\n");

    FILE    *f;
    int     i,n,l;
    char    *tmp;
    float   **data=&(m->data);

    f=fopen(path,"r");
    if(f==NULL)
        return 0;
    
    fread(&n,1,sizeof(int),f);
    if(endianness==kINTEL)
        swapint(&n);
    if(m->np==0)
        m->np=n;
    if(n!=m->np)
    {
        printf("ERROR: Annotation file corrupted. points:%i annotations:%i [FreeSurfer_load_annot]\n",m->np,n);
        return 1;
    }

    m->ddim=3;
    tmp=calloc(m->np,2*sizeof(int));
    *data=calloc(m->np,3*sizeof(float));
    if(tmp==NULL)
    {
        printf("ERROR: Cannot allocate memory [FreeSurfer_load_annot]\n");
        return 1;
    }

    fread(tmp,m->np,2*sizeof(int),f);
    for(i=0;i<min(m->np,n);i++)
    {
        l=((int*)tmp)[2*i+1];
        if(endianness==kINTEL)
            swapint(&l);
        (*data)[3*i+0]=(l&0xff);
        (*data)[3*i+1]=((l>>8)&0xff);
        (*data)[3*i+2]=((l>>16)&0xff);
    }
    free(tmp);
    fclose(f);

    return 0;
}
int FreeSurfer_load_mghdata(char *path, Mesh *m)
{
    // path: path to source thickness file in mgh format (non-compressed version of mgz)
    int     *np=&(m->np);
    float   **data=&(m->data);
    FILE    *f;
    int     v,ndim1,ndim2,ndim3,nframes,type,dof;
    int     i;

    f=fopen(path,"r");
    if(f==NULL)
        return 1;
    
    fread(&v,1,sizeof(int),f);          swapint(&v);
    fread(&ndim1,1,sizeof(int),f);      swapint(&ndim1);
    fread(&ndim2,1,sizeof(int),f);      swapint(&ndim2);
    fread(&ndim3,1,sizeof(int),f);      swapint(&ndim3);
    fread(&nframes,1,sizeof(int),f);    swapint(&nframes);
    fread(&type,1,sizeof(int),f);       swapint(&type);
    fread(&dof,1,sizeof(int),f);        swapint(&dof);
    
    if(verbose)
    {
        printf("version:%i\n",v);
        printf("ndim1:%i\n",ndim1);
        printf("ndim2:%i\n",ndim2);
        printf("ndim3:%i\n",ndim3);
        printf("nframes:%i\n",nframes);
        printf("type:%i\n",type);
        printf("dof:%i\n\n",dof);
    }
    
    *np=ndim1*ndim2*ndim3;
    *data=(float*)calloc(*np,sizeof(float));
    fseek(f,64*4,SEEK_CUR);
    for(i=0;i<*np;i++)
    {
        fread(&((*data)[i]),1,sizeof(float),f);
        swapfloat(&((*data)[i]));
    }
    fclose(f);
    
    return 0;
}

int FreeSurfer_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     id=16777214,a,b,c;
    int     NP,NT,i;
    char    date[6]="EMPTY",info[6]="EMPTY";
    float3D ftmp;
    int3D   itmp;

    f=fopen(path,"w");
    
    if(f==NULL)
        return 1;
    
    // write data identifier: 3 bytes
    a=id>>16;
    b=(id&0xff00)>>8;
    c=(id&0xff);
    fputc((char)a,f);
    fputc((char)b,f);
    fputc((char)c,f);
    
    // write date and info (EMPTY)
    date[5]=(char)10;
    info[5]=(char)10;
    fwrite(date,1,6,f);
    fwrite(info,1,6,f);

    // write number of vertices and triangles
    NP=*np;
    NT=*nt;
    if(endianness==kINTEL)
    {
        swapint(&NP);
        swapint(&NT);
    }
    fwrite(&NP,1,sizeof(int),f);
    fwrite(&NT,1,sizeof(int),f);

    // write vertices and triangles
    if(endianness==kINTEL)
    {
        for(i=0;i<*np;i++)
        {
            ftmp=p[i];
            swapfloat(&ftmp.x);
            swapfloat(&ftmp.y);
            swapfloat(&ftmp.z);
            fwrite(&ftmp,1,sizeof(float3D),f);
        }
        for(i=0;i<*nt;i++)
        {
            itmp=t[i];
            swapint(&itmp.a);
            swapint(&itmp.b);
            swapint(&itmp.c);
            fwrite(&itmp,1,sizeof(int3D),f);
        }
    }
    else
    {
        fwrite(p,*np,sizeof(float3D),f);
        fwrite(t,*nt,sizeof(int3D),f);
    }
    fclose(f);

    return 0;
}
int FreeSurfer_save_data(char *path, Mesh *m)
{
    int     *np=&(m->np);
    float   *data=m->data;
    FILE    *f;
    int     id=16777215,a,b,c;
    int     FaceCount=0,ValsPerVertex=1,i,n;
    float   x;

    f=fopen(path,"w");
    
    if(f==NULL)
        return 1;
    
    // write data identifier: 3 bytes
    a=id>>16;
    b=(id&0xff00)>>8;
    c=(id&0xff);
    fputc((char)a,f);
    fputc((char)b,f);
    fputc((char)c,f);
    
    n=*np;
    if(endianness==kINTEL)
    {
        swapint(&n);
        swapint(&FaceCount);
        swapint(&ValsPerVertex);
    }
    fwrite(&n,1,sizeof(int),f);
    fwrite(&FaceCount,1,sizeof(int),f);
    fwrite(&ValsPerVertex,1,sizeof(int),f);

    // write data
    if(endianness==kINTEL)
    {
        for(i=0;i<*np;i++)
        {
            x=data[i];
            swapfloat(&x);
            fwrite(&x,1,sizeof(float),f);
        }
    }
    else
        fwrite(data,*np,sizeof(float),f);
    fclose(f);

    return 0;
}
int FreeSurfer_save_mghdata(char *path, Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    FILE    *f;
    int     v,ndim1,ndim2,ndim3,nframes,type,dof;
    int     i;
    float   x;

    f=fopen(path,"w");
    if(f==NULL)
        return 1;
    
    v=1;
    ndim1=np;
    ndim2=1;
    ndim3=1;
    nframes=1;
    type=3;
    dof=0;
    
    swapint(&v);        fwrite(&v,1,sizeof(int),f);
    swapint(&ndim1);    fwrite(&ndim1,1,sizeof(int),f);
    swapint(&ndim2);    fwrite(&ndim2,1,sizeof(int),f);
    swapint(&ndim3);    fwrite(&ndim3,1,sizeof(int),f);
    swapint(&nframes);  fwrite(&nframes,1,sizeof(int),f);
    swapint(&type);     fwrite(&type,1,sizeof(int),f);
    swapint(&dof);      fwrite(&dof,1,sizeof(int),f);
    
    fseek(f,64*4,SEEK_CUR);
    for(i=0;i<np;i++)
    {
        x=data[i];
        swapfloat(&x);
        fwrite(&x,1,sizeof(float),f);
        
    }
    fclose(f);
    
    return 0;
}
int BrainVisa_load_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    char    tmp[6];
    int     i;
    int     endian,ignore;
    
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // READ HEADER
    // get format (ascii, binar)
    fread(tmp,5,sizeof(char),f); tmp[5]=(char)0;
    if(strcmp(tmp,"binar")==0)
    {
        for(i=0;i<4;i++) tmp[i]=fgetc(f); tmp[4]=(char)0;
        endian=-1;
        if(strcmp(tmp,"ABCD")==0)    endian=kMOTOROLA;    
        if(strcmp(tmp,"DCBA")==0)    endian=kINTEL;
        if(endian==-1){ printf("ERROR: Not ABCD nor DCBA order...exit.\n"); return 1;}
        fread(&ignore,4,sizeof(char),f);        // ignore "VOID" string length
        fread(&ignore,4,sizeof(char),f);        // ignore "VOID" string
        fread(&ignore,1,sizeof(int),f);         // verify number of vertices per polygon
        if(endian!=endianness)
            swapint(&ignore);
        if(ignore!=3){ printf("ERROR: Only able to read triangle meshes. This mesh has %i vertices per polygon.\n",ignore); return 1;}
        fread(&ignore,1,sizeof(int),f);         // ignore time steps
        fread(&ignore,1,sizeof(int),f);         // ignore time step index
        
        // READ VERTICES
        fread(np,1,sizeof(int),f);              // read number of vertices
        if(endian!=endianness)
            swapint(np);
        (*p) = (float3D*)calloc(*np,sizeof(float3D));
        if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
        fread((char*)(*p),*np,3*sizeof(float),f);  if(endian!=endianness) swapvertices(m);    
        if(verbose)
            printf("Read %i vertices\n",*np);
        
        // IGNORE NORMAL VECTORS
        fseek(f,sizeof(int),SEEK_CUR);          // ignore normal vectors
        fseek(f,*np*sizeof(float3D),SEEK_CUR);
        fread(&ignore,1,sizeof(int),f);         // ignore number of texture coordinates
        
        // READ TRIANGLES
        fread(nt,1,sizeof(int),f);              // read number of triangles
        if(endian!=endianness)
            swapint(nt);
        (*t) = (int3D*)calloc(*nt,sizeof(int3D));
        if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
        fread((char*)(*t),*nt,3*sizeof(int),f);    if(endian!=endianness) swaptriangles(m);
        if(verbose)
            printf("Read %i triangles\n",*nt);
    }
    else
    if(strcmp(tmp,"ascii")==0)
    {
        fscanf(f," %*s ");            // ignore VOID
        fscanf(f," %*i %*i %*i ");    // ignore 3 integers
        
        // READ 3-D COORDINATES
        fscanf(f," %i ",np);
        *p=(float3D*)calloc(*np,sizeof(float3D));
        for(i=0;i<*np;i++)
            fscanf(f," ( %f , %f , %f ) ", &((*p)[i].x),&((*p)[i].y),&((*p)[i].z));
        
        fscanf(f," %*i ");            // ignore number of normal vectors
        for(i=0;i<*np;i++)            // ignore normal vectors
            fscanf(f," ( %*f , %*f , %*f ) ");
        fscanf(f," %*i ");            // ignore an integer
        
        // READ TRIANGLES
        fscanf(f," %i ",nt);
        *t=(int3D*)calloc(*nt,sizeof(int3D));
        for(i=0;i<*nt;i++)
            fscanf(f," ( %i , %i , %i ) ", &((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
    }
    else
    {
        printf("ERROR: Cannot read '%s' format.\n",tmp);
        return 1;
    }
    
    fclose(f);
    
    return 0;
}
int BrainVisa_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
        fprintf(f,"ascii\n");
    
        fprintf(f,"VOID\n");    // ignore VOID
        fprintf(f,"3\n1\n0\n");    // ignore 3 integers
        
        // WRITE 3-D COORDINATES
        fprintf(f,"%i\n",*np);
        for(i=0;i<*np;i++)
            fprintf(f,"(%f,%f,%f) ",p[i].x,p[i].y,p[i].z);
        fprintf(f,"\n");
        fprintf(f,"0\n");            // ignore number of normal vectors
        fprintf(f,"0\n");            // ignore an integer
        
        // WRITE TRIANGLES
        fprintf(f,"%i\n",*nt);
        for(i=0;i<*nt;i++)
            fprintf(f,"(%i,%i,%i) ",t[i].a,t[i].b,t[i].c);
        fprintf(f,"\n");
    fclose(f);
    
    return 0;
}
int Text_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt),nt_tmp;
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    float   **data=&(m->data);
    FILE    *f;
    int     i;
    char    str[512];
    
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // READ HEADER
    fgets(str,511,f);
    sscanf(str," %i %i ",np,&nt_tmp);
    
    if(nt_tmp==1)    // mesh data file, dimension 1
    {
        *data=(float*)calloc(*np,sizeof(float));
        if(data==NULL){printf("ERROR: Not enough memory for mesh data\n");return 1;}
        for(i=0;i<*np;i++)
            fscanf(f," %f ",&((*data)[i]));    
        if(verbose)
            printf("Read %i data values\n",*np);

    }
    else
    {   // mesh file
        // READ VERTICES
        *p=(float3D*)calloc(*np,sizeof(float3D));
        if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
        for(i=0;i<*np;i++)
            fscanf(f," %f %f %f ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));    
        if(verbose)
            printf("Read %i vertices\n",*np);
    
        // READ TRIANGLES
        *nt=nt_tmp;
        *t = (int3D*)calloc(*nt,sizeof(int3D));
        if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
        for(i=0;i<*nt;i++)
            fscanf(f," %i %i %i ",&((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
        if(verbose)
            printf("Read %i triangles\n",*nt);
    }

    fclose(f);
    
    return 0;
}
int Text_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // WRITE HEADER
    fprintf(f,"%i %i\n",*np,*nt);

    // WRITE VERTICES
    for(i=0;i<*np;i++)
        fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);    

    // WRITE TRIANGLES
    for(i=0;i<*nt;i++)
        fprintf(f,"%i %i %i\n",t[i].a,t[i].b,t[i].c);

    fclose(f);
    
    return 0;
}
int Asc_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt),nt_tmp;
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     i;
    char    str[512];
    
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // READ HEADER
    fgets(str,511,f); // ignore identification line
    fgets(str,511,f);
    sscanf(str," %i %i ",np,&nt_tmp);
    
    // READ VERTICES
    *p=(float3D*)calloc(*np,sizeof(float3D));
    if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
    for(i=0;i<*np;i++)
        fscanf(f," %f %f %f %*i ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));    
    if(verbose)
        printf("Read %i vertices\n",*np);

    // READ TRIANGLES
    *nt=nt_tmp;
    *t = (int3D*)calloc(*nt,sizeof(int3D));
    if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
    for(i=0;i<*nt;i++)
        fscanf(f," %i %i %i %*i ",&((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
    if(verbose)
        printf("Read %i triangles\n",*nt);

    fclose(f);
    
    return 0;
}
int Asc_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // WRITE HEADER
    fprintf(f,"#!ascii by meshgeometry with love\n"); // identification line
    fprintf(f,"%i %i\n",*np,*nt);

    // WRITE VERTICES
    for(i=0;i<*np;i++)
        fprintf(f,"%f %f %f 0\n",p[i].x,p[i].y,p[i].z);

    // WRITE TRIANGLES
    for(i=0;i<*nt;i++)
        fprintf(f,"%i %i %i 0\n",t[i].a,t[i].b,t[i].c);

    fclose(f);

    return 0;
}
int Text_save_data(char *path, Mesh *m)
{
    if(verbose)
        printf("* Text_save_data\n");

    int     *np=&(m->np);
    float   *data=m->data;
    FILE    *f;
    int     i,j;
    
    f=fopen(path,"w");
    if(f==NULL)
    {
        printf("ERROR: Cannot open file\n");
        return 1;
    }
    
    // WRITE HEADER
    fprintf(f,"%i %i 3\n",*np,m->ddim);

    // WRITE DATA
    for(i=0;i<*np;i++)
    {
        for(j=0;j<m->ddim;j++)
            if(j<m->ddim-1)
                fprintf(f,"%f ",data[m->ddim*i+j]);
            else
                fprintf(f,"%f\n",data[m->ddim*i+j]);
    }

    fclose(f);
    
    return 0;
}
int VRML_load_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     i,loop;
    char    str[256];

    f=fopen(path,"r");

    *np=0;
    *nt=0;
    
    loop=1;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"point"))
            loop=0;
    }
    loop=1;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"]"))
            loop=0;
        else
            (*np)++;
    }
    loop=1;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"coordIndex"))
            loop=0;
    }
    loop=1;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"]"))
            loop=0;
        else
            (*nt)++;
    }

    *p = (float3D*)calloc(*np,sizeof(float3D));
    *t = (int3D*)calloc(*nt,sizeof(int3D));
    fseek(f,0,SEEK_SET);
   
    loop=1;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"point"))
            loop=0;
    }
    loop=1;
    i=0;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"]"))
            loop=0;
        else
        {
            sscanf(str," %f %f %f ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));
            i++;
        }
    }
    loop=1;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"coordIndex"))
            loop=0;
    }
    loop=1;
    i=0;
    while(loop)
    {
        fgets(str,255,f);
        if(strstr(str,"]"))
            loop=0;
        else
        {
            sscanf(str," %i , %i , %i ",&((*t)[i].a),&((*t)[i].c),&((*t)[i].b));
            i++;
        }
    }
    fclose(f);

    return 0;
}
int VRML_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    fprintf(f,"#VRML V1.0 ascii\n");
    fprintf(f,"Separator {\n");
    fprintf(f,"ShapeHints {\n");
    fprintf(f,"vertexOrdering COUNTERCLOCKWISE\n");
    fprintf(f,"faceType CONVEX\n");
    fprintf(f,"}\n");
    fprintf(f,"Coordinate3 {\n");
    fprintf(f,"point [\n");
    for(i=0;i<*np;i++)
        fprintf(f,"%f %f %f,\n",p[i].x,p[i].y,p[i].z);
    fprintf(f,"]\n");
    fprintf(f,"}\n");
    fprintf(f,"IndexedFaceSet {\n");
    fprintf(f,"coordIndex [\n");
    for(i=0;i<*nt;i++)
        fprintf(f,"%i,%i,%i,-1\n",t[i].a,t[i].b,t[i].c);
    fprintf(f,"]\n");
    fprintf(f,"}\n");
    fprintf(f,"}\n");
    fclose(f);

    return 0;
}
int CivetObj_load_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     i;
    char    str[512];
    
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // READ HEADER
    fgets(str,511,f);
    sscanf(str," %*c %*f %*f %*f %*i %*i %i ",np);
    
    // READ VERTICES
    *p = (float3D*)calloc(*np,sizeof(float3D));
    if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
    for(i=0;i<*np;i++)
        fscanf(f," %f %f %f ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));    
    if(verbose)
        printf("Read %i vertices\n",*np);
    
    // IGNORE NORMALS
//    fgets(str,511,f); // skip empty line
    for(i=0;i<*np;i++)
        fscanf(f," %*f %*f %*f ");    

    // READ TRIANGLES
//    fgets(str,511,f);     // skip empty line
    fscanf(f," %i ",nt);
//    fgets(str,511,f);     // skip empty line
    for(i=0;i<5+*nt;i++)    // skip nt+5 integers (the first 5, all integers, then nt multiples of 3)
        fscanf(f," %*i ");
    *t = (int3D*)calloc(*nt,sizeof(int3D));
    if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
    for(i=0;i<*nt;i++)
        fscanf(f," %i %i %i ",&((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
    if(verbose)
        printf("Read %i triangles\n",*nt);

    fclose(f);
    
    return 0;
}
int CivetObj_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // WRITE HEADER
    fprintf(f,"P 0.3 0.3 0.4 10 1 %i\n",*np);

    // WRITE VERTICES
    for(i=0;i<*np;i++)
        fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    fprintf(f,"\n");
    
    // WRITE DUMMY NORMALS
    for(i=0;i<*np;i++)
        fprintf(f,"0 0 0\n");
    fprintf(f,"\n");
    
    // WRITE NUMBER OF TRIANGLES
    fprintf(f,"%i\n",*nt);
    
    // WRITE 5 DUMMY NUMBERS
    fprintf(f,"0 1 1 1 1 1\n");
    fprintf(f,"\n");
    
    // WRITE nt MULTIPLES OF 3, IN ROWS OF EIGHT
    for(i=0;i<*nt;i++)
    {
        fprintf(f,"%i ",(i+1)*3);
        if(i>0 && i%8==0)
            fprintf(f,"\n");
    }
    fprintf(f,"\n");

    // WRITE TRIANGLES
    for(i=0;i<*nt;i++)
        fprintf(f,"%i %i %i\n",t[i].a,t[i].b,t[i].c);

    fclose(f);
    
    return 0;
}
int Obj_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    char    str[1024],s[16];
    int     n;
    
    f=fopen(path,"r");
    *np=*nt=0;
    while(!feof(f))
    {
        fgets(str,1024,f);
        sscanf(str," %s ",s);
        if(strcmp(s,"v")==0)
            (*np)++;
        if(strcmp(s,"f")==0)
            (*nt)++;
    }
    fclose(f);
    
    *p = (float3D*)calloc(*np,sizeof(float3D));
    if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
    *t = (int3D*)calloc(*nt,sizeof(int3D));
    if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n");return 1;}

    f=fopen(path,"r"); 
    *np=*nt=0;
    while(!feof(f))
    {
        fgets(str,1024,f);
        sscanf(str," %s ",s);
        if(strcmp(s,"v")==0)
        {
            sscanf(str," %*s %f %f %f ",&((*p)[*np].x),&((*p)[*np].y),&((*p)[*np].z));
            (*np)++;
        }
        if(strcmp(s,"f")==0)
        {
            n=sscanf(str,"f %i %i %i ",&((*t)[*nt].a),&((*t)[*nt].b),&((*t)[*nt].c));
            if(n<3)
                n=sscanf(str,"f %i/%*i %i/%*i %i/%*i ",&((*t)[*nt].a),&((*t)[*nt].b),&((*t)[*nt].c));
            if(n<3)
                n=sscanf(str,"f %i//%*i %i//%*i %i//%*i ",&((*t)[*nt].a),&((*t)[*nt].b),&((*t)[*nt].c));
            (*t)[*nt].a--;
            (*t)[*nt].b--;
            (*t)[*nt].c--;
            (*nt)++;
        }
    }
    if(verbose)
        printf("Read %i vertices and %i triangles\n",*np,*nt);

    return 0;
}
int Obj_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    for(i=0;i<*np;i++)
        fprintf(f,"v %f %f %f\n",p[i].x,p[i].y,p[i].z);
    
    for(i=0;i<*nt;i++)
        fprintf(f,"f %i %i %i\n",t[i].a+1,t[i].b+1,t[i].c+1);

    fclose(f);
    
    return 0;
}
int Ply_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     i,x;
    char    str[512],str1[256],str2[256];
        
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}

    // READ HEADER
    *np=*nt=0;
    do
    {
        fgets(str,511,f);
        sscanf(str," %s %s %i ",str1,str2,&x);
        if(strcmp(str1,"element")==0&&strcmp(str2,"vertex")==0)
            *np=x;
        else
        if(strcmp(str1,"element")==0&&strcmp(str2,"face")==0)
            *nt=x;
    }
    while(strcmp(str1,"end_header")!=0 && !feof(f));
    if((*np)*(*nt)==0)
    {
        printf("ERROR: Bad Ply file header format\n");
        return 1;
    }
    // READ VERTICES
    *p = (float3D*)calloc(*np,sizeof(float3D));
    if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
    for(i=0;i<*np;i++)
        fscanf(f," %f %f %f ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));    
    if(verbose)
        printf("Read %i vertices\n",*np);

    // READ TRIANGLES
    *t = (int3D*)calloc(*nt,sizeof(int3D));
    if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
    for(i=0;i<*nt;i++)
        fscanf(f," 3 %i %i %i ",&((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
    if(verbose) {
        printf("Read %i triangles\n",*nt);
        printf("%i %i %i\n",(*t)[0].a,(*t)[0].b,(*t)[0].c);
    }

    fclose(f);
    
    return 0;
}
int Ply_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // WRITE HEADER
    fprintf(f,"ply\n");
    fprintf(f,"format ascii 1.0\n");
    fprintf(f,"comment meshconvert, R. Toro 2010\n");
    fprintf(f,"element vertex %i\n",*np);
    fprintf(f,"property float x\n");
    fprintf(f,"property float y\n");
    fprintf(f,"property float z\n");
    fprintf(f,"element face %i\n",*nt);
    fprintf(f,"property list uchar int vertex_indices\n");
    fprintf(f,"end_header\n");

    // WRITE VERTICES
    for(i=0;i<*np;i++)
        fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);    

    // WRITE TRIANGLES
    for(i=0;i<*nt;i++)
        fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);

    fclose(f);
    
    return 0;
}
int STL_load(char *path, Mesh *m)
{
    printf("ERROR: meshconvert DOES NOT LOAD STL DATA FOR THE MOMENT, ONLY SAVES IT\n");
    return 1;
    /*
    FILE    *f;
    int        i;
    char    str[512];
    
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // SKIP HEADER
    fgets(str,511,f);
    
        p = (float3D*)calloc(np,sizeof(float3D));
        if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
        for(i=0;i<np;i++)
            fscanf(f," %f %f %f ",&p[i].x,&p[i].y,&p[i].z);    
        printf("Read %i vertices\n",np);
    
        // READ TRIANGLES
        t = (int3D*)calloc(nt,sizeof(int3D));
        if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
        for(i=0;i<nt;i++)
            fscanf(f," %i %i %i ",&t[i].a,&t[i].b,&t[i].c);
        printf("Read %i triangles\n",nt);

    fclose(f);
    */
    
    return 0;
}
int STL_save_mesh(char *path, Mesh *m)
{
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    float3D n;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // WRITE HEADER
    fprintf(f,"solid mySolid\n");

    // WRITE FACES
    for(i=0;i<*nt;i++)
    {
        n=normal3D(i,m);
        fprintf(f,"facet normal %e %e %e\n",n.x,n.y,n.z);
        fprintf(f,"outer loop\n");
        fprintf(f,"vertex %e %e %e\n",p[t[i].a].x,p[t[i].a].y,p[t[i].a].z);
        fprintf(f,"vertex %e %e %e\n",p[t[i].b].x,p[t[i].b].y,p[t[i].b].z);
        fprintf(f,"vertex %e %e %e\n",p[t[i].c].x,p[t[i].c].y,p[t[i].c].z);
        fprintf(f,"endLoop\n");
        fprintf(f,"endFacet\n");
    }
    fprintf(f,"endSolid mySolid\n");
    fclose(f);
    
    return 0;
}
int Smesh_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     i;
    char    str[512];
    
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // READ POINTS HEADER
    fgets(str,511,f);
    sscanf(str," %i ",np);
    // READ VERTICES
    *p = (float3D*)calloc(*np,sizeof(float3D));
    if(*p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
    for(i=0;i<*np;i++)
        fscanf(f," %*i %f %f %f ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));    
    if(verbose)
        printf("Read %i vertices\n",*np);
    
    // READ TRIANGLES HEADER
    // READ TRIANGLES
    fgets(str,511,f);
    sscanf(str," %i ",nt);
    *t = (int3D*)calloc(*nt,sizeof(int3D));
    if(*t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
    for(i=0;i<*nt;i++)
        fscanf(f," %*i %i %i %i ",&((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
    if(verbose)
        printf("Read %i triangles\n",*nt);

    fclose(f);
    return 0;
}
int Smesh_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    // WRITE VERTICES HEADER
    fprintf(f,"%i 3 0 0\n",*np);
    // WRITE VERTICES
    for(i=0;i<*np;i++)
        fprintf(f,"%i %f %f %f 0\n",i,p[i].x,p[i].y,p[i].z);    
    // WRITE TRIANGLES HEADER
    fprintf(f,"%i 0\n",*nt);
    // WRITE TRIANGLES
    for(i=0;i<*nt;i++)
        fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);
    fprintf(f,"0\n0\n");
    fclose(f);   
    return 0;
}
int Bin_load(char *path, Mesh *m)
{
    printf("ERROR: Load is not yet implemented for Bin filetype\n");
    return 1;
}
int Bin_save_mesh(char *path, Mesh *m)
{
    // Bin filetype is the binary format used by n-e-r-v-o-u-s system
    // to display meshes in the web
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    int     itmp;
    short   stmp;
    float   ftmp;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    // WRITE NUMBER OF VERTICES
    itmp=*np;
    swapint(&itmp);
    fwrite(&itmp,1,sizeof(int),f);
    // WRITE NUMBER OF TRIANGLES
    itmp=*nt;
    swapint(&itmp);
    fwrite(&itmp,1,sizeof(int),f);
    // WRITE VERTICES
    for(i=0;i<*np;i++)
    {
        ftmp=p[i].x;
        //swapfloat(&ftmp);
        fwrite(&ftmp,1,sizeof(float),f);
        ftmp=p[i].y;
        //swapfloat(&ftmp);
        fwrite(&ftmp,1,sizeof(float),f);
        ftmp=p[i].z;
        //swapfloat(&ftmp);
        fwrite(&ftmp,1,sizeof(float),f);
    }
    // WRITE TRIANGLES
    for(i=0;i<*nt;i++)
    {
        stmp=t[i].a;
        //swapshort(&stmp);
        fwrite(&stmp,1,sizeof(short),f);
        stmp=t[i].b;
        //swapshort(&stmp);
        fwrite(&stmp,1,sizeof(short),f);
        stmp=t[i].c;
        //swapshort(&stmp);
        fwrite(&stmp,1,sizeof(short),f);
    }
    fclose(f);   
    return 0;
}
int Off_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    int     i;
    FILE    *f;
    char    str[512];
    
    f=fopen(path,"r");

    fgets(str,512,f);    // skip OFF
    fscanf(f," %i %i %*i ",np,nt);
    *p=(float3D*)calloc(*np,sizeof(float3D));
    *t=(int3D*)calloc(*nt,sizeof(int3D));
    for(i=0;i<*np;i++)
        fscanf(f," %f %f %f ",&((*p)[i].x),&((*p)[i].y),&((*p)[i].z));
    if(verbose)
        printf("Read %i vertices\n",*np);
    for(i=0;i<*nt;i++)
        fscanf(f," %*i %i %i %i ",&((*t)[i].a),&((*t)[i].b),&((*t)[i].c));
    if(verbose)
        printf("Read %i triangles\n",*nt);
    fclose(f);
    return 0;
}
int Off_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;

    f=fopen(path,"w");
    if(f==NULL)
    {
        printf("ERROR: Cannot write to file %s\n",path);
        return 1;
    }
    
    fprintf(f,"OFF\n");
    fprintf(f,"%i %i 0\n",*np,*nt);
    for(i=0;i<*np;i++)
        fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    for(i=0;i<*nt;i++)
        fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);
    fclose(f);
    return 0;
}
int FloatData_save_data(char *path, Mesh *m)
{
    if(verbose)
        printf("* FloatData_save_data\n");

    int     *np=&(m->np);
    float   *data=m->data;
    FILE    *f;

    f=fopen(path,"w");
    
    if(f==NULL)
        return 1;
    
    fprintf(f,"%i %i 3\n",*np,m->ddim);
    fwrite(data,(m->np)*m->ddim,sizeof(float),f);
    fclose(f);
    
    return 0;
}
int RawFloatData_save_data(char *path, Mesh *m)
{
    if(verbose)
        printf("* RawFloatData_save_data\n");

    float   *data=m->data;
    FILE    *f;

    f=fopen(path,"w");
    
    if(f==NULL)
        return 1;
    
    fwrite(data,(m->np)*m->ddim,sizeof(float),f);
    fclose(f);
    
    return 0;
}
int VTK_load_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    int     ip,it,j,k,nval,x;
    char    str[512],str1[256],str2[256];
        
    f=fopen(path,"r");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}

    // READ HEADER
    *np=*nt=ip=it=0;
    do
    {
        fgets(str,511,f);
        sscanf(str," %s %i %*s ",str1,&x);

        if(strcmp(str1,"POINTS")==0)
        {
            *np=x;
            *p = (float3D*)calloc(*np,sizeof(float3D));
        }
        else
        if(*np>0 && ip<*np)
        {
            j=0;
            nval=0;            
            do
            {
                while(str[j]==' '||str[j]=='\t')
                    j++;
                k=0;
                while(str[j]!=' '&&str[j]!='\t'&&str[j]!='\r'&&str[j]!='\n')
                {
                    str2[k]=str[j];
                    j++;
                    k++;
                }
                str2[k]=(char)0;

                if(nval==0)
                    (*p)[ip].x=atof(str2);
                else
                if(nval==1)
                    (*p)[ip].y=atof(str2);
                else
                if(nval==2)
                    (*p)[ip].z=atof(str2);
                nval++;
                if(nval==3)
                {
                    nval=0;
                    ip++;
                }
            }
            while(str[j]!='\r'&&str[j]!='\n');
        }
        else
        if(ip==*np && strcmp(str1,"POLYGONS")==0)
        {
            *nt=x;
            *t = (int3D*)calloc(*nt,sizeof(int3D));
        }
        else
        if(*nt>0 && it<*nt)
        {
            sscanf(str," 3 %i %i %i ",&((*t)[it].a),&((*t)[it].b),&((*t)[it].c));
            it++;
        }
        else
        if(*nt>0 && it==*nt)
            break;
    }
    while(!feof(f));
    fclose(f);

    if(verbose)
        printf("Read %i vertices and %i triangles\n",ip,it);

    return 0;
}
int VTK_save_mesh(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    FILE    *f;
    int     i;
    
    f=fopen(path,"w");
    if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
    
    // WRITE HEADER
    fprintf(f,"# vtk DataFile Version 3.0\n");
    fprintf(f,"vtk output\n");
    fprintf(f,"ASCII\n");
    fprintf(f,"DATASET POLYDATA\n");

    // WRITE VERTICES
    fprintf(f,"POINTS %i float\n",*np);
    for(i=0;i<*np;i++)
    {
        fprintf(f,"%f %f %f ", p[i].x,p[i].y,p[i].z);
        if(i%3==0)
            fprintf(f,"\n");
    }
    fprintf(f,"\n");
    
    // WRITE TRIANGLES
    fprintf(f,"POLYGONS %i %i\n",*nt,*nt*4);
    for(i=0;i<*nt;i++)
        fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);

    fclose(f);
    
    return 0;
}
int DPV_load_data(char *path, Mesh *m)
{
    FILE    *f;
    char    str[255];
    int     *np=&(m->np);
    float   **data=&(m->data);
    int     i,n;

    if(verbose)
        printf("* DPV_load_data\n");

    f=fopen(path,"r");
    if(f==NULL)
        return 1;

    n=0;
    while(!feof(f))
    {
        fgets(str,255,f);
        n++;
    }
    fclose(f);
    
    f=fopen(path,"r");
    *np=n;
    *data=(float*)calloc(*np,sizeof(float));
    for(i=0;i<*np;i++)
    {
        fgets(str,255,f);
        sscanf(str," %*i %*f %*f %*f %f ",&((*data)[i]));
    }
    fclose(f);
    
    return 0;
}
void read_giiElement(char *path, char *el, char **data, int *n, int *d)
{
    FILE        *f;
    long        sz,expsz;
    char        *gzdata;
    z_stream    strm;
    char        str[256],cmd[2048];
    int         nd;

    *d = 0;
    nd = 0;

    // Read number of values
    sprintf(cmd,"awk 'BEGIN{s=0}/%s/{s=1}/Dimensionality/{if(s==1){split($0,a,\"\\\"\");d=a[2];s=2}}/Dim0/{if(s==2){split($0,a,\"\\\"\");n=a[2];s=0}}END{print d,n}' %s",el,path);
    f=popen(cmd,"r");
    fgets(str,256,f);
    sscanf(str," %i %i ", d, n);
    pclose(f);

    // We currently handle only two cases. In the first one, Dimensionality=1 and Dim0
    // gives the number of values to read. In the second case, Dimensionality=2, Dim0 is
    // a number of vectors, each with 3 values. More generally, this value could be
    // obtained from the Dim1 field
    if(*d==1)
    {
        nd = 1;
    }
    else if(*d==2)
    {
        nd = 3;
    }
    printf("dim: %i\n",nd);
    printf("n: %i\n",*n);

    // Calculate length of gzip compressed data
    sprintf(cmd,"awk 'BEGIN{s=0}/%s/{s=1}/<Data>/{if(s==1){split($0,a,\"[<>]\");d=a[3];s=0}}END{print d}' %s|base64 --decode",el,path);
    f=popen(cmd,"r");
    sz=0;
    while(!feof(f))
    {
        fread(str,256,sizeof(char),f);
        sz+=256;
    }
    pclose(f);
    
    // Allocate memory for gzip data
    gzdata=(char*)calloc(sz,sizeof(char));
    
    // Read gzip data
    f=popen(cmd,"r");
    sz=0;
    while(!feof(f))
    {
        fread(gzdata+sz,256,sizeof(char),f);
        sz+=256;
    }
    pclose(f);
    
    // Inflate gzip data
    expsz=(*n)*nd*4; // because 4: float or int 4 byte values
    *data=calloc(expsz,sizeof(char));
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = (uInt)sz; // size of input
    strm.next_in = (Bytef *)gzdata; // input char array
    strm.avail_out = (uInt)expsz; // size of output
    strm.next_out = (Bytef *)(*data); // output char array
    inflateInit(&strm);
    inflate(&strm, Z_NO_FLUSH);
    inflateEnd(&strm);

    // Free memory
    free(gzdata);
}
int Gii_load(char *path, Mesh *m)
{
    int     d = 0;
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);

    read_giiElement(path,"NIFTI_INTENT_POINTSET",(char**)p,np,&d);
    read_giiElement(path,"NIFTI_INTENT_TRIANGLE",(char**)t,nt,&d);
    if(verbose)
    {
        printf("Read %i vertices\n",*np);
        printf("Read %i triangles\n",*nt);
    }
    
    return 0;    
}
int Gii_load_data(char *path, Mesh *m)
{
    printf("Gii_load_data\n");

    int     d = 0;
    int     *np=&(m->np);
    float   **data=&(m->data);

    read_giiElement(path,"NIFTI_INTENT_NONE",(char**)data,np,&d);
    if(verbose)
    {
        printf("Read %i values\n",*np);
    }

    return 0;
}


#pragma mark -
void printVertex(float3D p)
{
    printf("x: %g, y:%g, z:%g\n",p.x,p.y,p.z);
}
void printTriangle(int3D t)
{
    printf("a: %i, b:%i, c:%i\n",t.a,t.b,t.c);
}
void printTriangleAndVertices(Mesh *m, int i)
{
    float3D *p=m->p;
    int3D   *t=m->t;
    printf("t %i, a: %i (%g, %g, %g), b:%i (%g, %g, %g), c:%i (%g, %g, %g)\n",
        i,
        t[i].a, p[t[i].a].x,p[t[i].a].y,p[t[i].a].z,
        t[i].b, p[t[i].b].x,p[t[i].b].y,p[t[i].b].z,
        t[i].c, p[t[i].c].x,p[t[i].c].y,p[t[i].c].z);
}
int freeMesh(Mesh *m)
{
    if(verbose) printf("* freeMesh\n");

    free(m->p);
    free(m->t);
    if(m->data)
        free(m->data);
    if(m->selection)
        free(m->selection);
    if(m->NT)
        free(m->NT);

    return 0;
}
int loadMesh(char *path, Mesh *m,int iformat)
{
    if(verbose) printf("* imesh: %s\n",path);

    int        err,format;

    if(iformat==0)
        format=getformatindex(path);
    else
        format=iformat;

    switch(format)
    {
        case kFreeSurferMesh:
            err=FreeSurfer_load_mesh(path,m);
            break;
        case kFreeSurferData:
            err=FreeSurfer_load_data(path,m);
            break;
        case kFreeSurferAnnot:
            err=FreeSurfer_load_annot(path,m);
            break;
        case kMGHData:
            err=FreeSurfer_load_mghdata(path,m);
            break;
        case kBrainVisaMesh:
            err=BrainVisa_load_mesh(path,m);
            break;
        case kText:
            err=Text_load(path,m);
            break;
        case kVRMLMesh:
            err=VRML_load_mesh(path,m);
            break;
        case kObjMesh:
            err=Obj_load(path,m);
            break;
        case kPlyMesh:
            err=Ply_load(path,m);
            break;
        case kSTLMesh:
            err=STL_load(path,m);
            break;
        case kSmeshMesh:
            err=Smesh_load(path,m);
            break;
        case kBinMesh:
            err=Bin_load(path,m);
            break;
        case kOffMesh:
            err=Off_load(path,m);
            break;
        case kVTKMesh:
            err=VTK_load_mesh(path,m);
            break;
        case kDPVData:
            err=DPV_load_data(path,m);
            break;
        case kCivetObjMesh:
            err=CivetObj_load_mesh(path,m);
            break;
        case kGiiMesh:
            err=Gii_load(path,m);
            break;
        case kGiiData:
            err=Gii_load_data(path,m);
            break;
        case kAscMesh:
            err=Asc_load(path,m);
            break;
        default:
            printf("ERROR: Input mesh format not recognised\n");
            return 1;
    }

    if(err==0)
    {
        // Mesh read without error
        // initialise selection
        m->selection = (char*)calloc(m->np,sizeof(char));
    }
    else
    {
        printf("ERROR: cannot read file: %s\n",path);
        char cwd[1024];
        if(getcwd(cwd,sizeof(cwd))!=NULL)
            printf("Current working directory: %s\n", cwd);

        return 1;
    }

    return 0;
}
int addMesh(char *path, Mesh *m0,int iformat)
{
    if(verbose) printf("* addMesh: %s\n",path);

    // m0: old mesh
    // m1: new mesh
    Mesh    amesh;
    Mesh    *m1;
    int     i;

    amesh.p=NULL;
    amesh.t=NULL;
    amesh.data=NULL;
    amesh.NT=NULL;

    loadMesh(path,&amesh,iformat);

    m1=(Mesh*)calloc(1,sizeof(Mesh));
    m1->np=m0->np+amesh.np;
    m1->nt=m0->nt+amesh.nt;
    m1->p=(float3D*)calloc(m0->np+amesh.np,sizeof(float3D));
    m1->t=(int3D*)calloc(m0->nt+amesh.nt,sizeof(int3D));
    m1->selection = (char*)calloc(m0->nt+amesh.nt,sizeof(char));

    // add points
    for(i=0;i<m0->np;i++)
        m1->p[i]=m0->p[i];
    for(i=0;i<amesh.np;i++)
        m1->p[m0->np+i]=amesh.p[i];
        
    // add triangles
    for(i=0;i<m0->nt;i++)
        m1->t[i]=m0->t[i];
    for(i=0;i<amesh.nt;i++)
        m1->t[m0->nt+i]=(int3D){amesh.t[i].a+m0->np,amesh.t[i].b+m0->np,amesh.t[i].c+m0->np};

    // free tmp mesh
    freeMesh(&amesh);

    // configure new mesh
    freeMesh(m0);
    m0->np=m1->np;
    m0->nt=m1->nt;
    m0->p=m1->p;
    m0->t=m1->t;
    m0->selection=m1->selection;

    return 0;
}
int saveMesh(char *path, Mesh *m, int oformat)
{
    if(verbose) printf("* omesh: %s\n",path);

    int    err=0,format;
    
    if(oformat==0)
        format=getformatindex(path);
    else
        format=oformat;

    switch(format)
    {
        case kFreeSurferMesh:
            err=FreeSurfer_save_mesh(path,m);
            break;
        case kFreeSurferData:
            err=FreeSurfer_save_data(path,m);
            break;
        case kMGHData:
            err=FreeSurfer_save_mghdata(path,m);
            break;
        case kBrainVisaMesh:
            err=BrainVisa_save_mesh(path,m);
            break;
        case kText:
            err=Text_save_mesh(path,m);
            break;
        case kTextData:
            err=Text_save_data(path,m);
            break;
        case kVRMLMesh:
            err=VRML_save_mesh(path,m);
            break;
        case kObjMesh:
            err=Obj_save_mesh(path,m);
            break;
        case kPlyMesh:
            err=Ply_save_mesh(path,m);
            break;
        case kSTLMesh:
            err=STL_save_mesh(path,m);
            break;
        case kSmeshMesh:
            err=Smesh_save_mesh(path,m);
            break;
        case kBinMesh:
            err=Bin_save_mesh(path,m);
            break;
        case kOffMesh:
            err=Off_save_mesh(path,m);
            break;
        case kFloatData:
            err=FloatData_save_data(path,m);
            break;
        case kRawFloatData:
            err=RawFloatData_save_data(path,m);
            break;
        case kVTKMesh:
            err=VTK_save_mesh(path,m);
            break;
        case kCivetObjMesh:
            err=CivetObj_save_mesh(path,m);
            break;
        case kAscMesh:
            err=Asc_save_mesh(path,m);
            break;
        default:
            printf("ERROR: Output data format not recognised\n");
            err=1;
            break;
    }
    if(err!=0)
    {
        printf("ERROR: cannot write to file: %s\n",path);
        return 1;
    }
    
    return 0;
}

#pragma mark -
#pragma mark [ Save TIFF ]
void WriteHexString(FILE *f, char *str)
{
    int		i,j,len=strlen(str);
    int		a;
    short	b;
    char	c[5];
    
    for(i=0;i<len;i+=4)
    {
        for(j=0;j<4;j++)
            c[j]=str[i+j];
        c[4]=(char)0;
        sscanf(c,"%x",&a);
        b=(short)a;
        fwrite(&((char*)&b)[1],1,1,f);
        fwrite(&((char*)&b)[0],1,1,f);
    }
}
void writeTIFF(char *path, char *addr, int nx, int ny)
{
    FILE	*fptr;
    int		offset;
    int		i,j;
    char	red,green,blue;
    
    fptr=fopen(path,"w");
    
    /* Write the header */
    WriteHexString(fptr,"4d4d002a");    /* Little endian & TIFF identifier */
    offset = nx * ny * 3 + 8;
    putc((offset & 0xff000000) / 16777216,fptr);
    putc((offset & 0x00ff0000) / 65536,fptr);
    putc((offset & 0x0000ff00) / 256,fptr);
    putc((offset & 0x000000ff),fptr);

    /* Write the binary data */
    for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
    
         red=addr[4*(j*nx+i)+0];
         green=addr[4*(j*nx+i)+1];
         blue=addr[4*(j*nx+i)+2];
         
         fputc(red,fptr);
         fputc(green,fptr);
         fputc(blue,fptr);
      }
    }
   
    WriteHexString(fptr,"000e");  						/* Write the footer */ /* The number of directory entries (14) */
    WriteHexString(fptr,"0100000300000001");			/* Width tag, short int */
    fputc((nx & 0xff00) / 256,fptr);    /* Image width */
    fputc((nx & 0x00ff),fptr);
    WriteHexString(fptr,"0000");
    WriteHexString(fptr,"0101000300000001");			/* Height tag, short int */
    fputc((ny & 0xff00) / 256,fptr);    /* Image height */
    fputc((ny & 0x00ff),fptr);
    WriteHexString(fptr,"0000");
    WriteHexString(fptr,"0102000300000003");			/* Bits per sample tag, short int */
    offset = nx * ny * 3 + 182;
    putc((offset & 0xff000000) / 16777216,fptr);
    putc((offset & 0x00ff0000) / 65536,fptr);
    putc((offset & 0x0000ff00) / 256,fptr);
    putc((offset & 0x000000ff),fptr);
    WriteHexString(fptr,"010300030000000100010000");	/* Compression flag, short int */
    WriteHexString(fptr,"010600030000000100020000");	/* Photometric interpolation tag, short int */
    WriteHexString(fptr,"011100040000000100000008");	/* Strip offset tag, long int */
    WriteHexString(fptr,"011200030000000100010000");	/* Orientation flag, short int */
    WriteHexString(fptr,"011500030000000100030000");	/* Sample per pixel tag, short int */
    WriteHexString(fptr,"0116000300000001");			/* Rows per strip tag, short int */
    fputc((ny & 0xff00) / 256,fptr);
    fputc((ny & 0x00ff),fptr);
    WriteHexString(fptr,"0000");
    WriteHexString(fptr,"0117000400000001");			/* Strip byte count flag, long int */
    offset = nx * ny * 3;
    putc((offset & 0xff000000) / 16777216,fptr);
    putc((offset & 0x00ff0000) / 65536,fptr);
    putc((offset & 0x0000ff00) / 256,fptr);
    putc((offset & 0x000000ff),fptr);
    WriteHexString(fptr,"0118000300000003");			/* Minimum sample value flag, short int */
    offset = nx * ny * 3 + 188;
    putc((offset & 0xff000000) / 16777216,fptr);
    putc((offset & 0x00ff0000) / 65536,fptr);
    putc((offset & 0x0000ff00) / 256,fptr);
    putc((offset & 0x000000ff),fptr);
    WriteHexString(fptr,"0119000300000003");			/* Maximum sample value tag, short int */
    offset = nx * ny * 3 + 194;
    putc((offset & 0xff000000) / 16777216,fptr);
    putc((offset & 0x00ff0000) / 65536,fptr);
    putc((offset & 0x0000ff00) / 256,fptr);
    putc((offset & 0x000000ff),fptr);
    WriteHexString(fptr,"011c00030000000100010000");	/* Planar configuration tag, short int */
    WriteHexString(fptr,"0153000300000003");			/* Sample format tag, short int */
    offset = nx * ny * 3 + 200;
    putc((offset & 0xff000000) / 16777216,fptr);
    putc((offset & 0x00ff0000) / 65536,fptr);
    putc((offset & 0x0000ff00) / 256,fptr);
    putc((offset & 0x000000ff),fptr);
    WriteHexString(fptr,"00000000");					/* End of the directory entry */
    WriteHexString(fptr,"000800080008");				/* Bits for each colour channel */
    WriteHexString(fptr,"000000000000");				/* Minimum value for each component */
    WriteHexString(fptr,"00ff00ff00ff");				/* Maximum value per channel */
    WriteHexString(fptr,"000100010001");				/* Samples per pixel for each channel */
    fclose(fptr);
}
// plot long rainbow RGB from https://www.particleincell.com/2014/colormap/
int greyscale(float val, float *r, float *g, float *b)
{
    *r=val;
    *g=val;
    *b=val;
    return 1;
}
int rainbow(float val, float *r, float *g, float *b)
{
    float   a=(1-val)/0.2;	//invert and group
    int     X=floor(a);	//this is the integer part
    float   Y=(a-X); //fractional part from 0 to 255
    switch(X)
    {
        case 0: *r=1;*g=Y;*b=0;break;
        case 1: *r=1-Y;*g=1;*b=0;break;
        case 2: *r=0;*g=1;*b=Y;break;
        case 3: *r=0;*g=1-Y;*b=1;break;
        case 4: *r=Y;*g=0;*b=1;break;
        case 5: *r=1;*g=0;*b=1;break;
    }
    return 1;
}
#pragma mark -
#pragma mark [ Geometry functions ]
void absgi(Mesh *m)
{
    float   S,V,logAbsGI;

    S=area(m);
    V=volume(m);
    
    // log(absGI)    = log(Sx)-2log(Vx)/3-log(36π)/3
    // absGI        = Sx/(Vx^(2/3)(36π)^(1/3))
    logAbsGI=log(S)-2*log(V)/3.0-log(36*M_PI)/3.0;
    
    printf("absgi: %f\n",exp(logAbsGI));
}
void align(Mesh *m, char *path)
{
/*
    Align (translate/rotate) the mesh m to the mesh pointed by path. Both have the same topology,
    different geometry.
*/
    Mesh    target;
    int     np=m->np;
    float3D *p=m->p;
    int     npt;    // number of vertices in target
    float3D *pt;    // vertices in target
    float3D pp; // moving vertex
    int     i,j,niter;
    float3D c0={0,0,0},c1={0,0,0};
    float   x0,y0,z0,x,y,z,M[9];
    float   err,minerr;
    int     iminerr;
    
    loadMesh(path, &target,0);
    npt=target.np;
    pt=target.p;
    if(np!=npt)
    {
        printf("ERROR: original and target meshes have to have the same number of points\n");
        return;
    }

    // Move barycentre of m to zero
    //-------------------------------
    for(i=0;i<np;i++)
    {
        c0=add3D(c0,p[i]);
        c1=add3D(c1,pt[i]);
    }
    c0=sca3D(c0,1/(float)np);
    c1=sca3D(c1,1/(float)np);
    for(i=0;i<np;i++)
    {
        p[i]=sub3D(p[i],c0);
        pt[i]=sub3D(pt[i],c1);
    }
    
    // Rotate m to minimise (m-p)^2
    //-------------------------------
    // Guess an initial rotation
    niter=1000;
    iminerr=0;
    srand(0); // seed the random number generator
    for(j=0;j<niter;j++)
    {
        x=2*M_PI*rand()/(float)RAND_MAX;
        y=2*M_PI*rand()/(float)RAND_MAX;
        z=2*M_PI*rand()/(float)RAND_MAX;

        M[0]=cos(z)*cos(y);
        M[1]=-sin(z)*cos(x)+cos(z)*sin(y)*sin(x);
        M[2]=sin(z)*sin(x)+cos(z)*sin(y)*cos(x);
        M[3]=sin(z)*cos(y);
        M[4]=cos(z)*cos(x)+sin(z)*sin(y)*sin(x);
        M[5]=-cos(z)*sin(x)+sin(z)*sin(y)*cos(x);
        M[6]=-sin(y);
        M[7]=cos(y)*sin(x);
        M[8]=cos(y)*cos(x);

        err=0;
        for(i=0;i<np;i++)
        {
            pp.x=M[0]*p[i].x+M[1]*p[i].y+M[2]*p[i].z;
            pp.y=M[3]*p[i].x+M[4]*p[i].y+M[5]*p[i].z;
            pp.z=M[6]*p[i].x+M[7]*p[i].y+M[8]*p[i].z;
            err+=norm3D(sub3D(pt[i],pp));
        }
        if(j==0)
            minerr=err;
        else if(err<minerr)
        {
            minerr=err;
            iminerr=j;
        }
    }
    printf("Minimum err2 for 1st guess: %f\n",minerr);
    srand(0);
    for(j=0;j<iminerr+1;j++)
    {
        x0=2*M_PI*rand()/(float)RAND_MAX;
        y0=2*M_PI*rand()/(float)RAND_MAX;
        z0=2*M_PI*rand()/(float)RAND_MAX;
    }

    iminerr=0;
    srand(0); // seed the random number generator
    for(j=0;j<niter;j++)
    {
        x=x0+(10/M_PI)-(20/M_PI)*rand()/(float)RAND_MAX;
        y=y0+(10/M_PI)-(20/M_PI)*rand()/(float)RAND_MAX;
        z=z0+(10/M_PI)-(20/M_PI)*rand()/(float)RAND_MAX;

        M[0]=cos(z)*cos(y);
        M[1]=-sin(z)*cos(x)+cos(z)*sin(y)*sin(x);
        M[2]=sin(z)*sin(x)+cos(z)*sin(y)*cos(x);
        M[3]=sin(z)*cos(y);
        M[4]=cos(z)*cos(x)+sin(z)*sin(y)*sin(x);
        M[5]=-cos(z)*sin(x)+sin(z)*sin(y)*cos(x);
        M[6]=-sin(y);
        M[7]=cos(y)*sin(x);
        M[8]=cos(y)*cos(x);

        err=0;
        for(i=0;i<np;i++)
        {
            pp.x=M[0]*p[i].x+M[1]*p[i].y+M[2]*p[i].z;
            pp.y=M[3]*p[i].x+M[4]*p[i].y+M[5]*p[i].z;
            pp.z=M[6]*p[i].x+M[7]*p[i].y+M[8]*p[i].z;
            err+=norm3D(sub3D(pt[i],pp));
        }
        if(j==0)
            minerr=err;
        else if(err<minerr)
        {
            minerr=err;
            iminerr=j;
        }
    }
    printf("Minimum err2 for 2nd guess: %f\n",minerr);
    srand(0);
    for(j=0;j<iminerr+1;j++)
    {
        x=x0+(10/M_PI)-(20/M_PI)*rand()/(float)RAND_MAX;
        y=y0+(10/M_PI)-(20/M_PI)*rand()/(float)RAND_MAX;
        z=z0+(10/M_PI)-(20/M_PI)*rand()/(float)RAND_MAX;
    }

    printf("x: %g, y: %g, z: %g\n",x*180/M_PI,y*180/M_PI,z*180/M_PI);

    M[0]=cos(z)*cos(y);
    M[1]=-sin(z)*cos(x)+cos(z)*sin(y)*sin(x);
    M[2]=sin(z)*sin(x)+cos(z)*sin(y)*cos(x);
    M[3]=sin(z)*cos(y);
    M[4]=cos(z)*cos(x)+sin(z)*sin(y)*sin(x);
    M[5]=-cos(z)*sin(x)+sin(z)*sin(y)*cos(x);
    M[6]=-sin(y);
    M[7]=cos(y)*sin(x);
    M[8]=cos(y)*cos(x);
    for(i=0;i<np;i++)
    {
        pp.x=M[0]*p[i].x+M[1]*p[i].y+M[2]*p[i].z;
        pp.y=M[3]*p[i].x+M[4]*p[i].y+M[5]*p[i].z;
        pp.z=M[6]*p[i].x+M[7]*p[i].y+M[8]*p[i].z;
        p[i]=pp;
    }

    // move mesh m to target barycentre
    /*
    for(i=0;i<np;i++)
        p[i]=add3D(p[i],c1);
    */
}
float area(Mesh *m)
{
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    int     i;
    float   area=0;
    
    for(i=0;i<*nt;i++)
        area+=triangle_area(p[t[i].a],p[t[i].b],p[t[i].c]);
    printf("area: %f\n",area);
    
    return area;
}
int areaMap(float *C, Mesh *m)
{
/*
    Computes a vector with the area around each vertex.
    The area of each triangle is distributed among its
    three vertices.
*/
    if(verbose) printf("* areaMap\n");
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    int     i;
    float   a;

    for(i=0;i<*nt;i++)
    {
        a=triangle_area(p[t[i].a],p[t[i].b],p[t[i].c]);
        C[t[i].a]+=a/3;
        C[t[i].b]+=a/3;
        C[t[i].c]+=a/3;
    }    
    return 0;
}
int average(int N, char *paths[], Mesh *m)
{
    if(verbose) printf("* average\n");
    
    int     i,j;
    int     np;
    float3D *p;
    Mesh    m1;
    
    loadMesh(paths[0],m,0);
    np=m->np;
    p=m->p;
    
    for(j=1;j<N;j++)
    {
        loadMesh(paths[j],&m1,0);
        for(i=0;i<np;i++)
            p[i]=add3D(p[i],m1.p[i]);
        free(m1.p);
        free(m1.t);
    }
    for(i=0;i<np;i++)
        p[i]=sca3D(p[i],1/(float)N);

    return 0;
}
int applyMatrix(float *M, Mesh *m)
{
    /*
    Multiply all mesh vertices by the 4x4 matrix m
    */
    int i;
    for(i=0;i<m->np;i++)
    {
        multMatVec(M,m->p[i],&(m->p[i]));
    }
    
    return 0;
}
int barycentricProjection(char *path_rm, Mesh *m)
{
    /*
    If the mesh is a smooth mesh (for example, spherical), and given
    a smooth reference mesh (rm) print for each vertex in rm its
    barycentric coordinates within the triangles of the mesh. This
    permits, for example, to map data defined over the vertices of
    the mesh over rm.
    */
    Mesh    rm;
    int     nt;
    int3D   *t;
    int     np_rm;
    float3D *p,*p_rm;
    float   c0,c1;
    int     i,j,result;
    float3D n;
    float   flipTest;
    
    loadMesh(path_rm,&rm,0);
    
    // Check whether the meshes are properly oriented
    n=normal3D(0,m);
    flipTest=dot3D(m->p[m->t[0].a],n);
    if(flipTest<0)
    {
        printf("ERROR: the mesh is mis-oriented\n");
        return 1;
    }
    n=normal3D(0,&rm);
    flipTest=dot3D(rm.p[rm.t[0].a],n);
    if(flipTest<0)
    {
        printf("ERROR: rm is mis-oriented\n");
        return 1;
    }
    
    // Actual mesh (smooth)
    p=m->p;
    nt=m->nt;
    t=m->t;

    // Reference mesh (smooth)
    np_rm=rm.np;
    p_rm=rm.p;
    
    for(i=0;i<np_rm;i++)
    {
        for(j=0;j<nt;j++)
        {
            result=intersect_VectorTriangle(p_rm[i],j,&c0,&c1,m);
            if(result==1)
            {
                printf("%i %f %f\n",j,c0,c1);
                break;
            }
        }
        if(j==nt)
        {
            int k,imin=0;
            float   d,dmin=1000;
            printf("ERROR: could not resample point %i of the reference mesh (%f %f %f)\n",i,p_rm[i].x,p_rm[i].y,p_rm[i].z);
            for(k=0;k<m->np;k++)
            {
                d=norm3D(sub3D(m->p[k],p_rm[i]));
                if(d<dmin)
                {
                    dmin=d;
                    imin=k;
                }
            }
            printf("closest vertex %i (%f,%f,%f), dist: %f\n",imin,p[imin].x,p[imin].y,p[imin].z,dmin);

            result=intersect_VectorTriangle(p_rm[i],1984,&c0,&c1,m);
            printf(">> t[1984] c0,c1: %f, %f\n",c0,c1);
            result=intersect_VectorTriangle(p_rm[i],1920,&c0,&c1,m);
            printf(">> t[1920] c0,c1: %f, %f\n",c0,c1);

            printf("triangles including vertex %i\n",imin);
            for(k=0;k<nt;k++)
            {
                if(t[k].a==imin||t[k].b==imin||t[k].c==imin)
                {
                    printf("%i. %i,%i,%i\n",k,t[k].a,t[k].b,t[k].c);
                    printf("   %i: %f,%f,%f\n",t[k].a,p[t[k].a].x,p[t[k].a].y,p[t[k].a].z);
                    printf("   %i: %f,%f,%f\n",t[k].b,p[t[k].b].x,p[t[k].b].y,p[t[k].b].z);
                    printf("   %i: %f,%f,%f\n",t[k].c,p[t[k].c].x,p[t[k].c].y,p[t[k].c].z);
                }
            }
            return 1;
        }
    }
    
    return 0;
}
void barycentre(Mesh *m)
{
    int     *np=&(m->np);
    float3D *p=m->p;
    int     i;
    float3D centre={0,0,0};
    
    for(i=0;i<*np;i++)
        centre=add3D(centre,p[i]);
    centre=sca3D(centre,1/(float)*np);
    if(verbose)
    {
        printf("centre %g,%g,%g\n",centre.x,centre.y,centre.z);
    }
    for(i=0;i<*np;i++)
        p[i]=sub3D(p[i],centre);
}
int boundingBox(Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p;
    float3D min,max;
    int     i;
    
    min=max=p[0];
    for(i=0;i<np;i++)
    {
        min.x=(p[i].x<min.x)?p[i].x:min.x;
        min.y=(p[i].y<min.y)?p[i].y:min.y;
        min.z=(p[i].z<min.z)?p[i].z:min.z;
        max.x=(p[i].x>max.x)?p[i].x:max.x;
        max.y=(p[i].y>max.y)?p[i].y:max.y;
        max.z=(p[i].z>max.z)?p[i].z:max.z;
    }
    printf("boundingBox: %f,%f,%f,%f,%f,%f\n",min.x,min.y,min.z,max.x,max.y,max.z);

    return 0;
}
void translate(float a,float b,float c,Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p;
    int     i;
    float3D new_centre;
    new_centre=(float3D){a,b,c};
    for(i=0;i<np;i++)
        p[i]=add3D(p[i],new_centre);
}
void centre(Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p;
    int     i;
    float3D mi,ma,centre;
    
    mi=ma=p[0];
    for(i=0;i<np;i++)
    {
        mi.x=(mi.x>p[i].x)?p[i].x:mi.x;
        mi.y=(mi.y>p[i].y)?p[i].y:mi.y;
        mi.z=(mi.z>p[i].z)?p[i].z:mi.z;
        ma.x=(ma.x<p[i].x)?p[i].x:ma.x;
        ma.y=(ma.y<p[i].y)?p[i].y:ma.y;
        ma.z=(ma.z<p[i].z)?p[i].z:ma.z;
    }
    centre=(float3D){(mi.x+ma.x)/2.0,(mi.y+ma.y)/2.0,(mi.z+ma.z)/2.0};
    if(verbose)
        printf("centre %g,%g,%g\n",centre.x,centre.y,centre.z);
    for(i=0;i<np;i++)
        p[i]=sub3D(p[i],centre);
}
void checkOrientation(Mesh *m)
{
    float3D n=normal3D(0,m);
    float   flipTest=dot3D(m->p[m->t[0].a],n);
    
    printf("orientation: %c\n",(flipTest>0)?'+':'-');

}
int clip(Mesh *m, float min, float max)
{
    float   *data=m->data;
    int     i;
    int     np=m->np;
    
    if(data==NULL)
    {
        printf("ERROR: In clip, no data available\n");
        return 0;
    }
    
    if(verbose) 
        printf("clipping to [%f,%f]\n",min,max);
        
    for(i=0;i<np;i++)
    {
        if(data[i]>max)
            data[i]=max;
        if(data[i]<min)
            data[i]=min;
    }
    
    return 1;
}
double  sum;
int     *tmark,icmax,ncverts;
void cluster(int ip, float *thrsrc, float thr, Mesh *m)
{
    int3D   *t=m->t;
    NTriRec *NT=m->NT;
    int     i,j,it;
    int     *tt;
    
    tmark[ip]=1;
    ncverts++;    
    for(i=0;i<=NT[ip].n;i++)
    {
        it=NT[ip].t[i];        
        tt=(int*)&(t[it]);
        for(j=0;j<3;j++)
            if(thrsrc[tt[j]]>=thr && tmark[tt[j]]==0)
            {
                cluster(tt[j],thrsrc,thr,m);
                if(thrsrc[tt[j]]>thrsrc[icmax])
                    icmax=tt[j];
            }
    }
}
float cot(float3D a, float3D b)
{
    // cot(a,b)=cos(a,b)/sin(a,b)=a*b/sqrt(|a|^2*|b|^2-(a*b)^2)
    float   ab=dot3D(a,b);
    float   na2=norm3Dsqr(a);
    float   nb2=norm3Dsqr(b);
    
    return ab/sqrt(na2*nb2-ab*ab);
}
void countClusters(float thr, Mesh *m)
{
    int     *np=&(m->np);
    float3D *p=m->p;
    float   *data=m->data;
    int     n;
    int     i;
    
    neighbours(m);
    
    n=1;
    tmark=(int*)calloc(*np,sizeof(int));
    if(verbose)
        printf("number\tnVertices\timax\tmax\tXmax\tYmax\tZmax\n");
    for(i=0;i<*np;i++)
    if(data[i]>=thr && tmark[i]==0)
    {
        icmax=i;
        ncverts=0;
        cluster(i,data,thr,m);
        if(verbose)
            printf("%i\t%i\t%i\t%f\t%f\t%f\t%f\n",n,ncverts,icmax,data[icmax],p[icmax].x,p[icmax].y,p[icmax].z);
        n++;
    }
    if(verbose)
        printf("\n");
    else
        printf("countClusters: %i\n",n-1);
        
    free(tmark);
}
int curvature(float *C, Mesh *m)
{
    if(verbose>1)
        printf("mean curvature\n");
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp,*tmp1;
    int     *n;
    float3D nn;
    float   absmax;
    int     i;

    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    // compute smoothing direction as the vector to the average of neighbour vertices
    for(i=0;i<*nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
        tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    for(i=0;i<*np;i++)
        tmp[i]=sub3D(sca3D(tmp[i],1/(float)n[i]),p[i]);
        
    tmp1=(float3D*)calloc(*np,sizeof(float3D));
    // compute normal direction as the average of neighbour triangle normals
    for(i=0;i<*nt;i++)
    {
        nn=cross3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
        nn=sca3D(nn,1/norm3D(nn));
        tmp1[t[i].a]=add3D(tmp1[t[i].a],nn);
        tmp1[t[i].b]=add3D(tmp1[t[i].b],nn);
        tmp1[t[i].c]=add3D(tmp1[t[i].c],nn);
    }
    for(i=0;i<*np;i++)
        tmp1[i]=sca3D(tmp1[i],1/(float)n[i]);
    free(n);
    
    for(i=0;i<*np;i++)
        C[i]=-dot3D(tmp1[i],tmp[i]);
    free(tmp);
    free(tmp1);
    
    absmax=-1;
    for(i=0;i<*np;i++)
        absmax=(fabs(C[i])>absmax)?fabs(C[i]):absmax;
    absmax*=0.95;
    for(i=0;i<*np;i++)
    {
        C[i]/=absmax;
        if(C[i]>1)    C[i]=1;
        if(C[i]<-1)   C[i]=-1;
    }
    
    return 0;
}
int curvature_exact(float *C, Mesh *m)
{
    // Adapted from Sullivan (2008) Curvatures of smooth and discrete surfaces
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp,*tmp1;
    int     *n;
    float3D a,b,c,nn;
    float   cotab,cotbc,cotca;
    int     i;

    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    for(i=0;i<*nt;i++)
    {
        a=p[t[i].a];
        b=p[t[i].b];
        c=p[t[i].c];
        cotab=cot(sub3D(a,c),sub3D(b,c));
        cotbc=cot(sub3D(b,a),sub3D(c,a));
        cotca=cot(sub3D(c,b),sub3D(a,b));
        tmp[t[i].a]=add3D(tmp[t[i].a], sca3D(sub3D(a,b),cotab));
        tmp[t[i].a]=add3D(tmp[t[i].a], sca3D(sub3D(a,c),cotca));
        tmp[t[i].b]=add3D(tmp[t[i].b], sca3D(sub3D(b,c),cotbc));
        tmp[t[i].b]=add3D(tmp[t[i].b], sca3D(sub3D(b,a),cotab));
        tmp[t[i].c]=add3D(tmp[t[i].c], sca3D(sub3D(c,a),cotca));
        tmp[t[i].c]=add3D(tmp[t[i].c], sca3D(sub3D(c,b),cotbc));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    for(i=0;i<*np;i++)
        tmp[i]=sca3D(tmp[i],1/(float)n[i]);
        
    tmp1=(float3D*)calloc(*np,sizeof(float3D));
    // compute normal direction as the average of neighbour triangle normals
    for(i=0;i<*nt;i++)
    {
        nn=cross3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
        nn=sca3D(nn,1/norm3D(nn));
        tmp1[t[i].a]=add3D(tmp1[t[i].a],nn);
        tmp1[t[i].b]=add3D(tmp1[t[i].b],nn);
        tmp1[t[i].c]=add3D(tmp1[t[i].c],nn);
    }
    for(i=0;i<*np;i++)
        tmp1[i]=sca3D(tmp1[i],1/(float)n[i]);
    free(n);
    
    for(i=0;i<*np;i++)
        C[i]=dot3D(tmp1[i],tmp[i]);
    free(tmp);
    free(tmp1);

    /*
    float   min,max;
    min=max=0;
    for(i=0;i<*np;i++)
    {
        min=(C[i]<min)?C[i]:min;
        max=(C[i]>max)?C[i]:max;
    }
    printf("min,max: %f,%f\n",min,max);
    */
    
    return 0;
}
void depth(float *C, Mesh *m)
{
    if(verbose)
        printf("depth\n");
    int			i;
    float		n,max;
    float3D		ce={0,0,0},ide,siz;
    int         np=m->np;
    float3D     *p=m->p;
    
    // compute sulcal depth
    for(i=0;i<np;i++)
    {
        ce=(float3D){ce.x+p[i].x,ce.y+p[i].y,ce.z+p[i].z};
        
        if(i==0) ide=siz=p[i];
        
        if(ide.x<p[i].x) ide.x=p[i].x;
        if(ide.y<p[i].y) ide.y=p[i].y;
        if(ide.z<p[i].z) ide.z=p[i].z;
        
        if(siz.x>p[i].x) siz.x=p[i].x;
        if(siz.y>p[i].y) siz.y=p[i].y;
        if(siz.z>p[i].z) siz.z=p[i].z;
    }
    ce=(float3D){ce.x/(float)np,ce.y/(float)np,ce.z/(float)np};

    max=0;
    for(i=0;i<np;i++)
    {
        n=	pow(2*(p[i].x-ce.x)/(ide.x-siz.x),2) +
            pow(2*(p[i].y-ce.y)/(ide.y-siz.y),2) +
            pow(2*(p[i].z-ce.z)/(ide.z-siz.z),2);

        C[i] = sqrt(n);
        if(C[i]>max)	max=C[i];
    }
    for(i=0;i<np;i++)
        C[i]=C[i]/max;
}
int drawSurface(Mesh *m,char *cmap,char *tiff_path, int toonFlag)
{
    int		i;
    char	*addr;      // memory for tiff image
    int     width=512;  // tiff width
    int     height=512; // tiff height
    float	zoom=1/4.0;
    float3D	a,b,c;
    float3D back={0xff,0xff,0xff};  // background colour
    int     np=m->np;
    int     nt=m->nt;
    int3D   *t=m->t;
    float3D *p=m->p;
    float   *data=m->data;
    float   R,G,B;
    float3D *color;
    int     argc = 1;
    char    *argv[1] = {(char*)"Something"};
    float   min,max,val;
    
    // configure data
    min=minData(m);
    max=maxData(m);
    if(min==max)
    {
        printf("ERROR: In drawSurface, min and max values are the same\n");
        return 0;
    }
    color=(float3D*)calloc(np,sizeof(float3D));
    for(i=0;i<np;i++)
    {
        val=(data[i]-min)/(max-min);
        if(strcmp(cmap,"rainbow")==0)
            rainbow(val,&R,&G,&B);
        else
        if(strcmp(cmap,"grey")==0 || strcmp(cmap,"level2")==0 || strcmp(cmap,"level4")==0)
            greyscale(val,&R,&G,&B);
        else
        {
            printf("ERROR: Unknown colour map %s\n",cmap);
            return 0;
        }
        color[i]=(float3D){R,G,B};
    }
    
    // draw
    if(g_gluInitFlag==0)
    {
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
        g_gluInitFlag=1;
    }
    glutInitWindowSize(width,height);
    glutCreateWindow("meshgeometry");
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_SMOOTH);
    glClearColor(back.x,back.y,back.z,1);

    // init projection
        glViewport(0, 0, (GLsizei)width, (GLsizei)height);
        glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT+GL_STENCIL_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(zoom*width/2,-zoom*width/2,-zoom*height/2,zoom*height/2, -100000.0, 100000.0);

    // prepare drawing
        glMatrixMode (GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt (0,0,-10, 0,0,0, 0,1,0); // eye,center,updir

    // draw
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3,GL_FLOAT,0,(GLfloat*)p);
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(3,GL_FLOAT,0,(GLfloat*)color);
        glDrawElements(GL_TRIANGLES,nt*3,GL_UNSIGNED_INT,(GLuint*)t);

    // toon shading
        if(toonFlag)
        {
            glEnable( GL_CULL_FACE );
            glPolygonMode( GL_BACK, GL_FILL );
            glCullFace( GL_FRONT );

            glPolygonMode(GL_FRONT, GL_LINE);
            glLineWidth(3.0);
            glCullFace(GL_BACK);
            glDepthFunc(GL_LESS);
            glColor3f(0,0,0);
            glBegin(GL_TRIANGLES);
            for(i=0;i<nt;i++)
            {
                a=p[t[i].a];
                b=p[t[i].b];
                c=p[t[i].c];

                glVertex3fv((float*)&a);
                glVertex3fv((float*)&b);
                glVertex3fv((float*)&c);
            }
            glEnd();
            glDisable( GL_CULL_FACE );
        }
    
    // Write image in TIFF format
    addr=(char*)calloc(width*height,sizeof(char)*4);
    glReadPixels(0,0,width,height,GL_RGBA,GL_UNSIGNED_BYTE,addr);
    
    if(strcmp(cmap,"level2")==0)
    for(i=0;i<width*height*4;i++)
        addr[i]=(char)((addr[i]%128>=120 && addr[i]%128<128)?0:255);

    if(strcmp(cmap,"level4")==0)
    for(i=0;i<width*height*4;i++)
        addr[i]=(char)((addr[i]%64>=60 && addr[i]%64<64)?0:255);

    writeTIFF(tiff_path,addr,width,height);
    
    free(color);
    free(addr);
    
    return 0;
}
int edgeLength(Mesh *m)
{
    if(verbose)
        printf("edgeLength\n");
    int         i;
    int         nt=m->nt;
    float3D     *p=m->p;
    int3D       *t=m->t;
    float       length;
    float       s, ss;
    float       avr, std;

    s=ss=0;
    for(i=0;i<nt;i++)
    {
        length=(norm3D(sub3D(p[t[i].a],p[t[i].b]))
               +norm3D(sub3D(p[t[i].b],p[t[i].c]))
               +norm3D(sub3D(p[t[i].c],p[t[i].a])))/3.0;
        s+=length;
        ss+=length*length;
    }
    avr=s/(float)nt;
    std=sqrt(fabs(ss/(float)nt-pow(avr,2)));
    printf("edgeLength: %f±%f\n",avr,std);

    return 0;
}
int edgeLengthMinMax(Mesh *m)
{
    if(verbose)
        printf("edgeLengthMinMax\n");
    int         i;
    int         nt=m->nt;
    float3D     *p=m->p;
    int3D       *t=m->t;
    float       l0,l1,l2;
    float       min, max;

    for(i=0;i<nt;i++)
    {
        l0=norm3D(sub3D(p[t[i].a],p[t[i].b]));
        l1=norm3D(sub3D(p[t[i].b],p[t[i].c]));
        l2=norm3D(sub3D(p[t[i].c],p[t[i].a]));
        if(i==0) { min=max=l0; }
        if(l0<min) min=l0;
        if(l1<min) min=l1;
        if(l2<min) min=l2;
        if(l0>max) max=l0;
        if(l1>max) max=l1;
        if(l2>max) max=l2;
    }
    printf("edgeLengthMinMax: %f %f\n",min,max);
    return 1;
}
int fixflip(Mesh *m)
{
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    NTriRec *NT;
    int     i,j,pos,nflipped=0;
    float3D *nn=(float3D*)calloc(nt,sizeof(float3D));
    float3D tmp;
    
    // compute all triangle normals
    for(i=0;i<nt;i++)
        nn[i]=normal3D(i,m);
    
    // find neighbouring triangles for every vertex
    neighbours(m);
    NT=m->NT;
    
    // find vertices with 1 inverted triangle
    for(i=0;i<np;i++)
    {
        tmp=nn[NT[i].t[0]];    // take the 1st triangle as reference
        pos=0;
        for(j=0;j<NT[i].n;j++)
            if(dot3D(nn[NT[i].t[j]],tmp)>0)
                pos++;
        if(pos==NT[i].n)    // no problem: all triangles have the same orientation
            continue;

        // flip detected: move the vertex to the barycentre of its neighbours to fix it
        nflipped++;
        if(verbose>1)
            printf("Vertex %i is in a flipped triangle. Fixing it.\n",i);
        tmp=(float3D){0,0,0};
        for(j=0;j<NT[i].n;j++)
        {
            if(t[NT[i].t[j]].a!=i)
                tmp=add3D(tmp,p[t[NT[i].t[j]].a]);
            if(t[NT[i].t[j]].b!=i)
                tmp=add3D(tmp,p[t[NT[i].t[j]].b]);
            if(t[NT[i].t[j]].c!=i)
                tmp=add3D(tmp,p[t[NT[i].t[j]].c]);
        }
        tmp=sca3D(tmp,1/(float)(NT[i].n*2));
        p[i]=tmp;
    }
    if(verbose)
        printf("Detected %i flipped triangles\n",nflipped);

    return 0;
}
int fixflipSphere(Mesh *m)
{
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    NTriRec *NT;
    int     i,j,nflipped=0;
    float3D *nn=(float3D*)calloc(nt,sizeof(float3D));
    float3D tmp;
    
    // compute all triangle normals
    for(i=0;i<nt;i++)
        nn[i]=normal3D(i,m);
    
    // find neighbouring triangles for every vertex
    neighbours(m);
    NT=m->NT;
    
    // find vertices with 1 inverted triangle
    for(i=0;i<np;i++)
    {
        tmp=sca3D(p[i], 1/norm3D(p[i]));
        if(dot3D(nn[NT[i].t[j]],tmp)>0)
            continue;

        // flip detected: move the vertex to the barycentre of its neighbours to fix it
        nflipped++;
        if(verbose>1)
            printf("Vertex %i is in a flipped triangle. Fixing it.\n",i);
        tmp=(float3D){0,0,0};
        for(j=0;j<NT[i].n;j++)
        {
            if(t[NT[i].t[j]].a!=i)
                tmp=add3D(tmp,p[t[NT[i].t[j]].a]);
            if(t[NT[i].t[j]].b!=i)
                tmp=add3D(tmp,p[t[NT[i].t[j]].b]);
            if(t[NT[i].t[j]].c!=i)
                tmp=add3D(tmp,p[t[NT[i].t[j]].c]);
        }
        tmp=sca3D(tmp,1/(float)(NT[i].n*2));
        p[i]=tmp;
    }
    if(verbose)
        printf("Detected %i flipped triangles.\n",nflipped);

    return 0;
}
int fixNonmanifold_verts(Mesh *mesh)
{
    NTriRec *ne;
    int3D   *e;
    int e_length;
    int i,j,k,l,found,loop;
    int np=mesh->np;
    float3D *p=mesh->p;
    int3D *t=mesh->t;
    float3D *p1;
    int *t1,t1_length;
    int *i1,i1_length;
    
    neighbours(mesh);
    ne=mesh->NT;

    for(i=0;i<np;i++)
    {
        e=(int3D*)calloc(ne[i].n*3,sizeof(int3D));
        e_length=0;
        for(j=0;j<ne[i].n;j++)
        {
            if(t[ne[i].t[j]].a==i)
                e[e_length++]=(int3D){t[ne[i].t[j]].b,t[ne[i].t[j]].c,ne[i].t[j]};
            else
            if(t[ne[i].t[j]].b==i)
                e[e_length++]=(int3D){t[ne[i].t[j]].c,t[ne[i].t[j]].a,ne[i].t[j]};
            else
                e[e_length++]=(int3D){t[ne[i].t[j]].a,t[ne[i].t[j]].b,ne[i].t[j]};
        }
        
        //printf("p[%i]: ",i); for(j=0;j<e_length;j++) printf("(%i,%i,[%i]) ",e[j].a,e[j].b,e[j].c); printf("\n");
        j=0;
        i1=(int*)calloc(ne[i].n,sizeof(int));
        t1=(int*)calloc(ne[i].n*3,sizeof(int));
        i1_length=0;
        t1_length=0;
        i1[i1_length++]=t1_length;
        t1[t1_length++]=e[j].c; //printf("  t[%i] ",e[j].c);
        while(j<e_length-1)
        {
            loop=0;
            k=j+1;
            while(k<e_length)
            {
                found=1;
                if(e[j].a==e[k].a)
                    e[j]=(int3D){e[j].b,e[k].b,e[j].c};
                else if(e[j].a==e[k].b)
                    e[j]=(int3D){e[j].b,e[k].a,e[j].c};
                else if(e[j].b==e[k].a)
                    e[j]=(int3D){e[j].a,e[k].b,e[j].c};
                else if(e[j].b==e[k].b)
                    e[j]=(int3D){e[j].a,e[k].a,e[j].c};
                else
                {
                    //printf("j%i k%i. *\n",j,k);
                    found=0;
                }
                if(found)
                {
                    //printf("j%i k%i. t[%i], t[%i]\n",j,k,e[j].c,e[k].c);
                    t1[t1_length++]=e[k].c; //printf("t[%i] ",e[k].c);
                    e[k]=e[--e_length];
    
                    //printf("p[%i]: ",i); for(m=0;m<e_length;m++) printf("(%i,%i,[%i]) ",e[m].a,e[m].b,e[m].c); printf("\n");
        
                    if(e[j].a==e[j].b)
                    {
                        //printf("\n");
                        loop=1;
                        j++;
                        if(j<e_length)
                        {
                            i1[i1_length++]=t1_length;
                            t1[t1_length++]=e[j].c; //printf("  t[%i] ",e[j].c);
                        }
                        break;
                    }
                }
                else
                    k++;
            }
        }
        free(e);

        if(ne[i].n==0 && verbose)
            printf("WARNING, %i is isolated: remove it\n",i);
        else
        if(ne[i].n==1 && verbose)
            printf("WARNING, %i is dangling: remove the the vertex and its triangle\n",i);
        else
        if(loop==0 && verbose)
            printf("\nWARNING, %i is in a degenerate region: examine more in detail\n",i);
        else
        if(e_length>1)
        {
            if(verbose)
                printf("WARNING, %i has %i loops: split the vertex into %i vertices and remesh\n",i,e_length,e_length);
            
            p1=(float3D*)calloc(np+e_length-1,sizeof(float3D));
            for(l=0;l<np;l++)
                p1[l]=p[l];     // copy the original vertices
            for(l=0;l<e_length-1;l++)
                p1[np+l]=p[i]; // make e_length-1 copies of vertex i at the end of the vertex vector
            
            k=1;
            for(j=0;j<t1_length;j++)
            {
                if(i1[k]==j)
                {
                    if(verbose)
                        printf(" | ");
                    k++;
                }
                if(verbose)
                    printf("%i ",t1[j]);
                if(k>1)
                {
                    if(t[t1[j]].a==i)   t[t1[j]].a=np+(k-2);
                    if(t[t1[j]].b==i)   t[t1[j]].b=np+(k-2);
                    if(t[t1[j]].c==i)   t[t1[j]].c=np+(k-2);
                }
            }
            if(verbose)
                printf("\n");
            free(mesh->p);
            mesh->p=p1;
            mesh->np=np+e_length-1;
            return 1;
        }
        else
        {
            // printf("%i: ok\n",i);
        }
        free(t1);
        free(i1);
    }
    
    return 0;
}
int fixnonmanifold_tris(Mesh *mesh)
{
	int     i,j,found;
	int3D   *t=mesh->t,t1;
	int     nt=mesh->nt;
	int     np=mesh->np;
	NTriRec *ne;
	float3D *p1,*p=mesh->p;
	
	found=nonmanifold_tris(mesh);
	
	if(found==0)
	{
	    printf("no nonmanifold triangles found\n");
	    return 0;
	}
	
	p1=(float3D*)calloc(np+found*3,sizeof(float3D));
	for(i=0;i<np;i++)
	    p1[i]=p[i];

    neighbours(mesh);
    ne=mesh->NT;
    
    found=0;
    for(i=0;i<nt;i++)
    {
        for(j=0;j<ne[t[i].a].n;j++)
        if(ne[t[i].a].t[j]!=i)
        {
            t1=t[ne[t[i].a].t[j]];
            if( (t1.a==t[i].a||t1.a==t[i].b||t1.a==t[i].c)&&
                (t1.b==t[i].a||t1.b==t[i].b||t1.b==t[i].c)&&
                (t1.c==t[i].a||t1.c==t[i].b||t1.c==t[i].c))
            {
                found++;
                printf("triangle %i doubles triangle %i. Fixing it\n",i,ne[t[i].a].t[j]);
                break;
            }
        }
    }
    printf("%i double triangles found\n",found/2);
    return found/2;
}
int fixSmall(Mesh *m)
{
    int     nt=m->nt;
    int     *T;
    float3D *p=m->p;
    int3D   *t=m->t;
    NTriRec *NT;
    int     i,j,k,l;
    float3D x,y;
    float   angle;
    float3D tmp;
    int     didFix;
    
    // find neighbouring triangles for every vertex
    neighbours(m);
    NT=m->NT;

    // compute all triangle areas
    didFix=1;
    for(l=0;l<5 && didFix;l++)
    {
        didFix=0;
        for(i=0;i<nt;i++)
        {
            T=(int*)&(t[i]);
            for(k=0;k<3;k++)
            {
                x=sub3D(p[T[(k+1)%3]],p[T[k]]); // vector b-a
                x=sca3D(x,1/norm3D(x));
                y=sub3D(p[T[(k+2)%3]],p[T[k]]); // vector c-a
                y=sca3D(y,1/norm3D(y));
                angle=acos(dot3D(x,y))*180/M_PI;
                if(angle>150)
                {
                    printf("Vertex %i makes a triangle with an angle of %g degrees. Fixing it.\n",T[k],angle);
                    // triangle is fixed by moving the offending vertex to the barycentre of its neighbours
                    tmp=(float3D){0,0,0};
                    for(j=0;j<NT[T[k]].n;j++)
                    {
                        if(t[NT[T[k]].t[j]].a!=T[k])
                            tmp=add3D(tmp,p[t[NT[T[k]].t[j]].a]);
                        if(t[NT[T[k]].t[j]].b!=T[k])
                            tmp=add3D(tmp,p[t[NT[T[k]].t[j]].b]);
                        if(t[NT[T[k]].t[j]].c!=T[k])
                            tmp=add3D(tmp,p[t[NT[T[k]].t[j]].c]);
                    }
                    tmp=sca3D(tmp,1/(float)(NT[T[k]].n*2));
                    p[T[k]]=tmp;
                    didFix=1;
                }
            }
        }
    }
    if(l==5)
        printf("WARNING: fixSmall: There may still be small triangles\n");
    
    return 0;
}
int flip(Mesh *m)
{
    if(verbose) printf("* flip\n");

    int     nt=m->nt;
    int3D   *t=m->t;
    int     i;    

    for(i=0;i<nt;i++)
    {
       // printf("%i %i %i\n",t[i].a,t[i].b,t[i].c);
        t[i]=(int3D){t[i].a,t[i].c,t[i].b};
    }
    return 0;
}
int foldLength(Mesh *m)
{
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    float   *data=m->data;
    int     i,j;
    float   length=0,a,x;
    float3D p0[3];
    
    if(data==NULL)
    {
        printf("ERROR: there is no data to calculate foldLength (use -curv before)\n");
        return 1;
    }

    for(i=0;i<nt;i++)
    {
        j=0;
        if(data[t[i].a]*data[t[i].b]<0)
        {
            a=fabs(data[t[i].a]);
            x=a/(a+fabs(data[t[i].b]));
            p0[j++]=add3D(sca3D(p[t[i].a],1-x),sca3D(p[t[i].b],x));
        }
        if(data[t[i].b]*data[t[i].c]<0)
        {
            a=fabs(data[t[i].b]);
            x=a/(a+fabs(data[t[i].c]));
            p0[j++]=add3D(sca3D(p[t[i].b],1-x),sca3D(p[t[i].c],x));
        }
        if(data[t[i].c]*data[t[i].a]<0)
        {
            a=fabs(data[t[i].c]);
            x=a/(a+fabs(data[t[i].a]));
            p0[j++]=add3D(sca3D(p[t[i].c],1-x),sca3D(p[t[i].a],x));
        }
        if(j==2)
            length+=norm3D(sub3D(p0[0],p0[1]));
    }
    printf("foldLength: %f\n",length/2.0);
    return 0;
}
int icurvature(int iter, Mesh *m)
{
    if(verbose)
        printf("icurvature\n");
    int     *np=&(m->np);
    float   *data=m->data;
    float   *tmp;
    float   absmax;
    int     i,k;

    tmp=(float*)calloc(*np,sizeof(float));
    for(k=0;k<iter;k++)
    {
        curvature(tmp,m);
        for(i=0;i<*np;i++)
            data[i]+=tmp[i]*(1+k/(float)iter);
        
        smooth(m);
        smooth(m);
    }
    absmax=-1;
    for(i=0;i<*np;i++)
        absmax = (fabs(data[i])>absmax)?fabs(data[i]):absmax;
    for(i=0;i<*np;i++)
        data[i]/=absmax;
    
    return 0;
}
/**
  * @function invert
  * @desc Invert a specific coordinate of the mesh vertices
  * @param axis string The axis, x, y or z
  * @param m pointer Pointer to the mesh to invert
  */
int invert(char *axis, Mesh *m)
{
    int x, y, z;
    int i;
    int *np=&(m->np);
    float3D *p=m->p;

    x = y = z = 1;

    if(strcmp(axis, "x") == 0)
    {
        x = -1;
    }
    else
    if(strcmp(axis, "y") == 0)
    {
        y = -1;
    }
    else
    if(strcmp(axis, "z") == 0)
    {
        z = -1;
    }

    for(i=0;i<*np;i++)
    {
        p[i].x *= x;
        p[i].y *= y;
        p[i].z *= z;
    }

    return 0;
}
int isolatedVerts(Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    int3D   *t=m->t;
    int     *n;
    int     i,sum;
    
    n=(int*)calloc(*np,sizeof(int));
    for(i=0;i<*nt;i++)
    {
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    sum=0;
    for(i=0;i<*np;i++)
    {
        if(n[i]==0)
            sum++;
    }
    free(n);

    printf("isolatedVerts: %i\n",sum);

    return 0;
}
int removeIsolatedVerts(Mesh *m)
{
    if(verbose)
        printf("* removeIsolatedVerts\n");

    int     np0,*np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    int     *n,*ip;
    int     i,j;
    
    ip=(int*)calloc(*np,sizeof(int));
    n=(int*)calloc(*np,sizeof(int));
    
    // count neighbours
    for(i=0;i<*nt;i++)
    {
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    
    // make a lookup table for the new vertex indices
    // and re-index the vertices
    j=0;
    for(i=0;i<*np;i++)
    {
        if(n[i]>0)
        {
            ip[i]=j;    // lookup table
            p[j]=p[i];  // re-indexed vertices
            j++;
        }
        else
            ip[i]=-1;   // this doesn't really matter (it'll never be read)
    }
    np0=*np;
    *np=j;  // j is the new number of vertices
    
    // re-index triangles
    for(i=0;i<*nt;i++)
        t[i]=(int3D){ip[t[i].a],ip[t[i].b],ip[t[i].c]};
    free(n);
    
    if(verbose)
        printf("%i vertices were removed\n",np0-j);
    
    return 0;
}
int removeVerts(Mesh *m)
{
    if(verbose)
        printf("* removeVerts\n");

    int     np0,*np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    int     *n,*ip;
    int     i,j,sum;
    float   *data=m->data;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data to removeVerts (use -curv, for example)\n");
        return 1;
    }
    
    // remove triangles that contain vertices with negative data values,
    // vertices with negative data will be then isolated
    sum=0;
    for(i=0;i<*nt;i++)
    {
        if(data[t[i].a]<0 || data[t[i].b]<0 || data[t[i].c]<0)
        {
            (*nt)--;
            for(j=i;j<*nt;j++)
                t[j]=t[j+1];
            sum++;
            i--;
        }
    }
    if(verbose)
        printf("%i triangles with negative data vertices removed\n",sum);
    
    // remove isolated vertices
    ip=(int*)calloc(*np,sizeof(int));
    n=(int*)calloc(*np,sizeof(int));
    
    // count neighbours
    for(i=0;i<*nt;i++)
    {
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    
    // make a lookup table for the new vertex indices
    // and re-index the vertices
    j=0;
    for(i=0;i<*np;i++)
    {
        if(n[i]>0)
        {
            ip[i]=j;    // lookup table
            p[j]=p[i];  // re-indexed vertices
            data[j]=data[i]; // re-index data
            j++;
        }
        else
            ip[i]=-1;   // this doesn't really matter (it'll never be read)
    }
    np0=*np;
    *np=j;  // j is the new number of vertices
    
    // re-index triangles
    for(i=0;i<*nt;i++)
        t[i]=(int3D){ip[t[i].a],ip[t[i].b],ip[t[i].c]};
    free(n);
    
    if(verbose)
        printf("%i vertices were removed\n",np0-j);
    
    return 0;
}
int laplace(float lambda, Mesh *m)
{
    if(verbose)
        printf("laplace smooth\n");
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp,x,dx;
    int     *n;
    int     i;
    
    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    for(i=0;i<*nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
        tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    for(i=0;i<*np;i++)
    {
        if(n[i]==0)
        {
            printf("WARNING: isolated vertex %i\n",i);
            x=tmp[i];
        }
        else
            x=sca3D(tmp[i],1/(float)n[i]);
        dx=sub3D(x,p[i]);
        p[i]=add3D(p[i],sca3D(dx,lambda));    // p=p+l(x-p)
    }
    free(tmp);
    free(n);
    return 0;
}
int laplaceSelection(float lambda, Mesh *m)
{
    if(verbose)
        printf("laplace smooth on selected vertices\n");
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp,x,dx;
    int     *n;
    int     i, nprocessed=0;
    char    *selection=m->selection;

    tmp=(float3D*)calloc(np,sizeof(float3D));
    n=(int*)calloc(np,sizeof(int));
    for(i=0;i<nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
        tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }

    for(i=0;i<np;i++)
        if(selection[i]>0)
        {
            if(n[i]==0)
            {
                printf("WARNING: isolated vertex %i\n",i);
                x=tmp[i];
            }
            else
                x=sca3D(tmp[i],1/(float)n[i]);
            dx=sub3D(x,p[i]);
            p[i]=add3D(p[i],sca3D(dx,lambda));    // p=p+l(x-p)
            nprocessed++;
        }
    free(tmp);
    free(n);
    if(verbose)
        printf("%i vertices laplace smoothed\n",nprocessed);

    return 0;
}
int level(float v, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    float   **data=&(m->data);
    float   da,db,a,x;
    int     i,k,l,v1,v2,flag;
    float3D p0;
    int     max,i0,*T,*n;
    int2D   *newplut;
    int     nnp,nnt,ip,it;
    float3D *newp;
    int3D   *newt;
    float   *newdata;
    float   datamax,datamin;
    
    if((*data)==NULL)
    {
        printf("ERROR: there is no data to calculate level (use -curv, for example)\n");
        return 1;
    }
    
    // Zeroth pass
    datamax=datamin=(*data)[0]-v;
    for(i=0;i<*np;i++)
    {
        if((*data)[i]-v>datamax)
            datamax=(*data)[i]-v;
        if((*data)[i]-v<datamin)
            datamin=(*data)[i]-v;
    }
    for(i=0;i<*np;i++)
        if(fabs((*data)[i]-v)/(datamax-datamin)<0.05)
            (*data)[i]=v;

    // First pass
    n=(int*)calloc(*np,sizeof(int));    // number of times that a given vertex is in an edge that will be cut
    nnp=0;  // total number of new vertices that will be created
    nnt=0;  // total number of new triangles that will be created
    for(i=0;i<(*nt);i++)
    {
        T=(int*)&((*t)[i]);
        for(k=0;k<3;k++)
        {
            da=(*data)[T[k]]-v;
            db=(*data)[T[(k+1)%3]]-v;
            if(da*db<0)
            {
                n[T[k]]++;
                n[T[(k+1)%3]]++;
                nnp++;
                nnt++;
                if(verbose)
                    printf("1st pass. tri %i (%i,%i,%i), edge %i\n",i,T[0],T[1],T[2],k);
            }
        }
    }
    nnp/=2;
    max=0;
    for(i=0;i<(*np);i++)
        if(n[i]>max)
        {
            max=n[i];
            i0=i;
        }
    free(n);
    max=max/2;
    
    if(verbose)
    {
        printf("MSG: level will add %i vertices and %i triangles.\n",nnp,nnt);
        printf("     The maximum number of times a vertex (%i) will be involved\n",i0);
        printf("     in a cut is %i\n",max);
    }

    // Second pass
    newplut=(int2D*)calloc(((*np)+nnp)*max,sizeof(int2D));
    newp=(float3D*)calloc((*np)+nnp,sizeof(float3D));
    for(i=0;i<(*np);i++)
        newp[i]=(*p)[i];
    newt=(int3D*)calloc((*nt)+nnt,sizeof(int3D));
    for(i=0;i<(*nt);i++)
        newt[i]=(*t)[i];
    newdata=(float*)calloc((*np)+nnp,sizeof(float));
    for(i=0;i<(*np);i++)
        newdata[i]=(*data)[i];
    ip=0;
    it=0;
    for(i=0;i<(*nt)+nnt;i++)
    {
        T=(int*)&(newt[i]);
        for(k=0;k<3;k++)
        {
            // check if the edge (a,b) has to be cut
            da=newdata[T[k]]-v;
            db=newdata[T[(k+1)%3]]-v;
            if(da*db<0)
            {   // yes, it has to be cut
                if(T[k]<T[(k+1)%3])
                {
                    v1=T[k];
                    v2=T[(k+1)%3];
                }
                else
                {
                    v1=T[(k+1)%3];
                    v2=T[k];
                }
                if(verbose)
                    printf("%i. should cut vtx %i,%i\n",i,v1,v2);
                // check if it was already cut (flag==0 means no)
                flag=0;
                for(l=0;l<max;l++)
                {
                    if(newplut[v1*max+l].a==0)      // 0 means empty
                        break;
                    if(newplut[v1*max+l].a==v2+1)
                    {   // yes, already cut, then get the vertex index
                        i0=newplut[v1*max+l].b;
                        if(verbose)
                            printf("%i. recycled vtx %i\n",i,i0);
                        flag=1;
                        break;
                    }
                }
                if(flag==0)
                {   // no, never cut, then make it and add it to the list
                    da=newdata[v1]-v;
                    db=newdata[v2]-v;
                    a=fabs(da);
                    x=a/(a+fabs(db));
                    p0=add3D(sca3D(newp[v1],1-x),sca3D(newp[v2],x));
                    i0=(*np)+ip;
                    if(i0>=(*np)+nnp)
                    {
                        printf("ERROR: at triangle %i, vtx out of bounds\n",i);
                        break;
                    }
                    if(verbose)
                        printf("%i. added vtx %i (%g)\n",i,i0,x);
                    newp[i0]=p0;
                    for(l=0;l<max;l++)
                        if(newplut[v1*max+l].a==0)
                            break;
                    if(l>=max)
                        printf("ERROR: at triangle %i, vtx %i included too many times\n",i,v1);
                    newplut[v1*max+l]=(int2D){v2+1,i0};
                    ip++;
                }
                if(i0>=(*np)+nnp)
                {
                    printf("ERROR: at triangle %i, vtx out of bounds\n",i);
                    break;
                }
                newt[(*nt)+it]=(int3D){T[k],i0,T[(k+2)%3]};  // add new triangle
                if(verbose)
                {
                    printf("%i. added tri %i,%i,%i\n",i,T[k],i0,T[(k+2)%3]);
                    printf("%i. changed tri %i,%i,%i, to ",i,T[k],T[(k+1)%3],T[(k+2)%3]);
                }
                T[k]=i0;                                           // change current triangle
                if(verbose)
                    printf("%i,%i,%i\n",T[k],T[(k+1)%3],T[(k+2)%3]);
                it++;
            }
        }
    }
    
    // replace old mesh with new mesh
    free(*p);
    free(*t);
    free(*data);
    if(m->NT)
        free(m->NT);
    *np+=nnp;
    *nt+=nnt;
    *p=newp;
    *t=newt;
    *data=newdata;
    free(newplut);
    
    return 0;
}
int lissencephalic(int iter, Mesh *m)
/*
 given a mesh (a brain) with data initialised to mean curvature (-curv),
 1st add vertices to the regions where the mean curvature =0, and then
 smooth out everything but those vertices. The idea is to try to find
 the surface of the coast, smoothing out hills and valleys. It's not
 quite working, though.
*/
{
    int     np;
    float3D *tmp,*p;
    float   *data=m->data;
    int     i,j,k;
    
    level(0,m);
    np=m->np;
    p=m->p;
    data=m->data;
    tmp=(float3D*)calloc(m->np,sizeof(float3D));
    for(i=0;i<np;i++)
        tmp[i]=p[i];
    for(j=0;j<iter;j++)
    {
        smooth(m);
        
        for(k=0;k<10;k++)
        {
            for(i=0;i<np;i++)
                if(data[i]==0)
                    p[i]=tmp[i];
            taubin(0.5,-0.53,1,m);
        }
    }
    free(tmp);
    
    return 0;
}
float maxData(Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    float   max;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    max=data[0];
    for(i=1;i<np;i++)
        max=(data[i]>max)?(data[i]):max;
    printf("maxData: %f\n",max);
    return max;
}
float meanData(Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    float   mean;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    mean=0;
    for(i=0;i<np;i++)
        mean+=data[i];
    mean/=(float)np;
    printf("meanData: %f\n",mean);
    return mean;
}
float minData(Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    float   min;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    min=data[0];
    for(i=1;i<np;i++)
        min=(data[i]<min)?(data[i]):min;
    printf("minData: %f\n",min);
    return min;
}
float mirror(Mesh *m, char *coord)
{
    int np=m->np;
    float3D *p=m->p;
    int i;
    float3D c={0,0,0};
    
    // compute barycentre
    for(i=0;i<np;i++)
        c=add3D(c,p[i]);
    c=sca3D(c,1/(float)np);
    
    // mirror
    switch((char)coord[0])
    {
        case 'x':
            for(i=0;i<np;i++)
                p[i]=(float3D){2*c.x-p[i].x,p[i].y,p[i].z};
            break;
        case 'y':
            for(i=0;i<np;i++)
                p[i]=(float3D){p[i].x,2*c.y-p[i].y,p[i].z};
            break;
        case 'z':
            for(i=0;i<np;i++)
                p[i]=(float3D){p[i].x,p[i].y,2*c.z-p[i].z};
            break;
    }
    
    return 0;
}
float stdData(Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    float   s,ss,std;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    s=ss=0;
    for(i=0;i<np;i++)
    {
        s+=data[i];
        ss+=data[i]*data[i];
    }
    std=sqrt(fabs(ss/(float)np-pow(s/(float)np,2)));
    printf("stdData: %f\n",std);
    return std;
}
/*
3 5 4 5		c1: e1.a<e2.a					r -1
4 5 3 5		c2: e1.a>e2.a					r  1
3 5 3 7		c3: e1.a==e2.a, e1.b<e2.b		r -1
3 7 3 5		c4: e1.a==e2.a, e1.b>e2.b		r  1
            c5: e1.a==e2.a, e1.b==e2.b		r  0
*/
int compareEdges (const void *a, const void *b)
{
    int2D	e1=*(int2D*)a;
    int2D	e2=*(int2D*)b;

    if(e1.a==e2.a)
    {
        if(e1.b==e2.b)
            return 0;
        else
        if(e1.b<e2.b)
            return -1;
        else
            return 1;
    }
    else
    {
        if(e1.a<e2.a)
            return -1;
        else
            return	1;
    }
}
int nonmanifold_verts(Mesh *mesh)
{
    NTriRec *ne;
    int3D   *e;
    int e_length;
    int i,j,k,found,loop;
    int np=mesh->np;
    int3D *t=mesh->t;
    int *t1,t1_length;
    int *i1,i1_length;
    int sum=0;

    neighbours(mesh);
    ne=mesh->NT;

    printf("non manifold vertices\n");
    for(i=0;i<np;i++)
    {
        // store all the edges in the triangles connected to vertex p[i]
        // that do not contain vertex p[i]
        e=(int3D*)calloc(ne[i].n*3,sizeof(int3D));
        e_length=0;
        for(j=0;j<ne[i].n;j++)
        {
            if(t[ne[i].t[j]].a==i)
                e[e_length++]=(int3D){t[ne[i].t[j]].b,t[ne[i].t[j]].c,ne[i].t[j]};
            else
            if(t[ne[i].t[j]].b==i)
                e[e_length++]=(int3D){t[ne[i].t[j]].c,t[ne[i].t[j]].a,ne[i].t[j]};
            else
                e[e_length++]=(int3D){t[ne[i].t[j]].a,t[ne[i].t[j]].b,ne[i].t[j]};
        }
/*        
        printf("p[%i]: ",i); for(j=0;j<e_length;j++) printf("(%i,%i) ",e[j].a,e[j].b); printf("\n");
*/       
        // scan the list of edges, if 2 edges share a vertex,
        // delete the vertex and connect the points directly.
        // at the end, there should be only one edge remaining
        // connecting one vertex to itself. All remaining a-a
        // edges represent supplementary loops, then, non-manifoldness
        j=0;
        i1=(int*)calloc(ne[i].n,sizeof(int));
        t1=(int*)calloc(ne[i].n*3,sizeof(int));
        i1_length=0;
        t1_length=0;
        i1[i1_length++]=t1_length;
        t1[t1_length++]=e[j].c;
        while(j<e_length-1)
        {
            loop=0;
            k=j+1;
            while(k<e_length)
            {
                found=1;
                if(e[j].a==e[k].a)
                    e[j]=(int3D){e[j].b,e[k].b,e[j].c};
                else if(e[j].a==e[k].b)
                    e[j]=(int3D){e[j].b,e[k].a,e[j].c};
                else if(e[j].b==e[k].a)
                    e[j]=(int3D){e[j].a,e[k].b,e[j].c};
                else if(e[j].b==e[k].b)
                    e[j]=(int3D){e[j].a,e[k].a,e[j].c};
                else
                    found=0;
                if(found)
                {
                    t1[t1_length++]=e[k].c; //printf("t[%i] ",e[k].c);
                    e[k]=e[--e_length];
                    if(e[j].a==e[j].b)
                    {
                        //printf("\n");
                        loop=1;
                        j++;
                        if(j<e_length)
                        {
                            i1[i1_length++]=t1_length;
                            t1[t1_length++]=e[j].c; //printf("  t[%i] ",e[j].c);
                        }
                        break;
                    }
                }
                else
                    k++;
            }
        }
        // printf("p[%i]: ",i); for(j=0;j<e_length;j++) printf("(%i,%i) ",e[j].a,e[j].b); printf("\n");
        free(e);
        free(t1);
        free(i1);

        if(ne[i].n==0 && verbose)      printf("WARNING, %i is isolated: remove it\n",i);
        else if(ne[i].n==1 && verbose) printf("WARNING, %i is dangling: remove the the vertex and its triangle\n",i);
        else if(loop==0 && verbose)    printf("WARNING, %i is in a degenerate region: examine more in detail\n",i);
        else if(e_length>1)
        {
            if(verbose)
                printf("WARNING, %i has %i loops: split the vertex into %i vertices and remesh\n",i,e_length,e_length);
            sum++;
        }
        else
        {
            // printf("%i: ok\n",i);
        }
    }
    return sum;
}
void nonmanifold_eds(Mesh *mesh)
{
    int     i,j,k,equal;
    int     n;  // # manifold edges
    int3D	*e;
    int		*t;	
    int3D   *tris=mesh->t;
    int     nt=mesh->nt;

    // make a list of all edges
    e=(int3D*)calloc(nt*3,sizeof(int3D));
    for(i=0;i<nt;i++)
    {
        t=(int*)&(tris[i]);
        for(j=0;j<3;j++)
        {
            if(t[j]<t[(j+1)%3])
            {
                e[3*i+j].a=t[j];        // 1st vertex
                e[3*i+j].b=t[(j+1)%3];  // 2nd vertex
                e[3*i+j].c=i;           // # triangle
            }
            else
            {
                e[3*i+j].a=t[(j+1)%3];  // 1st vertex
                e[3*i+j].b=t[j];        // 2nd vertex
                e[3*i+j].c=i;           // # triangle
            }
        }
    }
    
    // sort edges
    qsort(e,nt*3,sizeof(int3D),compareEdges);

    //for(i=0;i<nt*3;i++) printf("e[%i]=(%i,%i), t[%i]\n",i,e[i].a,e[i].b,e[i].c);
    
    printf("non manifold edges\n");
    // count nonmanifold edges
    j=0;
    k=0;
    n=0;
    do
    {
        k=1;
        do
        {
    		equal=0;
            if(e[j].a==e[j+k].a && e[j].b==e[j+k].b)
            {
                k++;
                equal=1;
            }
            else
            {
                if(k!=2) {
                    printf("edge (%i,%i) belongs to %i triangles\n",e[j].a,e[j].b,k);
                    n++;
                }
                j+=k;
                k=1;
            }
        }
        while(equal);
    }
    while(j<nt*3);
    printf("nonmanifold edges:%i\n",n);
}
int nonmanifold_tris(Mesh *mesh)
{
    int     i,j,found;
    int3D   *t=mesh->t,t1;
    int     nt=mesh->nt;
    NTriRec *ne;

    neighbours(mesh);
    ne=mesh->NT;
    
    found=0;
    for(i=0;i<nt;i++)
    {
        for(j=0;j<ne[t[i].a].n;j++)
        if(ne[t[i].a].t[j]!=i)
        {
            t1=t[ne[t[i].a].t[j]];
            if( (t1.a==t[i].a||t1.a==t[i].b||t1.a==t[i].c)&&
                (t1.b==t[i].a||t1.b==t[i].b||t1.b==t[i].c)&&
                (t1.c==t[i].a||t1.c==t[i].b||t1.c==t[i].c))
            {
                found++;
                printf("triangle %i doubles triangle %i\n",i,ne[t[i].a].t[j]);
                break;
            }
        }
    }
    printf("%i double triangles found\n",found/2);
    return found/2;
}
void normalise(Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p,c={0,0,0},x;
    int     i;
    
    for(i=0;i<np;i++)
        c=add3D(c,p[i]);
    c=sca3D(c,1/(float)np);
    for(i=0;i<np;i++)
    {
        x=sub3D(p[i],c);
        x=sca3D(x,100/norm3D(x));
        p[i]=x;
    }
}
void printBarycentre(Mesh *m)
{
    int     *np=&(m->np);
    float3D *p=m->p;
    int     i;
    float3D centre={0,0,0};
    
    for(i=0;i<*np;i++)
        centre=add3D(centre,p[i]);
    centre=sca3D(centre,1/(float)*np);
    printf("barycentre %g,%g,%g\n",centre.x,centre.y,centre.z);
}
void printCentre(Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p;
    int     i;
    float3D mi,ma,centre;
    
    mi=ma=p[0];
    for(i=0;i<np;i++)
    {
        mi.x=(mi.x>p[i].x)?p[i].x:mi.x;
        mi.y=(mi.y>p[i].y)?p[i].y:mi.y;
        mi.z=(mi.z>p[i].z)?p[i].z:mi.z;
        ma.x=(ma.x<p[i].x)?p[i].x:ma.x;
        ma.y=(ma.y<p[i].y)?p[i].y:ma.y;
        ma.z=(ma.z<p[i].z)?p[i].z:ma.z;
    }
    centre=(float3D){(mi.x+ma.x)/2.0,(mi.y+ma.y)/2.0,(mi.z+ma.z)/2.0};
    printf("centre %g,%g,%g\n",centre.x,centre.y,centre.z);
}
int relax(char *path, Mesh *m0,int iformat)
{
    // m0: actual mesh
    // m1: target mesh
    Mesh    *m1;
    int     i,j;
    float3D a,b,c,o,a1,b1,c1;
    int     *n;
    int     niter;
    float   J;
    float3D *f,*g,*p0,*p1,zero={0,0,0};
    int3D   *t;
    float   area0,area1;
    float   alpha,beta;
    float   dot0,dot1;
    float3D *s0,*s1,*nn0,*nn1,nn;
    
    niter=400000;

    alpha=0.5;
    beta=0.5;
    
    m1=(Mesh*)calloc(1,sizeof(Mesh));
    loadMesh(path,m1,iformat);
    
    p0=m0->p;
    p1=m1->p;
    t=m0->t;
    
    printf("area source: %g\n",area(m0));
    printf("area target: %g\n",area(m1));
    
    f=(float3D*)calloc(m0->np,sizeof(float3D));
    g=(float3D*)calloc(m0->np,sizeof(float3D));
    n=(int*)calloc(m0->np,sizeof(int));
    
    s0=(float3D*)calloc(m0->np,sizeof(float3D));
    s1=(float3D*)calloc(m0->np,sizeof(float3D));
    nn0=(float3D*)calloc(m0->np,sizeof(float3D));
    nn1=(float3D*)calloc(m0->np,sizeof(float3D));
    
    for(i=0;i<m0->nt;i++)
    {
        n[t[i].a]++;
        n[t[i].b]++;
        n[t[i].c]++;
    }
    
    // compute target smoothing direction
    for(i=0;i<m0->nt;i++)
    {
        s1[t[i].a]=add3D(s1[t[i].a],add3D(p1[t[i].b],p1[t[i].c]));
        s1[t[i].b]=add3D(s1[t[i].b],add3D(p1[t[i].c],p1[t[i].a]));
        s1[t[i].c]=add3D(s1[t[i].c],add3D(p1[t[i].a],p1[t[i].b]));
    }
    for(i=0;i<m0->np;i++)
        s1[i]=sub3D(p1[i],sca3D(s1[i],1/(float)(2*n[i])));
    // compute target normal direction
    for(i=0;i<m0->nt;i++)
    {
        nn=cross3D(sub3D(p1[t[i].b],p1[t[i].a]),sub3D(p1[t[i].c],p1[t[i].a]));
        nn=sca3D(nn,1/norm3D(nn));
        nn1[t[i].a]=add3D(nn1[t[i].a],nn);
        nn1[t[i].b]=add3D(nn1[t[i].b],nn);
        nn1[t[i].c]=add3D(nn1[t[i].c],nn);
    }
    for(i=0;i<m0->np;i++)
        nn1[i]=sca3D(nn1[i],1/(float)n[i]);

    // relax source mesh to target
    for(j=0;j<niter;j++)
    {
        if(j*1000/niter>(j-1)*1000/niter)
        {
            //laplace(0.1,m0);
            printf("%i%%\n",j*1000/niter);
        }
        // 1. Relax curvature
            // compute smoothing direction
        for(i=0;i<m0->nt;i++)
        {
            s0[t[i].a]=add3D(s0[t[i].a],add3D(p0[t[i].b],p0[t[i].c]));
            s0[t[i].b]=add3D(s0[t[i].b],add3D(p0[t[i].c],p0[t[i].a]));
            s0[t[i].c]=add3D(s0[t[i].c],add3D(p0[t[i].a],p0[t[i].b]));
        }
        for(i=0;i<m0->np;i++)
            s0[i]=sub3D(p0[i],sca3D(s0[i],1/(float)(2*n[i])));
            // compute normal direction
        for(i=0;i<m0->nt;i++)
        {
            nn=cross3D(sub3D(p0[t[i].b],p0[t[i].a]),sub3D(p0[t[i].c],p0[t[i].a]));
            nn=sca3D(nn,1/norm3D(nn));
            nn0[t[i].a]=add3D(nn0[t[i].a],nn);
            nn0[t[i].b]=add3D(nn0[t[i].b],nn);
            nn0[t[i].c]=add3D(nn0[t[i].c],nn);
        }
        for(i=0;i<m0->np;i++)
            nn0[i]=sca3D(nn0[i],1/(float)n[i]);
            // compute dot product
        for(i=0;i<m0->np;i++)
        {
            dot0=dot3D(s0[i],nn0[i]);
            dot1=dot3D(s1[i],nn1[i]);
            g[i]=sub3D(sca3D(nn0[i],dot1),sca3D(nn0[i],dot0));
            //printf("d0:%g d1:%g\n",dot0,dot1);
        }

        // 2. Relax surface area (only along the normal direction)
        for(i=0;i<m0->nt;i++)
        {
            a=p0[t[i].a];
            b=p0[t[i].b];
            c=p0[t[i].c];
            o=sca3D(add3D(a,add3D(b,c)),1/3.0);
            area0=triangle_area(a,b,c);
            area1=triangle_area(p1[t[i].a],p1[t[i].b],p1[t[i].c]);
            if(area0==0)
            {
                printf("denominator==0 at triangle %i\n",i);
                continue;
            }
            J=area1/area0;
            a1=add3D(sca3D(sub3D(a,o),pow(J,0.5)),o);
            b1=add3D(sca3D(sub3D(b,o),pow(J,0.5)),o);
            c1=add3D(sca3D(sub3D(c,o),pow(J,0.5)),o);
            
            // printf("%g %g %g\n",area0,area1,triangle_area(a1,b1,c1));
        
            f[t[i].a]=add3D(f[t[i].a],sub3D(a1,a));
            f[t[i].b]=add3D(f[t[i].b],sub3D(b1,b));
            f[t[i].c]=add3D(f[t[i].c],sub3D(c1,c));
        }
        for(i=0;i<m0->np;i++)
        {
            f[i]=sca3D(f[i],1/(float)n[i]);
            //f[i]=sca3D(nn0[i],dot3D(f[i],nn0[i]));
        }
    
        // 3. Apply forces
        for(i=0;i<m0->np;i++)
        {
            p0[i]=add3D(p0[i],sca3D(f[i],alpha));
            p0[i]=add3D(p0[i],sca3D(g[i],beta));
        }
        
        // 4. Reinitialise
        for(i=0;i<m0->np;i++)
        {
            s0[i]=zero;
            nn0[i]=zero;
            f[i]=zero;
            g[i]=zero;          
            //printf("d0:%g d1:%g\n",dot0,dot1);
        }
    }
 
     printf("area source: %g\n",area(m0));
   
    // free m1
    freeMesh(m1);
    free(f);
    //free(g);
    free(n);
    return 0;
}
int repulse(Mesh *m)
{
    if(verbose)
        printf("repulse vertices\n");
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *f,v;
    float3D a,b,c;
    float3D at,bt,ct;
    float3D ab, bc, ca;
    float3D fa,fb,fc;
    int     *n;
    float   alpha=0.3;
    float   *beta;
    float   tau=0.2;
    float   *deg;
    float   *min;
    float   A,s,ar;
    float   dab,dbc,dca;
    float   ha,hb,hc;
    float   max;
    int     i;
    float   totDisp=0;

    f=(float3D*)calloc(np,sizeof(float3D));
    n=(int*)calloc(np,sizeof(int));
    beta=(float*)calloc(np,sizeof(float));
    deg=(float*)calloc(np,sizeof(float));
    min=(float*)calloc(nt,sizeof(float));
    m->data=(float*)calloc(np,sizeof(float));

    for(i=0;i<np;i++)
        beta[i]=100;
    for(i=0;i<nt;i++)
    {
        // triangle vertices
        a=p[t[i].a];
        b=p[t[i].b];
        c=p[t[i].c];

        // triangle area
        A=triangle_area(a, b, c);

        // edges and edge lengths
        ab=sub3D(a,b);
        bc=sub3D(b,c);
        ca=sub3D(c,a);
        dab=norm3D(ab);
        dbc=norm3D(bc);
        dca=norm3D(ca);

        // triangle aspect ratio
        s=(dab+dbc+dca)/2.0;
        ar=dab*dbc*dca/(8*(s-dab)*(s-dbc)*(s-dca));

        // mark degenerate triangles
        if(ar>25)
            deg[i]=ar;

        // triangle altitudes
        ha=A/dbc;
        hb=A/dca;
        hc=A/dab;

        // compute the minimum triangle dimension
        min[i]=ha;
        if(min[i]>hb) min[i] = hb;
        if(min[i]>hc) min[i] = hc;
        if(min[i]>dab) min[i] = dab;
        if(min[i]>dbc) min[i] = dbc;
        if(min[i]>dca) min[i] = dca;

        // repulsion forces
        fa=sca3D(cross3D(bc,cross3D(ab,bc)),1/dbc/dbc);
        fb=sca3D(cross3D(ca,cross3D(bc,ca)),1/dca/dca);
        fc=sca3D(cross3D(ab,cross3D(ca,ab)),1/dab/dab);
        fa=sca3D(fa,1/pow(ha,3));
        fb=sca3D(fb,1/pow(hb,3));
        fc=sca3D(fc,1/pow(hc,3));

        // add forces
        f[t[i].a]=add3D(f[t[i].a], fa);
        f[t[i].b]=add3D(f[t[i].b], fb);
        f[t[i].c]=add3D(f[t[i].c], fc);
        n[t[i].a]++;
        n[t[i].b]++;
        n[t[i].c]++;
    }
    for(i=0;i<np;i++)
        f[i]=sca3D(f[i],1/(float)n[i]);

    // adapt step
    for(i=0;i<nt;i++)
    {
        if(deg[i]>0)
            continue;
        max=norm3D(f[t[i].a]);
        if(norm3D(f[t[i].b])>max) max=norm3D(f[t[i].b]);
        if(norm3D(f[t[i].c])>max) max=norm3D(f[t[i].c]);
        beta[t[i].a]=(beta[t[i].a]<alpha*min[i]/max)?(beta[t[i].a]):(alpha*min[i]/max);
        beta[t[i].b]=(beta[t[i].b]<alpha*min[i]/max)?(beta[t[i].b]):(alpha*min[i]/max);
        beta[t[i].c]=(beta[t[i].c]<alpha*min[i]/max)?(beta[t[i].c]):(alpha*min[i]/max);
    }
    for(i=0;i<np;i++)
        if(beta[i]>1)
            beta[i]=1;
    for(i=0;i<np;i++)
        f[i]=sca3D(f[i],beta[i]);

    // prevent triangles from flipping, and degenerate triangles to worsen
    for(i=0;i<nt;i++)
    {
        a=p[t[i].a];
        b=p[t[i].b];
        c=p[t[i].c];
        at=add3D(a,f[t[i].a]);
        bt=add3D(b,f[t[i].b]);
        ct=add3D(c,f[t[i].c]);
        at=sca3D(at,1/norm3D(at));
        bt=sca3D(bt,1/norm3D(bt));
        ct=sca3D(ct,1/norm3D(ct));
        ab=sub3D(at,bt);
        bc=sub3D(bt,ct);
        ca=sub3D(ct,at);
        if(deg[i]>0)
        {
            dab=norm3D(ab);
            dbc=norm3D(bc);
            dca=norm3D(ca);
            s=(dab+dbc+dca)/2.0;
            ar=dab*dbc*dca/(8.0*(s-dab)*(s-dbc)*(s-dca));
            // worsening degenerate triangle
            if(ar>deg[i])
            {
                f[t[i].a]=(float3D){0,0,0};
                f[t[i].b]=(float3D){0,0,0};
                f[t[i].c]=(float3D){0,0,0};
            }
        }
        // inverting triangle
        if(dot3D(cross3D(ab,bc),a)<0)
        {
            f[t[i].a]=(float3D){0,0,0};
            f[t[i].b]=(float3D){0,0,0};
            f[t[i].c]=(float3D){0,0,0};
        }
    }

    // apply displacement to mesh
    for(i=0;i<np;i++)
    {
        v=add3D(p[i],f[i]);
        v=sca3D(v,1/norm3D(v));
        p[i]=add3D(sca3D(p[i],(1-tau)),sca3D(v,tau));
        totDisp+=norm3D(sub3D(v,p[i]));
    }
    printf("totDisp: %g\n",totDisp);
    free(f);
    free(n);
    free(deg);
    free(min);
    free(beta);

    return 0;
}

int rotate(Mesh *m, float x, float y, float z)
{
    if(verbose)
        printf("rotate %f %f %f\n",x,y,z);
    int i;
    float   M[9];
    float3D *p=m->p,pp;
    
    x*=M_PI/180.0;
    y*=M_PI/180.0;
    z*=M_PI/180.0;
    
    M[0]=cos(z)*cos(y);
    M[1]=-sin(z)*cos(x)+cos(z)*sin(y)*sin(x);
    M[2]=sin(z)*sin(x)+cos(z)*sin(y)*cos(x);
    
    M[3]=sin(z)*cos(y);
    M[4]=cos(z)*cos(x)+sin(z)*sin(y)*sin(x);
    M[5]=-cos(z)*sin(x)+sin(z)*sin(y)*cos(x);
    
    M[6]=-sin(y);
    M[7]=cos(y)*sin(x);
    M[8]=cos(y)*cos(x);
    
    for(i=0;i<m->np;i++)
    {
        pp.x=M[0]*p[i].x+M[1]*p[i].y+M[2]*p[i].z;
        pp.y=M[3]*p[i].x+M[4]*p[i].y+M[5]*p[i].z;
        pp.z=M[6]*p[i].x+M[7]*p[i].y+M[8]*p[i].z;
        p[i]=pp;
    }
    return 0;
}
int scale(float t, Mesh *m)
{
    int     *np=&(m->np);
    float3D *p=m->p;
    int     i;
    
    for(i=0;i<*np;i++)
        p[i]=sca3D(p[i],t);
    
    return 0;
}
int scale3(float x, float y, float z, Mesh *m)
{
    int     *np=&(m->np);
    float3D *p=m->p;
    int     i;
    
    for(i=0;i<*np;i++)
    {
        p[i].x=p[i].x*x;
        p[i].y=p[i].y*y;
        p[i].z=p[i].z*z;
    }

    return 0;
}
int selectionAll(Mesh *m)
{
    int i;
    char *selection = m->selection;
    for(i=0;i<m->np;i++)
        selection[i]=1;

    return 0;
}
int selectionNone(Mesh *m)
{
    int i;
    char *selection = m->selection;
    for(i=0;i<m->np;i++)
        selection[i]=0;

    return 0;
}
int selectionInvert(Mesh *m)
{
    int i;
    char *selection = m->selection;
    for(i=0;i<m->np;i++)
    {
        if(selection[i])
            selection[i]=0;
        else
            selection[i]=1;
    }

    return 0;
}
int selectionErode(Mesh *m)
{
    int     np=m->np;
    NTriRec *NT;
    int     i,j;
    char    *selection = m->selection;
    char    *tmp = (char*)calloc(m->np,sizeof(char));

    // find neighbouring triangles for every vertex
    neighbours(m);
    NT=m->NT;
    // copy the selection into the temporary selection
    for(i=0;i<np;i++)
        tmp[i]=selection[i];

    // look for selected vertices
    for(i=0;i<np;i++)
    {
        if(selection[i])
        {
            // if it has a non selected neighbour, unselect it
            for(j=0;j<NT[i].n;j++)
                if(NT[i].t[j]!=i && selection[NT[i].t[j]]==0)
                {
                    tmp[i]=0;
                    break;
                }
        }
    }

    // copy the temporary selection into the selection
    for(i=0;i<np;i++)
        selection[i]=tmp[i];

    // free the temporary selection
    free(tmp);

    return 0;
}
int selectionDilate(Mesh *m)
{
    int     np=m->np;
    NTriRec *NT;
    int     i,j;
    char    *selection = m->selection;
    char    *tmp = (char*)calloc(m->np,sizeof(char));
    int     nprocessed;

    // find neighbouring triangles for every vertex
    neighbours(m);
    NT=m->NT;
    // copy the selection into the temporary selection
    for(i=0;i<np;i++)
        tmp[i]=selection[i];

    // look for unselected vertices
    for(i=0;i<np;i++)
    {
        if(selection[i]==0)
        {
            // if it has a selected neighbour, select it
            for(j=0;j<NT[i].n;j++)
                if(NT[i].t[j]!=i && selection[NT[i].t[j]])
                {
                    tmp[i]=1;
                    break;
                }
        }
    }

    // copy the temporary selection into the selection
    nprocessed=0;
    for(i=0;i<np;i++) {
        selection[i]=tmp[i];
        if(selection[i]>0)
            nprocessed++;
    }
    if(verbose)
        printf("%i vertices are now selected\n", nprocessed);

    // free the temporary selection
    free(tmp);

    return 0;
}
int selection(char *option, Mesh *m)
{
    if(strcmp(option,"all")==0)
        selectionAll(m);
    else
    if(strcmp(option,"invert")==0)
        selectionInvert(m);
    else
    if(strcmp(option,"none")==0)
        selectionNone(m);
    else
    if(strcmp(option,"erode")==0)
        selectionErode(m);
    else
    if(strcmp(option,"dilate")==0)
        selectionDilate(m);
    else
        printf("ERROR: select option '%s' is not known\n", option);

    return 0;
}
int selectFlipTriangle(Mesh *m)
{
    int     np=m->np;
    int     nt=m->nt;
    NTriRec *NT;
    int     i,j,pos,nflipped=0;
    float3D *nn=(float3D*)calloc(nt,sizeof(float3D));
    float3D tmp;
    char    *selection = m->selection;
    
    // compute all triangle normals
    for(i=0;i<nt;i++)
        nn[i]=normal3D(i,m);

    // find neighbouring triangles for every vertex
    neighbours(m);
    NT=m->NT;

    // find vertices with 1 inverted triangle
    for(i=0;i<np;i++)
    {
        tmp=nn[NT[i].t[0]];    // take the 1st triangle as reference
        pos=0;
        for(j=0;j<NT[i].n;j++)
            if(dot3D(nn[NT[i].t[j]],tmp)>0)
                pos++;
        if(pos==NT[i].n)    // no problem: all triangles have the same orientation
            continue;

        // flip detected
        selection[i]=1;
        nflipped++;
    }
    if(verbose)
        printf("Selected %i flipped triangles\n",nflipped);

    return 0;
}
int selectFlipTriangleSphere(Mesh *m)
{
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    int     i,nflipped;
    float3D *nn=(float3D*)calloc(nt,sizeof(float3D));
    float3D tmp;
    char    *selection = m->selection;

    // compute all triangle normals
    for(i=0;i<nt;i++)
        nn[i]=normal3D(i,m);

    // find vertices with 1 inverted triangle
    for(i=0;i<nt;i++)
    {
        tmp=sca3D(p[t[i].a], 1/norm3D(p[t[i].a]));
        if(dot3D(nn[i],tmp)<0)
        {
            selection[t[i].a]=1;
            selection[t[i].b]=1;
            selection[t[i].c]=1;
        }
    }
    free(nn);

    nflipped=0;
    for(i=0;i<np;i++)
        if(selection[i])
            nflipped++;
    if(verbose)
        printf("Selected %i vertices on flipped triangles on the sphere\n",nflipped);

    return 0;
}

int size(Mesh *m)
{
    int     np=m->np;
    float3D *p=m->p;
    float3D min,max;
    int     i;
    
    min=max=p[0];
    for(i=0;i<np;i++)
    {
        min.x=(p[i].x<min.x)?p[i].x:min.x;
        min.y=(p[i].y<min.y)?p[i].y:min.y;
        min.z=(p[i].z<min.z)?p[i].z:min.z;
        max.x=(p[i].x>max.x)?p[i].x:max.x;
        max.y=(p[i].y>max.y)?p[i].y:max.y;
        max.z=(p[i].z>max.z)?p[i].z:max.z;
    }
    printf("size: %f,%f,%f\n",fabs(max.x-min.x),fabs(max.y-min.y),fabs(max.z-min.z));

    return 0;
}
int smooth(Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp;
    int     *n;
    int     i;
    
    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    for(i=0;i<*nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
        tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    for(i=0;i<*np;i++)
        p[i]=sca3D(tmp[i],1/(float)n[i]);
    free(tmp);
    free(n);

    return 0;
}
int smoothData(Mesh *m,float l,int niter)
{
    if(verbose) printf("* smoothData lambda:%f N:%i\n",l,niter);

    int     np=m->np;
    int     nt=m->nt;
    int3D   *t=m->t;
    float   *data=m->data;
    float   *tmp;
    int     *ntmp;
    int     i,k;

    ntmp=(int*)calloc(np,sizeof(int));
    for(i=0;i<nt;i++)
    {
        ntmp[t[i].a]+=2;
        ntmp[t[i].b]+=2;
        ntmp[t[i].c]+=2;
    }
    
    tmp=(float*)calloc(np,sizeof(float));
    for(k=0;k<niter;k++)
    {
        for(i=0;i<nt;i++)
        {
            tmp[t[i].a]+=data[t[i].b]+data[t[i].c];
            tmp[t[i].b]+=data[t[i].c]+data[t[i].a];
            tmp[t[i].c]+=data[t[i].a]+data[t[i].b];
        }
        for(i=0;i<np;i++)
            if(ntmp[i]>0)
            {
                //printf("%i. %i\n",i,ntmp[i]);
                tmp[i]/=(float)ntmp[i];
                data[i]=data[i]*(1-l)+tmp[i]*l;
                tmp[i]=0;
            }
            else
                data[i]=0;
    }
    free(tmp);
    free(ntmp);

    return 0;
}
int sphereLaplace(float lambda, Mesh *m)
{
    if(verbose)
        printf("sphere laplace smooth\n");
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp,x,dx;
    float   length1, length2;
    int     *n;
    int     i;

    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    for(i=0;i<*nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
        tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    for(i=0;i<*np;i++)
    {
        if(n[i]==0)
        {
            printf("WARNING: isolated vertex %i\n",i);
            x=tmp[i];
        }
        else
            x=sca3D(tmp[i],1/(float)n[i]);
        dx=sub3D(x,p[i]);
        length1=norm3D(p[i]);
        p[i]=add3D(p[i],sca3D(dx,lambda));    // p=p+l(x-p)
        length2=norm3D(p[i]);
        p[i]=sca3D(p[i],length1/length2); // conserve the lenght
    }
    free(tmp);
    free(n);
    return 0;
}

int stereographic(Mesh *m)
{
    int     i,nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    int3D   *t1;    
    float	a,b;

    
    for(i=0;i<m->np;i++)
    {    
        a=atan2(p[i].y,p[i].x);
        b=acos(p[i].z/norm3D(p[i]));
        p[i].x=b*cos(a);
        p[i].y=b*sin(a);
        p[i].z=0;
    }
    
    // delete triangles close to the border
    t1=(int3D*)calloc(m->nt,sizeof(int3D));
    nt=0;
    for(i=0;i<m->nt;i++)
    {
        if( norm3D(sub3D(p[t[i].a],p[t[i].b]))<1.5 &&
            norm3D(sub3D(p[t[i].b],p[t[i].c]))<1.5 &&
            norm3D(sub3D(p[t[i].c],p[t[i].a]))<1.5)
            t1[nt++]=t[i];
    }
    free(m->t);
    m->t=t1;
    m->nt=nt;
    printf("new nt: %i\n",nt);
    
    return 0;
}
int subdivide(Mesh *m)
{
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *newp;
    int3D   *newt;
    int     newnp=np+nt;
    int     newnt=nt*3;
    int     i,j,k,found,*tt;
    int     p1,p2,s;
    int     i1,j1,i2,j2;
    NTriRec *T,*I;
    float3D x;
    float3D *sump;
    int     *n;
    float   beta;
    
    sump=(float3D*)calloc(np,sizeof(float3D));
    n=(int*)calloc(np,sizeof(int));

    // allocate memory for new face barycentre vertices
    newp=(float3D*)calloc(newnp,sizeof(float3D));

    // add new face vertices
    for(i=0;i<nt;i++)
    {
        x=add3D(p[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        x=sca3D(x,1/3.0);
        newp[np+i]=x;
        
        sump[t[i].a]=add3D(sump[t[i].a],x);
        sump[t[i].b]=add3D(sump[t[i].b],x);
        sump[t[i].c]=add3D(sump[t[i].c],x);
        n[t[i].a]++;
        n[t[i].b]++;
        n[t[i].c]++;
    }
    
    // update triangles
    newt=(int3D*)calloc(newnt,sizeof(int3D));
    for(i=0;i<nt;i++)
    {
        newt[i]=(int3D){t[i].a,t[i].b,np+i};
        newt[nt+i]=(int3D){t[i].b,t[i].c,np+i};
        newt[2*nt+i]=(int3D){t[i].c,t[i].a,np+i};
    }
    
    // flip old edges
    T=(NTriRec*)calloc(np,sizeof(NTriRec));
    I=(NTriRec*)calloc(np,sizeof(NTriRec));
    for(i=0;i<nt;i++)
    for(j=0;j<3;j++)
    {
        p1=((int*)&(t[i].a))[j];
        p2=((int*)&(t[i].a))[(j+1)%3];
        if(p1>p2)
        {
            s=p1;
            p1=p2;
            p2=s;
        }
        found=0;
        for(k=0;k<T[p1].n;k++) {
            tt=(int*)&(t[T[p1].t[k]]);
            // check whether the edge is already present
            // (either inverted or non inverted)
            if(tt[I[p1].t[k]]==p2 || tt[(I[p1].t[k]+1)%3]==p2)
            {
                found=1;
                break;
            }
        }
        if(found==1)
        {
            // all edge data found, that is:
            // index of triangle i1=i, index of starting vertex of the edge
            // within triangle i1, j1=j
            i1=i;
            j1=j;
            // same for triangle i2=T[p1].t[k], index j2=I[p1].t[k]
            i2=T[p1].t[k];
            j2=I[p1].t[k];
            // the edge can be processed
            // (eventually removed from T and I, to shorten future searches)
            // to each old tri t correspond 3 new tris newt[0*nt+i], newt[1*nt+i] and newt[2*nt+i]
            // the 1st tri affected is t1=newt[j1*nt+i1],
            // the 2nd is t2=newt[j2*nt+i2]
            // the first has to be changed to
            newt[j1*nt+i1].b=np+i2; // np+i2 is the new vertex added inside triangle i2
            newt[j1*nt+i1].c=np+i1;
            // the second to
            newt[j2*nt+i2].b=np+i1;
            newt[j2*nt+i2].c=np+i2;
        }
        else
        {
            // 1st reference to edge, store data
            T[p1].t[T[p1].n++]=i;
            I[p1].t[I[p1].n++]=j;
        }
    }
    free(T);
    free(I);
    
    // update position of old vertices
    for(i=0;i<np;i++)
    {
        // beta=(4-2cos(2M_PI/n))/(9n)
        beta=(4-2*cos(2*M_PI/n[i]))/(9*n[i]);
        
        // p(k+1)=(1-n*beta)p(k) + beta*sum(neighbours)
        newp[i]=add3D(sca3D(p[i],1-n[i]*beta),sca3D(sump[i],beta));
    }
    free(sump);
    free(n);
    
    m->np=newnp;
    m->nt=newnt;
    m->p=newp;
    m->t=newt;
    
    return 1;
}
int tangentLaplace(float lambda, Mesh *m)
{

    if(verbose>1)
        printf("tangent laplace smooth\n");
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float3D *tmp,dx,*tmp1,nn;
    int     *n;
    int     i;
    
    normalise(m);
    
    // compute barycentre
    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    for(i=0;i<*nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
        tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
        tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
        n[t[i].a]+=2;
        n[t[i].b]+=2;
        n[t[i].c]+=2;
    }
    for(i=0;i<*np;i++)
    {
        if(n[i]==0)
        {
            printf("WARNING: isolated vertex %i\n",i);
        }
        else
            tmp[i]=sca3D(tmp[i],1/(float)n[i]);
    }
 
    // compute normal direction as the average of neighbour triangle normals
     tmp1=(float3D*)calloc(*np,sizeof(float3D));
    for(i=0;i<*nt;i++)
    {
        nn=cross3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
        nn=sca3D(nn,1/norm3D(nn));
        tmp1[t[i].a]=add3D(tmp1[t[i].a],nn);
        tmp1[t[i].b]=add3D(tmp1[t[i].b],nn);
        tmp1[t[i].c]=add3D(tmp1[t[i].c],nn);
    }
    for(i=0;i<*np;i++)
        tmp1[i]=sca3D(tmp1[i],1/(float)n[i]);

    // apply only the smoothing component orthogonal to the normal direction
    for(i=0;i<*np;i++)
    {
        dx=sub3D(tmp[i],p[i]);  // displacement vector
        dx=sub3D(dx,sca3D(tmp1[i],dot3D(dx,tmp1[i]))); // tangential displacement vector
        p[i]=add3D(p[i],sca3D(dx,lambda));    // tangential displacement weighted by lambda
    }
    free(tmp);
    free(tmp1);
    free(n);
    return 0;   
}
int taubin(float lambda, float mu, int N, Mesh *m)
{
    int j;
    
    if(verbose)
        printf("* taubinSmooth %f %f %i\n",lambda,mu,N);
        
    for(j=0;j<2*N;j++)
        if(j%2==0)
            laplace(lambda,m);
        else
            laplace(mu,m);
    return 0;
}
void threshold(float thr, int direction, Mesh *m)
{
    int     *np=&(m->np);
    float   *data=m->data;
    int     i;

    if(direction==0)
    {
        for(i=0;i<*np;i++)
            if(data[i]<=thr)
                data[i]=1;
            else
                data[i]=-1;
    }
    else
    {
        for(i=0;i<*np;i++)
            if(data[i]>=thr)
                data[i]=1;
            else
                data[i]=-1;
    }
}
// Code By S Melax from http://www.melax.com/volint/
float volume(Mesh *m)
{
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    float   vol=0;
    int     i;
    
    for(i=0;i<*nt;i++)
        vol += determinant(p[t[i].a],p[t[i].b],p[t[i].c]); //divide by 6 later for efficiency
    vol/=6.0;// since the determinant give 6 times tetra volume
    printf("volume: %f\n",vol);
    return vol;
}
int normal(Mesh *m)
{
    // IMPORTANT NOTE: everywhere C is just np scalars, but here is np 3d vectors.
    // There is no way, for the moment, of saving this normal.
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    int3D   *t=m->t;
    float   *data=m->data;
    float3D *tmp;
    int     *n;
    int     i;
    
    printf("WARNING: \"normal\" has not been tested\n");

    tmp=(float3D*)calloc(*np,sizeof(float3D));
    n=(int*)calloc(*np,sizeof(int));
    // All normals are weighted the same, but more
    // correctly, the values should be weighted by the angle
    // at the given vertex.
    for(i=0;i<*nt;i++)
    {
        tmp[t[i].a]=add3D(tmp[t[i].a],normal3D(i,m));
        tmp[t[i].b]=add3D(tmp[t[i].b],normal3D(i,m));
        tmp[t[i].c]=add3D(tmp[t[i].c],normal3D(i,m));
        n[t[i].a]++;
        n[t[i].b]++;
        n[t[i].c]++;
    }
    for(i=0;i<*np;i++)
        ((float3D*)data)[i]=sca3D(tmp[i],1/(float)n[i]);
    free(tmp);
    
    return 0;
}
int subVal(float val,Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    for(i=0;i<np;i++)
        data[i]-=val;    
    return 0;
}
int addVal(float val,Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    for(i=0;i<np;i++)
        data[i]+=val;
    return 0;
}
int multVal(float val,Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    for(i=0;i<np;i++)
        data[i]*=val;
    return 0;
}
int divVal(float val,Mesh *m)
{
    int     np=m->np;
    float   *data=m->data;
    int     i;
    
    if(data==NULL)
    {
        printf("ERROR: there is no data\n");
        return 1;
    }

    for(i=0;i<np;i++)
        data[i]/=val;
    return 0;
}
void randverts(int nrv, Mesh *m)
{
    // This function is meant to be used with qhull, like this:
    // meshgeometry -imesh $mesh -centre -normalise -randverts 100 | qhull o
    int     *nt=&(m->nt);
    float3D *p=m->p;
    int3D   *t=m->t;
    int     i,j;
    float   s0,s1;
    float3D v;
                
    srand(time(NULL)+getpid());
    
    printf("3 meshgeometry randverts\n");
    printf("%i\n",nrv);
    // generate random vertices
    for(i=0;i<nrv;i++)
    {
        j=(*nt-0.01)*(rand()/(float)RAND_MAX);    // pick a random triangle

        // pick a random point within the triangle
        s0=rand()/(float)RAND_MAX;
        s1=rand()/(float)RAND_MAX;
        
        // only values such that s0+s1<=0 are in the triangle, if
        // that is not the case, reflect the point on the s1=1-s0 axis
        if(s0+s1>1)
        {
            s0=(1-s1);
            s1=(1-s0);
        }
        v=add3D(p[t[j].a],add3D(sca3D(sub3D(p[t[j].b],p[t[j].a]),s0),sca3D(sub3D(p[t[j].c],p[t[j].a]),s1)));

        printf("%f %f %f\n",v.x,v.y,v.z);
    }
}
int resample(char *path_m1, char *path_rm, Mesh *m)
{
    /*
    m:          original mesh
    path_m1:    path to spherical version of m
    path_rm:    path to the spherical version of the target mesh (reference mesh)

    The function produces a mesh with the same geometry of m, but resampled with
    the topology of the reference mesh rm.
    rm is a spherical mesh, normally associated to a native mesh.
    The landmarks over m1 (spherical version of m) and rm (spherical version of the
    reference mesh) should have been previously made to correspond.
    */
    Mesh    m1; // spherical version of mesh m
    Mesh    rm; // spherical version of the target mesh 
    int     nt;
    int3D   *t; // original mesh topology
    int     np_rm;
    float3D *p;     // original mesh points
    float3D *p_m1;  // spherical version of the original mesh
    float3D *p_rm;  // spherical target (reference) mesh
    float3D *tmp;
    float   c0,c1;
    int     i,j,k,imindist,result;
    float3D n;
    float   flipTest,mindist;
    int     non_mapped=0;
    int case_deg,case_parl,case_disj;
    case_deg=case_parl=case_disj=0;
    
    loadMesh(path_m1,&m1,0);    // load spherical version of the original mesh
    loadMesh(path_rm,&rm,0);    // load spherical target mesh
    
    // Check that m and m1 have the same number of vertices
    if(m->np!=m1.np)
    {
        printf("ERROR: m1 does not have the same number of vertices as m\n");
        return 1;
    }
    
    // Centre spherical meshes m1 and rm
    centre(&m1);
    centre(&rm);
    
    // Check whether the meshes are properly oriented
    n=normal3D(0,&m1);
    flipTest=dot3D(m1.p[m1.t[0].a],n);
    if(flipTest<0)
    {
        printf("ERROR: m1 is mis-oriented. Path: %s\n",path_m1);
        return 1;
    }
    n=normal3D(0,&rm);
    flipTest=dot3D(rm.p[rm.t[0].a],n);
    if(flipTest<0)
    {
        printf("ERROR: rm is mis-oriented. Path: %s\n",path_rm);
        return 1;
    }
    
    // Points and triangles of the original mesh
    p=m->p;
    nt=m->nt;
    t=m->t;

    // Points of the spherical version of the original mesh
    p_m1=m1.p;

    // Points of the target spherical mesh (reference)
    np_rm=rm.np;
    p_rm=rm.p;
    
    // Interpolate coordinates on reference mesh
    tmp=(float3D*)calloc(np_rm,sizeof(float3D));    // the new points are stored in tmp
    for(i=0;i<np_rm;i++)
    {
        // display progress
        if(verbose)
        if((int)(i*100/np_rm)>(int)((i-1)*100/np_rm))
        {
            printf("%i%% ",i*100/np_rm);
            fflush(stdout);
        }
        
        // look for a triangle in the spherical original mesh containing
        // point i of the target mesh
        for(j=0;j<nt;j++)
        {
            c0=c1=0;
            result=intersect_VectorTriangle(p_rm[i],j,&c0,&c1,&m1);

            if(result==1)
            {
                tmp[i]=sca3D(p[t[j].a],1-c0-c1);
                tmp[i]=add3D(tmp[i],sca3D(p[t[j].b],c0));
                tmp[i]=add3D(tmp[i],sca3D(p[t[j].c],c1));
                break;
            }
            
            if(result==-1) case_deg++;
            if(result==0) case_disj++;
            if(result==2) case_parl++;
        }
        
        // if j==nt no triangle was found in the spherical original mesh that contained
        // point i. As an approximation, pick the closest point in the spherical original
        // mesh.
        if(j==nt)
        {
            //printVertex(p_rm[i]);
            mindist=norm3D(sub3D(p_rm[i],p_m1[0]));
            for(k=0;k<m->np;k++)
                if(norm3D(sub3D(p_rm[i],p_m1[k]))<mindist)
                {
                    mindist=norm3D(sub3D(p_rm[i],p_m1[k]));
                    imindist=k;
                }
            tmp[i]=p[imindist];
            non_mapped++;
            //printf("\nWARNING: could not resample point %i, mapped it to closest vertex, %i\n",i,imindist);
            //printVertex(p_m1[imindist]);
        }
    }
    if(verbose)
        printf("\n");
    if(non_mapped) printf("\nWARNING: %i vertices could not be resampled and were mapped to the closest vertex. Reference mesh: %s\n",non_mapped,path_rm);
    if(verbose) printf("degenerate: %i\ndisjoint: %i\nparallel: %i\n",case_deg,case_disj,case_parl);
    
    // Free data in original mesh (m), i.e., vertices, triangles, etc.
    free(m->p);
    free(m->t);
    if(m->data)
        free(m->data);
    if(m->NT)
        free(m->NT);
    
    // Free spherical version of the actual mesh (m1)
    free(m1.p);
    free(m1.t);
    
    // Reconfigure original mesh with resampled data
    m->p=tmp;
    m->t=rm.t;
    m->np=np_rm;
    m->nt=rm.nt;
    
    if(verbose)
        printf("resampled np: %i, nt: %i\n",m->np,m->nt);
    
    return 0;
}

/*
    re-index the mesh's triangles so that they are sorted along the x axis
*/
int sortTrianglesFunction(const void *a, const void *b)
{
    float3D	v1=*(float3D*)a;
    float3D	v2=*(float3D*)b;

    if(v1.x==v2.x)
    {
        return 0;
    }
    else
    {
        if(v1.x<v2.x)
            return -1;
        else
            return	1;
    }
}
void sortTriangles(Mesh *m)
{
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t,*t2;
    float3D   *x;
    int     i;
    
    x=(float3D*)calloc(nt,sizeof(float3D));
    for(i=0;i<nt;i++)
    {
        x[i].x=(p[t[i].a].x+p[t[i].b].x+p[t[i].c].x)/3.0;
        x[i].y=i;
    }
    
    // sort vertices per cos(angle)
    qsort(x,nt,sizeof(float3D),sortTrianglesFunction);
    
    t2=(int3D*)calloc(nt,sizeof(int3D));

    // redistribute uniformly the triangles along the x axis
    for(i=0;i<nt;i++)
        t2[i]=t[(int)round(x[i].y)];
    for(i=0;i<nt;i++)
        t[i]=t2[i];

    free(x);
    free(t2);
}
int sortVerticesFunction(const void *a, const void *b)
{
    float3D v1=*(float3D*)a;
    float3D v2=*(float3D*)b;

    if(v1.x==v2.x)
    {
        return 0;
    }
    else
    {
        if(v1.x<v2.x)
            return -1;
        else
            return	1;
    }
}
void uniform(Mesh *m)
{
    float3D *p=m->p,p0;
    int     np=m->np;
    float3D *a,d,v,ma,mi,c;
    int     niter=2000,maxiter;
    int     i,j;
    float   x,s,ss,std,maxstd,r1,r2,n;
    
    mi=ma=p[0];
    for(i=0;i<np;i++)
    {
        mi.x=(mi.x>p[i].x)?p[i].x:mi.x;
        mi.y=(mi.y>p[i].y)?p[i].y:mi.y;
        mi.z=(mi.z>p[i].z)?p[i].z:mi.z;
        ma.x=(ma.x<p[i].x)?p[i].x:ma.x;
        ma.y=(ma.y<p[i].y)?p[i].y:ma.y;
        ma.z=(ma.z<p[i].z)?p[i].z:ma.z;
    }
    c=(float3D){(mi.x+ma.x)/2.0,(mi.y+ma.y)/2.0,(mi.z+ma.z)/2.0};
    
    srand(0);
    for(j=0;j<niter;j++)
    {
        d.x=rand()/(float)RAND_MAX;
        d.y=rand()/(float)RAND_MAX;
        d.z=rand()/(float)RAND_MAX;
        d=sca3D(d,1/norm3D(d));
        s=ss=0;
        for(i=0;i<np;i++)
        {
            p0=sub3D(p[i],c);
            x=dot3D(p0,d);
            s+=x;
            ss+=x*x;
        }
        std=sqrt(fabs(ss/(float)np-pow(s/(float)np,2)));
        //printf("%f\n",std);
        if(j==0)
        {
            maxstd=std;
            maxiter=j;
        }
        else if(maxstd<std)
        {
            maxstd=std;
            maxiter=j;
        }   
    }

    srand(0);
    for(j=0;j<=maxiter;j++)
    {
        d.x=rand()/(float)RAND_MAX;
        d.y=rand()/(float)RAND_MAX;
        d.z=rand()/(float)RAND_MAX;
        d=sca3D(d,1/norm3D(d));
    }
    printf("d=(%g,%g,%g)\n",d.x,d.y,d.z);
    s=ss=0;
    a=(float3D*)calloc(np,sizeof(float3D));
    for(i=0;i<np;i++)
    {
        p0=sub3D(p[i],c);
        x=dot3D(sca3D(p0,1/norm3D(p0)),d);
        a[i].x=x;
        a[i].y=i;
        s+=x;
        ss+=x*x;
    }
    printf("angle mean=%g, angle s.d.=%g\n",s/(float)np,std);
    
    // sort vertices per cos(angle)
    qsort(a,np,sizeof(float3D),sortVerticesFunction);

    // redistribute uniformly the vertices along the axis d
    for(i=0;i<np;i++)
    {
        j=(int)(a[i].y+0.1);
        p0=sub3D(p[j],c);
        n=norm3D(p0);
        r1=dot3D(p0,d);
        v=sub3D(p0,sca3D(d,r1));
        v=sca3D(v,1/norm3D(v));
        r2=-(0.999-2*0.999*i/(float)np);
        p0=add3D(sca3D(d,r2),sca3D(v,sqrt(1-r2*r2)));
        p[j]=add3D(sca3D(p0,n/norm3D(p0)),c);
    }
    free(a);
}
void printHelp(void)
{
     printf("\
 Commands\n\
\n\
   Input/output\n\
    -iformat format_name                             Force input format (needs to precede imesh)\n\
    -oformat format_name                             Force output format (needs to precede omesh)\n\
    -i filename                                      Input file (also accepts -imesh)\n\
    -o filename                                      Output file (also accepts -omesh)\n\
    -odata filename                                  Output data\n\
\n\
   Mesh measurements\n\
    -absgi                                           Print absolute gyrification index\n\
    -area                                            Print surface area\n\
    -boundingBox                                     Print mesh bounding box minx, miny, minz, maxx, maxy, maxz.\n\
    -checkOrientation                                Check that normals point outside\n\
    -edgeLength                                      Print average edge length and standard deviation\n\
    -edgeLengthMinMax                                Print min and max edge length\n\
    -euler                                           Print Euler characteristic\n\
    -foldLength                                      Print total fold length\n\
    -printCentre                                     Print coordinates of the centre of the mesh\n\
    -printBarycentre                                 Print coordinates of the barycentre of the mesh\n\
    -printBarycentre                                 Print coordinates of the average of all mesh vertices\n\
    -printCentre                                     Print coordinates of the point at half width, length and height of the mesh\n\
    -size                                            Print mesh dimensions\n\
    -tris                                            Print number of triangles\n\
    -verts                                           Print number of vertices\n\
    -volume                                          Print mesh volume\n\
    -nonmanifold                                     Print number of vertices in edges with more than 2 triangles\n\
    -isolatedVerts                                   Print number of isolated vertices\n\
\n\
   Mesh modification\n\
    -add filename                                    Add mesh at filename to the current mesh\n\
    -align filename                                  Align mesh to the mesh pointed by filename, which has the same topology but different geometry\n\
    -applyMatrix m11,m12,...,m34                     Transform the vertex coordinates by multiplying them by a 4x4 matrix (the last row is set to 0,0,0,1)\n\
    -average n_meshes path1 path2 ... pathn          Compute an average of n_meshes all\n\
                                                       of the same topology\n\
    -barycentre                                      Put the mesh origin at the average of all its vertices\n\
    -barycentricProjection reference_mesh            Print barycentric coordinates for each vertex in reference_mesh\n\
    -centre                                          Put the mesh origin at half its width, length and height\n\
    -flip                                            Flip normals\n\
                                                       degrees and fix them\n\
    -fixFlip                                         Detect flipped triangles and fix them\n\
    -fixFlipSphere                                   Detect flipped triangles and fix them, for spheres\n\
    -fixSmall                                        Detect triangles with an angle >160\n\
    -invert axis                                     Invert the sign of the mesh vertices along the specified axis\n\
    -laplaceSmooth lambda num_iter                   Laplace smoothing\n\
    -laplaceSmoothSelection lambda num_iter          Laplace smoothing only of the selected vertices\n\
    -lissencephalic                                  Smooth valleys and hills, not the coast\n\
    -level level_value                               Adds new vertices (and triangles) to the\n\
                                                       edges that cross level_value in the\n\
                                                       vertex data (f.ex., mean curvature)\n\
    -mirror coord                                    Mirror vertices in the coordinate 'coord' relative to the mesh's barycentre\n\
    -normalise                                       Place all vertices at distance 100 from\n\
                                                       the origin\n\
    -randverts number_of_vertices                    Generate homogeneously distributed\n\
                                                       random vertices over the mesh\n\
    -relax filename                                  Relax current mesh to mesh at filename\n\
                                                        (both meshes have the same topology)\n\
    -removeIsolatedVerts                             Remove isolated vertices\n\
    -removeVerts                                     Remove vertices with negative vertex data values\n\
    -repulse num_iter                                Repulse vertices with a force f=1/dist for num_iter iterations\n\
    -resample smooth_mesh reference_mesh             Resample the mesh to match the vertices\n\
                                                       and the topology of the argument mesh\n\
    -rotate x y z                                    Rotate with angles x, y and z in degrees\n\
    -scale scale_value(s)                            Multiply each vertex by \"scale\", where scale can be one number, or 3 numbers separated by commas: x,y,z\n\
    -sphereLaplaceSmooth lambda num_iter             Laplace smoothing that conserves the norms of the mesh vertices, i.e., spheres stay spheres\n\
    -sortTriangles                                   Sort the triangles in the mesh file along the x axis\n\
    -stereographic                                   Stereographic projection\n\
    -subdivide                                       Subdivide the mesh using 1 iteration of Kobbelt's sqrt(3) algorithm\n\
    -tangentLaplace lambda number_of_iterations      Laplace smoothing tangential to the mesh surface\n\
    -taubinSmooth lambda mu number_of_iterations     Taubin Smoothing\n\
    -translate x y z                                 Translatory motion x,y,z\n\
\n\
   Vertex value modification\n\
    -addVal                                          Add value data\n\
    -subVal                                          Subtract value from data\n\
    -multVal                                         Multiply data time value\n\
    -divVal                                          Divide data by value\n\
    -divVal                                          Divide data by value\n\
    -clip min max                                    Clip data values to the interval [min,max]\n\
    -areaMap                                         Compute surface area per vertex\n\
    -countClusters  value                            Count clusters in texture data\n\
    -curv                                            Compute curvature\n\
    -depth                                           Compute sulcal depth\n\
    -icurv number_of_iterations                      Integrated curvature\n\
    -smoothData lambda number_of_iterations          Laplace smoothing of data, lambda=0 -> no smoothing, lambda=1 -> each vertex value to neighbour's average\n\
    -threshold value 0:down/1:up                     Threshold texture data\n\
    -normal                                          Store mesh normal vectors as vertex values\n\
\n\
   Vertex value measurements\n\
    -mean                                            Print mean data value\n\
    -min                                             Print minimum data value\n\
    -max                                             Print maximum data value\n\
\n\
   Selections\n\
    -selection option                                Select vertices. Options are: none, all, invert, erode, dilate\n\
    -selectFlip                                      Select flipped triangles\n\
    -selectFlipSphere                                Select flipped triangles, for spherical meshes\n\
\n\
   General\n\
    -drawSurface colourmap path                      draw surface in tiff format, colourmap is grey, rainbow, level2 or level4\n\
    -drawSurfaceToon colourmap path                  draw surface using toon rendering in tiff format, colourmap is grey, rainbow, level2 or level4\n\
    -h                                               Help\n\
    -v                                               Verbose mode\n\
\n\
   File formats\n\
    Meshgeometry can read and write several formats, based on the file extension:\n\
    .orig, .pial, .white, .inflated, .sphere, .reg   Freesurfer meshes\n\
    .curv, .sulc, .sratio                            Freesurfer data\n\
    .mesh                                            BrainVisa meshe\n\
    .txt                                             RT's mesh plain text format\n\
    .float                                           Raw float data\n\
    .txt1                                            RT's data format\n\
    .bin                                             n-e-r-v-o-u-s system web binary mesh\n\
    .obj                                             Civet's .obj format. Has to be used with -iformat civet_obj\n\
    .asc                                             Freesurfer's ascii format\n\
    .gii                                             Gifti format\n\
    .wrl, .obj, .ply, .stl, .smesh, .off, .vtk       Other mesh formats\n\
");
}
int main(int argc, char *argv[])
{
    checkEndianness();
    
    int    i;
    int    iformat=0;
    int    oformat=0;
    
    mesh.p=NULL;
    mesh.t=NULL;
    mesh.data=NULL;
    mesh.selection=NULL;
    mesh.ddim=1;
    mesh.NT=NULL;
    
    i=1;
    while(i<argc)
    {
        if(strcmp(argv[i],"-iformat")==0)
        {
            char    str[256];
            sprintf(str," .%s",argv[++i]);
            printf("iformat: %s\n",str);
            iformat=getformatindex(str);
        }
        else
        if(strcmp(argv[i],"-oformat")==0)
        {
            char    str[256];
            sprintf(str," .%s",argv[++i]);
            oformat=getformatindex(str);
        }
        else
        if(strcmp(argv[i],"-i")==0)
        {
            loadMesh(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-o")==0)
        {
            saveMesh(argv[++i],&mesh,oformat);
        }
        else
        if(strcmp(argv[i],"-imesh")==0)
        {
            printf("WARNING: -imesh still works, but better change to -i\n");
            loadMesh(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-omesh")==0)
        {
            printf("WARNING: -omesh still works, but better change to -o\n");
            saveMesh(argv[++i],&mesh,oformat);
        }
        else
        if(strcmp(argv[i],"-max")==0)
        {
            maxData(&mesh);
        }
        else
        if(strcmp(argv[i],"-mean")==0)
        {
            meanData(&mesh);
        }
        else
        if(strcmp(argv[i],"-min")==0)
        {
            minData(&mesh);
        }
        else
        if(strcmp(argv[i],"-odata")==0)
        {
            Text_save_data(argv[++i],&mesh);
        }
        else
        if(strcmp(argv[i],"-addVal")==0)
        {
            addVal(atof(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-subVal")==0)
        {
            subVal(atof(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-multVal")==0)
        {
            multVal(atof(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-divVal")==0)
        {
            divVal(atof(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-clip")==0)
        {
            float min=atof(argv[++i]);
            float max=atof(argv[++i]);
            clip(&mesh,min,max);
        }
        else
        if(strcmp(argv[i],"-add")==0)
        {
            addMesh(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-applyMatrix")==0)
        {
            char *str=argv[++i];
            float m[16]={0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1};
            sscanf(str," %f,%f,%f,%f, %f,%f,%f,%f, %f,%f,%f,%f ",
                &(m[0]),&(m[1]),&(m[2]),&(m[3]),
                &(m[4]),&(m[5]),&(m[6]),&(m[7]),
                &(m[8]),&(m[9]),&(m[10]),&(m[11]));
            applyMatrix(m,&mesh);
        }
        else
        if(strcmp(argv[i],"-align")==0)
        {
            align(&mesh,argv[++i]);
        }
        else
        if(strcmp(argv[i],"-areaMap")==0)
        {
            if(mesh.data==NULL)
                mesh.data=(float*)calloc(mesh.np,sizeof(float));
            areaMap(mesh.data,&mesh);
        }
        else
        if(strcmp(argv[i],"-relax")==0)
        {
            relax(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-barycentricProjection")==0)    // print barycentric coordinates for each vertex of reference mesh
        {
            barycentricProjection(argv[i+1],&mesh);
            i+=1;
        }
        else
        if(strcmp(argv[i],"-boundingBox")==0)    // mesh bounding box
        {
            boundingBox(&mesh);
        }
        else
        if(strcmp(argv[i],"-checkOrientation")==0)
        {
            checkOrientation(&mesh);
        }
        else
        if(strcmp(argv[i],"-curv")==0)
        {
            if(mesh.data==NULL)
                mesh.data=(float*)calloc(mesh.np,sizeof(float));
            curvature(mesh.data,&mesh);
        }
        else
        if(strcmp(argv[i],"-depth")==0)
        {
            if(mesh.data==NULL)
                mesh.data=(float*)calloc(mesh.np,sizeof(float));
            depth(mesh.data,&mesh);
        }
        else
        if(strcmp(argv[i],"-edgeLength")==0)    // print average edge length and std
        {
            edgeLength(&mesh);
        }
        else
        if(strcmp(argv[i],"-edgeLengthMinMax")==0)    // print min and max edge length
        {
            edgeLengthMinMax(&mesh);
        }
        else
        if(strcmp(argv[i],"-icurv")==0)
        {
            if(mesh.data==NULL)
                mesh.data=(float*)calloc(mesh.np,sizeof(float));
            icurvature(atoi(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-invert")==0)
        {
            char *axis = argv[++i];
            invert(axis,&mesh);
        }
        else
        if(strcmp(argv[i],"-isolatedVerts")==0)
        {
            isolatedVerts(&mesh);
        }
        else
        if(strcmp(argv[i],"-removeIsolatedVerts")==0)
        {
            removeIsolatedVerts(&mesh);
        }
        else
        if(strcmp(argv[i],"-removeVerts")==0)
        {
            removeVerts(&mesh);
        }
        else
        if(strcmp(argv[i],"-laplaceSmooth")==0)
        {
            int     j,N;
            float   l;
            
            l=atof(argv[++i]);
            N=atoi(argv[++i]);
            for(j=0;j<N;j++)
                laplace(l,&mesh);
        }
        else
        if(strcmp(argv[i],"-repulse")==0)
        {
            int j,niter;
            niter=atoi(argv[++i]);
            for(j=0;j<niter;j++)
                repulse(&mesh);
        }
        else
        if(strcmp(argv[i],"-sphereLaplaceSmooth")==0)
        {
            int     j,N;
            float   l;
            
            l=atof(argv[++i]);
            N=atoi(argv[++i]);
            for(j=0;j<N;j++)
                sphereLaplace(l,&mesh);
        }
        else
        if(strcmp(argv[i],"-laplaceSmoothSelection")==0)
        {
            int     j,N;
            float   l;
            
            l=atof(argv[++i]);
            N=atoi(argv[++i]);
            for(j=0;j<N;j++)
                laplaceSelection(l,&mesh);
        }
        else
        if(strcmp(argv[i],"-nonmanifold")==0)
        {
            nonmanifold_verts(&mesh);
            nonmanifold_eds(&mesh);
            nonmanifold_tris(&mesh);
        }
        else
        if(strcmp(argv[i],"-taubinSmooth")==0)
        {
            int     N;
            float   lambda,mu;
            
            lambda=atof(argv[++i]);
            mu=atof(argv[++i]);
            N=atoi(argv[++i]);
            taubin(lambda,mu,N,&mesh);
        }
        else
        if(strcmp(argv[i],"-smoothData")==0)
        {
            int     N;
            float   l;
            
            l=atof(argv[++i]);
            N=atoi(argv[++i]);
            smoothData(&mesh,l,N);
        }
        else
        if(strcmp(argv[i],"-printCentre")==0)
        {
            printCentre(&mesh);
        }
        else
        if(strcmp(argv[i],"-printBarycentre")==0)
        {
            printBarycentre(&mesh);
        }
        else
        if(strcmp(argv[i],"-euler")==0)
        {
            printf("euler: %i\n",mesh.np-mesh.nt/2);
        }
        else
        if(strcmp(argv[i],"-average")==0)
        {
            int N=atoi(argv[i+1]);
            average(N,&(argv[i+2]),&mesh);
            i+=N+1;
        }
        else
        if(strcmp(argv[i],"-flip")==0)
        {
            flip(&mesh);
        }
        else
        if(strcmp(argv[i],"-fixFlip")==0)
        {
            fixflip(&mesh);
        }
        else
        if(strcmp(argv[i],"-fixFlipSphere")==0)
        {
            fixflipSphere(&mesh);
        }
        else
        if(strcmp(argv[i],"-fixNonmanifold")==0)
        {
            int i,nm=nonmanifold_verts(&mesh);
            for(i=0;i<nm;i++)
                fixNonmanifold_verts(&mesh);
        }
        else
        if(strcmp(argv[i],"-fixSmall")==0)
        {
            fixSmall(&mesh);
        }
        else
        if(strcmp(argv[i],"-scale")==0)
        {
            char *vals = argv[++i];
            float x, y, z;
            int n;
            n = sscanf(vals, " %f,%f,%f ", &x, &y, &z);
            if(n==1)
            {
                scale(x,&mesh);
            }
            else if(n==3)
            {
                scale3(x,y,z,&mesh);
            }
            else
            {
                printf("ERROR: wrong arguments in -scale switch\n");
            }
        }
        else
        if(strcmp(argv[i],"-selection")==0)
        {
            char *option = argv[++i];
            selection(option, &mesh);
        }
        else
        if(strcmp(argv[i],"-selectFlip")==0)
        {
            selectFlipTriangle(&mesh);
        }
        else
        if(strcmp(argv[i],"-selectFlipSphere")==0)
        {
            selectFlipTriangleSphere(&mesh);
        }
        else
        if(strcmp(argv[i],"-sortTriangles")==0)
        {
            sortTriangles(&mesh);
        }
        else
        if(strcmp(argv[i],"-lissencephalic")==0)
        {
            lissencephalic(atof(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-level")==0)
        {
            level(atof(argv[++i]),&mesh);
        }
        else
        if(strcmp(argv[i],"-area")==0)
        {
            area(&mesh);
        }
        else
        if(strcmp(argv[i],"-volume")==0)
        {
            volume(&mesh);
        }
        else
        if(strcmp(argv[i],"-absgi")==0)
        {
            absgi(&mesh);
        }
        else
        if(strcmp(argv[i],"-tangentLaplace")==0)
        {
            int     j,N;
            float   l;
            
            l=atof(argv[++i]);
            N=atoi(argv[++i]);
            for(j=0;j<N;j++)
                tangentLaplace(l,&mesh);
        }
        else
        if(strcmp(argv[i],"-threshold")==0)    // threshold value 0=down/1=up
        {
            float   value=atof(argv[++i]);
            int     direction=atoi(argv[++i]);
            threshold(value,direction,&mesh);
        }
        else
        if(strcmp(argv[i],"-foldLength")==0)
        {
            foldLength(&mesh);
        }
        else
        if(strcmp(argv[i],"-barycentre")==0)
        {
            barycentre(&mesh);
        }
        else
        if(strcmp(argv[i],"-translate")==0)    // translate x, y, z
        {
            translate(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]),&mesh);
            i+=3;
        }
        else
        if(strcmp(argv[i],"-centre")==0)
        {
            centre(&mesh);
        }
        else
        if(strcmp(argv[i],"-drawSurface")==0)
        {
            char   *cmap=argv[++i];
            char   *tiff_path=argv[++i];
            drawSurface(&mesh,cmap,tiff_path,0);
        }
        else
        if(strcmp(argv[i],"-drawSurfaceToon")==0)
        {
            char   *cmap=argv[++i];
            char   *tiff_path=argv[++i];
            drawSurface(&mesh,cmap,tiff_path,1);
        }
        else
        if(strcmp(argv[i],"-mirror")==0)
        {
            char   *coord=argv[++i];
            mirror(&mesh,coord);
        }
        else
        if(strcmp(argv[i],"-normalise")==0)
        {
            normalise(&mesh);
        }
        else
        if(strcmp(argv[i],"-countClusters")==0)
        {
            float    thr=atof(argv[++i]);
            if(mesh.data==NULL)
                printf("ERROR: no texture data available\n");
            countClusters(thr,&mesh);
        }
        else
        if(strcmp(argv[i],"-normal")==0)    // surface normal vectors
        {
            if(mesh.data==NULL)
                mesh.data=(float*)calloc(3*mesh.np,sizeof(float));
            normal(&mesh);
        }
        else
        if(strcmp(argv[i],"-randverts")==0)    // randverts n
        {
            int    n=atoi(argv[++i]);
            randverts(n,&mesh);
        }
        else
        if(strcmp(argv[i],"-resample")==0)    // resample actual mesh using smooth mesh as reference mesh
        {
            resample(argv[i+1],argv[i+2],&mesh);
            i+=2;
        }
        else
        if(strcmp(argv[i],"-rotate")==0)    // rotate with angles x, y and z
        {
            rotate(&mesh,atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]));
            i+=3;
        }
        else
        if(strcmp(argv[i],"-size")==0)    // mesh dimensions
        {
            size(&mesh);
        }
        else
        if(strcmp(argv[i],"-stereographic")==0)    // stereographic projection
        {
            stereographic(&mesh);
        }
        else
        if(strcmp(argv[i],"-subdivide")==0)    // stereographic projection
        {
            subdivide(&mesh);
        }
        else
        if(strcmp(argv[i],"-uniform")==0)    // make distribution of vertices uniform along the axis of maximum variability (for spherical meshes)
        {
            uniform(&mesh);
        }
        else
        if(strcmp(argv[i],"-verts")==0)    // display number of vertices
        {
            printf("verts: %i\n",mesh.np);
        }
        else
        if(strcmp(argv[i],"-tris")==0)    // display number of triangles
        {
            printf("tris: %i\n",mesh.nt);
        }
        else
        if(strcmp(argv[i],"-h")==0)    // help
        {
            printHelp();
            return 0;
        }
        else
        if(strcmp(argv[i],"-v")==0)    // turn on verbose mode
        {
            verbose+=1;
        
            // print some information
            printf("%s\n",version);
            printf("CPU: %s\n",(endianness==kMOTOROLA)?"Motorola":"Intel");
        }
        else
        {
            printf("ERROR: Unknown argument '%s'.\n",argv[i]);
            return 1;
        }
        i++;
    }
    
    freeMesh(&mesh);

    if(verbose)
        printf("Done.\n");
    return 0;
}
