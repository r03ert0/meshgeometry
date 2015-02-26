//char version[]="meshgeometry, version 1, roberto toro, 12 August 2009";
//char version[]="meshgeometry, version 2, roberto toro, 28 April 2010";
//char version[]="meshgeometry, version 3, roberto toro, 1 May 2010";     // added -laplace and -taubinLM
//char version[]="meshgeometry, version 4, roberto toro, 21 May 2010";    // commands are processed as a chain
//char version[]="meshgeometry, version 5, roberto toro, 28 May 2012";    // added several commands: foldLength, volume, absgi, texture threshold, countClusters, and includes meshconvert v8
//char version[]="meshgeometry, version 6, roberto toro, 10 November 2012"; // added randomverts, help, centre, normalise, normal, verbose, off mesh format (load and save), added to github
//char version[]="meshgeometry, version 7, roberto toro, 17 Decembre 2014"; // vtk support
char version[]="meshgeometry, version 8, roberto toro, 26 Decembre 2015";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define pi 3.14159265358979323846264338327950288419716939927510
#define EPSILON  0.000001 // small enough to avoid division overflow
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
    int     np;     // number of vertices
    int     nt;     // number of triangles
    int     ddim;   // data dimensions (default: 1)
    float3D *p;     // vertices
    int3D   *t;     // triangles
    float   *data;  // data
    NTriRec *NT;    // neighbouring triangles
}Mesh;

float area(Mesh *m);
float volume(Mesh *m);
int smooth(Mesh *m);
int taubin(float lambda, float mu, int N, Mesh *m);

Mesh    mesh;
float   R;
int     verbose=0;

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
int intersect_VectorTriangle(float3D x, int i, float *c0, float *c1, Mesh *m)
{
    float3D *p=m->p;
    int3D   *t=m->t;
    int3D   T=t[i];
    float3D xx;
    float3D u, v, n;             // triangle vectors
    float3D dir,w0, w;               // ray vectors
    float   r, a, b;             // params to calc ray-plane intersect

    u=sub3D(p[T.b],p[T.a]);
    v=sub3D(p[T.c],p[T.a]);
    n=cross3D(u,v);
    if(norm3D(n)<=EPSILON)         // triangle is degenerate, do not deal with this case
        return -1;

    dir=x;
    w0 = sca3D(p[T.a],-1);
    a = dot3D(n,w0);
    b = dot3D(n,dir);
    if (fabs(b) < EPSILON) {        // ray is parallel to triangle plane
        if (a == 0)                 // ray lies in triangle plane
            return 2;
        else
            return 0;              // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = -a/b;
    if (r < 0.0)                    // ray goes away from triangle
        return 0;                   // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    xx = sca3D(dir,r);    // intersect point of ray and plane

    // is I inside T?
    float    uu, uv, vv, wu, wv, D;
    uu = dot3D(u,u);
    uv = dot3D(u,v);
    vv = dot3D(v,v);
    w = sub3D(xx,p[T.a]);
    wu = dot3D(w,u);
    wv = dot3D(w,v);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    *c0 = (uv * wv - vv * wu) / D;
    if(fabs(*c0)<EPSILON)
        *c0=0;
    if(fabs(1-*c0)<EPSILON)
        *c0=1;
    if (*c0 < 0.0 || *c0 > 1.0)        // I is outside T
        return 0;

    *c1 = (uv * wu - uu * wv) / D;
    if(fabs(*c1)<EPSILON)
        *c1=0;
    if(fabs(1-*c1)<EPSILON)
        *c1=1;
    if (*c1 < 0.0 || (*c0 + *c1) > 1.0)  // I is outside T
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
#pragma mark -
#pragma mark [ Format conversion ]
int getformatindex(char *path)
{
    char    *formats[]={"orig","pial","white","mesh","sratio","float","curv","txt","inflated","sphere","sulc","reg","txt1","wrl","obj","ply","stl","smesh","off","bin","mgh","annot","raw","vtk"};
    int     i,n=24; // number of recognised formats
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

	FILE	*f;
	int		i,n,l;
	char	*tmp;
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
    fread(&ndim1,1,sizeof(int),f);		swapint(&ndim1);
    fread(&ndim2,1,sizeof(int),f);		swapint(&ndim2);
    fread(&ndim3,1,sizeof(int),f);		swapint(&ndim3);
    fread(&nframes,1,sizeof(int),f);	swapint(&nframes);
    fread(&type,1,sizeof(int),f);		swapint(&type);
    fread(&dof,1,sizeof(int),f);		swapint(&dof);
    
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
    int     *nt=&(m->nt);
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
    sscanf(str," %i %i ",np,nt);
    
    if(*nt==1)    // mesh data file, dimension 1
    {
        *data=(float*)calloc(*np,sizeof(float));
        if(data==NULL){printf("ERROR: Not enough memory for mesh data\n");return 1;}
        for(i=0;i<*np;i++)
            fscanf(f," %f ",&((*data)[i]));    
        if(verbose)
            printf("Read %i values\n",*np);
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
/* What is this obj format??
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
*/
/*
What is this obj format??
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
*/
int Obj_load(char *path, Mesh *m)
{
    int     *np=&(m->np);
    int     *nt=&(m->nt);
    float3D **p=&(m->p);
    int3D   **t=&(m->t);
    FILE    *f;
    char    str[1024],s[16];
    
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
            sscanf(str,"f %i %i %i ",&((*t)[*nt].a),&((*t)[*nt].b),&((*t)[*nt].c));
            (*t)[*nt].a--;
            (*t)[*nt].b--;
            (*t)[*nt].c--;
            (*nt)++;
        }
    }
    if(verbose)
        printf("Read %i vertices and %i triangles\n",*np,*nt);

    // READ TRIANGLES
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
    if(verbose)
        printf("Read %i triangles\n",*nt);

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
    int		i;

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

#pragma mark -
int freeMesh(Mesh *m)
{
    if(verbose) printf("* freeMesh\n");
    
    free(m->p);
    free(m->t);
    if(m->data)
        free(m->data);
    if(m->NT)
        free(m->NT);
    return 0;
}
int loadMesh(char *path, Mesh *m,int iformat)
{
    if(verbose) printf("* imesh\n");

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
        default:
            printf("ERROR: Input mesh format not recognised\n");
            return 1;
    }
    if(err!=0)
    {
        printf("ERROR: cannot read file: %s\n",path);
        return 1;
    }
    
    return 0;
}
int addMesh(char *path, Mesh *m0,int iformat)
{
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
        
    return 0;
}
int saveMesh(char *path, Mesh *m, int oformat)
{
    if(verbose) printf("* omesh\n");

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
#pragma mark [ Geometry functions ]
void absgi(Mesh *m)
{
    float   S,V,logAbsGI;

    S=area(m);
    V=volume(m);
    
    // log(absGI)    = log(Sx)-2log(Vx)/3-log(36π)/3
    // absGI        = Sx/(Vx^(2/3)(36π)^(1/3))
    logAbsGI=log(S)-2*log(V)/3.0-log(36*pi)/3.0;
    
    printf("S=%f, V=%f, gi=%f, log(gi)=%f\n",S,V,exp(logAbsGI),logAbsGI);
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
            p[i]=add3D(p[i],(m->p)[i]);
        free(m->p);
        free(m->t);
    }
    for(i=0;i<np;i++)
        p[i]=sca3D(p[i],1/(float)N);

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
            printf("closest vertex %i (%f,%f,%f), dist=%f\n",imin,p[imin].x,p[imin].y,p[imin].z,dmin);

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
void centre(Mesh *m)
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
void checkOrientation(Mesh *m)
{
    float3D n=normal3D(0,m);
    float   flipTest=dot3D(m->p[m->t[0].a],n);
    
    printf("orientation: %c\n",(flipTest>0)?'+':'-');

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
    printf("min,max=%f,%f\n",min,max);
    */
    
    return 0;
}
int fixflip(Mesh *m)
{
    int     np=m->np;
    int     nt=m->nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    NTriRec *NT;
    int     i,j,pos;
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
    return 0;
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
                angle=acos(dot3D(x,y))*180/pi;
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
int laplace(float lambda, Mesh *m)
{
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
        x=sca3D(tmp[i],1/(float)n[i]);
        dx=sub3D(x,p[i]);
        p[i]=add3D(p[i],sca3D(dx,lambda));    // p=p+l(x-p)
    }
    free(tmp);
    free(n);
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
    int     np,np0;
    float3D *tmp,*p;
    float   *data=m->data;
    int     i,j,k;
    
    np0=m->np;
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
    std=ss/(float)np+pow(s/(float)np,2);
    printf("stdData: %f\n",std);
    return std;
}
void normalise(Mesh *m)
{
    int     *np=&(m->np);
    float3D *p=m->p;
    int     i;
    
    for(i=0;i<*np;i++)
        p[i]=sca3D(p[i],1/norm3D(p[i]));
}
int relaxMesh(char *path, Mesh *m0,int iformat)
{
/*
    // m0: actual mesh
    // m1: reference mesh
    Mesh    *m1;
    int     i,a,b;
    float   sum,l1,l0,k=1;
    float3D f,*p0,*p;
    int     *t0;

    m1=(Mesh*)calloc(1,sizeof(Mesh));
    loadMesh(path,m1,iformat);
    
    p0=m0->p;
    p1=m1->p;
    
    // compute total force
    sum=0;
    for(i=0;i<m0->nt;i++)
    {
        t0=&(m0->t[i]);
        for(j=0;j<3;j++)
        {
            a=t[j];
            b=t[(j+1)%3];
            if(a<b])
            {
                l0=norm3D(sub3D(p0[a],p0[b]));
                l1=norm3D(sub3D(p1[a],p1[b]));
                f=sca3D(sub3D(p0[b],p0[a]),k*(1-l1/l0));
                sum+=norm3D(f);//silly...
            }
        }
    }
    f=sca3D(f,m0->nt*3/2);
    printf
    
    // free m1
    freeMesh(m1);
*/
    return 0;
}
int rotate(Mesh *m, float x, float y, float z)
{
    printf("rotate %f %f %f\n",x,y,z);
    int i;
    float   M[9];
    float3D *p=m->p,pp;
    
   /*
    M[0]=cos(z)*cos(y)*cos(x)-sin(z)*sin(x);
    M[1]=cos(z)*cos(y)*sin(x)+sin(z)*cos(x);
    M[2]=-cos(z)*sin(y);
    
    M[3]=-sin(z)*cos(y)*cos(x)-cos(z)*sin(x);
    M[4]=-sin(z)*cos(y)*sin(x)+cos(z)*cos(x);
    M[5]=sin(z)*sin(y);
    
    M[6]=sin(y)*cos(x);
    M[7]=sin(y)*sin(x);
    M[8]=cos(y);
    */
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
int stereographic(Mesh *m)
{
    int     i,nt;
    float3D *p=m->p;
    int3D   *t=m->t;
    int3D   *t1;    
    float	x,y,z,n;
    float	h,v,delta;
    float   min,max;
    
    // find max and min height
    for(i=0;i<m->np;i++)
    {    
        n=sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
        if(i==0)
            min=max=n;
        if(n<min)
            min=n;
        if(n>max)
            max=n;
    }   
    
    // universal polar stereographic projection 
    for(i=0;i<m->np;i++)
    {    
        n=sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
        x = acos(p[i].x/n);
        y = acos(p[i].y/n);
        z = acos(p[i].z/n);

        if(z*z<0.000001)	delta=1;
        else		        delta = cos(y)/sin(z);
        if(delta<-1) delta=-1;
        if(delta>1) delta=1;
    
        h = z*sin(acos(delta));
        v = z*delta;
        if(x>pi/2.0)
            p[i]=(float3D){-h,v,(max-min)?((n-min)/(max-min)):0};
        else
            p[i]=(float3D){ h,v,(max-min)?((n-min)/(max-min)):0};
    }
    
    // delete triangles close to the border
    t1=(int3D*)calloc(m->nt,sizeof(int3D));
    nt=0;
    for(i=0;i<m->nt;i++)
    {
        if( norm3D(sub3D(p[t[i].a],p[t[i].b]))<pi/4.0 &&
            norm3D(sub3D(p[t[i].b],p[t[i].c]))<pi/4.0 &&
            norm3D(sub3D(p[t[i].c],p[t[i].a]))<pi/4.0)
            t1[nt++]=t[i];
    }
    m->t=t1;
    m->nt=nt;
    printf("new nt: %i\n",nt);
    return 0;
}
int taubin(float lambda, float mu, int N, Mesh *m)
{
    int j;
    
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
    ...I would like to resample the actual mesh to match a reference mesh,
    for example, to average folds, or to average thickness data over many
    meshes.
    1. to average folds, i would need to find for each vertex of the reference
    mesh, the corresponding position in a smooth version of the actual mesh,
    next i would get the coordinates of that point in the folded version of
    the actual mesh, and assign it to the reference.
    2. same thing for scalar data, but instead of 3D coordinates, there would
    be just one value per vertex.
    so, i need, the actual mesh m, the smooth version of the actual mesh m1,
    and the reference mesh rm. I'll have to find for each vertex in rm the
    corresponding triangle and c0, c1 coordinates in m1, and store the coordinates
    of that point in m

    m:          the mesh
    path_m1:    path to smoothed version of the mesh (for example, spherical)
    path_rm:    path to the smoothed version of the reference mesh (for example, spherical)
    */
    Mesh    m1,rm;
    int     nt;
    int3D   *t;
    int     np_rm;
    float3D *p,*p_m1,*p_rm,*tmp;
    float   c0,c1;
    int     i,j,result;
    float3D n;
    float   flipTest;
    
    loadMesh(path_m1,&m1,0);
    loadMesh(path_rm,&rm,0);
    
    // Check whether the meshes are properly oriented
    n=normal3D(0,&m1);
    flipTest=dot3D(m1.p[m1.t[0].a],n);
    if(flipTest<0)
    {
        printf("ERROR: m1 is mis-oriented\n");
        return 1;
    }
    n=normal3D(0,&rm);
    flipTest=dot3D(rm.p[rm.t[0].a],n);
    if(flipTest<0)
    {
        printf("ERROR: rm is mis-oriented\n");
        return 1;
    }
    
    // Actual mesh
    p=m->p;
    nt=m->nt;
    t=m->t;

    // Smoothed version of the actual mesh
    p_m1=m1.p;

    // Reference mesh
    np_rm=rm.np;
    p_rm=rm.p;
    
    tmp=(float3D*)calloc(np_rm,sizeof(float3D));
    for(i=0;i<np_rm;i++)
    {
        if(verbose)
            if(i%1000==0)
                printf("%i\n",i);
        for(j=0;j<nt;j++)
        {
            result=intersect_VectorTriangle(p_rm[i],j,&c0,&c1,&m1);
            if(result==1)
            {
                tmp[i]=p[t[j].a];
                tmp[i]=add3D(tmp[i],sca3D(sub3D(p[t[j].b],p[t[j].a]),c0));
                tmp[i]=add3D(tmp[i],sca3D(sub3D(p[t[j].c],p[t[j].a]),c1));
                break;
            }
        }
        if(j==nt)
        {
            printf("ERROR: could not resample point %i on mesh %s\n",i,path_m1);
            return 1;
        }
    }
    
    // Free actual mesh (m)
    free(m->p);
    free(m->t);
    if(m->data)
        free(m->data);
    if(m->NT)
        free(m->NT);
    
    // Free smooth version of the actual mesh (m1)
    free(m1.p);
    free(m1.t);
    
    // Reconfigure mesh with resampled data
    m->p=tmp;
    m->t=rm.t;
    m->np=rm.np;
    m->nt=rm.nt;
    
    return 0;
}
void printHelp(void)
{
     printf("\
 Commands\n\
    -iformat format name                             Force input format (needs to precede imesh)\n\
    -oformat format name                             Force output format (needs to precede omesh)\n\
    -i filename                                      Input file (also accepts -imesh)\n\
    -o filename                                      Output file (also accepts -omesh)\n\
    -odata filename                                  Output data\n\
    \n\
    -absgi                                           Compute absolute gyrification index\n\
    -add filename                                    Add mesh at filename to the current mesh\n\
    -area                                            Surface area\n\
    -average n_meshes path1 path2 ... pathn          Compute an average of n_meshes all\n\
                                                       of the same topology\n\
    -barycentricProjection reference_mesh            Print barycentric coordinates for each vertex in reference_mesh\n\
    -checkOrientation                                Check that normals point outside\n\
    -centre                                          Move the mesh's barycentre to the origin\n\
    -countClusters  value                            Count clusters in texture data\n\
    -curv                                            Compute curvature\n\
    -euler                                           Print Euler characteristic\n\
    -fixFlip                                         Detect flipped triangles and fix them\n\
    -fixSmall                                        Detect triangles with an angle >160\n\
    -flip                                            Flip normals\n\
                                                       degrees and fix them\n\
    -foldLength                                      Compute total fold length\n\
    -h                                               Help\n\
    -icurv number_of_iterations                      Integrated curvature\n\
    -laplaceSmooth lambda number_of_iterations       Laplace smoothing\n\
    -lissencephalic                                  Smooth valleys and hills, not the coast\n\
    -level level_value                               Adds new vertices (and triangles) to the\n\
                                                       edges that cross level_value in the\n\
                                                       vertex data (f.ex., mean curvature)\n\
    -addVal                                          Add value data\n\
    -subVal                                          Subtract value from data\n\
    -multVal                                         Multiply data time value\n\
    -divVal                                          Divide data by value\n\
    -max                                             Maximum data value\n\
    -mean                                            Mean data value\n\
    -min                                             Minimum data value\n\
    -normal                                          Mesh normal vectors\n\
    -normalise                                       Place all vertices at distance 1 from\n\
                                                       the origin\n\
    -randverts number_of_vertices                    Generate homogeneously distributed\n\
                                                       random vertices over the mesh\n\
    -relax filename                                  Relax current mesh to mesh at filename\n\
                                                        (both meshes have the same topology)\n\
    -resample smooth_mesh reference_mesh             Resample the mesh to match the vertices\n\
                                                       and the topology of the argument mesh\n\
    -rotate x y z                                    Rotate with angles x, y and z\n\
    -scale scale_value                               Multiply each vertex by \"scale\"\n\
    -size                                            Display mesh dimensions\n\
    -stereographic                                   Stereographic projection\n\
    -taubinSmooth lambda mu number_of_iterations     Taubin Smoothing\n\
    -threshold value 0:down/1:up                     Threshold texture data\n\
    -tris                                            Display number of triangles\n\
    -v                                               Verbose mode\n\
    -verts                                           Display number of vertices\n\
    -volume                                          Compute mesh volume\n\
    \n\
    Meshgeometry can read and write several formats, based on the file extension:\n\
    .orig, .pial, .white, .inflated, .sphere, .reg   Freesurfer meshes\n\
    .curv, .sulc, .sratio                            Freesurfer data\n\
    .mesh                                            BrainVisa meshe\n\
    .txt                                             RT's mesh plain text format\n\
    .float                                           Raw float data\n\
    .txt1                                            RT's data format\n\
    .bin                                             n-e-r-v-o-u-s system web binary mesh\n\
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
        if(strcmp(argv[i],"-add")==0)
        {
            addMesh(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-relax")==0)
        {
            relaxMesh(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-barycentricProjection")==0)    // print barycentric coordinates for each vertex of reference mesh
        {
            barycentricProjection(argv[i+1],&mesh);
            i+=1;
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
        if(strcmp(argv[i],"-icurv")==0)
        {
            if(mesh.data==NULL)
                mesh.data=(float*)calloc(mesh.np,sizeof(float));
            icurvature(atoi(argv[++i]),&mesh);
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
        if(strcmp(argv[i],"-euler")==0)
        {
            printf("euler=%i\n",mesh.np-mesh.nt/2);
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
        if(strcmp(argv[i],"-fixSmall")==0)
        {
            fixSmall(&mesh);
        }
        else
        if(strcmp(argv[i],"-scale")==0)
        {
            scale(atof(argv[++i]),&mesh);
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
        if(strcmp(argv[i],"-centre")==0)
        {
            centre(&mesh);
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
            verbose=1;
        
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