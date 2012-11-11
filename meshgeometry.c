//char version[]="meshgeometry, version 1, roberto toro, 12 August 2009";
//char version[]="meshgeometry, version 2, roberto toro, 28 April 2010";
//char version[]="meshgeometry, version 3, roberto toro, 1 May 2010"; // added -laplace and -taubinLM
//char version[]="meshgeometry, version 4, roberto toro, 21 May 2010"; // commands are processed as a chain
//char version[]="meshgeometry, version 5, roberto toro, 28 May 2012"; // added several commands: foldLength, volume, absgi, texture threshold, countClusters, and includes meshconvert v8
char version[]="meshgeometry, version 6, roberto toro, 10 November 2012"; // added randomverts, help, centre, normalise, normal, verbose, off mesh format (load and save)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define pi 3.14159265358979323846264338327950288419716939927510
#define min(a,b) (((a)<(b))?(a):(b))

#define kMAXNETRIS			100
#define kFreeSurferMesh		1
#define kBrainVisaMesh		2
#define kFreeSurferData		3
#define kSRatioFloatData	4
#define kText				5
#define kTextData			6
#define kVRMLMesh			7
#define kObjMesh			8
#define kPlyMesh			9
#define kSTLMesh			10
#define kSmeshMesh			11
#define kBinMesh			12
#define kOffMesh			13

typedef struct
{
	float	x,y,z;
}float3D;
typedef struct
{
	int	a,b,c;
}int3D;

int		np;
int		nt;
float3D	*p;
int3D	*t;
float	*data;
float	R;
int		verbose=0;

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
float3D	sca3D(float3D a, float t)
{
	return (float3D){a.x*t,a.y*t,a.z*t};
}
float norm3D(float3D a)
{
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}
float3D normal3D(int i)
{
	float3D	N;
	N=cross3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
	return sca3D(N,1/norm3D(N));
}

float determinant(float3D a, float3D b, float3D c)
{
	float	D=	a.x*(b.y*c.z-c.y*b.z)+
				a.y*(b.z*c.x-c.z*b.x)+
				a.z*(b.x*c.y-c.x*b.y);
	return D;
}

#pragma mark -
#pragma mark [ Utilities ]
int		endianness;
#define kMOTOROLA	1
#define kINTEL		2
void checkEndianness(void)
{
	char	b[]={1,0,0,0};
	int		num=*(int*)b;
	
	if(num==16777216)
		endianness=kMOTOROLA;
	else
		endianness=kINTEL;
}
void swapint(int *n)
{
	char 	*by=(char*)n;
	char	sw[4]={by[3],by[2],by[1],by[0]};
	
	*n=*(int*)sw;
}
void swapfloat(float *n)
{
	char 	*by=(char*)n;
	char	sw[4]={by[3],by[2],by[1],by[0]};
	
	*n=*(float*)sw;
}
void swaptriangles()
{
	int		i;
	for(i=0;i<nt;i++)
	{
		swapint(&t[i].a);
		swapint(&t[i].b);
		swapint(&t[i].c);
	}
}
void swapvertices()
{
	int		i;
	for(i=0;i<np;i++)
	{
		swapfloat(&p[i].x);
		swapfloat(&p[i].y);
		swapfloat(&p[i].z);
	}
}
#pragma mark -
#pragma mark [ Format conversion ]
int getformatindex(char *path)
{
	char	*formats[]={"orig","pial","white","mesh","sratio","sratiofloat","curv","txt","inflated","sphere","sulc","reg","txt1","wrl","obj","ply","stl","smesh","off","bin"};
	int		i,j,n=20; // number of recognised formats
	int		found,index;
	char	*extension;
	
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
		j=strlen(formats[i]);
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
		index=kSRatioFloatData;
		if(verbose)
			printf("Format: SRatioFloat Data\n");
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
		
	return index;
}
#pragma mark -
int FreeSurfer_load_mesh(char *path)
{
    FILE	*f;
    int		id,a,b,c,i;
    char	date[256],info[256];

    f=fopen(path,"r");
	
	if(f==NULL)
		return 1;

	// read triangle/quad identifier: 3 bytes
    a=((int)(u_int8_t)fgetc(f))<<16;
    b=((int)(u_int8_t)fgetc(f))<<8;
    c=(u_int8_t)fgetc(f);
    id=a+b+c;
    if(id==16777214)	// triangle mesh
    {
		fgets(date,256,f);
        fgets(info,256,f);
		fread(&np,1,sizeof(int),f); if(endianness==kINTEL) swapint(&np);
		fread(&nt,1,sizeof(int),f);	if(endianness==kINTEL) swapint(&nt);
        // read vertices
        p=(float3D*)calloc(np,sizeof(float3D));
			if(p==NULL) printf("ERROR: Cannot allocate memory for points [FreeSurfer_load_mesh]\n");
			else
			{
		fread((char*)p,np,3*sizeof(float),f);	if(endianness==kINTEL) swapvertices();
        // read triangles
        t=(int3D*)calloc(nt,sizeof(int3D));
			if(t==NULL) printf("ERROR: Cannot allocate memory for triangles [FreeSurfer_load_mesh]\n");
			else
			{
		fread((char*)t,nt,3*sizeof(int),f);		if(endianness==kINTEL) swaptriangles();
			}
			}
    }
	fclose(f);

	for(i=0;i<np;i++)
	{
		p[i].x+=128;
		p[i].y+=128;
		p[i].z+=128;
	}
	
	return 0;
}
int FreeSurfer_load_data(char *path)
{
    FILE	*f;
	int		i,j;
    int		id,a,b,c;
	char	byte4[4];

	char testInt[]={1,0,0,0};
	int endianness=*(int*)testInt;

    f=fopen(path,"r");
	if(f==NULL)
		return 1;

	// read identifier: 3 bytes
    a=((int)(u_int8_t)fgetc(f))<<16;
    b=((int)(u_int8_t)fgetc(f))<<8;
    c=(u_int8_t)fgetc(f);
    id=a+b+c;
    if(id==16777215)	// triangle mesh
    {
		if(endianness==1)	// Intel
			for(i=0;i<4;i++) byte4[3-i]=fgetc(f);
		else
			fread(byte4,4,sizeof(char),f);
		np=*(int*)byte4;
		if(verbose)
			printf("FS #vertex_data %i\n",np);
        
		data=(float*)calloc(np,sizeof(float));
		
		// disregard FaceCount and ValsPerVertex
		fgetc(f);fgetc(f);fgetc(f);fgetc(f);
		fgetc(f);fgetc(f);fgetc(f);fgetc(f);
		
        // read vertex data
        for(j=0;j<np;j++)
		{
			if(endianness==1)	// Intel
				for(i=0;i<4;i++) byte4[3-i]=fgetc(f);
			else
				fread(byte4,4,sizeof(char),f);
            data[j]=*(float*)byte4;
		}
    }
	if(verbose)
		printf("FSData finished\n");
	
	fclose(f);
	
	return 0;
}

int FreeSurfer_save_mesh(char *path)
{
	FILE	*f;
    int		id=16777214,a,b,c;
    int		NP,NT,i;
    char	date[6]="EMPTY",info[6]="EMPTY";
    float3D	ftmp;
    int3D	itmp;

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
	NP=np;
	NT=nt;
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
		for(i=0;i<np;i++)
		{
			ftmp=p[i];
			swapfloat(&ftmp.x);
			swapfloat(&ftmp.y);
			swapfloat(&ftmp.z);
			fwrite(&ftmp,1,sizeof(float3D),f);
		}
		for(i=0;i<nt;i++)
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
		fwrite(p,np,sizeof(float3D),f);
		fwrite(t,nt,sizeof(int3D),f);
	}
	fclose(f);

	return 0;
}
int FreeSurfer_save_data(char *path)
{
	FILE	*f;
    int		id=16777215,a,b,c;
    int		FaceCount=0,ValsPerVertex=1,i,n;
    float	x;

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
    
	n=np;
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
		for(i=0;i<np;i++)
		{
			x=data[i];
			swapfloat(&x);
			fwrite(&x,1,sizeof(float),f);
		}
	}
	else
		fwrite(data,np,sizeof(float),f);
	fclose(f);

	return 0;
}
int BrainVisa_load_mesh(char *path)
{
    FILE	*f;
	char	tmp[6];
    int		i;
    int		endian,ignore;
    
    f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// READ HEADER
	// get format (ascii, binar)
    fread(tmp,5,sizeof(char),f); tmp[5]=(char)0;
    if(strcmp(tmp,"binar")==0)
    {
    	for(i=0;i<4;i++) tmp[i]=fgetc(f); tmp[4]=(char)0;
		endian=-1;
		if(strcmp(tmp,"ABCD")==0)	endian=kMOTOROLA;	
		if(strcmp(tmp,"DCBA")==0)	endian=kINTEL;
		if(endian==-1){ printf("ERROR: Not ABCD nor DCBA order...exit.\n"); return 1;}
		fread(&ignore,4,sizeof(char),f);		// ignore "VOID" string length
		fread(&ignore,4,sizeof(char),f);		// ignore "VOID" string
		fread(&ignore,1,sizeof(int),f);		// verify number of vertices per polygon
		if(endian!=endianness)
			swapint(&ignore);
		if(ignore!=3){ printf("ERROR: Only able to read triangle meshes. This mesh has %i vertices per polygon.\n",ignore); return 1;}
		fread(&ignore,1,sizeof(int),f);		// ignore time steps
		fread(&ignore,1,sizeof(int),f);		// ignore time step index
		
		// READ VERTICES
		fread(&np,1,sizeof(int),f);			// read number of vertices
		if(endian!=endianness)
			swapint(&np);
		p = (float3D*)calloc(np,sizeof(float3D));
		if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
		fread((char*)p,np,3*sizeof(float),f);	if(endian!=endianness) swapvertices();	
		if(verbose)
			printf("Read %i vertices\n",np);
		
		// IGNORE NORMAL VECTORS
		fseek(f,sizeof(int),SEEK_CUR);		// ignore normal vectors
		fseek(f,np*sizeof(float3D),SEEK_CUR);
		fread(&ignore,1,sizeof(int),f);		// ignore number of texture coordinates
		
		// READ TRIANGLES
		fread(&nt,1,sizeof(int),f);			// read number of triangles
		if(endian!=endianness)
			swapint(&nt);
		t = (int3D*)calloc(nt,sizeof(int3D));
		if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
		fread((char*)t,nt,3*sizeof(int),f);		if(endian!=endianness) swaptriangles();
		if(verbose)
			printf("Read %i triangles\n",nt);
	}
	else
    if(strcmp(tmp,"ascii")==0)
    {
    	fscanf(f," %*s ");			// ignore VOID
    	fscanf(f," %*i %*i %*i ");	// ignore 3 integers
    	
    	// READ 3-D COORDINATES
    	fscanf(f," %i ",&np);
		p = (float3D*)calloc(np,sizeof(float3D));
    	for(i=0;i<np;i++)
    		fscanf(f," ( %f , %f , %f ) ", &p[i].x,&p[i].y,&p[i].z);
    	
    	fscanf(f," %*i ");			// ignore number of normal vectors
    	for(i=0;i<np;i++)			// ignore normal vectors
    		fscanf(f," ( %*f , %*f , %*f ) ");
    	fscanf(f," %*i ");			// ignore an integer
    	
    	// READ TRIANGLES
    	fscanf(f," %i ",&nt);
		t = (int3D*)calloc(nt,sizeof(int3D));
    	for(i=0;i<nt;i++)
    		fscanf(f," ( %i , %i , %i ) ", &t[i].a,&t[i].b,&t[i].c);
    }
    else
    {
    	printf("ERROR: Cannot read '%s' format.\n",tmp);
    	return 1;
    }
	
    
	fclose(f);
    
    return 0;
}
int BrainVisa_save_mesh(char *path)
{
    FILE	*f;
    int		i;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
    	fprintf(f,"ascii\n");
    
    	fprintf(f,"VOID\n");	// ignore VOID
    	fprintf(f,"3\n1\n0\n");	// ignore 3 integers
    	
    	// WRITE 3-D COORDINATES
    	fprintf(f,"%i\n",np);
    	for(i=0;i<np;i++)
    		fprintf(f,"(%f,%f,%f) ",p[i].x,p[i].y,p[i].z);
    	fprintf(f,"\n");
    	fprintf(f,"0\n");			// ignore number of normal vectors
    	fprintf(f,"0\n");			// ignore an integer
    	
    	// WRITE TRIANGLES
    	fprintf(f,"%i\n",nt);
    	for(i=0;i<nt;i++)
    		fprintf(f,"(%i,%i,%i) ",t[i].a,t[i].b,t[i].c);
    	fprintf(f,"\n");
	fclose(f);
    
    return 0;
}
int Text_load(char *path)
{
    FILE	*f;
    int		i;
    char	str[512];
    
    f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// READ HEADER
	fgets(str,511,f);
	sscanf(str," %i %i ",&np,&nt);
	
	if(nt==1)	// mesh data file, dimension 1
	{
		data=(float*)calloc(np,sizeof(float));
		if(data==NULL){printf("ERROR: Not enough memory for mesh data\n");return 1;}
		for(i=0;i<np;i++)
			fscanf(f," %f ",&data[i]);	
		if(verbose)
			printf("Read %i values\n",np);
	}
	else
	{			// mesh file
		// READ VERTICES
		p = (float3D*)calloc(np,sizeof(float3D));
		if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
		for(i=0;i<np;i++)
			fscanf(f," %f %f %f ",&p[i].x,&p[i].y,&p[i].z);	
		if(verbose)
			printf("Read %i vertices\n",np);
	
		// READ TRIANGLES
		t = (int3D*)calloc(nt,sizeof(int3D));
		if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
		for(i=0;i<nt;i++)
			fscanf(f," %i %i %i ",&t[i].a,&t[i].b,&t[i].c);
		if(verbose)
			printf("Read %i triangles\n",nt);
	}

	fclose(f);
    
    return 0;
}
int Text_save_mesh(char *path)
{
    FILE	*f;
    int		i;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// WRITE HEADER
	fprintf(f,"%i %i\n",np,nt);

	// WRITE VERTICES
	for(i=0;i<np;i++)
		fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);	

	// WRITE TRIANGLES
	for(i=0;i<nt;i++)
		fprintf(f,"%i %i %i\n",t[i].a,t[i].b,t[i].c);

	fclose(f);
    
    return 0;
}
int Text_save_data(char *path)
{
    FILE	*f;
    int		i;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// WRITE HEADER
	fprintf(f,"%i 1\n",np);

	// WRITE DATA
	for(i=0;i<np;i++)
		fprintf(f,"%f\n",data[i]);	

	fclose(f);
    
    return 0;
}
int VRML_load_mesh(char *path)
{
    FILE	*f;
	int		i,loop;
	char	str[256];

    f=fopen(path,"r");

	np=0;
	nt=0;
	
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
			np++;
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
			nt++;
	}

	p = (float3D*)calloc(np,sizeof(float3D));
	t = (int3D*)calloc(nt,sizeof(int3D));
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
			sscanf(str," %f %f %f ",&p[i].x,&p[i].y,&p[i].z);
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
			sscanf(str," %i , %i , %i ",&t[i].a,&t[i].c,&t[i].b);
			i++;
		}
	}
	fclose(f);

	return 0;
}
int VRML_save_mesh(char *path)
{
	FILE	*f;
    int		i;
    
    f=fopen(path,"w");
    fprintf(f,"#VRML V1.0 ascii\n");
    fprintf(f,"Separator {\n");
    fprintf(f,"ShapeHints {\n");
    fprintf(f,"vertexOrdering COUNTERCLOCKWISE\n");
    fprintf(f,"faceType CONVEX\n");
    fprintf(f,"}\n");
    fprintf(f,"Coordinate3 {\n");
    fprintf(f,"point [\n");
    for(i=0;i<np;i++)
        fprintf(f,"%f %f %f,\n",p[i].x,p[i].y,p[i].z);
	fprintf(f,"]\n");
	fprintf(f,"}\n");
	fprintf(f,"IndexedFaceSet {\n");
	fprintf(f,"coordIndex [\n");
    for(i=0;i<nt;i++)
        fprintf(f,"%i,%i,%i,-1\n",t[i].a,t[i].b,t[i].c);
	fprintf(f,"]\n");
	fprintf(f,"}\n");
	fprintf(f,"}\n");
	fclose(f);

    return 0;
}
int Obj_load(char *path)
{
    FILE	*f;
    int		i;
    char	str[512];
    
    f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// READ HEADER
	fgets(str,511,f);
	sscanf(str," %*c %*f %*f %*f %*i %*i %i ",&np);
	
	// READ VERTICES
	p = (float3D*)calloc(np,sizeof(float3D));
	if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
	for(i=0;i<np;i++)
		fscanf(f," %f %f %f ",&p[i].x,&p[i].y,&p[i].z);	
	if(verbose)
		printf("Read %i vertices\n",np);
	
	// IGNORE NORMALS
//	fgets(str,511,f); // skip empty line
	for(i=0;i<np;i++)
		fscanf(f," %*f %*f %*f ");	

	// READ TRIANGLES
//	fgets(str,511,f); 	// skip empty line
	fscanf(f," %i ",&nt);
//	fgets(str,511,f); 	// skip empty line
	for(i=0;i<5+nt;i++)	// skip nt+5 integers (the first 5, all integers, then nt multiples of 3)
		fscanf(f," %*i ");
	t = (int3D*)calloc(nt,sizeof(int3D));
	if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
	for(i=0;i<nt;i++)
		fscanf(f," %i %i %i ",&t[i].a,&t[i].b,&t[i].c);
	if(verbose)
		printf("Read %i triangles\n",nt);

	fclose(f);
    
    return 0;
}
int Obj_save_mesh(char *path)
{
    FILE	*f;
    int		i;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// WRITE HEADER
	fprintf(f,"P 0.3 0.3 0.4 10 1 %i\n",np);

	// WRITE VERTICES
	for(i=0;i<np;i++)
		fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
	fprintf(f,"\n");
	
	// WRITE DUMMY NORMALS
	for(i=0;i<np;i++)
		fprintf(f,"0 0 0\n");
	fprintf(f,"\n");
	
	// WRITE NUMBER OF TRIANGLES
	fprintf(f,"%i\n",nt);
	
	// WRITE 5 DUMMY NUMBERS
	fprintf(f,"0 1 1 1 1 1\n");
	fprintf(f,"\n");
	
	// WRITE nt MULTIPLES OF 3, IN ROWS OF EIGHT
	for(i=0;i<nt;i++)
	{
		fprintf(f,"%i ",(i+1)*3);
		if(i>0 && i%8==0)
			fprintf(f,"\n");
	}
	fprintf(f,"\n");

	// WRITE TRIANGLES
	for(i=0;i<nt;i++)
		fprintf(f,"%i %i %i\n",t[i].a,t[i].b,t[i].c);

	fclose(f);
    
    return 0;
}
int Ply_load(char *path)
{
    FILE	*f;
    int		i,x;
    char	str[512],str1[256],str2[256];
    
    f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// READ HEADER
	np=nt=0;
	do
	{
		fgets(str,511,f);
		sscanf(str," %s %s %i ",str1,str2,&x);
		if(strcmp(str1,"element")==0&&strcmp(str2,"vertex")==0)
			np=x;
		else
		if(strcmp(str1,"element")==0&&strcmp(str2,"face")==0)
			nt=x;
	}
	while(strcmp(str1,"end_header")!=0 && !feof(f));
	if(np*nt==0)
	{
		printf("ERROR: Bad Ply file header format\n");
		return 1;
	}
	// READ VERTICES
	p = (float3D*)calloc(np,sizeof(float3D));
	if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
	for(i=0;i<np;i++)
		fscanf(f," %f %f %f ",&p[i].x,&p[i].y,&p[i].z);	
	if(verbose)
		printf("Read %i vertices\n",np);

	// READ TRIANGLES
	t = (int3D*)calloc(nt,sizeof(int3D));
	if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
	for(i=0;i<nt;i++)
		fscanf(f," 3 %i %i %i ",&t[i].a,&t[i].b,&t[i].c);
	if(verbose)
		printf("Read %i triangles\n",nt);

	fclose(f);
    
    return 0;
}
int Ply_save_mesh(char *path)
{
    FILE	*f;
    int		i;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// WRITE HEADER
	fprintf(f,"ply\n");
	fprintf(f,"format ascii 1.0\n");
	fprintf(f,"comment meshconvert, R. Toro 2010\n");
	fprintf(f,"element vertex %i\n",np);
	fprintf(f,"property float x\n");
	fprintf(f,"property float y\n");
	fprintf(f,"property float z\n");
	fprintf(f,"element face %i\n",nt);
	fprintf(f,"property list uchar int vertex_indices\n");
	fprintf(f,"end_header\n");

	// WRITE VERTICES
	for(i=0;i<np;i++)
		fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);	

	// WRITE TRIANGLES
	for(i=0;i<nt;i++)
		fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);

	fclose(f);
    
    return 0;
}
int STL_load(char *path)
{
    printf("ERROR: meshconvert DOES NOT LOAD STL DATA FOR THE MOMENT, ONLY SAVES IT\n");
    return 1;
    /*
    FILE	*f;
    int		i;
    char	str[512];
    
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
int STL_save_mesh(char *path)
{
    FILE	*f;
    int		i;
    float3D	n;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// WRITE HEADER
	fprintf(f,"solid mySolid\n");

	// WRITE FACES
	for(i=0;i<nt;i++)
	{
		n=normal3D(i);
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
int Smesh_load(char *path)
{
    FILE	*f;
    int		i;
    char	str[512];
    
    f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// READ POINTS HEADER
	fgets(str,511,f);
	sscanf(str," %i ",&np);
	// READ VERTICES
	p = (float3D*)calloc(np,sizeof(float3D));
	if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
	for(i=0;i<np;i++)
		fscanf(f," %*i %f %f %f ",&p[i].x,&p[i].y,&p[i].z);	
	if(verbose)
		printf("Read %i vertices\n",np);
	
	// READ TRIANGLES HEADER
	// READ TRIANGLES
	fgets(str,511,f);
	sscanf(str," %i ",&nt);
	t = (int3D*)calloc(nt,sizeof(int3D));
	if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
	for(i=0;i<nt;i++)
		fscanf(f," %*i %i %i %i ",&t[i].a,&t[i].b,&t[i].c);
	if(verbose)
		printf("Read %i triangles\n",nt);

	fclose(f);
    return 0;
}
int Smesh_save_mesh(char *path)
{
    FILE	*f;
    int		i;
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	// WRITE VERTICES HEADER
	fprintf(f,"%i 3 0 0\n",np);
	// WRITE VERTICES
	for(i=0;i<np;i++)
		fprintf(f,"%i %f %f %f 0\n",i,p[i].x,p[i].y,p[i].z);	
	// WRITE TRIANGLES HEADER
	fprintf(f,"%i 0\n",nt);
	// WRITE TRIANGLES
	for(i=0;i<nt;i++)
		fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);
	fprintf(f,"0\n0\n");
	fclose(f);   
    return 0;
}
int Bin_load(char *path)
{
    printf("ERROR: Load is not yet implemented for Bin filetype\n");
    return 1;
}
int Bin_save_mesh(char *path)
{
	// Bin filetype is the binary format used by n-e-r-v-o-u-s system
	// to display meshes in the web
    FILE	*f;
    int		i;
    int		itmp;
    short	stmp;
    float	ftmp;
    
    f=fopen(path,"w");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	// WRITE NUMBER OF VERTICES
	itmp=np;
	swapint(&itmp);
	fwrite(&itmp,1,sizeof(int),f);
	// WRITE NUMBER OF TRIANGLES
	itmp=nt;
	swapint(&itmp);
	fwrite(&itmp,1,sizeof(int),f);
	// WRITE VERTICES
	for(i=0;i<np;i++)
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
	for(i=0;i<nt;i++)
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
int Off_load(char *path)
{
	int			i;
	FILE		*f;
	char		str[512];
	
	f=fopen(path,"r");

	fgets(str,512,f);	// skip OFF
	fscanf(f," %i %i %*i ",&np,&nt);
	p=(float3D*)calloc(np,sizeof(float3D));
	t=(int3D*)calloc(nt,sizeof(int3D));
	for(i=0;i<np;i++)
		fscanf(f," %f %f %f ",&p[i].x,&p[i].y,&p[i].z);
	if(verbose)
		printf("Read %i vertices\n",np);
	for(i=0;i<nt;i++)
		fscanf(f," %*i %i %i %i ",&t[i].a,&t[i].b,&t[i].c);
	if(verbose)
		printf("Read %i triangles\n",nt);
	fclose(f);
	return 0;
}
int Off_save_mesh(char *path)
{
    int		i;
    FILE	*f;

    f=fopen(path,"w");
    
    fprintf(f,"OFF\n");
    fprintf(f,"%i %i 0\n",np,nt);
    for(i=0;i<np;i++)
        fprintf(f,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    for(i=0;i<nt;i++)
        fprintf(f,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);
    fclose(f);
    return 0;
}
#pragma mark -
int loadMesh(char *path)
{
	int		err,format;
	
	format=getformatindex(path);
	
	switch(format)
	{
		case kFreeSurferMesh:
			err=FreeSurfer_load_mesh(path);
			break;
		case kFreeSurferData:
			err=FreeSurfer_load_data(path);
			break;
		case kBrainVisaMesh:
			err=BrainVisa_load_mesh(path);
			break;
		case kText:
			err=Text_load(path);
			break;
		case kVRMLMesh:
			err=VRML_load_mesh(path);
			break;
		case kObjMesh:
			err=Obj_load(path);
			break;
		case kPlyMesh:
			err=Ply_load(path);
			break;
		case kSTLMesh:
			err=STL_load(path);
			break;
		case kSmeshMesh:
			err=Smesh_load(path);
			break;
		case kBinMesh:
			err=Bin_load(path);
			break;
		case kOffMesh:
			err=Off_load(path);
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
int saveMesh(char *path)
{
	int	err,format;
	
	format=getformatindex(path);

	switch(format)
	{
		case kFreeSurferMesh:
			FreeSurfer_save_mesh(path);
			break;
		case kFreeSurferData:
			FreeSurfer_save_data(path);
			break;
		case kBrainVisaMesh:
			BrainVisa_save_mesh(path);
			break;
		case kText:
			Text_save_mesh(path);
			break;
		case kTextData:
			Text_save_data(path);
			break;
		case kVRMLMesh:
			VRML_save_mesh(path);
			break;
		case kObjMesh:
			Obj_save_mesh(path);
			break;
		case kPlyMesh:
			Ply_save_mesh(path);
			break;
		case kSTLMesh:
			STL_save_mesh(path);
			break;
		case kSmeshMesh:
			Smesh_save_mesh(path);
			break;
		case kBinMesh:
			Bin_save_mesh(path);
			break;
		case kOffMesh:
			Off_save_mesh(path);
			break;
		default:
			printf("ERROR: Output data format not recognised\n");
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
int smooth(void)
{
	float3D	*tmp;
	int		*n;
	int		i;
	
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
    	p[i]=sca3D(tmp[i],1/(float)n[i]);
    free(tmp);
    free(n);
	return 0;
}
int laplace(float lambda)
{
	float3D	*tmp,x,dx;
	int		*n;
	int		i;
	
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
    {
    	x=sca3D(tmp[i],1/(float)n[i]);
    	dx=sub3D(x,p[i]);
    	p[i]=add3D(p[i],sca3D(dx,lambda));	// p=p+l(x-p)
    }
    free(tmp);
    free(n);
	return 0;
}
int curvature(float	*C)
{
    float3D		*tmp,*tmp1;
    int			*n;
    float3D		nn;
    float		absmax;
    int			i;

    tmp=(float3D*)calloc(np,sizeof(float3D));
    n=(int*)calloc(np,sizeof(int));
    // compute smoothing direction as the vector to the average of neighbour vertices
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
    	tmp[i]=sub3D(sca3D(tmp[i],1/(float)n[i]),p[i]);
    	
    tmp1=(float3D*)calloc(np,sizeof(float3D));
    // compute normal direction as the average of neighbour triangle normals
    for(i=0;i<nt;i++)
    {
    	nn=cross3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
    	nn=sca3D(nn,1/norm3D(nn));
    	tmp1[t[i].a]=add3D(tmp1[t[i].a],nn);
    	tmp1[t[i].b]=add3D(tmp1[t[i].b],nn);
    	tmp1[t[i].c]=add3D(tmp1[t[i].c],nn);
    }
    for(i=0;i<np;i++)
    	tmp1[i]=sca3D(tmp1[i],1/(float)n[i]);
    free(n);
    
    for(i=0;i<np;i++)
		C[i]=-dot3D(tmp1[i],tmp[i]);
	free(tmp);
	free(tmp1);
    absmax=-1;
    for(i=0;i<np;i++)
        absmax=(fabs(C[i])>absmax)?fabs(C[i]):absmax;
    absmax*=0.95;
    for(i=0;i<np;i++)
    {
        C[i]/=absmax;
        if(C[i]>1)	C[i]=1;
        if(C[i]<-1)	C[i]=-1;
    }
    
    return 0;
}
int icurvature(float *C, int iter)
{
    float		*tmp;
    float		absmax;
    int			i,k;

    tmp=(float*)calloc(np,sizeof(float));
    for(k=0;k<iter;k++)
    {
        curvature(tmp);
        for(i=0;i<np;i++)
            C[i] += tmp[i]*(1+k/(float)iter);
        
        smooth();
        smooth();
    }
    absmax=-1;
    for(i=0;i<np;i++)
        absmax = (fabs(C[i])>absmax)?fabs(C[i]):absmax;
    for(i=0;i<np;i++)
        C[i]/=absmax;
    
    return 0;
}
int scale(float t)
{
	int	i;
	
	for(i=0;i<np;i++)
		p[i]=sca3D(p[i],t);
	
	return 0;
}
int flip(void)
{
	int	i;
	for(i=0;i<nt;i++)
		t[i]=(int3D){t[i].a,t[i].c,t[i].b};
	return 0;
}
double	sum;
int		*tmark,icmax,ncverts;
#define SIZESTACK	64
typedef struct {
	int		n;
	int		t[SIZESTACK];
}NTriRec, *NtriPtr;
void cluster(NTriRec *NT, int ip, float *thrsrc, float thr)
{
	int		i,j,it;
	int		inside;
	int		*tt;
	
	tmark[ip]=1;
	ncverts++;	
	for(i=0;i<=NT[ip].n;i++)
	{
		it=NT[ip].t[i];		
		tt=(int*)&(t[it]);
		inside=0;
		for(j=0;j<3;j++)
			if(thrsrc[tt[j]]>=thr && tmark[tt[j]]==0)
			{
				cluster(NT,tt[j],thrsrc,thr);
				if(thrsrc[tt[j]]>thrsrc[icmax])
					icmax=tt[j];
			}
	}
}
void countClusters(float *thrsrc, float thr)
{
	int		n;
	int		i;
	
	NTriRec	*NT;
	NT=(NTriRec*)calloc(np,sizeof(NTriRec));		
	for(i=0;i<nt;i++)
	{
		NT[t[i].a].t[NT[t[i].a].n++] = i;
		NT[t[i].b].t[NT[t[i].b].n++] = i;
		NT[t[i].c].t[NT[t[i].c].n++] = i;
	}

	n=1;
	tmark=(int*)calloc(np,sizeof(int));
	printf("number\tnVertices\timax\tmax\tXmax\tYmax\tZmax\n");
	for(i=0;i<np;i++)
	if(thrsrc[i]>=thr && tmark[i]==0)
	{
		icmax=i;
		ncverts=0;
		cluster(NT,i,thrsrc,thr);
		printf("%i\t%i\t%i\t%f\t%f\t%f\t%f\n",n,ncverts,icmax,thrsrc[icmax],p[icmax].x,p[icmax].y,p[icmax].z);
		n++;
	}
	printf("\n");
		
	free(tmark);
}
void foldLength(void)
{
	int		i,j;
	float	length=0,a,x;
	float3D	p0[3];
	
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
	printf("foldLength: %g\n",length/2.0);
}
void centre(void)
{
	int		i;
	float3D	centre={0,0,0};
	
	for(i=0;i<np;i++)
		centre=add3D(centre,p[i]);
	centre=sca3D(centre,1/(float)np);
	for(i=0;i<np;i++)
		p[i]=sub3D(p[i],centre);
}
void normalise(void)
{
	int		i;
	
	for(i=0;i<np;i++)
		p[i]=sca3D(p[i],1/norm3D(p[i]));
}
// triangle area using Heron's formula
float triangle_area(float3D p0, float3D p1, float3D p2)
{
	float	a,b,c;	// side lengths
	float	s;		// semiperimeter
	float	area;
	
	a=norm3D(sub3D(p0,p1));
	b=norm3D(sub3D(p1,p2));
	c=norm3D(sub3D(p2,p0));
	s=(a+b+c)/2.0;
	
	area=sqrt(s*(s-a)*(s-b)*(s-c));
	
	return area;
}
float area(void)
{
	int		i;
	float	area=0;
	
	for(i=0;i<nt;i++)
		area+=triangle_area(p[t[i].a],p[t[i].b],p[t[i].c]);
	printf("area: %g\n",area);
	
	return area;
}
float volume(void) // Code By S Melax from http://www.melax.com/volint/
{
	float	vol=0;
	int		i;
	
	for(i=0;i<nt;i++)
		vol += determinant(p[t[i].a],p[t[i].b],p[t[i].c]); //divide by 6 later for efficiency
	vol/=6.0;// since the determinant give 6 times tetra volume
	printf("volume: %g\n",vol);
	return vol;
}
void absgi(void)
{
	float	S,V,logAbsGI;

	S=area();
	V=volume();
	
	// log(absGI)	= log(Sx)-2log(Vx)/3-log(36π)/3
	// absGI		= Sx/(Vx^(2/3)(36π)^(1/3))
	logAbsGI=log(S)-2*log(V)/3.0-log(36*pi)/3.0;
	
	printf("S=%f, V=%f, gi=%f, log(gi)=%f\n",S,V,exp(logAbsGI),logAbsGI);
}
void threshold(float thr, int direction)
{
	int		i;

	if(direction==0)
	{
		for(i=0;i<np;i++)
			if(data[i]<=thr)
				data[i]=1;
			else
				data[i]=-1;
	}
	else
	{
		for(i=0;i<np;i++)
			if(data[i]>=thr)
				data[i]=1;
			else
				data[i]=-1;
	}
}
int normal(float *C)
{
	// IMPORTANT NOTE: everywhere C is just np scalars, but here is np 3d vectors.
	// There is no way, for the moment, of saving this normal.
	
	printf("WARNING: \"normal\" has not been tested\n");
    
    float3D		*tmp;
    int			*n;
    int			i;

    tmp=(float3D*)calloc(np,sizeof(float3D));
    n=(int*)calloc(np,sizeof(int));
    // All normals are weighted the same, but more
    // correctly, the values should be weighted by the angle
    // at the given vertex.
    for(i=0;i<nt;i++)
    {
    	tmp[t[i].a]=add3D(tmp[t[i].a],normal3D(i));
    	tmp[t[i].b]=add3D(tmp[t[i].b],normal3D(i));
    	tmp[t[i].c]=add3D(tmp[t[i].c],normal3D(i));
    	n[t[i].a]++;
    	n[t[i].b]++;
    	n[t[i].c]++;
    }
    for(i=0;i<np;i++)
    	((float3D*)C)[i]=sca3D(tmp[i],1/(float)n[i]);
	free(tmp);
    
    return 0;
}
void randverts(int nrv)
{
	// This function is meant to be used with qhull, like this:
	// meshgeometry -imesh $mesh -centre -normalise -randverts 100 | qhull o
	int			i,j;
	float		s0,s1;
	float3D		v;
				
	printf("3 meshgeometry randverts\n");
	printf("%i\n",nrv);
	// generate random vertices
	for(i=0;i<nrv;i++)
	{
		j=(nt-0.01)*(rand()/(float)RAND_MAX);	// pick a random triangle

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
void printHelp(void)
{
	 printf("\
 Commands\n\
	imesh filename                                  Input mesh\n\
	omesh filename                                  Output mesh\n\
	odata filename                                  Output data\n\
 	curv                                            Compute curvature\n\
 	icurv number_of_iterations                      Integrated curvature\n\
	laplaceSmooth lambda number_of_iterations       Laplace smoothing\n\
	taubinSmooth lambda mu number_of_iterations     Taubin Smoothing\n\
	scale scale                                     Scale multiplying each point by \"scale\"\n\
	euler                                           Print Euler characteristic\n\
	flip                                            Flip normals\n\
	foldLength                                      Compute total fold length\n\
	centre                                          Move the mesh's barycentre to the origin\n\
	normalise                                       Place all vertices at distance 1 from the origin\n\
	volume                                          Compute mesh volume\n\
	absgi                                           Compute absolute gyrification index\n\
	threshold                                       Threshold texture data\n\
	countClusters                                   Count clusters in texture data\n\
	normal                                          Mesh normal\n\
	randverts number_of_vertices                    Generate homogeneously\n\
	                                                 distributed random vertices\n\
	                                                 over the mesh\n\
    v                                               Verbose mode\n\
	h                                               Help\n\
    \n\
	Meshgeometry can read and write several formats, guessed from the file\n\
    extension. These are:\n\
    Freesurfer meshes with extension orig, pial, white, inflated, sphere,reg;\n\
    Freesurfer data with extension curv, sulc, sratio;\n\
    BrainVisa meshes with extension mesh;\n\
    RT's mesh format txt, and data format sratiofloat, txt1, bin (for use in the web);\n\
    and mesh formats wrl, obj, ply, stl, smesh and off.\n\
");
}
int main(int argc, char *argv[])
{
	checkEndianness();
	srand (time(NULL));
	
	int	i;
	
	data=NULL;
	
	i=1;
	while(i<argc)
	{
		if(strcmp(argv[i],"-imesh")==0)
		{
			loadMesh(argv[++i]);
		}
		else
		if(strcmp(argv[i],"-omesh")==0)
		{
			saveMesh(argv[++i]);
		}
		else
		if(strcmp(argv[i],"-odata")==0)
		{
			Text_save_data(argv[++i]);
		}
		else
		if(strcmp(argv[i],"-curv")==0)
		{
			if(data==NULL)
				data=(float*)calloc(np,sizeof(float));
			curvature(data);
		}
		else
		if(strcmp(argv[i],"-icurv")==0)
		{
			if(data==NULL)
				data=(float*)calloc(np,sizeof(float));
			icurvature(data,atoi(argv[++i]));
		}
		else
		if(strcmp(argv[i],"-laplaceSmooth")==0)
		{
			int		j,N;
			float	l;
			
			l=atof(argv[++i]);
			N=atoi(argv[++i]);
			for(j=0;j<N;j++)
				laplace(l);
		}
		else
		if(strcmp(argv[i],"-taubinSmooth")==0)
		{
			int		j,N;
			float	l,m;
			
			l=atof(argv[++i]);
			m=atof(argv[++i]);
			N=atoi(argv[++i]);
			for(j=0;j<2*N;j++)
				if(j%2==0)
					laplace(l);
				else
					laplace(m);
		}
		else
		if(strcmp(argv[i],"-euler")==0)
		{
			printf("euler=%i\n",np-nt/2);
		}
		else
		if(strcmp(argv[i],"-flip")==0)
		{
			flip();
		}
		else
		if(strcmp(argv[i],"-scale")==0)
		{
			scale(atof(argv[++i]));
		}
		else
		if(strcmp(argv[i],"-area")==0)
		{
			area();
		}
		else
		if(strcmp(argv[i],"-volume")==0)
		{
			volume();
		}
		else
		if(strcmp(argv[i],"-absgi")==0)
		{
			absgi();
		}
		else
		if(strcmp(argv[i],"-threshold")==0)	// threshold value 0=down/1=up
		{
			float	value=atof(argv[++i]);
			int		direction=atoi(argv[++i]);
			threshold(value,direction);
		}
		else
		if(strcmp(argv[i],"-foldLength")==0)
		{
			foldLength();
		}
		else
		if(strcmp(argv[i],"-centre")==0)
		{
			centre();
		}
		else
		if(strcmp(argv[i],"-normalise")==0)
		{
			normalise();
		}
		else
		if(strcmp(argv[i],"-countClusters")==0)
		{
			float	thr=atof(argv[++i]);
			if(data==NULL)
				printf("ERROR: no texture data available\n");
			countClusters(data,thr);
		}
		else
		if(strcmp(argv[i],"-normal")==0)	// surface normal vectors
		{
			if(data==NULL)
				data=(float*)calloc(3*np,sizeof(float));
			normal(data);
		}
		else
		if(strcmp(argv[i],"-randverts")==0)	// randverts n
		{
			int	n=atoi(argv[++i]);
			randverts(n);
		}
		else
		if(strcmp(argv[i],"-h")==0)	// randverts n
		{
			printHelp();
			return 0;
		}
		else
		if(strcmp(argv[i],"-v")==0)	// turn on verbose mode
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
	
	free(p);
	free(t);
	if(data)
		free(data);

	return 0;
}