#include "meshgeometry.h"

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
            puts("WARNING: -imesh still works, but better change to -i");
            loadMesh(argv[++i],&mesh,iformat);
        }
        else
        if(strcmp(argv[i],"-omesh")==0)
        {
            puts("WARNING: -omesh still works, but better change to -o");
            saveMesh(argv[++i],&mesh,oformat);
        }
        else
        if(strcmp(argv[i],"-max")==0)
        {
            float val = maxData(&mesh);
            printf("maxData: %f\n", val);
        }
        else
        if(strcmp(argv[i],"-mean")==0)
        {
            meanData(&mesh);
        }
        else
        if(strcmp(argv[i],"-min")==0)
        {
            float val = minData(&mesh);
            printf("minData: %f\n", val);
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
        if(strcmp(argv[i],"-extrude")==0)
        {
            float   d;

            d=atof(argv[++i]);
            extrude(d,&mesh);
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
            char *path = argv[i+1];
            int optional_path=0;
            if(path[0]!='-')
            {
                optional_path=1;
                i+=1;
            }
            removeVerts(&mesh,optional_path,path);
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
                puts("ERROR: wrong arguments in -scale switch");
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
                puts("ERROR: no texture data available");
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
            checkVersion(argv[0]);
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
        puts("Done.");
    return 0;
}
