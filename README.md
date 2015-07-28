meshgeometry
============

Mesh geometry tools. Several commands can be chained together, with the same operation
appearing many times in the same command line. They are processed sequentially.

I/O commands         |Description
---------------------|----------------
-iformat format_name |Force input format (needs to precede imesh)
-oformat format_name |Force output format (needs to precede omesh)
-i filename          |Input file (also accepts -imesh)
-o filename          |Output file (also accepts -omesh)
-odata filename      |Output data

Geometry commands                                |Description
-------------------------------------------------|----------------
-absgi                                           |Compute absolute gyrification index
-add filename                                    |Add mesh at filename to the current mesh
-area                                            |Surface area
-areaMap                                         |Compute surface area per vertex
-average n_meshes path1 path2 ... pathn          |Compute an average of n_meshes all of the same topology
-barycentricProjection reference_mesh            |Print barycentric coordinates for each vertex in reference_mesh
-checkOrientation                                |Check that normals point outside
-centre                                          |Move the mesh's barycentre to the origin
-countClusters  value                            |Count clusters in texture data
-curv                                            |Compute curvature
-depth                                           |Compute sulcal depth
-drawSurface path orientation                    |draw surface in tiff format, orientation is 'lat' or 'med'
-euler                                           |Print Euler characteristic
-fixFlip                                         |Detect flipped triangles and fix them
-fixSmall                                        |Detect triangles with an angle >160
-flip                                            |Flip normals degrees and fix them
-foldLength                                      |Compute total fold length
-h                                               |Help
-icurv number_of_iterations                      |Integrated curvature
-laplaceSmooth lambda number_of_iterations       |Laplace smoothing
-lissencephalic                                  |Smooth valleys and hills, not the coast
-level level_value                               |Adds new vertices (and triangles) to the edges that cross level_value in the vertex data (f.ex., mean curvature)
-addVal                                          |Add value data
-subVal                                          |Subtract value from data
-multVal                                         |Multiply data time value
-divVal                                          |Divide data by value
-clip min max                                    |Clip data values to the interval [min,max]\n\
-max                                             |Maximum data value
-mean                                            |Mean data value
-min                                             |Minimum data value
-normal                                          |Mesh normal vectors
-normalise                                       |Place all vertices at distance 1 from the origin
-randverts number_of_vertices                    |Generate homogeneously distributedrandom vertices over the mesh
-relax filename                                  |Relax current mesh to mesh at filename (both meshes have the same topology)
-resample smooth_mesh reference_mesh             |Resample the mesh to match the vertices and the topology of the argument mesh
-rotate x y z                                    |Rotate with angles x, y and z (in degrees)
-scale scale_value                               |Multiply each vertex by 'scale'
-size                                            |Display mesh dimensions
-stereographic                                   |Stereographic projection
-subdivide                                       |Subdivide the mesh with one iteration of Kobbelt's sqrt(3) algorithm
-taubinSmooth lambda mu number_of_iterations     |Taubin Smoothing
-smoothData lambda number_of_iterations          |Laplace smoothing of data, lambda=0 -> no smoothing, lambda=1 -> each vertex value to neighbour's average
-threshold value 0:down/1:up                     |Threshold texture data
-tris                                            |Display number of triangles
-v                                               |Verbose mode
-verts                                           |Display number of vertices
-volume                                          |Compute mesh volume

Meshgeometry can read and write several formats, based on the file extension:

Extension                                        |Description
-------------------------------------------------|---------------------------
.orig, .pial, .white, .inflated, .sphere, .reg   |Freesurfer meshes
.curv, .sulc, .sratio                            |Freesurfer data
.mesh                                            |BrainVisa meshe
.txt                                             |RT's mesh plain text format
.float                                           |Raw float data
.txt1                                            |RT's data format
.bin                                             |n-e-r-v-o-u-s system web binary mesh
.wrl, .obj, .ply, .stl, .smesh, .off, .vtk       |Other mesh formats
