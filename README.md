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

Geometry measurement                             |Description
-------------------------------------------------|----------------
-absgi                                           |Compute absolute gyrification index
-area                                            |Surface area
-checkOrientation                                |Check that normals point outside
-centre                                          |Move the mesh's barycentre to the origin
-euler                                           |Print Euler characteristic
-foldLength                                      |Compute total fold length
-size                                            |Display mesh dimensions
-tris                                            |Display number of triangles
-verts                                           |Display number of vertices
-volume                                          |Compute mesh volume
-curv                                            |Compute curvature
-depth                                           |Compute sulcal depth
-areaMap                                         |Compute surface area per vertex
-icurv number_of_iterations                      |Integrated curvature (warning: icurv changes the geometry of the mesh)
-isolatedVerts                                   |Count the number of isolated vertices in the mesh

Geometry modification                            |Description
-------------------------------------------------|----------------
-add filename                                    |Add mesh at filename to the current mesh
-average n_meshes path1 path2 ... pathn          |Compute an average of n_meshes all of the same topology
-barycentricProjection reference_mesh            |Print barycentric coordinates for each vertex in reference_mesh
-fixFlip                                         |Detect flipped triangles and fix them
-fixSmall                                        |Detect triangles with an angle >160
-flip                                            |Flip normals degrees and fix them
-laplaceSmooth lambda number_of_iterations       |Laplace smoothing
-level level_value                               |Adds new vertices (and triangles) to the edges that cross level_value in the vertex data (f.ex., mean curvature)
-normalise                                       |Place all vertices at distance 1 from the origin
-randverts number_of_vertices                    |Generate homogeneously distributedrandom vertices over the mesh
-relax filename                                  |Relax current mesh to mesh at filename (both meshes have the same topology)
-resample smooth_mesh reference_mesh             |Resample the mesh to match the vertices and the topology of the argument mesh
-removeIsolatedVerts                             |Removes isolated vertices in the mesh (if present)
-rotate x y z                                    |Rotate with angles x, y and z (in degrees)
-scale scale_value                               |Multiply each vertex by 'scale'
-stereographic                                   |Stereographic projection
-subdivide                                       |Subdivide the mesh with one iteration of Kobbelt's sqrt(3) algorithm
-taubinSmooth lambda mu number_of_iterations     |Taubin Smoothing
-translate x y z                                 |Translate mesh by x, y, and z.
-lissencephalic                                  |Smooth valleys and hills, not the coast
-normal                                          |Mesh normal vectors

Data measurement                                 |Description
-------------------------------------------------|----------------
-countClusters  value                            |Count clusters in texture data
-max                                             |Maximum data value
-mean                                            |Mean data value
-min                                             |Minimum data value

Data modification                                |Description
-------------------------------------------------|----------------
-addVal                                          |Add value data
-subVal                                          |Subtract value from data
-multVal                                         |Multiply data time value
-divVal                                          |Divide data by value
-clip min max                                    |Clip data values to the interval [min,max]\n\
-smoothData lambda number_of_iterations          |Laplace smoothing of data, lambda=0 -> no smoothing, lambda=1 -> each vertex value to neighbour's average
-threshold value 0:down/1:up                     |Threshold texture data

Other                                            |Description
-------------------------------------------------|----------------
-drawSurface colourmap path                      |draw surface in tiff format, colormap is 'grey' or 'rainbow'
-v                                               |Verbose mode
-h                                               |Help

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
.wrl, .obj, .ply, .stl, .smesh, .off, .vtk, .gii |Other mesh formats
