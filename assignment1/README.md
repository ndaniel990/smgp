# Assignment 1

## Required results

### Tasks
3) Add a text dump of the content of the two data structures for the provided mesh “plane.off”.

Vertex 0 connected to: faces 5 21

Vertex 1 connected to: faces 31

Vertex 2 connected to: faces 10 26

Vertex 3 connected to: faces 0

Vertex 4 connected to: faces 14 27 30

Vertex 5 connected to: faces 13 23 29

Vertex 6 connected to: faces 1 4 17

Vertex 7 connected to: faces 2 8 18

Vertex 8 connected to: faces 6 9 12 19 22 25

Vertex 9 connected to: faces 15 30 31

Vertex 10 connected to: faces 11 26 27

Vertex 11 connected to: faces 7 21 23

Vertex 12 connected to: faces 15 29 31

Vertex 13 connected to: faces 4 5 20

Vertex 14 connected to: faces 0 1 16

Vertex 15 connected to: faces 8 10 24

Vertex 16 connected to: faces 0 2 16

Vertex 17 connected to: faces 3 8 9 18 19 24

Vertex 18 connected to: faces 3 4 6 17 19 20

Vertex 19 connected to: faces 7 12 13 22 23 28

Vertex 20 connected to: faces 11 12 14 25 27 28

Vertex 21 connected to: faces 1 2 3 16 17 18

Vertex 22 connected to: faces 5 6 7 20 21 22

Vertex 23 connected to: faces 9 10 11 24 25 26

Vertex 24 connected to: faces 13 14 15 28 29 30


Vertex 0 connected to: vertices 11 13 22

Vertex 1 connected to: vertices 9 12

Vertex 2 connected to: vertices 10 15 23

Vertex 3 connected to: vertices 14 16

Vertex 4 connected to: vertices 9 10 20 24

Vertex 5 connected to: vertices 11 12 19 24

Vertex 6 connected to: vertices 13 14 18 21

Vertex 7 connected to: vertices 15 16 17 21

Vertex 8 connected to: vertices 17 18 19 20 22 23

Vertex 9 connected to: vertices 1 4 12 24

Vertex 10 connected to: vertices 2 4 20 23

Vertex 11 connected to: vertices 0 5 19 22

Vertex 12 connected to: vertices 1 5 9 24

Vertex 13 connected to: vertices 0 6 18 22

Vertex 14 connected to: vertices 3 6 16 21

Vertex 15 connected to: vertices 2 7 17 23

Vertex 16 connected to: vertices 3 7 14 21

Vertex 17 connected to: vertices 7 8 15 18 21 23

Vertex 18 connected to: vertices 6 8 13 17 21 22

Vertex 19 connected to: vertices 5 8 11 20 22 24

Vertex 20 connected to: vertices 4 8 10 19 23 24

Vertex 21 connected to: vertices 6 7 14 16 17 18

Vertex 22 connected to: vertices 0 8 11 13 18 19

Vertex 23 connected to: vertices 2 8 10 15 17 20

Vertex 24 connected to: vertices 4 5 9 12 19 20

4) Show three screenshots of the 'fandisk.off' model using 'per-face shading', 'per-vertex shading' and 'per-corner shading'. For the per-corner shading a threshold of 80° was used.

![alt text][face]

![alt text][vertex]

![alt text][corner]

[face]: results/per_face_shading.png "per-face shading"
[vertex]: results/per_vertex_shading.png "per-vertex shading"
[corner]: results/per_corner_shading.png "per-corner shading"

5) Show screenshots of the provided meshes with each connected component colored differently. Show the number of connected components and the size of each component (measured in number
of faces) for all the provided models.

Number of components: 1

Component 0 has 2496 faces.

![alt text][comp_bumpy_cube]

Number of components: 1

Component 0 has 27864 faces.

![alt text][comp_bunny]

Number of components: 1

Component 0 has 3072 faces.

![alt text][comp_cube]

Number of components: 1

Component 0 has 12946 faces.

![alt text][comp_sandisk]

Number of components: 1

Component 0 has 13500 faces.

![alt text][comp_gargo]

Number of components: 1

Component 0 has 32 faces.

![alt text][comp_plane]

Number of components: 1

Component 0 has 1280 faces.

![alt text][comp_sphere]

Number of components: 2

Component 0 has 3360 faces.

Component 1 has 2304 faces.

![alt text][coffeecup]

Number of components: 11

Component 0 has 90 faces.

Component 1 has 192 faces.

Component 2 has 192 faces.

Component 3 has 13216 faces.

Component 4 has 704 faces.

Component 5 has 1088 faces.

Component 6 has 1088 faces.

Component 7 has 1088 faces.

Component 8 has 1088 faces.

Component 9 has 736 faces.

Component 10 has 736 faces.

![alt text][honda]

[comp_bumpy_cube]: results/comp_bumpy_cube.png "bumpy_cube components"
[comp_bunny]: results/comp_bunny.png "comp_bunny components"
[comp_cube]: results/comp_cube.png "comp_cube components"
[comp_sandisk]: results/comp_sandisk.png "comp_sandisk components"
[comp_gargo]: results/comp_gargo.png "comp_gargo components"
[comp_plane]: results/comp_plane.png "comp_plane components"
[comp_sphere]: results/comp_sphere.png "comp_sphere components"
[coffeecup]: results/coffeecup.png "coffeecup components"
[honda]: results/honda.png "honda components"

6) Show screenshots of the subdivided meshes.

![alt text][bumpy_cube_sub]

![alt text][bunny_sub]

![alt text][coffeecup_sub]

![alt text][cube_sub]

![alt text][fandisk_sub]

![alt text][gargo_sub]

![alt text][honda_sub]

![alt text][plane_sub]

![alt text][sphere_lo_norm_sub]

[bumpy_cube_sub]: results/bumpy_cube_sub.png "bumpy_cube_sub"
[bunny_sub]: results/bunny_sub.png "bunny_sub"
[coffeecup_sub]: results/coffeecup_sub.png "coffeecup_sub"
[cube_sub]: results/cube_sub.png "cube_sub"
[fandisk_sub]: results/fandisk_sub.png "fandisk_sub"
[gargo_sub]: results/gargo_sub.png "gargo_sub"
[honda_sub]: results/honda_sub.png "honda_sub"
[plane_sub]: results/plane_sub.png "plane_sub"
[sphere_lo_norm_sub]: results/sphere_lo_norm_sub.png "sphere_lo_norm_sub"
