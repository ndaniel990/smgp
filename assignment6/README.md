# Assignment 6: Final Project

## Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow.

Computation is split into 2 parts. When the mesh is loaded we call the computePrefactor() function which computes all the information which is not heat source dependent,
as soon as the user clicks on the mesh the solve() function is called which solves for \phi and sets the UV coordinates. We will now discuss both functions in more detail.


### computePrefactor():

We first compute the mass and cotangent matrix followed by the timestep which is suggested from the paper to be the square of the average edge length. Following the computation of the paper we compute P = (A - t&ast;L)
where A is the mass matrix, t is the square of the average edge length and L is the cotangent matrix. P is then decomposed using the simplicial cholesky solver from Eigen.

Furthermore we decompose the cotangent matrix using again the simplicial cholesky solver, we compute the gradient matrix and its inverse which is the same as the divergence as it will be used in the solve() function. We also compute the area of the triangles as it will be needed to obtain the correct scaling of \phi as well as a texture vector and a scaling factor for the texture.


### solve():

This method is called when the user clicks somewhere on the mesh.

We first construct the vector U(0) which is filled with zeros except at the vertex position where the user clicked on, which is set to 1. Using the cholesky solver from P we solve U(dt) = P\\U(0),
we then multiply the gradient matrix with U(dt), we call the result X which is a vector of size 3&ast;F.rows(). We then reshape X to F.rows() by 3 such that we can normalize for every face and then convert it back to a vector of size 3&ast;F.rows() and invert the sign as stated in the paper in Algorithm 1, step II.

This normalized vector \hat{X} is then multiplied with the divergence and the triangle areas to obtain a vector t = div &ast; triA &ast; \hat{X}. This vector t is then used to solve for \phi = L\\t using the second cholesky solver computed in the prefactor step. As mentioned in the paper \phi can contain negative values up to a constant so we shift it such that the smallest value is zero. This shifted \phi is then used to color the mesh using igl::jet.

For the isolines we use the UV coordinates and set them to be UV = (\phi, \phi) this will result in the textures being mapped in the correct way spreading out from the heat source. The UV coordinates are then multiplied by a scaling factor also determined in the prefactor step which relies on the length of the bounding box diagonal to ensure a correct scaling for the different kinds of mesh sizes.

### Timestep parameter discussion:

One would assume that as the timestep dt goes towards zero the accuracy would increase and converge to the true geodesic distance, this is true for the continuous case but not for discrete meshes on which we are operating. A too small timestep would result in the combinatorial distance and not the geodesic distance which we are looking for. Choosing a too big timestep however leads to a smoothed geodesic distance which can be visible at points where there is no unique shortest path to the heat source. To obtain a good balance it turns out using a constant times the average edge length squared as the timestep is a good heuristic.

### GUI:

The implemented GUI functionalities are as follows:

timestep: set the timestep manually, by default when the mesh is loaded the timestep is computed as the square of the average edge length.

isoline scaling: set the isoline scale manually, by default when the mesh is loaded the value 12 / bbDiagNorm is used, where 12 was estimated empirically.

cont. update: if this is true, the solve() function is called whenever the mouse is moved and hovering above the mesh, this can be used for realtime geodesic distance drawing following the mouse.

multiple points: if this is true, U(0) will not be cleared at the beginning of the solve() function and thus compute the geodesic distance for multiple points simultaneously.

![alt text][3holes]
![alt text][cow]
![alt text][hand]
![alt text][hand2]
![alt text][bumpyplane]
![alt text][cylinder]

cylinder using two heat sources on opposite sides:

![alt text][cylinder_2points]
![alt text][armadillo]

armadillo using 2 heat sources:

![alt text][armadillo_2points]

armadillo with average edge length squared timestep:

![alt text][armadillo_correct]

armadillo with 10 times bigger timestep, clearly visible smoothing of geodesic distance as discussed in timestep section:

![alt text][armadillo_10xbiggerstepsize]

![alt text][bunny]

bunny using 3 heat sources:

![alt text][bunny_3points]
![alt text][cactus]

cactus using 1/6 of default scale for isolines:

![alt text][cactus_02xisoscale]

cactus using 10 times the default scale for isolines:

![alt text][cactus_10xisoscale]
![alt text][woody-hi]

woody-hi using much smaller timestep:

![alt text][woody-hi_reduceddt]

woody-hi using greater timestep:

![alt text][woody-hi_increaseddt]

[3holes]: results/3holes.png "3holes"
[cow]: results/cow.png "cow"
[hand]: results/hand.png "hand"
[hand2]: results/hand2.png "hand2"
[bumpyplane]: results/bumpyplane.png "bumpyplane"
[cylinder]: results/cylinder.png "cylinder"
[cylinder_2points]: results/cylinder_2points.png "cylinder_2points"
[cylinder_2points]: results/cylinder_2points.png "cylinder_2points"
[armadillo]: results/armadillo.png "armadillo"
[armadillo_2points]: results/armadillo_2points.png "armadillo_2points"
[armadillo_correct]: results/armadillo_correct.png "armadillo_correct"
[armadillo_10xbiggerstepsize]: results/armadillo_10xbiggerstepsize.png "armadillo_10xbiggerstepsize"
[bunny]: results/bunny.png "bunny"
[bunny_3points]: results/bunny_3points.png "bunny_3points"
[cactus]: results/cactus.png "cactus"
[cactus_02xisoscale]: results/cactus_02xisoscale.png "cactus_02xisoscale"
[cactus_10xisoscale]: results/cactus_10xisoscale.png "cactus_10xisoscale"
[woody-hi]: results/woody-hi.png "woody-hi"
[woody-hi_reduceddt]: results/woody-hi_reduceddt.png "woody-hi_reduceddt"
[woody-hi_increaseddt]: results/woody-hi_increaseddt.png "woody-hi_increaseddt"
