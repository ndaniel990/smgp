# Assignment 2

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.

![alt text][cat_constraints]

[cat_constraints]: results/cat_constraints.png "cat_constraints"

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).

![alt text][cat_grid]

![alt text][luigi_grid]

[cat_grid]: results/cat_grid.png "cat_grid"

[luigi_grid]: results/luigi_grid.png "luigi_grid"

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results).

The meshes and their respective settings can be seen in the following images. We noted that the wendland radius is heavily mesh dependent and that the polynomial degree 1 and especially 2 often give rise to artifacts due to overfitting.

![alt text][sphere_mesh]

![alt text][bunny_500_mesh]

![alt text][bunny_1000_mesh]

![alt text][cat_mesh]

![alt text][horse_mesh]

![alt text][dog_mesh]

![alt text][luigi_mesh]

[sphere_mesh]: results/sphere_mesh.png "sphere_mesh"
[bunny_500_mesh]: results/bunny_500_mesh.png "bunny_500_mesh"
[bunny_1000_mesh]: results/bunny_1000_mesh.png "bunny_1000_mesh"
[cat_mesh]: results/cat_mesh.png "cat_mesh"
[horse_mesh]: results/horse_mesh.png "horse_mesh"
[dog_mesh]: results/dog_mesh.png "dog_mesh"
[luigi_mesh]: results/luigi_mesh.png "luigi_mesh"

4) Theory question: Save your notes to assignment2/results and add a link to this page.

![alt text][theory]

[theory]: results/theory.png "theory_proof"

### Optional Tasks

1) Save your notes and add a link to this page.

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

3) Compare your MLS reconstruction results to the surfaces obtained with this method, and try to understand the differences. Report your findings.
