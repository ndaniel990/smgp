#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

using namespace std;
using Viewer = igl::viewer::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
unsigned int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
unsigned int resolutionX = 20;
unsigned int resolutionY = 20;
unsigned int resolutionZ = 20;

// Intermediate result: grid points, at which the implicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

// additional vector
Eigen::VectorXd constrained_values_nonZero;

// accelerated datastructure
vector<vector<vector<vector<int> > > > accelGrid;
vector<vector<vector<vector<int> > > > accelGridconstraints;

// eps
double eps;

// use object or axis aligned bounding box
bool PCAbool = true;

// mesh enlargment
double enlargement = 1.03;

// resolution x,y,z
int resX, resY, resZ;
int resXcon, resYcon, resZcon;

// offset
Eigen::RowVectorXd offset;

// name of mesh
string strName;

// matrix containing neighbors of grid point
Eigen::MatrixXd neighbours;

// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines. resize(0, 6);
    grid_values.resize(0);
    V. resize(0, 3);
    F. resize(0, 3);
    FN.resize(0, 3);

	Eigen::MatrixXd tempRotP = P;
	Eigen::MatrixXd transform;
	Eigen::RowVectorXd mean;

	if (PCAbool) {
		//PCA
		mean = P.colwise().mean();
		Eigen::MatrixXd centerP = P.rowwise() - mean;
		Eigen::MatrixXd covariance = centerP.adjoint() * centerP;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(covariance);
		Eigen::MatrixXd eigenVectors = eigen.eigenvectors();
		if (eigenVectors.determinant() < 0)
			eigenVectors = eigenVectors * (-1);
		transform = eigenVectors.rightCols(3);
		tempRotP = centerP * transform;
		//end PCA
	}

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;

	bb_min = tempRotP.colwise().minCoeff();
	bb_max = tempRotP.colwise().maxCoeff();
	// Bounding box dimensions
	Eigen::RowVector3d dim = bb_max - bb_min;
	dim *= enlargement;

    // Grid spacing
    const double dx = dim[0] / (double)(resolutionX - 1);
    const double dy = dim[1] / (double)(resolutionY - 1);
    const double dz = dim[2] / (double)(resolutionZ - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolutionX * resolutionY * resolutionZ, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min*enlargement + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }

	if (PCAbool) {
		grid_points = grid_points * transform.adjoint();
		grid_points = grid_points.rowwise() + mean;
	}

}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x<resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionX - 1) {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1) {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1) {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

void createAccelGrid() {
	//accelGrid.resize(0, 3);
	/*
	Eigen::RowVectorXd mean = P.colwise().mean();
	Eigen::MatrixXd centerP = P.rowwise() - mean;
	Eigen::MatrixXd covariance = centerP.adjoint() * centerP;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(covariance);
	Eigen::MatrixXd eigenVectors = eigen.eigenvectors();
	Eigen::MatrixXd transform = eigenVectors.rightCols(3);
	Eigen::MatrixXd tempRotP = centerP * transform;
	*/
	Eigen::RowVector3d bb_min, bb_max;


	bb_min = P.colwise().minCoeff();// tempRotP.colwise().minCoeff();
	bb_max = P.colwise().maxCoeff();// tempRotP.colwise().maxCoeff();

	Eigen::RowVector3d dim = bb_max - bb_min;
	/*
	const double dx = dim[0] / (double)(eps - 1);
	const double dy = dim[1] / (double)(eps - 1);
	const double dz = dim[2] / (double)(eps - 1);
	*/
	resX = dim[0] / eps + 1;
	resY = dim[1] / eps + 1;
	resZ = dim[2] / eps + 1;
	//accelGrid.resize(resX * resY * resZ, 1); //only one value per cell possible??
	//accelGrid = accelGrid.setOnes() * (-1);
	//accelGrid = vector<vector<int> >(resX*resY*resZ, vector<int>(1, -1));
	accelGrid = vector<vector<vector<vector<int> > > >(resX, vector<vector<vector<int> > >(resY, vector<vector<int> >(resZ, vector<int>(1,-1))));

	// Create each gridpoint
	/*
	for (unsigned int x = 0; x < eps; ++x) {
		for (unsigned int y = 0; y < eps; ++y) {
			for (unsigned int z = 0; z < eps; ++z) {
				// Linear index of the point at (x,y,z)
				int index = x + eps * (y + eps * z);
				// 3D point at (x,y,z)
				accelGrid.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
			}
		}
	}
	*/
	/*
	offset = bb_min.cwiseAbs();
	for (int i = 0; i < P.rows(); i++) {
		Eigen::VectorXd temp = P.row(i) + offset;
		temp /= eps;
		int index = temp.x() + resX * (temp.y() + resY * temp.z());
		if (index >= accelGrid.size())
			index = 1;
		if (accelGrid[index][0] == -1)
			accelGrid[index][0] = i;
		else
			accelGrid[index].push_back(i);
	}
	*/
	offset = bb_min;
	for (int i = 0; i < P.rows(); i++) {
		Eigen::Vector3d t = P.row(i) - offset;
		t /= eps;
		if (accelGrid[t.x()][t.y()][t.z()][0] == -1)
			accelGrid[t.x()][t.y()][t.z()][0] = i;
		else
			accelGrid[t.x()][t.y()][t.z()].push_back(i);
	}

}

void addToAccelGridconstraints(int x, int y, int z, int i) {
	if (x < resXcon && y < resYcon && z < resZcon) {
		if (accelGridconstraints[x][y][z][0] == -1)
			accelGridconstraints[x][y][z][0] = i;
		else
			accelGridconstraints[x][y][z].push_back(i);
	}

}

void createAccelGridconstraints() {
	Eigen::RowVector3d bb_min, bb_max;


	bb_min = constrained_points.colwise().minCoeff();
	bb_max = constrained_points.colwise().maxCoeff();

	Eigen::RowVector3d dim = bb_max - bb_min;

	resXcon = dim[0] / eps + 1;
	resYcon = dim[1] / eps + 1;
	resZcon = dim[2] / eps + 1;
	//cout << resXcon << " " << resYcon << " " << resZcon << endl;

	accelGridconstraints = vector<vector<vector<vector<int> > > >(resXcon, vector<vector<vector<int> > >(resYcon, vector<vector<int> >(resZcon, vector<int>(1, -1))));
	for (int i = 0; i < constrained_points.rows(); i++) {
		Eigen::Vector3d t = constrained_points.row(i) - offset;
		t /= eps;
		//add point to every corner of cube
		addToAccelGridconstraints(t.x(), t.y(), t.z(), i);
		/*
		addToAccelGridconstraints(t.x() + 1, t.y(), t.z(), i);
		addToAccelGridconstraints(t.x() + 1, t.y() + 1, t.z(), i);
		addToAccelGridconstraints(t.x() + 1, t.y() + 1, t.z() + 1, i);
		addToAccelGridconstraints(t.x() + 1, t.y(), t.z() + 1, i);
		addToAccelGridconstraints(t.x(), t.y() + 1, t.z(), i);
		addToAccelGridconstraints(t.x(), t.y() + 1, t.z() + 1, i);
		addToAccelGridconstraints(t.x(), t.y(), t.z() + 1, i);
		*/
	}
}

bool getProximity(int p, int q) {
	//fast search
	/*
	double distSq = (constrained_points.row(p) - constrained_points.row(q)).squaredNorm();
	Eigen::RowVectorXd pointofinterest = constrained_points.row(q) + offset;
	//int index = pointofinterest.x + eps * (pointofinterest.y + eps * pointofinterest.z);
	for (int i = pointofinterest.x()-1; i < pointofinterest.x()+1; i++){
		i = (i < 0) ? 0 : i;
		i = (i >= resX) ? resX - 1 : i;
		for (int j = pointofinterest.y()-1; j < pointofinterest.y()+1; j++){
			j = (j < 0) ? 0 : j;
			j = (j >= resY) ? resY - 1 : j;
			for (int k = pointofinterest.z()-1; k < pointofinterest.z()+1; k++){
				k = (k < 0) ? 0 : k;
				k = (k >= resZ) ? resZ - 1 : k;
				int lookup = i + resX * (j + resY * k);
				if (accelGrid[lookup][0] != -1) {
					for (int p = 0; p < accelGrid[lookup].size(); p++) {
						//cout << P.row(accelGrid[lookup][p]) << endl;
						//cout << pointofinterest << endl;
						if ((pointofinterest - P.row(accelGrid[lookup][p])).squaredNorm() < distSq)
							return false;
					}
				}
			}
		}
	}
	return true;
	*/

	double distSq = (constrained_points.row(p) - constrained_points.row(q)).squaredNorm();
	Eigen::RowVectorXd pointofinterest = constrained_points.row(q) - offset;
	pointofinterest /= eps;
	for (int x = pointofinterest.x() - 1; x <= pointofinterest.x() + 1; x++) {
		if (x < 0 || x >= resX) break;
		for (int y = pointofinterest.y() - 1; y <= pointofinterest.y() + 1; y++) {
			if (y < 0 || y >= resY) break;
			for (int z = pointofinterest.z() - 1; z <= pointofinterest.z() + 1; z++) {
				if (z < 0 || z >= resZ) break;
				if (accelGrid[x][y][z][0] != -1) {
					for (int p = 0; p < accelGrid[x][y][z].size(); p++) {
						//cout << P.row(accelGrid[lookup][p]) << endl;
						//cout << pointofinterest << endl;
						if ((constrained_points.row(q) - P.row(accelGrid[x][y][z][p])).squaredNorm() < distSq) {
							//cout << "false" << endl;
							return false;
						}
					}
				}
			}
		}
	}
	return true;

	/*
	//naive for loop
	double distSq = (constrained_points.row(p) - constrained_points.row(q)).squaredNorm();
	for (int i = 0; i < P.rows(); i++) {
		if ((P.row(i)-constrained_points.row(q)).squaredNorm() < distSq){
			return false;
		}

	}
	return true;
	*/
}

void getProximityList(int i) {
	Eigen::RowVectorXd data = grid_points.row(i);
	//fast search
	/*
	Eigen::MatrixXd result;
	result.setZero(constrained_points.rows(), 3);
	constrained_values_nonZero.setZero(constrained_points.rows(), 1);
	int c = 0;
	double o = 0;
	for (int i = x.x() - o; i < x.x() + o; i++) {
		i = (i < 0) ? 0 : i;
		i = (i >= resX) ? resX - 1 : i;
		for (int j = x.y() - o; j < x.y() + o; j++) {
			j = (j < 0) ? 0 : j;
			j = (j >= resY) ? resY - 1 : j;
			for (int k = x.z() - o; k < x.z() + o; k++) {
				k = (k < 0) ? 0 : k;
				k = (k >= resZ) ? resZ - 1 : k;
				int lookup = i + resX * (j + resY * k);
				if (accelGrid[lookup][0] != -1) {
					for (int p = 0; p < accelGrid[lookup].size(); p++) {
						result.row(c) = c
					}
				}
			}
		}
	}
	constrained_values_nonZero = constrained_values_nonZero.block(0, 0, c, 1);
	return result.block(0, 0, c, 3);
	*/

	Eigen::MatrixXd result;
	result.setZero(constrained_points.rows(), 3);
	constrained_values_nonZero.setZero(constrained_points.rows(), 1);
	int c = 0;

	Eigen::RowVectorXd pointofinterest = data - offset;
	pointofinterest /= eps;
	double searchRad = wendlandRadius / eps + 1;
	//to ensure we visit all cells which could contain points <= wendlandRadius
	//searchRad *= 2;
	for (int x = pointofinterest.x() - searchRad; x <= pointofinterest.x() + searchRad; x++) {
		if (x < 0 || x >= resXcon) continue;
		for (int y = pointofinterest.y() - searchRad; y <= pointofinterest.y() + searchRad; y++) {
			if (y < 0 || y >= resYcon) continue;
			for (int z = pointofinterest.z() - searchRad; z <= pointofinterest.z() + searchRad; z++) {
				if (z < 0 || z >= resZcon) continue;
				if (accelGridconstraints[x][y][z][0] != -1) {
					for (int p = 0; p < accelGridconstraints[x][y][z].size(); p++) {
						//cout << P.row(accelGrid[lookup][p]) << endl;
						//cout << pointofinterest << endl;
						if ((data - constrained_points.row(accelGridconstraints[x][y][z][p])).squaredNorm() <= powf(wendlandRadius, 2)) {
							result.row(c) = constrained_points.row(accelGridconstraints[x][y][z][p]);
							constrained_values_nonZero.row(c) = constrained_values.row(accelGridconstraints[x][y][z][p]);
							c++;
						}
					}
				}
			}
		}
	}
	constrained_values_nonZero.conservativeResize(c, 1);
	result.conservativeResize(c,3);
	//cout << result << endl;
	neighbours = result;

	//naive 
	/*
	Eigen::MatrixXd result;
	result.setZero(constrained_points.rows(), 3);
	constrained_values_nonZero.setZero(constrained_points.rows(), 1);
	int c = 0;
	for (int i = 0; i < constrained_points.rows(); i++) {
		if ((data - constrained_points.row(i)).squaredNorm() <= powf(wendlandRadius, 2)) {
			result.row(c) = constrained_points.row(i);
			constrained_values_nonZero.row(c) = constrained_values.row(i);
			c++;
		}
	}
	constrained_values_nonZero = constrained_values_nonZero.block(0, 0, c, 1);
	return result.block(0, 0, c, 3);
	*/
}

double wendland(double r) {
	return powf(1 - r / wendlandRadius, 4) * (4 * r / wendlandRadius + 1);
}

int deg2coef() {
	switch (polyDegree)
	{
	case 0:
		return 1;
		break;
	case 1:
		return 4;
		break;
	case 2:
		return 10;
		break;
	default:
		//not allowed
		return -1;
		break;
	}

}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
        viewer.data.clear();
        viewer.core.align_camera_position(P);
        viewer.data.point_size = 11;
        viewer.data.add_points(P, Eigen::RowVector3d(0,0,0));
    }

    if (key == '2') {
        // Show all constraints
        viewer.data.clear();
        viewer.core.align_camera_position(P);

		// Add your code for computing auxiliary constraint points here
		//N.normalize(); //necessary? does this work? might take norm over wrong dimension
		constrained_points.setZero(P.rows() * 3, 3);
		constrained_values.setZero(P.rows() * 3, 1);

		eps = (P.colwise().maxCoeff() - P.colwise().minCoeff()).norm()*0.01;
		double varEps = eps;
		//create accelerated datastructure
		createAccelGrid();
		for (int i = 0; i < P.rows(); i++) {
			constrained_points.row(i) = P.row(i); //first entry is point on surface
			constrained_values(i) = 0;
			//need to check if p_i is closest point to p_i+N
			constrained_points.row(i + P.rows()) = P.row(i) + varEps*N.row(i); //second entry is point outside of surface
			while (!getProximity(i, i + P.rows())) {
				varEps *= 0.5;
				constrained_points.row(i + P.rows()) = P.row(i) + varEps*N.row(i);
			}
			constrained_values(i + P.rows()) = varEps;
			//reset eps
			varEps = eps;
			//need to check if p_i is closest point to p_i+2N
			constrained_points.row(i + 2 * P.rows()) = P.row(i) - varEps*N.row(i); //third entry is point inside of surface
			while (!getProximity(i, i + 2 * P.rows())) {
				varEps *= 0.5;
				constrained_points.row(i + 2 * P.rows()) = P.row(i) - varEps*N.row(i);
			}
			constrained_values(i + 2 * P.rows()) = -varEps;
			//reset eps
			varEps = eps;
		}

		cout << "Creating acceleration data structure..." << endl;
		createAccelGridconstraints();
		cout << "Done." << endl;

        // Add code for displaying all points, as above
		viewer.data.point_size = 11;

		viewer.data.add_points(constrained_points.block(0,0,P.rows(),3), Eigen::RowVector3d(0, 0, 1));
		viewer.data.add_points(constrained_points.block(P.rows(), 0, P.rows(), 3), Eigen::RowVector3d(1, 0, 0));
		viewer.data.add_points(constrained_points.block(2*P.rows(), 0, P.rows(), 3), Eigen::RowVector3d(0, 1, 0));
    }

    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
        viewer.data.clear();
        viewer.core.align_camera_position(P);

		createGrid();
        // Add your code for evaluating the implicit function at the grid points
		// Scalar values of the grid points (the implicit function values)
		grid_values.resize(resolutionX * resolutionY * resolutionZ);

		for (int i = 0; i < grid_points.rows(); i++) {
			//for every constrained point
			Eigen::MatrixXd basis, f;
			Eigen::VectorXd weights;
			//weights.setZero(constrained_points.rows());
			//basis.setZero(constrained_points.rows(),10);

			getProximityList(i);
			if (neighbours.rows() == 0)
				grid_values(i) = INT_MAX;
			else
			{

				weights.setZero(neighbours.rows(), 1);
				basis.setZero(neighbours.rows(), 10);
				for (int j = 0; j < neighbours.rows(); j++) {
					weights(j) = wendland((grid_points.row(i) - neighbours.row(j)).norm());

					basis.row(j) << 1, neighbours(j, 0), neighbours(j, 1), neighbours(j, 2),
						powf(neighbours(j, 0), 2), powf(neighbours(j, 1), 2), powf(neighbours(j, 2), 2),
						neighbours(j, 0)*neighbours(j, 1), neighbours(j, 1)*neighbours(j, 2), neighbours(j, 0)*neighbours(j, 2);
				}
				Eigen::MatrixXd A = weights.asDiagonal()*basis.block(0, 0, neighbours.rows(), deg2coef());
				//Eigen::VectorXd coef = (A.transpose() * A).ldlt().solve(A.transpose() * weights.asDiagonal()*constrained_values_nonZero);
				Eigen::VectorXd coef = A.colPivHouseholderQr().solve(weights.asDiagonal()*constrained_values_nonZero);
				//cout << coef << endl;
				Eigen::VectorXd basisBlock;
				basisBlock.setZero(10, 1);
				basisBlock << 1, grid_points(i, 0), grid_points(i, 1), grid_points(i, 2),
					powf(grid_points(i, 0), 2), powf(grid_points(i, 1), 2), powf(grid_points(i, 2), 2),
					grid_points(i, 0)*grid_points(i, 1), grid_points(i, 1)*grid_points(i, 2), grid_points(i, 0)*grid_points(i, 2);
				grid_values(i) = basisBlock.block(0, 0, 1, deg2coef()).dot(coef);
			}

		}
        // Add code for displaying points and lines
		getLines();
		grid_colors.setZero(grid_points.rows(), 3);

		// Build color map
		for (int i = 0; i < grid_points.rows(); ++i) {
			double value = grid_values(i);
			if (value < 0) {
				grid_colors(i, 1) = 1;
			}
			else {
				if (value > 0)
					grid_colors(i, 0) = 1;
			}
		}

		// Draw lines and points
		viewer.data.point_size = 8;
		viewer.data.add_points(grid_points, grid_colors);
		viewer.data.add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
			grid_lines.block(0, 3, grid_lines.rows(), 3),
			Eigen::RowVector3d(0.8, 0.8, 0.8));
		
    }

    if (key == '4') {
        // Show reconstructed mesh
        viewer.data.clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
        if (V.rows() == 0) {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data.set_mesh(V, F);
        viewer.data.show_lines = true;
        viewer.data.show_faces = true;
        viewer.data.set_normals(FN);

		string str = "../results/" + strName + ".off";
		igl::writeOFF(str, V, F);
    }

    return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
  igl::readOFF(filename,P,F,N);
  callback_key_down(viewer,'1',0);
  return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
      cout << "Usage ex2_bin <mesh.off>" << endl;
      igl::readOFF("../data/sphere.off",P,F,N);
	  strName = "sphere";
    }
	  else
	  {
		  // Read points and normals
		  igl::readOFF(argv[1],P,F,N);
	  }

    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    viewer.callback_load_mesh = callback_load_mesh;

    viewer.callback_init = [&](Viewer &v) {
        // Add widgets to the sidebar.
        v.ngui->addGroup("Reconstruction Options");
        v.ngui->addVariable("Resolution X", resolutionX);
		v.ngui->addVariable("Resolution Y", resolutionY);
		v.ngui->addVariable("Resolution Z", resolutionZ);
        v.ngui->addButton("Reset Grid", [&](){
            // Recreate the grid
            createGrid();
            // Switch view to show the grid
            callback_key_down(v, '3', 0);
        });

        // TODO: Add more parameters to tweak here...
		v.ngui->addVariable("Wendland Radius", wendlandRadius);
		v.ngui->addVariable("Mesh Enlargement", enlargement);
		v.ngui->addVariable("Poly. Deg.", polyDegree);
		v.ngui->addVariable("OABB", PCAbool);
		v.ngui->addVariable(".off Name", strName);


        v.screen->performLayout();
        return false;
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
