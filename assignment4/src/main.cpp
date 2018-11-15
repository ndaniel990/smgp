#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/dijkstra.h>
#include <igl/facet_components.h>
#include <igl/jet.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::viewer::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

// bool
bool first = true;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::viewer::ViewerCore temp3D;
igl::viewer::ViewerCore temp2D;

bool cotLap = false;
int anglePreserving = 1;

Eigen::MatrixXd color;
bool pressedFive;
bool considerOnlyBoundary;

void Redraw()
{
	viewer.data.clear();

	if (!showingUV)
	{
		viewer.data.set_mesh(V, F);
		viewer.data.set_face_based(false);

    if(UV.size() != 0)
    {
      viewer.data.set_uv(TextureResolution*UV);
	  viewer.data.show_texture = true;
    }
	}
	else
	{
		viewer.data.show_texture = false;
		viewer.data.set_mesh(UV, F);
	}

	if (pressedFive) {
		viewer.data.show_texture = false;
		viewer.data.set_colors(color);
	}
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::viewer::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	d.setZero(positions.rows() * 2);
	C.resize(positions.rows() * 2, V.rows() * 2);
	for (int i = 0; i < positions.rows(); i++) {
		d[i] = positions(i,0);
		d[i + positions.rows()] = positions(i,1);
	}
	std::vector<Eigen::Triplet<double> > c;
	for (int i = 0; i < indices.rows(); i++) {
		c.push_back(Eigen::Triplet<double>(i,indices[i],1));
		c.push_back(Eigen::Triplet<double>(indices.rows()+i,V.rows()+indices[i],1));
	}
	C.setFromTriplets(c.begin(), c.end());
}

void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;
	
	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary)
	{
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.
		igl::boundary_loop(F, fixed_UV_indices);
		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
	}
	else
	{
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
		//find most distant points
		vector<vector<int> > adjList;
		set<int> targets; //leave empty to find all
		VectorXd minDist;
		VectorXi previous;
		vector<double> maxDist;
		vector<int> maxDistIndex;
		double greatestDist = 0;
		int greatestDistInd1, greatestDistInd2;
		Eigen::VectorXd::Index index;
		if (considerOnlyBoundary) {
			vector<int> boundaryInd;
			igl::boundary_loop(F, boundaryInd);
			igl::adjacency_list(F, adjList);

			for (int i = 0; i < boundaryInd.size(); i++) {
				igl::dijkstra_compute_paths(boundaryInd[i], targets, adjList, minDist, previous);
				//for (int j = 0; j < minDist.rows(); j++)
				//	cout << minDist.row(i) << endl;
				maxDist.push_back(minDist.maxCoeff(&index));
				//cout << "greatest Dist " << minDist.maxCoeff() << endl;
				maxDistIndex.push_back(index);
			}
		}
		else { //consider all points
			igl::adjacency_list(F, adjList);

			for (int i = 0; i < V.rows(); i++) {
				igl::dijkstra_compute_paths(i, targets, adjList, minDist, previous);
				//for (int j = 0; j < minDist.rows(); j++)
				//	cout << minDist.row(i) << endl;
				maxDist.push_back(minDist.maxCoeff(&index));
				//cout << "greatest Dist " << minDist.maxCoeff() << endl;
				maxDistIndex.push_back(index);
			}
		}

		for (int i = 0; i < maxDist.size(); i++) {
			if (maxDist[i] > greatestDist) {
				greatestDist = maxDist[i];
				greatestDistInd1 = i;
				greatestDistInd2 = maxDistIndex[i];
			}
		}
		fixed_UV_indices.resize(2, 1);
		fixed_UV_indices << greatestDistInd1, greatestDistInd2;
		fixed_UV_positions.resize(2,2);
		fixed_UV_positions << -1, 0, 1, 0;
	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
		Eigen::SparseMatrix<double> tmp;
		igl::adjacency_matrix(F,tmp);
		SparseVector<double> tmpSum;
		igl::sum(tmp,1,tmpSum);
		SparseMatrix<double> tmpDiag;
		igl::diag(tmpSum,tmpDiag);
		SparseMatrix<double> U;
		U = tmp-tmpDiag;
		vector<Eigen::Triplet<double> > a;
		A.resize(2 * V.rows(), 2 * V.rows());
		for (int i = 0; i < U.outerSize(); i++) {
			for (SparseMatrix<double>::InnerIterator it(U, i); it; ++it) {
				a.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				a.push_back(Eigen::Triplet<double>(it.row() + V.rows(), it.col() + V.rows(), it.value()));
			}
		}
		A.setFromTriplets(a.begin(), a.end());
		b.setZero(2 * V.rows());
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~
		Eigen::SparseMatrix<double> L;
		if (cotLap) {
			Eigen::SparseMatrix<double> Dx, Dy;
			computeSurfaceGradientMatrix(Dx, Dy);
			VectorXd tmpVec;
			igl::doublearea(V, F, tmpVec);
			Eigen::SparseMatrix<double> Adouble;
			Adouble.resize(F.rows(), F.rows());
			vector<Eigen::Triplet<double> > tripletVec;
			for (int i = 0; i < F.rows(); i++)
				tripletVec.push_back(Eigen::Triplet<double>(i, i, tmpVec(i)));
			Adouble.setFromTriplets(tripletVec.begin(), tripletVec.end());
			Eigen::SparseMatrix<double> DxT, DyT;
			DxT = Dx.transpose();
			DyT = Dy.transpose();
			L = (DxT*Adouble*Dx + DyT*Adouble*Dy) * 0.5;
		}
		else {
			igl::cotmatrix(V, F, L);
		}

		vector<Eigen::Triplet<double> > a;
		A.resize(2 * V.rows(), 2 * V.rows());
		for (int i = 0; i < L.outerSize(); i++) {
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it) {
				a.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				a.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(a.begin(), a.end());
		b.setZero(2*V.rows());
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		Eigen::SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		VectorXd tmpVec;
		igl::doublearea(V, F, tmpVec);
		Eigen::SparseMatrix<double> Adouble;
		Adouble.resize(F.rows(),F.rows());
		vector<Eigen::Triplet<double> > tripletVec;
		for (int i = 0; i < F.rows(); i++)
			tripletVec.push_back(Eigen::Triplet<double>(i,i,tmpVec(i)));
		Adouble.setFromTriplets(tripletVec.begin(), tripletVec.end());
		A.resize(2 * V.rows(), 2 * V.rows());
		Eigen::SparseMatrix<double> tmpL, tmpR, DxT, DyT, tmpTL, tmpTR, tmpBL, tmpBR;
		DxT = Dx.transpose();
		DyT = Dy.transpose();
		tmpTL = (DxT*Adouble*Dx + DyT*Adouble*Dy);
		tmpBL = (-DyT*Adouble*Dx + DxT*Adouble*Dy);
		tmpTR = (-DxT*Adouble*Dy + DyT*Adouble*Dx);
		tmpBR = (DxT*Adouble*Dx + DyT*Adouble*Dy);
		igl::cat(1, tmpTL, tmpBL, tmpL);
		igl::cat(1, tmpTR, tmpBR, tmpR);
		igl::cat(2, tmpL, tmpR, A);
		b.setZero(2 * V.rows());
	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
		Eigen::SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		MatrixXd bigR;
		bigR.resize(F.rows(), 4);
		Eigen::VectorXd D1, D2, D3, D4;
		D1 = Dx * UV.col(0);
		D2 = Dy * UV.col(0);
		D3 = Dx * UV.col(1);
		D4 = Dy * UV.col(1);
		for (int i = 0; i < F.rows(); i++) {
			Eigen::Matrix2d J, U, S, V, R, tmpD, VTr;
			J.resize(2, 2);
			J << D1[i], D2[i], D3[i], D4[i];
			SSVD2x2(J, U, S, V);
			R.resize(2, 2);
			VTr = V.transpose();
			tmpD.resize(2, 2);
			double s = (signbit((U*VTr).determinant())) ? -1 : 1;
			tmpD << 1, 0, 0, s;
			R = U*tmpD*VTr;
			bigR.row(i) << R(0, 0), R(0, 1), R(1, 0), R(1, 1);
		}

		VectorXd tmpVec;
		igl::doublearea(V, F, tmpVec);
		Eigen::SparseMatrix<double> areaDouble;
		areaDouble.resize(F.rows(), F.rows());
		vector<Eigen::Triplet<double> > tripletVec;
		for (int i = 0; i < F.rows(); i++)
			tripletVec.push_back(Eigen::Triplet<double>(i, i, tmpVec(i)));
		areaDouble.setFromTriplets(tripletVec.begin(), tripletVec.end());
		A.resize(2 * V.rows(), 2 * V.rows());
		Eigen::SparseMatrix<double> tmpL, tmpR, tmpTL, tmpTR, tmpBL, tmpBR;
		tmpTL = (Dx.transpose()*areaDouble*Dx*0.5 + Dy.transpose()*areaDouble*Dy*0.5);
		tmpBL.resize(V.rows(),V.rows());
		tmpTR.resize(V.rows(),V.rows());
		tmpBR = (Dx.transpose()*areaDouble*Dx*0.5 + Dy.transpose()*areaDouble*Dy*0.5);
		igl::cat(1, tmpTL, tmpBL, tmpL);
		igl::cat(1, tmpTR, tmpBR, tmpR);
		igl::cat(2, tmpL, tmpR, A);
		VectorXd rhsTop, rhsBot;
		rhsTop = Dx.transpose()*areaDouble*bigR.col(0)*0.5 + Dy.transpose()*areaDouble*bigR.col(1)*0.5;
		rhsBot = Dx.transpose()*areaDouble*bigR.col(2)*0.5 + Dy.transpose()*areaDouble*bigR.col(3)*0.5;
		b.setZero(2 * V.rows());
		igl::cat(1, rhsTop, rhsBot, b);
	}

	// Solve the linear system. 
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail
	Eigen::SparseMatrix<double> tmpMatL, tmpMatR;
	Eigen::SparseMatrix<double, ColMajor> zeroMat, lhs;
	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
	zeroMat.resize(C.rows(), C.rows());
	igl::cat(1, A, C, tmpMatL);
	C = C.transpose();
	igl::cat(1, C, zeroMat, tmpMatR);
	igl::cat(2, tmpMatL, tmpMatR, lhs);
	VectorXd rhs, x;
	rhs.resize(b.rows() + d.rows(), 1);
	rhs << b, d;
	solver.analyzePattern(lhs);
	solver.factorize(lhs);
	x = solver.solve(rhs);

	// The solver will output a vector
	UV.resize(V.rows(), 2);
	UV.col(0) = x.block(0, 0, V.rows(), 1);
	UV.col(1) = x.block(V.rows(), 0, V.rows(), 1);
}

void computeDistortion() {
	double distortion;
	Eigen::SparseMatrix<double> Dx, Dy;
	computeSurfaceGradientMatrix(Dx, Dy);
	Eigen::Matrix2d I;
	I.setIdentity(); //does this work as intended?
	Eigen::VectorXd C;
	C.resize(F.rows());
	Eigen::VectorXd D1, D2, D3, D4;
	D1 = Dx * UV.col(0);
	D2 = Dy * UV.col(0);
	D3 = Dx * UV.col(1);
	D4 = Dy * UV.col(1);
	for (int i = 0; i < F.rows(); i++) {
		Eigen::Matrix2d J, Jtr, U, S, V, R, tmpD, VTr;
		J << D1[i], D2[i], D3[i], D4[i];
		//cout << J << endl << endl;
		if (anglePreserving == 3) {
			distortion = J.determinant()*J.determinant();
		}
		else {
			if (anglePreserving == 1) {
				Jtr = J.transpose();
				Eigen::Matrix2d tmp = I*J.trace();
				J = J + Jtr;
				J = J - tmp;
				//cout << J << endl << endl;
			}
			else if (anglePreserving == 2) {
				SSVD2x2(J, U, S, V);
				R.resize(2, 2);
				VTr = V.transpose();
				tmpD.resize(2, 2);
				double s = (signbit((U*VTr).determinant())) ? -1 : 1;
				tmpD << 1, 0, 0, s;
				R = U*tmpD*VTr;
				J = J - R;
			}

			distortion = J.squaredNorm();
		}

		//cout << distortion << endl;
		C.row(i) << distortion;
	}

	//find colors according to distortion
	color.resize(F.rows(), 3);
	double smallest = C.minCoeff();
	Eigen::VectorXd s;
	s.resize(F.rows());
	s.setConstant(smallest);
	C = C - s;
	double largest = C.maxCoeff();
	C = C / largest;
	Eigen::VectorXd ones;
	ones.resize(F.rows());
	ones.setOnes();
	color.col(0) << ones;
	color.col(1) << ones - C;
	color.col(2) << ones - C;
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
		computeParameterization(key);
		first = true;
		pressedFive = false;
		break;
	case '4':
		pressedFive = false;
		if (first) {
			computeParameterization('3');
			first = false;
		}
		else
			computeParameterization(key);
		break;
	case '5':
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
		computeDistortion();
		first = true;
		pressedFive = true;
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core;
      viewer.core = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core;
        viewer.core = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

void init_core_states()
{
  // save initial viewer core state
  temp3D = viewer.core;
  temp2D = viewer.core;
  temp2D.orthographic = true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core.align_camera_position(V);
  showingUV = false;

  return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
  load_mesh(filename);
  init_core_states();
  return true;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

  viewer.callback_init = [&](Viewer &v) {
    // Add widgets to the sidebar.
    v.ngui->addGroup("Parmaterization");
    v.ngui->addVariable("Free boundary",freeBoundary);
    // TODO: Add more parameters to tweak here...
	v.ngui->addVariable("consider only boundary",considerOnlyBoundary);
	v.ngui->addVariable("Cot. Laplacian", cotLap);
	v.ngui->addGroup("Preservation");
	v.ngui->addVariable("Angle(1)/Length(2)/Area(3)", anglePreserving);

    viewer.screen->performLayout();

    init_core_states();

    return false;
  };

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_load_mesh = callback_load_mesh;

  viewer.launch();
}
