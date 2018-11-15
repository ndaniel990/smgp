#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>

#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/slice.h>
#include <igl/per_vertex_normals.h>
#include <igl/invert_diag.h>
#include <igl/grad.h>
#include <igl/avg_edge_length.h>
#include <igl/face_areas.h>
#include <igl/jet.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/repdiag.h>
#include <igl/diag.h>

#include "Lasso.h"
#include "Colors.h"

//activate this for alternate UI (easier to debug)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::viewer::Viewer;

igl::viewer::Viewer viewer;

Eigen::VectorXi freeVertices;

//vertex array, #V x3
Eigen::MatrixXd V(0,3);
//face array, #F x3
Eigen::MatrixXi F(0,3);

//mouse interaction
enum MouseMode { SELECT, TRANSLATE, ROTATE, NONE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0,1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0,1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0,1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0,3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0,3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0,0,0);
Eigen::Vector4f rotation(0,0,0,1.);

//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

double dt, olddt, h;

//function declarations (see below for implementation)
bool solve(igl::viewer::Viewer& viewer);
void get_new_handle_locations();
Eigen::Vector3f computeTranslation (igl::viewer::Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
Eigen::Vector4f computeRotation(igl::viewer::Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
void compute_handle_centroids();
Eigen::MatrixXd readMatrix(const char *filename);

bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier);
bool callback_mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(igl::viewer::Viewer& viewer, int button, int modifier);
bool callback_pre_draw(igl::viewer::Viewer& viewer);
bool callback_key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers);
void onNewHandleID();
void applySelection();

Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > colpiv, colpivL;
Eigen::SparseMatrix<double> L, M, A, G, divM, diagA, tri_A;
Eigen::MatrixXd d;
Eigen::VectorXd T, Tdt, gradT, nsA;
Eigen::SparseVector<double> dA;
Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> gTexture;
double bbDiagLen, scalingFactor;
bool realtime = false, multiplePoints = false;

void computePrefactor(double v) {
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::cotmatrix(V, F, L);

	dt = v;
	A = M - dt*L;

	colpiv.compute(A);

	igl::doublearea(V, F, nsA);
	dA = nsA.sparseView();
	igl::diag(dA, diagA);
	igl::repdiag(diagA, 3, tri_A);

	igl::grad(V, F, G);
	divM = G.transpose();
	colpivL.compute(L);

	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture(10, 1);
	texture.setConstant(255);
	for (int i = 0; i < 2; i++)
		texture(i,0) = 0;

	gTexture = texture;
	bbDiagLen = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
	scalingFactor = 12 / bbDiagLen;

	viewer.data.set_texture(gTexture, gTexture, gTexture);
	T.setZero(V.rows());
	olddt = dt;
}

bool solve(igl::viewer::Viewer& viewer)
{
	if (olddt != dt)
		computePrefactor(dt);

	if (!multiplePoints)
		T.setZero(V.rows());

	for (int i = 0; i < handle_vertices.rows(); i++) {
		T[handle_vertices[i]] = 1;
	}

	Tdt = colpiv.solve(T);

	MatrixXd X = Map<const MatrixXd>((G*Tdt).eval().data(), F.rows(), 3);
	VectorXd Xstraight = Map<const MatrixXd>(X.rowwise().normalized().eval().data(), 3*F.rows(), 1);

	Eigen::VectorXd tmp = -divM*tri_A*Xstraight, phi;

	phi = colpivL.solve(tmp);
	double smallestVal = phi.minCoeff();
	Eigen::VectorXd addConst;
	addConst.setZero(phi.rows(), 1);
	addConst.setConstant(abs(smallestVal));
	phi = phi + addConst;

	Eigen::MatrixXd c;
	/*
	c.setZero(phi.rows(), 3);
	addConst.setConstant(1);
	phi = phi / phi.maxCoeff();
	c.col(0) = addConst - phi;
	c.col(1) = phi;
	*/
	igl::jet(phi, true, c);
	Eigen::MatrixXd UV(V.rows(), 2);
	UV << phi, phi;
	//cout << " bbVeclen: " << bbDiagLen << endl;
	viewer.data.set_uv(scalingFactor * UV);
	viewer.data.set_colors(c);

	return true;
};

void get_new_handle_locations()
{
  int count = 0;
  for (long vi = 0; vi < V.rows(); ++vi)
    if (handle_id[vi] >= 0)
    {
      Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
      if (handle_id[vi] == moving_handle) {
        if (mouse_mode == TRANSLATE)
          goalPosition += translation;
        else if (mouse_mode == ROTATE) {
          goalPosition -= handle_centroids.row(moving_handle).cast<float>();
          igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPosition.data());
          goalPosition += handle_centroids.row(moving_handle).cast<float>();
        }
      }
      handle_vertex_positions.row(count++) = goalPosition.cast<double>();
    }
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  viewer.data.clear();
  viewer.data.set_mesh(V, F);

  viewer.core.align_camera_position(V);
  handle_id.setConstant(V.rows(), 1, -1);

  // Initialize selector
  lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

  selected_v.resize(0,1);

  h = igl::avg_edge_length(V, F);
  computePrefactor(h*h);
  return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
  load_mesh(filename);
  return true;
}

int main(int argc, char *argv[])
{
  if(argc != 2) {
    cout << "Usage assignment5_bin mesh.off>" << endl;
    load_mesh("../data/woody-lo.off");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

  // Plot the mesh
  viewer.callback_key_down = callback_key_down;
  viewer.callback_init = [&](igl::viewer::Viewer& viewer)
  {

//      viewer.ngui->addGroup("Deformation Controls");

//      viewer.ngui->addVariable<MouseMode >("MouseMode",mouse_mode)->setItems({"SELECT", "TRANSLATE", "ROTATE", "NONE"});

//      viewer.ngui->addButton("ClearSelection",[](){ selected_v.resize(0,1); });
//      viewer.ngui->addButton("ApplySelection",[](){ applySelection(); });
//      viewer.ngui->addButton("ClearConstraints",[](){ handle_id.setConstant(V.rows(),1,-1); });

	  //viewer.ngui->addVariable("Smooth", boolSmooth);
	  //viewer.ngui->addVariable("DeformedSmooth", boolDeformedSmooth);
	  viewer.ngui->addVariable("timestep", dt);
	  viewer.ngui->addVariable("isoline scaling", scalingFactor);
	  viewer.ngui->addVariable("cont. update", realtime);
	  viewer.ngui->addVariable("multiple points", multiplePoints);

      viewer.screen->performLayout();
      return false;
  };

  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  //viewer.callback_pre_draw = callback_pre_draw;
  viewer.callback_load_mesh = callback_load_mesh;

  viewer.data.point_size = 10;
  viewer.core.set_rotation_type(igl::viewer::ViewerCore::ROTATION_TYPE_TRACKBALL);

  viewer.launch();
}


bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
  if (button == (int) igl::viewer::Viewer::MouseButton::Right)
    return false;

  down_mouse_x = viewer.current_mouse_x;
  down_mouse_y = viewer.current_mouse_y;

  int fid, vi;
  Vector3f bc;
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.data.model,
	  viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
  {
	  bc.maxCoeff(&vi);
	  handle_vertices.setZero(1);
	  handle_vertices[0] = F(fid, vi);
	  solve(viewer);
	  return true;
  }
  return false;

  /*
  if (mouse_mode == SELECT)
  {
    if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >=0)
      doit = true;
    else
      lasso->strokeReset();
  }
  else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
    if(vi>=0 && handle_id[vi]>=0)  //if a region was found, mark it for translation/rotation
    {
      moving_handle = handle_id[vi];
      get_new_handle_locations();
      doit = true;
    }
  }
  return doit;
  */
}

bool callback_mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y)
{
  if (!realtime)
    return false;

  int fid, vi;
  Vector3f bc;
  double x = mouse_x;
  double y = viewer.core.viewport(3) - mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.data.model,
	  viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
  {
	  bc.maxCoeff(&vi);
	  handle_vertices.setZero(1);
	  handle_vertices[0] = F(fid, vi);
	  solve(viewer);
	  return true;
  }
  return false;
  /*
  if (mouse_mode == SELECT)
  {
    lasso->strokeAdd(mouse_x, mouse_y);
    return true;
  }
  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    if (mouse_mode == TRANSLATE) {
      translation = computeTranslation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
    }
    else {
      rotation = computeRotation(viewer,
                                 mouse_x,
                                 down_mouse_x,
                                 mouse_y,
                                 down_mouse_y,
                                 handle_centroids.row(moving_handle));
    }
    get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
    solve(viewer);
    down_mouse_x = mouse_x;
    down_mouse_y = mouse_y;
#endif
    return true;

  }
  return false;
  */
}

bool callback_mouse_up(igl::viewer::Viewer& viewer, int button, int modifier)
{
	/*
	int fid, vi;
	Vector3f bc;
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.data.model,
		viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
	{
		bc.maxCoeff(&vi);
		handle_vertices.setZero(1);
		handle_vertices[0] = F(fid, vi);
		solve(viewer);
		return true;
	}
	return false;
	*/

	return false;
};


bool callback_pre_draw(igl::viewer::Viewer& viewer)
{
  // initialize vertex colors
  vertex_colors = Eigen::MatrixXd::Constant(V.rows(),3,.9);

  // first, color constraints
  int num = handle_id.maxCoeff();
  if (num == 0)
    num = 1;
  for (int i = 0; i<V.rows(); ++i)
    if (handle_id[i]!=-1)
    {
      int r = handle_id[i] % MAXNUMREGIONS;
      vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
    }
  // then, color selection
  for (int i = 0; i<selected_v.size(); ++i)
    vertex_colors.row(selected_v[i]) << 131./255, 131./255, 131./255.;
  
  //viewer.data.set_colors(vertex_colors);
  viewer.data.V_material_specular.fill(0);
  viewer.data.V_material_specular.col(3).fill(1);
  viewer.data.dirty |= viewer.data.DIRTY_DIFFUSE | viewer.data.DIRTY_SPECULAR;


  //clear points and lines
  viewer.data.set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
  viewer.data.set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));

  //draw the stroke of the selection
  for (unsigned int i = 0; i<lasso->strokePoints.size(); ++i)
  {
    viewer.data.add_points(lasso->strokePoints[i],Eigen::RowVector3d(0.4,0.4,0.4));
    if(i>1)
      viewer.data.add_edges(lasso->strokePoints[i-1], lasso->strokePoints[i], Eigen::RowVector3d(0.7,0.7,0.7));
  }

  // update the vertex position all the time
  viewer.data.V.resize(V.rows(),3);
  viewer.data.V << V;

  viewer.data.dirty |= viewer.data.DIRTY_POSITION;

#ifdef UPDATE_ONLY_ON_UP
  //draw only the moving parts with a white line
  if (moving_handle>=0)
  {
    Eigen::MatrixXd edges(3*F.rows(),6);
    int num_edges = 0;
    for (int fi = 0; fi<F.rows(); ++fi)
    {
      int firstPickedVertex = -1;
      for(int vi = 0; vi<3 ; ++vi)
        if (handle_id[F(fi,vi)] == moving_handle)
        {
          firstPickedVertex = vi;
          break;
        }
      if(firstPickedVertex==-1)
        continue;


      Eigen::Matrix3d points;
      for(int vi = 0; vi<3; ++vi)
      {
        int vertex_id = F(fi,vi);
        if (handle_id[vertex_id] == moving_handle)
        {
          int index = -1;
          // if face is already constrained, find index in the constraints
          (handle_vertices.array()-vertex_id).cwiseAbs().minCoeff(&index);
          points.row(vi) = handle_vertex_positions.row(index);
        }
        else
          points.row(vi) =  V.row(vertex_id);

      }
      edges.row(num_edges++) << points.row(0), points.row(1);
      edges.row(num_edges++) << points.row(1), points.row(2);
      edges.row(num_edges++) << points.row(2), points.row(0);
    }
    edges.conservativeResize(num_edges, Eigen::NoChange);
    viewer.data.add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9,0.9,0.9));

  }
#endif
  return false;

}

bool callback_key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
	/*
  bool handled = false;
  if (key == 'S')
  {
    mouse_mode = SELECT;
    handled = true;
  }

  if ((key == 'T') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = TRANSLATE;
    handled = true;
  }

  if ((key == 'R') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = ROTATE;
    handled = true;
  }
  if (key == 'A')
  {
    applySelection();
    callback_key_down(viewer, '1', 0);
    handled = true;
  }
  */
  viewer.ngui->refresh();
  return false;
}

void onNewHandleID()
{
  //store handle vertices too
  int numFree = (handle_id.array() == -1).cast<int>().sum();
  int num_handle_vertices = V.rows() - numFree;
  handle_vertices.setZero(num_handle_vertices);
  handle_vertex_positions.setZero(num_handle_vertices,3);
  freeVertices.setZero(numFree);

  int count = 0;
  int countFree = 0;
  for (long vi = 0; vi < V.rows(); ++vi) {
	  if (handle_id[vi] >= 0)
		  handle_vertices[count++] = vi;
	  else
		  freeVertices[countFree++] = vi;
  }

  compute_handle_centroids();
  h = igl::avg_edge_length(V, F);
  computePrefactor(h*h);
}

void applySelection()
{
  int index = handle_id.maxCoeff()+1;
  for (int i =0; i<selected_v.rows(); ++i)
  {
    const int selected_vertex = selected_v[i];
    if (handle_id[selected_vertex] == -1)
      handle_id[selected_vertex] = index;
  }
  selected_v.resize(0,1);

  onNewHandleID();
}

void compute_handle_centroids()
{
  //compute centroids of handles
  int num_handles = handle_id.maxCoeff()+1;
  handle_centroids.setZero(num_handles,3);

  Eigen::VectorXi num; num.setZero(num_handles,1);
  for (long vi = 0; vi<V.rows(); ++vi)
  {
    int r = handle_id[vi];
    if ( r!= -1)
    {
      handle_centroids.row(r) += V.row(vi);
      num[r]++;
    }
  }

  for (long i = 0; i<num_handles; ++i)
    handle_centroids.row(i) = handle_centroids.row(i).array()/num[i];

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f computeTranslation (igl::viewer::Viewer& viewer,
                                    int mouse_x,
                                    int from_x,
                                    int mouse_y,
                                    int from_y,
                                    Eigen::RowVector3d pt3D)
{
  Eigen::Matrix4f modelview = viewer.core.view * viewer.data.model;
  //project the given point (typically the handle centroid) to get a screen space depth
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core.proj,
                                      viewer.core.viewport);
  float depth = proj[2];

  double x, y;
  Eigen::Vector3f pos1, pos0;

  //unproject from- and to- points
  x = mouse_x;
  y = viewer.core.viewport(3) - mouse_y;
  pos1 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core.proj,
                        viewer.core.viewport);


  x = from_x;
  y = viewer.core.viewport(3) - from_y;
  pos0 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core.proj,
                        viewer.core.viewport);

  //translation is the vector connecting the two
  Eigen::Vector3f translation = pos1 - pos0;
  return translation;

}


//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f computeRotation(igl::viewer::Viewer& viewer,
                                int mouse_x,
                                int from_x,
                                int mouse_y,
                                int from_y,
                                Eigen::RowVector3d pt3D)
{

  Eigen::Vector4f rotation;
  rotation.setZero();
  rotation[3] = 1.;

  Eigen::Matrix4f modelview = viewer.core.view * viewer.data.model;

  //initialize a trackball around the handle that is being rotated
  //the trackball has (approximately) width w and height h
  double w = viewer.core.viewport[2]/8;
  double h = viewer.core.viewport[3]/8;

  //the mouse motion has to be expressed with respect to its center of mass
  //(i.e. it should approximately fall inside the region of the trackball)

  //project the given point on the handle(centroid)
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core.proj,
                                      viewer.core.viewport);
  proj[1] = viewer.core.viewport[3] - proj[1];

  //express the mouse points w.r.t the centroid
  from_x -= proj[0]; mouse_x -= proj[0];
  from_y -= proj[1]; mouse_y -= proj[1];

  //shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
  from_x += w/2; mouse_x += w/2;
  from_y += h/2; mouse_y += h/2;

  //get rotation from trackball
  Eigen::Vector4f drot = viewer.core.trackball_angle.coeffs();
  Eigen::Vector4f drot_conj;
  igl::quat_conjugate(drot.data(), drot_conj.data());
  igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());

  //account for the modelview rotation: prerotate by modelview (place model back to the original
  //unrotated frame), postrotate by inverse modelview
  Eigen::Vector4f out;
  igl::quat_mult(rotation.data(), drot.data(), out.data());
  igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
  return rotation;
}
