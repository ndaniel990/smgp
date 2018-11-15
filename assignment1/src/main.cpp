#include <iostream>
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
/*** insert any libigl headers here ***/
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/flip_edge.h>

using namespace std;
using Viewer = igl::viewer::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
		igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
		for (int i = 0; i < VF.size(); i++) {
			cout << "Vertex " << i << " connected to: faces ";
			for (int j = 0; j < VF[i].size(); j++) {
				cout << VF[i][j] << " ";
			}
			cout << endl;
		}
    }

    if (key == '2') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
		igl::adjacency_list(F, VV);
		for (int i = 0; i < VV.size(); i++) {
			cout << "Vertex " << i << " connected to: vertices ";
			for (int j = 0; j < VV[i].size(); j++) {
				cout << VV[i][j] << " ";
			}
			cout << endl;
		}
    }

    if (key == '3') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
		igl::per_face_normals(V, F, FN);
        // Set the viewer normals.
        viewer.data.set_normals(FN);
    }

    if (key == '4') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
		igl::per_vertex_normals(V, F, VN);
        // Set the viewer normals.
		viewer.data.set_normals(VN);
    }

    if (key == '5') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
		const double threshold = 80;
		igl::per_corner_normals(V, F, threshold, CN);
        //Set the viewer normals
		viewer.data.set_normals(CN);
    }

    if (key == '6') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
		igl::facet_components(F, cid);
		//compute number of components
		cout << "Number of components: " << cid.maxCoeff() + 1 << endl;
		//compute number of faces per component
		vector<int> facesPerComponent(cid.maxCoeff() + 1,0);
		for (int i = 0; i < F.rows(); i++) {
			facesPerComponent[cid(i)]++;
		}
		for (int i = 0; i < (cid.maxCoeff() + 1); i++) {
			cout << "Component " << i << " has " << facesPerComponent[i] << " faces." << endl;
		}
        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
		igl::jet(cid, true, component_colors_per_face);
        // Set the viewer colors
        viewer.data.set_colors(component_colors_per_face);
    }

    if (key == '7') {
        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;
        // Add your code for sqrt(3) subdivision here.
		Eigen::MatrixXd BC;
		igl::barycenter(V, F, BC);
		Vout.setZero(V.rows()+BC.rows(), V.cols());
		Vout << V, BC;
		Eigen::MatrixXi EV, FE, EF;
		igl::edge_topology(V, F, EV, FE, EF);
		Fout.setZero(F.rows()*3, F.cols());

		//Fout.resize(F.rows()*3,3);
		int c = 0, v1, v2, v, face1, face2;
		for (int i = 0; i < EV.rows(); i++) {
			if (EF(i,0) == -1 || EF(i,1) == -1) { //border edge
				//connect both vertices and barycenter
				if (EF(i, 0) == -1) {
					v1 = EV(i, 0);
					v2 = EV(i, 1);
					v = EF(i, 1) + V.rows();
					Fout.row(c) = Eigen::Vector3i(v, v2, v1);
				}
				else {
					v1 = EV(i, 0);
					v2 = EV(i, 1);
					v = EF(i, 0) + V.rows();
					Fout.row(c) = Eigen::Vector3i(v, v1, v2);
				}
				//cout << EV(i, 0) << " " << EV(i, 1) << " " << (j + V.rows()) << endl;
				c++;
			}
			else { //no border
				face1 = EF(i,0)+V.rows();
				face2 = EF(i,1)+V.rows();
				v1 = EV(i,0);
				v2 = EV(i,1);
				//cout << face1 << " " << face2 << " " << v1 << " " << v2 << endl;
				Fout.row(c) = Eigen::Vector3i(face2, face1, v1);
				//cout << Fout.row(c) << endl;
				//cout << Fout(c) << endl;
				Fout.row(c+1) = Eigen::Vector3i(face1,face2,v2);
				//cout << Fout.row(c + 1) << endl;
				c += 2;
			}

		}
		/*
		for(int i=0; i<F.rows(); i++){
		 Fout.row(3 * i) = Eigen::Vector3i(F(i, 0), F(i, 1), (i + V.rows()));
		 Fout.row(3 * i + 1) = Eigen::Vector3i(F(i, 1), F(i, 2), (i + V.rows()));
		 Fout.row(3 * i + 2) = Eigen::Vector3i(F(i, 2), F(i, 0), (i + V.rows()));
		}
		*/

		igl::adjacency_list(F, VV);
		Eigen::MatrixXd sumVec;
		sumVec.setZero(1, 3);
		for (int i = 0; i < V.rows(); i++) {
			int n = VV[i].size();
			sumVec.setZero();
			for (int j = 0; j < n; j++) {
				sumVec += V.row(VV[i][j]); //VV[i][j] is index of j-th vector adjacent to i-th vector
			}
			double a = (4 - 2 * cos((2 * M_PI) / n)) / 9;
			Vout.row(i) << (1 - a)*V.row(i) + a / n * sumVec;
		}
		// Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
    }

    return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
  igl::readOFF(filename,V,F);
  viewer.data.clear();
  viewer.data.set_mesh(V,F);
  viewer.data.compute_normals();
  viewer.core.align_camera_position(viewer.data);
  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    viewer.callback_load_mesh = callback_load_mesh;
    
    if (argc == 2)
    {
      // Read mesh
      igl::readOFF(argv[1],V,F);
      
    }
    else
    {
      // Read mesh
      igl::readOFF("../data/coffeecup.off",V,F); //bunny
    }

    callback_key_down(viewer, '6', 0);//1

    viewer.launch();
}
