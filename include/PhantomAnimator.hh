#ifndef PhantomAnimator_class
#define PhantomAnimator_class

#include <iostream>
#include <ctime>
#include <functions.h>
#include <functional>
#include <iomanip>

#include <igl/readTGF.h>
#include <igl/writeTGF.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/readMESH.h>
#include <igl/readPLY.h>
#include <igl/directed_edge_parents.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/boundary_conditions.h>
#include <igl/directed_edge_parents.h>
#include <igl/directed_edge_orientations.h>
#include <igl/deform_skeleton.h>
#include <igl/forward_kinematics.h>
#include <igl/doublearea.h>
#include <igl/dqs.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/Timer.h>
#include <igl/mat_max.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>


class PhantomAnimator
{
//functions
public:
    PhantomAnimator();
    PhantomAnimator(string prefix);
    ~PhantomAnimator();

    void ReadTetMesh(string prefix);
    void PreparePhantom(string prefix);
    // bool ReadFiles(string prefix);
    string CalibrateTo(string name);
    void Animate(RotationList vQ, Vector3d root);

    MatrixXd GetU(){ return U; }
    MatrixXi GetT(){ return T; }

public:
    bool ReadProfileData(string fileName);
    map<string, int> profileIDs;
    vector<map<int, double>> jointLengths;
    vector<Vector3d> eyeR_vec, eyeL_vec;

//variables
private:
    MatrixXd C, V, U, W, Wj;
    MatrixXi BE, T;
    VectorXi P;
    MatrixXd V_calib, C_calib;
    vector<int> eyeIDs;
    // vector<int> eye2ply;
    // RotationList alignRot;
    vector<map<int, double>> cleanWeights;
  
    //tmp
    map<int, double> lengths;
};

#endif
