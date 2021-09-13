#include "PhantomAnimator.hh"
#include "G4SystemOfUnits.hh"
PhantomAnimator::PhantomAnimator() {}
PhantomAnimator::~PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string prefix)
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE);
    ReadTetMesh(prefix);
    igl::readDMAT(prefix + ".W", W);
    igl::readDMAT(prefix + ".Wj", Wj);
    
    double epsilon(1e-2);
    for (int i = 0; i < W.rows(); i++)
    {
        double sum(0);
        map<int, double> vertexWeight;
        for (int j = 0; j < W.cols(); j++)
        {
            if (W(i, j) < epsilon)
                continue;
            vertexWeight[j] = W(i, j);
            sum += W(i, j);
        }
        for (auto &iter : vertexWeight)
            iter.second /= sum;
        cleanWeights.push_back(vertexWeight);
    }


    // ReadFiles(prefix);
    ReadProfileData("../phantoms/profile.txt");
}

void PhantomAnimator::ReadTetMesh(string prefix)
{
    cout << "Read " + prefix + ".ele/node" << endl;
    ifstream ifsNode(prefix + ".node");
    if (!ifsNode.is_open())
    {
        cout << "There is no " + prefix + ".node file!" << endl;
        exit(1);
    }
    int num, tmp, a, b, c, d, id;
    double x, y, z;
    ifsNode >> num >> tmp >> tmp >> tmp;
    V.resize(num, 3);
    for (int i = 0; i < num; i++)
    {
        ifsNode >> tmp >> x >> y >> z;
        V.row(i) = RowVector3d(x, y, z) * cm;
    }
    ifsNode.close();
    ifstream ifsEle(prefix + ".ele");
    if (!ifsEle.is_open())
    {
        cout << "There is no " + prefix + ".ele file!" << endl;
        exit(1);
    }
    ifsEle >> num >> tmp >> tmp;
    T.resize(num, 5); //last col. is for organ ID
    for (int i = 0; i < num; i++)
    {
        ifsEle >> tmp >> a >> b >> c >> d >> id;
        T(i, 0) = a;
        T(i, 1) = b;
        T(i, 2) = c;
        T(i, 3) = d;
        T(i, 4) = id;
    }
    ifsEle.close();
}

void PhantomAnimator::PreparePhantom(string prefix) //there should be prefix.ply file
{
    cout << "perform BBW with skin.ply" << endl;
    MatrixXd Vply;
    MatrixXi Fply;
    if (!igl::readPLY(prefix + ".ply", Vply, Fply)) //mesh that completely covers phantom
    {
        cout << "There is no " + prefix + ".ply file!" << endl;
        exit(1);
    }
    MatrixXd boneP = GenerateBonePoints(C, BE, 1.);
    Vply.conservativeResize(Vply.rows() + boneP.rows(), 3);
    Vply.bottomRows(boneP.rows()) = boneP;
    MatrixXd VT;
    MatrixXi FT, TT;
    igl::copyleft::tetgen::tetrahedralize(Vply, Fply, "pYq", VT, TT, FT);
    cout << "<Calculate Joint Weights>" << endl;
    MatrixXd C1 = C.block(0, 0, C.rows() - 1, 3);
    if (!CalculateScalingWeights(C1, VT, TT, Wj))
        exit(1);
    igl::normalize_row_sums(Wj, Wj);
    cout << "<Calculate Bone Weights>" << endl;
    MatrixXd bc;
    VectorXi B;
    igl::boundary_conditions(VT, TT, C, VectorXi(), BE, MatrixXi(), B, bc);
    cout << bc.rows() << " " << bc.cols() << endl;
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    if (!igl::bbw(VT, TT, B, bc, bbw_data, W))
        exit(1);
    igl::normalize_row_sums(W, W);

    cout<<"calulate barycentric coodinates"<<endl;
    auto baryMap = GenerateBarycentricCoord(VT, TT, V);
    bary = GenerateBarySparse(baryMap, VT.rows());

    W = bary*W;
    Wj = bary*Wj;
    cout<<"write "+prefix+".W"<<endl;
    igl::writeDMAT(prefix+".W", W, false);
    cout<<"write "+prefix+".Wj"<<endl;
    igl::writeDMAT(prefix+".Wj", Wj, false);
}

// bool PhantomAnimator::ReadFiles(string prefix)
// {
//     cout << "Read " + prefix + ".tgf" << endl;
//     igl::readTGF(prefix + ".tgf", C, BE);

//     //perform BBW
//     if ((!igl::readDMAT(prefix + ".W", W)) || (!igl::readDMAT(prefix + ".Wj", Wj)))
//     {
//         W.resize(V.rows(), BE.rows());
//         Wj.resize(V.rows(), C.rows() - 1);

//         MatrixXd boneP = GenerateBonePoints(C, BE, 1.);
//         MatrixXd V1(V.rows() + boneP.rows(), 3);
//         V1 << V, boneP;
//         MatrixXd VT, WT_j, WT;
//         MatrixXi TT, FT;
//         cout << "<Tetrahedralization>" << endl;
//         igl::copyleft::tetgen::tetrahedralize(V1, F, "pYq", VT, TT, FT);

//         cout << "<Calculate Joint Weights>" << endl;
//         MatrixXd C1 = C.block(0, 0, C.rows() - 1, 3);
//         if (!CalculateScalingWeights(C1, VT, TT, WT_j))
//             return EXIT_FAILURE;
//         igl::normalize_row_sums(WT_j, WT_j);

//         cout << "<Calculate Bone Weights>" << endl;
//         MatrixXd bc;
//         VectorXi b;
//         igl::boundary_conditions(VT, TT, C, VectorXi(), BE, MatrixXi(), b, bc);
//         cout << bc.rows() << " " << bc.cols() << endl;
//         igl::BBWData bbw_data;
//         bbw_data.active_set_params.max_iter = 10;
//         bbw_data.verbosity = 2;
//         if (!igl::bbw(VT, TT, b, bc, bbw_data, WT))
//             return EXIT_FAILURE;
//         igl::normalize_row_sums(WT, WT);

//         //matching between tetra & ply
//         cout << "matching between tetra & ply.." << flush;
//         auto grid = GenerateGrid(VT);
//         int count(0);
//         for (int i = 0; i < V.rows(); i++)
//         {
//             int x = floor(V(i, 0) + 0.5);
//             int y = floor(V(i, 1) + 0.5);
//             int z = floor(V(i, 2) + 0.5);
//             auto key = make_tuple(x, y, z);
//             for (int n : grid[key])
//             {
//                 if (fabs(V(n, 0) - VT(i, 0)) > 0.01)
//                     continue;
//                 if (fabs(V(n, 1) - VT(i, 1)) > 0.01)
//                     continue;
//                 if (fabs(V(n, 2) - VT(i, 2)) > 0.01)
//                     continue;
//                 W.row(i) = WT.row(n);
//                 Wj.row(i) = WT_j.row(n);
//                 cout << "\rmatching between tetra & ply.." << ++count << "/" << V.rows() << flush;
//                 break;
//             }
//         }
//         cout << endl;
//         igl::writeDMAT(prefix + ".W", W, false);
//         igl::writeDMAT(prefix + ".Wj", Wj, false);
//     }

//     //eye part
//     VectorXd data = VectorXd::Zero(Wj.rows());
//     for (int i = 0; i < Wj.rows(); i++)
//     {
//         // right eye
//         if ((C.row(23) - V.row(i)).squaredNorm() < 10)
//         {
//             eyeIDs.push_back(i);
//             data(i) = -(C.row(23) - V.row(i)).squaredNorm() + 10;
//         }
//         // left eye
//         else if ((C.row(22) - V.row(i)).squaredNorm() < 10)
//         {
//             eyeIDs.push_back(i);
//             data(i) = -(C.row(22) - V.row(i)).squaredNorm() + 10;
//         }
//     }

//     double epsilon(1e-2);
//     for (int i = 0; i < W.rows(); i++)
//     {
//         double sum(0);
//         map<int, double> vertexWeight;
//         for (int j = 0; j < W.cols(); j++)
//         {
//             if (W(i, j) < epsilon)
//                 continue;
//             vertexWeight[j] = W(i, j);
//             sum += W(i, j);
//         }
//         for (auto &iter : vertexWeight)
//             iter.second /= sum;
//         cleanWeights.push_back(vertexWeight);
//     }

//     return true;
// }

bool PhantomAnimator::Initialize()
{
    //additional variables
    igl::directed_edge_parents(BE, P);

    //distance to parent joint
    for (int i = 0; i < BE.rows(); i++)
        lengths[i] = (C.row(BE(i, 0)) - C.row(BE(i, 1))).norm();

    // joint number conv.
    vector<int> i2k = {0, 1, 2, 4, 5, 6, 7, 26, 11, 12, 13, 14, 18, 19, 20, 22, 23, 24};
    map<int, int> k2i;
    for (size_t i = 0; i < i2k.size(); i++)
        k2i[i2k[i]] = i;

    //preprocessing for KINECT draw_vidata
    map<int, Quaterniond> kinectDataRot;
    vector<int> groups;
    Matrix3d tf2;
    tf2 << 0, -1, 0, 0, 0, -1, 1, 0, 0;
    for (int id : groups = {0, 1, 2, 18, 19, 20, 26})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 0, 1, 0, 0, 0, 1, 1, 0, 0;
    for (int id : groups = {22, 23, 24})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 0, 1, 0, 0, 0, -1, -1, 0, 0;
    for (int id : groups = {5, 6})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 0, -1, 0, 0, 0, 1, -1, 0, 0;
    for (int id : groups = {12, 13})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 1, 0, 0, 0, 0, 1, 0, -1, 0;
    kinectDataRot[11] = Quaterniond(tf2);
    tf2 << 1, 0, 0, 0, 0, -1, 0, 1, 0;
    kinectDataRot[4] = Quaterniond(tf2);
    tf2 << 0, 1, 0, 1, 0, 0, 0, 0, -1;
    kinectDataRot[7] = Quaterniond(tf2);
    tf2 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    kinectDataRot[14] = Quaterniond(tf2);

    // Align Arm
    map<int, Vector3d> desiredOrt;
    groups = {12, 13, 5, 6, 22, 23, 18, 19};
    for (int id : groups)
        desiredOrt[id] = Vector3d(0, 1, 0);

    // Clavicle
    desiredOrt[4] = Vector3d(1, 0, 0);
    desiredOrt[11] = Vector3d(-1, 0, 0);

    for (int i = 0; i < BE.rows(); i++)
    {
        if (desiredOrt.find(i2k[BE(i, 0)]) == desiredOrt.end())
        {
            alignRot.push_back(kinectDataRot[i2k[BE(i, 0)]]);
        }

        else // Arm, Clavicle
        {
            Vector3d v = (C.row(k2i[i2k[BE(i, 0)] + 1]) - C.row(BE(i, 0))).transpose();
            alignRot.push_back(kinectDataRot[i2k[BE(i, 0)]] * GetRotMatrix(v, desiredOrt[i2k[BE(i, 0)]]));
        }
    }

    // cout << "generating F_area, F_incident, F_area_eye matrix..." << flush;
    // VectorXd F_area;
    // igl::doublearea(V, F, F_area);
    // F_area.cwiseSqrt();
    // VectorXd W_avgSkin = ArrayXd::Zero(V.rows());

    // for (int i = 0; i < F.rows(); i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //         W_avgSkin(F(i, j)) += F_area(i);
    // }
    // W_avgSkin = W_avgSkin / W_avgSkin.sum();

    // VectorXd W_Lens(eye2ply.size());
    // for (size_t i = 0; i < eye2ply.size(); i++)
    //     W_Lens(i) = W_avgSkin(eye2ply[i]);

    // W_avgSkin = W_avgSkin.array() / W_avgSkin.sum();
    // W_Lens = W_Lens.array() / W_Lens.sum();

    cout << "done" << endl;
    V_calib = V;
    C_calib = C;
    return true;
}

string PhantomAnimator::CalibrateTo(string name)
{
    int id = profileIDs[name];
    map<int, double> calibLengths = jointLengths[id];
    Vector3d eyeL_pos = eyeL_vec[id];
    Vector3d eyeR_pos = eyeR_vec[id];

    MatrixXd jointTrans = MatrixXd::Zero(C.rows(), 3);
    int headJ(24), eyeLJ(22), eyeRJ(23);
    stringstream ss;
    for (int i = 0; i < BE.rows(); i++)
    {
        if (calibLengths.find(i) == calibLengths.end())
        {
            // Bone Tip: Head, Hands, Feets
            calibLengths[i] = lengths[i];
            jointTrans.row(BE(i, 1)) = jointTrans.row(BE(P(i), 1));
            cout << setw(3) << i << " = " << lengths[i] << endl;
            continue;
        }
        // Bone Scaling
        double ratio = calibLengths[i] / lengths[i]; // kinedct measuring / tgf based phantom measuring
        cout << setw(3) << i << " = " << lengths[i] << " -> " << calibLengths[i] << " (" << ratio * 100 << " %)" << endl;
        ss << i << " : " << ratio * 100 << " %" << endl;

        // adjust to bone length
        jointTrans.row(BE(i, 1)) = (1 - ratio) * (C.row(BE(i, 0)) - C.row(BE(i, 1)));
        if (P(i) < 0)
            continue;
        jointTrans.row(BE(i, 1)) += jointTrans.row(BE(P(i), 1));
    }

    jointTrans.row(eyeLJ) = C.row(headJ) + jointTrans.row(headJ) + eyeL_pos.transpose() - C.row(eyeLJ);
    jointTrans.row(eyeRJ) = C.row(headJ) + jointTrans.row(headJ) + eyeR_pos.transpose() - C.row(eyeRJ);
    jointTrans.row(headJ) = MatrixXd::Zero(1, 3);
    jointTrans(headJ, 1) = (jointTrans(eyeLJ, 1) + jointTrans(eyeRJ, 1)) * 0.5;
    C_calib = C + jointTrans;

    cout << Wj.rows() << "*" << Wj.cols() << endl;
    cout << jointTrans.rows() << "*" << jointTrans.cols() << endl;
    V_calib = V + Wj * jointTrans.block(0, 0, C.rows() - 1, 3);
    //MatrixXd jt = jointTrans.block(0,0,C.rows()-1,3);

    return ss.str();
}

void PhantomAnimator::Animate(RotationList vQ, const MatrixXd &C_disp, MatrixXd &C_new, bool calibChk)
{
    vector<Vector3d> vT;

    if (calibChk)
    {
        C_new = C_calib;
        C_new.row(0) = C_disp.row(0); //set root

        for (int i = 0; i < BE.rows(); i++)
        {
            Affine3d a;
            a = Translation3d(Vector3d(C_new.row(BE(i, 0)).transpose())) * vQ[i].matrix() * Translation3d(Vector3d(-C_calib.row(BE(i, 0)).transpose()));
            vT.push_back(a.translation());
            C_new.row(BE(i, 1)) = a * Vector3d(C_new.row(BE(i, 1)));
        }
        myDqs(V_calib, cleanWeights, vQ, vT, U);
    }
    else
    {
        C_new = C;
        C_new.row(0) = C_disp.row(0); //set root
        for (int i = 0; i < BE.rows(); i++)
        {
            Affine3d a;
            a = Translation3d(Vector3d(C_new.row(BE(i, 0)).transpose())) * vQ[i].matrix() * Translation3d(Vector3d(-C.row(BE(i, 0)).transpose()));
            vT.push_back(a.translation());
            C_new.row(BE(i, 1)) = a * Vector3d(C_new.row(BE(i, 1)));
        }
        myDqs(V, cleanWeights, vQ, vT, U);
    }
}

void PhantomAnimator::Animate(RotationList vQ, MatrixXd &V_new)
{
    vector<Vector3d> vT;
    igl::forward_kinematics(C_calib, BE, P, vQ, vQ, vT);

    MatrixXd C_new = C_calib;
    for (int i = 0; i < BE.rows(); i++)
    {
        Affine3d a;
        a = Translation3d(Vector3d(C_new.row(BE(i, 0)).transpose())) * vQ[i].matrix() * Translation3d(Vector3d(-C_new.row(BE(i, 0)).transpose()));
        vT.push_back(a.translation());
        C_new.row(BE(i, 1)) = a * Vector3d(C_new.row(BE(i, 1)));
    }
    myDqs(V_calib, cleanWeights, vQ, vT, V_new);
}

bool PhantomAnimator::ReadProfileData(string fileName)
{
    ifstream ifs(fileName);
    if (!ifs.is_open())
    {
        cout << "There is no " + fileName << endl;
        return false;
    }
    profileIDs.clear();
    jointLengths.clear();
    eyeR_vec.clear();
    eyeL_vec.clear();

    int num;
    ifs >> num;
    string firstLine;
    getline(ifs, firstLine);
    jointLengths.resize(num);
    eyeR_vec.resize(num);
    eyeL_vec.resize(num);

    getline(ifs, firstLine);
    stringstream ss(firstLine);
    vector<int> boneIDs;
    for (int i = 0; i < 17; i++)
    {
        int boneID;
        ss >> boneID;
        boneIDs.push_back(boneID);
    }

    for (int i = 0; i < num; i++)
    {
        string aLine;
        getline(ifs, aLine);
        stringstream ss1(aLine);
        string name;
        double d;
        ss1 >> name;
        profileIDs[name] = i;
        for (int j = 0; j < 17; j++)
        {
            ss1 >> d;
            jointLengths[i][boneIDs[j]] = d;
        }
        for (int j = 0; j < 3; j++)
        {
            ss1 >> d;
            eyeL_vec[i](j) = d;
        }
        for (int j = 0; j < 3; j++)
        {
            ss1 >> d;
            eyeR_vec[i](j) = d;
        }
    }

    ifs.close();
    return true;
}

bool PhantomAnimator::WriteProfileData(string filename)
{
    if (!jointLengths.size())
        return false;
    ofstream ofs(filename);
    int num = jointLengths.size();
    ofs << num << endl
        << "\t";

    vector<int> boneIDs;
    for (auto iter : jointLengths[0])
        boneIDs.push_back(iter.first);

    for (int i : boneIDs)
    {
        ofs << i << "\t";
    }
    ofs << "eyeL_x\teyeL_y\teyeL_z\teyeR_x\teyeR_y\teyeR_z" << endl;

    vector<string> names;
    names.resize(num);
    for (auto iter : profileIDs)
        names[iter.second] = iter.first;

    for (int i = 0; i < num; i++)
    {
        ofs << names[i] << "\t";
        for (int j = 0; j < 17; j++)
            ofs << jointLengths[i][boneIDs[j]] << "\t";
        for (int j = 0; j < 3; j++)
            ofs << eyeL_vec[i](j) << "\t";
        for (int j = 0; j < 3; j++)
            ofs << eyeR_vec[i](j) << "\t";
        ofs << endl;
    }

    ofs.close();

    return true;
}
