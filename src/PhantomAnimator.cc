#include "PhantomAnimator.hh"
#include "G4SystemOfUnits.hh"
PhantomAnimator::PhantomAnimator() {}
PhantomAnimator::~PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string prefix)
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE);
    igl::directed_edge_parents(BE, P);

    //distance to parent joint
    for (int i = 0; i < BE.rows(); i++)
        lengths[i] = (C.row(BE(i, 0)) - C.row(BE(i, 1))).norm();

    ReadTetMesh(prefix);
    V_calib = V;
    C_calib = C;

    if (!igl::readDMAT(prefix + ".W", W) || !igl::readDMAT(prefix + ".Wj", Wj))
        PreparePhantom(prefix);

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
    cout << "perform BBW with "+prefix+".ply" << endl;
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

    cout << "calulate barycentric coodinates" << endl;
    auto baryMap = GenerateBarycentricCoord(VT, TT, V);
    SparseMatrix<double> bary = GenerateBarySparse(baryMap, VT.rows());

    W = bary * W;
    Wj = bary * Wj;
    cout << "write " + prefix + ".W" << endl;
    igl::writeDMAT(prefix + ".W", W, false);
    cout << "write " + prefix + ".Wj" << endl;
    igl::writeDMAT(prefix + ".Wj", Wj, false);
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
    U = V_calib;
    //MatrixXd jt = jointTrans.block(0,0,C.rows()-1,3);

    return ss.str();
}

void PhantomAnimator::Animate(RotationList vQ, Vector3d root)
{
    vector<Vector3d> vT;

    MatrixXd C_new = C_calib;
    C_new.row(0) = root; //set root

    for (int i = 0; i < BE.rows(); i++)
    {
        Affine3d a;
        a = Translation3d(Vector3d(C_new.row(BE(i, 0)).transpose())) * vQ[i].matrix() * Translation3d(Vector3d(-C_calib.row(BE(i, 0)).transpose()));
        vT.push_back(a.translation());
        C_new.row(BE(i, 1)) = a * Vector3d(C_new.row(BE(i, 1)));
    }
    myDqs(V_calib, cleanWeights, vQ, vT, U);
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

