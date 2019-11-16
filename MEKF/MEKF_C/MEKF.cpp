//#include "MEKF.h"
#include <iostream>
#include <fstream>
#include <string>
#include "../../eigen-git-mirror/Eigen/Dense"

//using Eigen::MatrixXd;
using namespace std;
using namespace Eigen;

MatrixXd triad_ad(MatrixXd M, MatrixXd V);

#define MAXBUFSIZE  ((int) 1e6)


int main() {

    // initialize W and V based on MEKF template data
    MatrixXd W(6,6);
    MatrixXd V(6,6);
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++){
            if (i==j){
                if (i < 3){
                    W(i,j) = 0.003;
                    V(i,j) = 0.003;
                }
                else{
                    W(i,j) = 0.0076;
                    V(i,j) = 0.0076;
                }
            }
            else {
                    W(i,j) = 0.0;
                    V(i,j) = 0.0;
                }
        }
    }

    cout << W << endl;
    cout << V << endl;
    /*MEKF_inputs = loadmat('mekf_inputs.mat')
    MEKF_truth  = loadmat('mekf_truth.mat')
    b_true = MEKF_truth['btrue']
    q_true = MEKF_truth['qtrue']
    rN1 = MEKF_inputs['rN1']
    rN2 = MEKF_inputs['rN2']
    rN = np.vstack((rN1,rN2))
    whist = MEKF_inputs['whist']
    rB2hist = MEKF_inputs['rB2hist']
    rB1hist = MEKF_inputs['rB1hist']
    W = MEKF_inputs['W']
    V = MEKF_inputs['V']
    dt = MEKF_inputs['dt']
    */


   return 0;
}

MatrixXd readMatrix(const char *filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
}

MatrixXd triad_ad(MatrixXd M, MatrixXd V) {
    /*
    Gives rotation matrix from inertial to body frame
    Inputs :
    M - Matrix where each column is a measurement vector in the body frame - Note, most accurate measurement should be in the first column
    V - Matrix where each column is a modeled vector of the corresponding measurement vector from M, in the inertial frame
    Outputs:
    R - rotation matrix from inertial to body frame
    */
    MatrixXd R(3,3);
    if (M.outerSize() == 2 && V.outerSize() == 2)
    {
        Vector3d m1 = M.col(0);
        Vector3d mtemp = M.col(1);
        Vector3d m2 = m1.cross(mtemp);
        Vector3d m3 = m1.cross(m2);

        Vector3d v1 = V.col(0);
        Vector3d vtemp = V.col(1);
        Vector3d v2 = v1.cross(vtemp);
        Vector3d v3 = v1.cross(v2);

        MatrixXd Rtemp(3,3);
        MatrixXd Vtemp(3, 3);
        Rtemp << m1, m2, m3;
        Vtemp << v1, v2, v3;
        R = Rtemp * Vtemp.inverse();
    }
    else
    {
        R = M * V.completeOrthogonalDecomposition().pseudoInverse();
    }
    return R;
}


