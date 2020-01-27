/*
 * Author: Nick Goodson
 * Jan 15th 2020
 *
*/

#include "iLQRsimple.h"
#include <sstream>
#include <fstream>

using namespace std;
using namespace Eigen;


/**
 * Converts primitive types to strings
 */
template<typename T>
string toString(T const& value) {
    ostringstream sstr;
    sstr << value;
    return sstr.str();
}


/**
 * Writes iLQR results to a csv file
 * Columns are (xtraj[0], xtraj[1], utraj, J)
 */
void writeToFile(const MatrixXd& xtraj, const MatrixXd& utraj, const vector<double>& Jhist) {
    auto N = static_cast<unsigned int>( xtraj.cols() );
    ofstream datafile;
    datafile.open("iLQR_pendulum_data.csv");
    for (int i = 0; i < N-1; ++i ) {
        string data_line = toString(xtraj(0, i)) + "," + toString(xtraj(1, i)) + ","
                           + toString(utraj(0, i)) + "," + toString(Jhist[i]) + "\n";
        datafile << data_line;
    }
    datafile.close();
    cout << "File written successfully" << endl;
}


/**
 * Prints an Eigen::Matrix to the console
 * @param mat
 */
void printMatrix(const MatrixXd& mat) {
    auto Nrow = static_cast<unsigned int>( mat.rows() );
    auto Ncol = static_cast<unsigned int>( mat.cols() );

    cout << "[";
    string rowString = "[";
    for( unsigned int i = 0; i < Nrow; ++i ) {
        for ( unsigned int j = 0; j < Ncol; ++j ) {
            if ( j == 0 ) {
                rowString += toString(mat(i, j));
            } else {
                rowString += ", " + toString(mat(i, j));
            }
        }
        rowString += "]";

        if ( i == Nrow -1 ) {
            // Add final bracket for last row
            cout << rowString << "]" << endl;
        } else {
            cout << rowString << endl;
        }
        rowString = "[";
    }
}

