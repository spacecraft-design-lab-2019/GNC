/*
 * Author: Nick Goodson
 * Jan 28th 2020
 *
*/

#include "iLQRsimple.h"
#include <sstream>

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
 * Columns are (xtraj-states, utraj-controls, J-cost)
 */
void writeToFile(const MatrixXd& xtraj, const MatrixXd& utraj, const vector<double>& Jhist) {
    auto N = static_cast<unsigned int>( xtraj.cols() );
    auto Nx = static_cast<unsigned int>( xtraj.rows() );
    auto Nu = static_cast<unsigned int>( utraj.rows() );

    ofstream datafile;
    datafile.open("iLQR_pendulum_data.csv");
    for (int i = 0; i < N-1; ++i) {
        string data_line;
        for (int state = 0; state < Nx; ++state) {
            data_line += toString(xtraj(state, i)) + ",";  // Add states
        }
        for (int contr = 0; contr < Nu; ++contr) {
            data_line += toString(utraj(contr, i)) + ",";  // Add controls
        }
        if (i < Jhist.size()) {
            data_line += toString(Jhist[i]) + "\n";  // Add cost
        } else {
            data_line += "0\n";  // Pad with zeros b/c cost is different length
        }
        datafile << data_line;
    }
    datafile.close();
    cout << "File written successfully" << endl;
}


/**
 * Prints an Eigen::Matrix to the console
 * @param mat, the matrix to print
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

