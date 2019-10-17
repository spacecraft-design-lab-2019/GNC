#include <iostream>
#include "eigen-git-mirror/Eigen/Dense"

using Eigen::MatrixXd;

int main() {
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0)=2.5;
    m(0,1)=-1;
    m(1,1)=m(1,0)+m(0,1);

    MatrixXd n(2,2);
    n(0,0) = 10;
    n(1,0)= 25;
    n(0,1)=-10;
    n(1,1)=m(1,0)+m(0,1);

    std::cout << "m: \n" << m << std::endl;
    std::cout << "n: \n" << n << std::endl;
    std::cout << "m*n: \n" << m*n << std::endl;
    return 0;
}