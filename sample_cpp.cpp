#include "sample_cpp.h"
#include <../../pybind11/include/pybind11/pybind11.h>
namespace py = pybind11;
using namespace std;

int main(){
	return 0;
}


float add(float x){
	return x + 1.0;
}

PYBIND11_MODULE(sample_cpp, m) {
    m.doc() = "Example for FSW"; // optional module docstring

    m.def("add", &add, "Adds 1 to the number");
    
}