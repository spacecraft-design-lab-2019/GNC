

#include "pointer_t.h"
#include <iostream>
#include "../../eigen-git-mirror/Eigen/Dense"
#include <../../pybind11/include/pybind11/pybind11.h>
#include <../../pybind11/include/pybind11/eigen.h>
using namespace Eigen;

using namespace std;

#include <memory>

namespace py = pybind11;

template <class T> class ptr_wrapper
{
    public:
        ptr_wrapper() : ptr(nullptr) {}
        ptr_wrapper(T* ptr) : ptr(ptr) {}
        ptr_wrapper(const ptr_wrapper& other) : ptr(other.ptr) {}
        T& operator* () const { return *ptr; }
        T* operator->() const { return  ptr; }
        T* get() const { return ptr; }
        void destroy() { delete ptr; }
        T& operator[](std::size_t idx) const { return ptr[idx]; }
    private:
        T* ptr;
};

ptr_wrapper<float> get_ptr(void);

int functest(float*);
float* get_ptr(float x);

int main(){
	float x;
	float * y = &x;
	x = 5.0;

	float array[8];
	for(int i = 0; i < 8; ++i) array[i] = i;
	float * pa = array;
	
	cout << *pa << endl;
	cout << array << endl;
	functest(pa);
	for (int n=0; n<8; n++)
    	cout << array[n] << ", ";

    return 0;
}

int functest(float* z){

	Map<MatrixXf> mf(z, 2, 4);
	MatrixXf R(4, 4);
	R << 0, 1, 5, 10,
		0, 1, 5, 10,
		0, 1, 5, 10,
		0, 1, 5, 10;

	mf = mf * R;
	// *z[0] = 100;
	return 0;
}



float* get_ptr(float x){
	return &x;
}

PYBIND11_MODULE(pointer_t_cpp, m) {
    m.doc() = "Pointer testing"; // optional module docstring

    m.def("functest", &functest, "");
    m.def("get_ptr", &get_ptr, "");
}