from GNC.cmake_build_debug import sun_utils_cpp
from GNC.cmake_build_debug import time_functions_cpp
x = sun_utils_cpp.sun_position(float(540000))
y = time_functions_cpp.date2MJD(float(1), float(1), float(2000), float(12), float(0), float(0))
print(type(x))
print(y)