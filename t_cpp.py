from GNC.cmake_build_debug import sun_utils_cpp
from GNC.cmake_build_debug import time_functions_cpp
x = sun_utils_cpp.sun_position(float(540000))
y = time_functions_cpp.date2MJD(1, 1, 2000, 12, 0, 0)
print(type(x))
print(y)