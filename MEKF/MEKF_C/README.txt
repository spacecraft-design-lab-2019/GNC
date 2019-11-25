The steps of the MEKF are encapsulated in the following functions:
  - predict.cpp
  - measurement.cpp
  - Innovation.cpp
  - update.cpp

and executed in MEKF.cpp with header file MEKF.hpp.

The current inputs to the MEKF (as of 11/17/2019) come from the MEKF starter code residing in a seperate folder within this github,
with rB1hist, rB2hist, and whist being imported as .txt files to the filter. 

The following Linux command line sequence is used to execute the filter:

g++ -std=c++11 -Wall  -Wextra -Wpedantic MEKF.cpp predict.cpp update.cpp measurement.cpp Innovation.cpp triad_ad.cpp DCM2q.cpp -o MEKF 
./MEKF

triad_ad.cpp and DCM2q.cpp are utility functions used in the MEKF to make an initial guess on the quaternion, which is then propagated in 
the MEKF.
