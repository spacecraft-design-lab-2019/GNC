/* foo.cpp */

#include "foo.h"

#include <iostream>

extern "C" {

void foo(int arg) {
  std::cout << arg << std::endl;
}

} /* extern "C" */