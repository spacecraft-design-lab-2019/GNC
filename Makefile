# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/cmake-3.16.2-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /opt/cmake-3.16.2-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nick/Documents/AA236/GNC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nick/Documents/AA236/GNC

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."

	/opt/cmake-3.16.2-Linux-x86_64/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/opt/cmake-3.16.2-Linux-x86_64/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/nick/Documents/AA236/GNC/CMakeFiles /home/nick/Documents/AA236/GNC/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/nick/Documents/AA236/GNC/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named sun_utils_cpp

# Build rule for target.
sun_utils_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sun_utils_cpp
.PHONY : sun_utils_cpp

# fast build rule for target.
sun_utils_cpp/fast:
	$(MAKE) -f CMakeFiles/sun_utils_cpp.dir/build.make CMakeFiles/sun_utils_cpp.dir/build
.PHONY : sun_utils_cpp/fast

#=============================================================================
# Target rules for targets named detumble_cpp

# Build rule for target.
detumble_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 detumble_cpp
.PHONY : detumble_cpp

# fast build rule for target.
detumble_cpp/fast:
	$(MAKE) -f CMakeFiles/detumble_cpp.dir/build.make CMakeFiles/detumble_cpp.dir/build
.PHONY : detumble_cpp/fast

#=============================================================================
# Target rules for targets named time_functions_cpp

# Build rule for target.
time_functions_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 time_functions_cpp
.PHONY : time_functions_cpp

# fast build rule for target.
time_functions_cpp/fast:
	$(MAKE) -f CMakeFiles/time_functions_cpp.dir/build.make CMakeFiles/time_functions_cpp.dir/build
.PHONY : time_functions_cpp/fast

#=============================================================================
# Target rules for targets named triad_cpp

# Build rule for target.
triad_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 triad_cpp
.PHONY : triad_cpp

# fast build rule for target.
triad_cpp/fast:
	$(MAKE) -f CMakeFiles/triad_cpp.dir/build.make CMakeFiles/triad_cpp.dir/build
.PHONY : triad_cpp/fast

#=============================================================================
# Target rules for targets named frame_conversions_cpp

# Build rule for target.
frame_conversions_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 frame_conversions_cpp
.PHONY : frame_conversions_cpp

# fast build rule for target.
frame_conversions_cpp/fast:
	$(MAKE) -f CMakeFiles/frame_conversions_cpp.dir/build.make CMakeFiles/frame_conversions_cpp.dir/build
.PHONY : frame_conversions_cpp/fast

#=============================================================================
# Target rules for targets named sample_cpp

# Build rule for target.
sample_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sample_cpp
.PHONY : sample_cpp

# fast build rule for target.
sample_cpp/fast:
	$(MAKE) -f CMakeFiles/sample_cpp.dir/build.make CMakeFiles/sample_cpp.dir/build
.PHONY : sample_cpp/fast

#=============================================================================
# Target rules for targets named euler_cpp

# Build rule for target.
euler_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 euler_cpp
.PHONY : euler_cpp

# fast build rule for target.
euler_cpp/fast:
	$(MAKE) -f CMakeFiles/euler_cpp.dir/build.make CMakeFiles/euler_cpp.dir/build
.PHONY : euler_cpp/fast

#=============================================================================
# Target rules for targets named magnetic_field_cpp

# Build rule for target.
MEKF_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 MEKF_cpp
.PHONY : MEKF_cpp

# fast build rule for target.
MEKF_cpp/fast:
	$(MAKE) -f CMakeFiles/MEKF_cpp.dir/build.make CMakeFiles/MEKF_cpp.dir/build
.PHONY : MEKF_cpp/fast

#=============================================================================
# Target rules for targets named MEKF_cpp

# Build rule for target.
MEKF_cpp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 MEKF_cpp
.PHONY : MEKF_cpp

# fast build rule for target.
MEKF_cpp/fast:
	$(MAKE) -f CMakeFiles/MEKF_cpp.dir/build.make CMakeFiles/MEKF_cpp.dir/build
.PHONY : MEKF_cpp/fast

MEKF/MEKF_cpp/MEKF_cpp.o: MEKF/MEKF_cpp/MEKF_cpp.cpp.o

.PHONY : MEKF/MEKF_cpp/MEKF_cpp.o

# target to build an object file
MEKF/MEKF_cpp/MEKF_cpp.cpp.o:
	$(MAKE) -f CMakeFiles/MEKF_cpp.dir/build.make CMakeFiles/MEKF_cpp.dir/MEKF/MEKF_cpp/MEKF_cpp.cpp.o
.PHONY : MEKF/MEKF_cpp/MEKF_cpp.cpp.o

MEKF/MEKF_cpp/MEKF_cpp.i: MEKF/MEKF_cpp/MEKF_cpp.cpp.i

.PHONY : MEKF/MEKF_cpp/MEKF_cpp.i

# target to preprocess a source file
MEKF/MEKF_cpp/MEKF_cpp.cpp.i:
	$(MAKE) -f CMakeFiles/MEKF_cpp.dir/build.make CMakeFiles/MEKF_cpp.dir/MEKF/MEKF_cpp/MEKF_cpp.cpp.i
.PHONY : MEKF/MEKF_cpp/MEKF_cpp.cpp.i

MEKF/MEKF_cpp/MEKF_cpp.s: MEKF/MEKF_cpp/MEKF_cpp.cpp.s

.PHONY : MEKF/MEKF_cpp/MEKF_cpp.s

# target to generate assembly for a file
MEKF/MEKF_cpp/MEKF_cpp.cpp.s:
	$(MAKE) -f CMakeFiles/MEKF_cpp.dir/build.make CMakeFiles/MEKF_cpp.dir/MEKF/MEKF_cpp/MEKF_cpp.cpp.s
.PHONY : MEKF/MEKF_cpp/MEKF_cpp.cpp.s

TRIAD/cpp/deterministic_ad.o: TRIAD/cpp/deterministic_ad.cpp.o

.PHONY : TRIAD/cpp/deterministic_ad.o

# target to build an object file
TRIAD/cpp/deterministic_ad.cpp.o:
	$(MAKE) -f CMakeFiles/triad_cpp.dir/build.make CMakeFiles/triad_cpp.dir/TRIAD/cpp/deterministic_ad.cpp.o
.PHONY : TRIAD/cpp/deterministic_ad.cpp.o

TRIAD/cpp/deterministic_ad.i: TRIAD/cpp/deterministic_ad.cpp.i

.PHONY : TRIAD/cpp/deterministic_ad.i

# target to preprocess a source file
TRIAD/cpp/deterministic_ad.cpp.i:
	$(MAKE) -f CMakeFiles/triad_cpp.dir/build.make CMakeFiles/triad_cpp.dir/TRIAD/cpp/deterministic_ad.cpp.i
.PHONY : TRIAD/cpp/deterministic_ad.cpp.i

TRIAD/cpp/deterministic_ad.s: TRIAD/cpp/deterministic_ad.cpp.s

.PHONY : TRIAD/cpp/deterministic_ad.s

# target to generate assembly for a file
TRIAD/cpp/deterministic_ad.cpp.s:
	$(MAKE) -f CMakeFiles/triad_cpp.dir/build.make CMakeFiles/triad_cpp.dir/TRIAD/cpp/deterministic_ad.cpp.s
.PHONY : TRIAD/cpp/deterministic_ad.cpp.s

util_funcs/cpp/sun_utils.o: util_funcs/cpp/sun_utils.cpp.o

.PHONY : util_funcs/cpp/sun_utils.o

# target to build an object file
util_funcs/cpp/sun_utils.cpp.o:
	$(MAKE) -f CMakeFiles/sun_utils_cpp.dir/build.make CMakeFiles/sun_utils_cpp.dir/util_funcs/cpp/sun_utils.cpp.o
.PHONY : util_funcs/cpp/sun_utils.cpp.o

util_funcs/cpp/sun_utils.i: util_funcs/cpp/sun_utils.cpp.i

.PHONY : util_funcs/cpp/sun_utils.i

# target to preprocess a source file
util_funcs/cpp/sun_utils.cpp.i:
	$(MAKE) -f CMakeFiles/sun_utils_cpp.dir/build.make CMakeFiles/sun_utils_cpp.dir/util_funcs/cpp/sun_utils.cpp.i
.PHONY : util_funcs/cpp/sun_utils.cpp.i

util_funcs/cpp/sun_utils.s: util_funcs/cpp/sun_utils.cpp.s

.PHONY : util_funcs/cpp/sun_utils.s

# target to generate assembly for a file
util_funcs/cpp/sun_utils.cpp.s:
	$(MAKE) -f CMakeFiles/sun_utils_cpp.dir/build.make CMakeFiles/sun_utils_cpp.dir/util_funcs/cpp/sun_utils.cpp.s
.PHONY : util_funcs/cpp/sun_utils.cpp.s

util_funcs/cpp/time_functions.o: util_funcs/cpp/time_functions.cpp.o

.PHONY : util_funcs/cpp/time_functions.o

# target to build an object file
util_funcs/cpp/time_functions.cpp.o:
	$(MAKE) -f CMakeFiles/time_functions_cpp.dir/build.make CMakeFiles/time_functions_cpp.dir/util_funcs/cpp/time_functions.cpp.o
.PHONY : util_funcs/cpp/time_functions.cpp.o

util_funcs/cpp/time_functions.i: util_funcs/cpp/time_functions.cpp.i

.PHONY : util_funcs/cpp/time_functions.i

# target to preprocess a source file
util_funcs/cpp/time_functions.cpp.i:
	$(MAKE) -f CMakeFiles/time_functions_cpp.dir/build.make CMakeFiles/time_functions_cpp.dir/util_funcs/cpp/time_functions.cpp.i
.PHONY : util_funcs/cpp/time_functions.cpp.i

util_funcs/cpp/time_functions.s: util_funcs/cpp/time_functions.cpp.s

.PHONY : util_funcs/cpp/time_functions.s

# target to generate assembly for a file
util_funcs/cpp/time_functions.cpp.s:
	$(MAKE) -f CMakeFiles/time_functions_cpp.dir/build.make CMakeFiles/time_functions_cpp.dir/util_funcs/cpp/time_functions.cpp.s
.PHONY : util_funcs/cpp/time_functions.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"

	@echo "... rebuild_cache"
	@echo "... sun_utils_cpp"
	@echo "... detumble_cpp"
	@echo "... time_functions_cpp"
	@echo "... triad_cpp"
	@echo "... edit_cache"
	@echo "... frame_conversions_cpp"
	@echo "... sample_cpp"
	@echo "... euler_cpp"
	@echo "... magnetic_field_cpp"
	@echo "... MEKF_cpp"
	@echo "... MEKF/MEKF_cpp/MEKF_cpp.o"
	@echo "... MEKF/MEKF_cpp/MEKF_cpp.i"
	@echo "... MEKF/MEKF_cpp/MEKF_cpp.s"
	@echo "... TRIAD/cpp/deterministic_ad.o"
	@echo "... TRIAD/cpp/deterministic_ad.i"
	@echo "... TRIAD/cpp/deterministic_ad.s"
	@echo "... util_funcs/cpp/sun_utils.o"
	@echo "... util_funcs/cpp/sun_utils.i"
	@echo "... util_funcs/cpp/sun_utils.s"
	@echo "... util_funcs/cpp/time_functions.o"
	@echo "... util_funcs/cpp/time_functions.i"
	@echo "... util_funcs/cpp/time_functions.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

