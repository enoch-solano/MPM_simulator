# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build

# Utility rule file for SuiteSparse.

# Include the progress variables for this target.
include CMakeFiles/SuiteSparse.dir/progress.make

CMakeFiles/SuiteSparse:
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/SuiteSparse && make library -j 8

SuiteSparse: CMakeFiles/SuiteSparse
SuiteSparse: CMakeFiles/SuiteSparse.dir/build.make

.PHONY : SuiteSparse

# Rule to build all files generated by this target.
CMakeFiles/SuiteSparse.dir/build: SuiteSparse

.PHONY : CMakeFiles/SuiteSparse.dir/build

CMakeFiles/SuiteSparse.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SuiteSparse.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SuiteSparse.dir/clean

CMakeFiles/SuiteSparse.dir/depend:
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/CMakeFiles/SuiteSparse.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SuiteSparse.dir/depend
