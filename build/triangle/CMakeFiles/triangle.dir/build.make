# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build

# Include any dependencies generated for this target.
include triangle/CMakeFiles/triangle.dir/depend.make

# Include the progress variables for this target.
include triangle/CMakeFiles/triangle.dir/progress.make

# Include the compile flags for this target's objects.
include triangle/CMakeFiles/triangle.dir/flags.make

triangle/CMakeFiles/triangle.dir/triangle.c.o: triangle/CMakeFiles/triangle.dir/flags.make
triangle/CMakeFiles/triangle.dir/triangle.c.o: ../Deps/libigl/external/triangle/triangle.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object triangle/CMakeFiles/triangle.dir/triangle.c.o"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/triangle.dir/triangle.c.o   -c /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/triangle/triangle.c

triangle/CMakeFiles/triangle.dir/triangle.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/triangle.dir/triangle.c.i"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/triangle/triangle.c > CMakeFiles/triangle.dir/triangle.c.i

triangle/CMakeFiles/triangle.dir/triangle.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/triangle.dir/triangle.c.s"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/triangle/triangle.c -o CMakeFiles/triangle.dir/triangle.c.s

# Object files for target triangle
triangle_OBJECTS = \
"CMakeFiles/triangle.dir/triangle.c.o"

# External object files for target triangle
triangle_EXTERNAL_OBJECTS =

triangle/libtriangle.a: triangle/CMakeFiles/triangle.dir/triangle.c.o
triangle/libtriangle.a: triangle/CMakeFiles/triangle.dir/build.make
triangle/libtriangle.a: triangle/CMakeFiles/triangle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libtriangle.a"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle && $(CMAKE_COMMAND) -P CMakeFiles/triangle.dir/cmake_clean_target.cmake
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/triangle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
triangle/CMakeFiles/triangle.dir/build: triangle/libtriangle.a

.PHONY : triangle/CMakeFiles/triangle.dir/build

triangle/CMakeFiles/triangle.dir/clean:
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle && $(CMAKE_COMMAND) -P CMakeFiles/triangle.dir/cmake_clean.cmake
.PHONY : triangle/CMakeFiles/triangle.dir/clean

triangle/CMakeFiles/triangle.dir/depend:
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/triangle /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/triangle/CMakeFiles/triangle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : triangle/CMakeFiles/triangle.dir/depend

