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
include stb_image/CMakeFiles/igl_stb_image.dir/depend.make

# Include the progress variables for this target.
include stb_image/CMakeFiles/igl_stb_image.dir/progress.make

# Include the compile flags for this target's objects.
include stb_image/CMakeFiles/igl_stb_image.dir/flags.make

stb_image/CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.o: stb_image/CMakeFiles/igl_stb_image.dir/flags.make
stb_image/CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.o: ../Deps/libigl/external/stb_image/igl_stb_image.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object stb_image/CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.o"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.o -c /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/stb_image/igl_stb_image.cpp

stb_image/CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.i"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/stb_image/igl_stb_image.cpp > CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.i

stb_image/CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.s"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/stb_image/igl_stb_image.cpp -o CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.s

# Object files for target igl_stb_image
igl_stb_image_OBJECTS = \
"CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.o"

# External object files for target igl_stb_image
igl_stb_image_EXTERNAL_OBJECTS =

stb_image/libigl_stb_image.a: stb_image/CMakeFiles/igl_stb_image.dir/igl_stb_image.cpp.o
stb_image/libigl_stb_image.a: stb_image/CMakeFiles/igl_stb_image.dir/build.make
stb_image/libigl_stb_image.a: stb_image/CMakeFiles/igl_stb_image.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libigl_stb_image.a"
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image && $(CMAKE_COMMAND) -P CMakeFiles/igl_stb_image.dir/cmake_clean_target.cmake
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/igl_stb_image.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
stb_image/CMakeFiles/igl_stb_image.dir/build: stb_image/libigl_stb_image.a

.PHONY : stb_image/CMakeFiles/igl_stb_image.dir/build

stb_image/CMakeFiles/igl_stb_image.dir/clean:
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image && $(CMAKE_COMMAND) -P CMakeFiles/igl_stb_image.dir/cmake_clean.cmake
.PHONY : stb_image/CMakeFiles/igl_stb_image.dir/clean

stb_image/CMakeFiles/igl_stb_image.dir/depend:
	cd /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/Deps/libigl/external/stb_image /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image /home/enoch/Desktop/cis563/projects/project03/cis563-2019-assignment/build/stb_image/CMakeFiles/igl_stb_image.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : stb_image/CMakeFiles/igl_stb_image.dir/depend

