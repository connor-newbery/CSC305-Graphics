# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build

# Include any dependencies generated for this target.
include CMakeFiles/assignment1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/assignment1.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/assignment1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assignment1.dir/flags.make

CMakeFiles/assignment1.dir/src/main.cpp.o: CMakeFiles/assignment1.dir/flags.make
CMakeFiles/assignment1.dir/src/main.cpp.o: /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/src/main.cpp
CMakeFiles/assignment1.dir/src/main.cpp.o: CMakeFiles/assignment1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/assignment1.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/assignment1.dir/src/main.cpp.o -MF CMakeFiles/assignment1.dir/src/main.cpp.o.d -o CMakeFiles/assignment1.dir/src/main.cpp.o -c /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/src/main.cpp

CMakeFiles/assignment1.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assignment1.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/src/main.cpp > CMakeFiles/assignment1.dir/src/main.cpp.i

CMakeFiles/assignment1.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assignment1.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/src/main.cpp -o CMakeFiles/assignment1.dir/src/main.cpp.s

# Object files for target assignment1
assignment1_OBJECTS = \
"CMakeFiles/assignment1.dir/src/main.cpp.o"

# External object files for target assignment1
assignment1_EXTERNAL_OBJECTS =

assignment1: CMakeFiles/assignment1.dir/src/main.cpp.o
assignment1: CMakeFiles/assignment1.dir/build.make
assignment1: CMakeFiles/assignment1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable assignment1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assignment1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assignment1.dir/build: assignment1
.PHONY : CMakeFiles/assignment1.dir/build

CMakeFiles/assignment1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assignment1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assignment1.dir/clean

CMakeFiles/assignment1.dir/depend:
	cd /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1 /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1 /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build /Users/connornewbery/uvic/summer-2023/csc305/Assignment_1/build/CMakeFiles/assignment1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/assignment1.dir/depend

