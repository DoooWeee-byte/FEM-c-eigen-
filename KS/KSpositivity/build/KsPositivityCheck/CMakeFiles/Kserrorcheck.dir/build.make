# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build

# Include any dependencies generated for this target.
include KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/depend.make

# Include the progress variables for this target.
include KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/progress.make

# Include the compile flags for this target's objects.
include KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/flags.make

KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.o: KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/flags.make
KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.o: ../KsPositivityCheck/Kserrorcheck.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.o"
	cd /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.o -c /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/KsPositivityCheck/Kserrorcheck.cpp

KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.i"
	cd /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/KsPositivityCheck/Kserrorcheck.cpp > CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.i

KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.s"
	cd /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/KsPositivityCheck/Kserrorcheck.cpp -o CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.s

# Object files for target Kserrorcheck
Kserrorcheck_OBJECTS = \
"CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.o"

# External object files for target Kserrorcheck
Kserrorcheck_EXTERNAL_OBJECTS =

KsPositivityCheck/Kserrorcheck: KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/Kserrorcheck.cpp.o
KsPositivityCheck/Kserrorcheck: KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/build.make
KsPositivityCheck/Kserrorcheck: /usr/lib/x86_64-linux-gnu/libpython3.8.so
KsPositivityCheck/Kserrorcheck: KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Kserrorcheck"
	cd /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Kserrorcheck.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/build: KsPositivityCheck/Kserrorcheck

.PHONY : KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/build

KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/clean:
	cd /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck && $(CMAKE_COMMAND) -P CMakeFiles/Kserrorcheck.dir/cmake_clean.cmake
.PHONY : KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/clean

KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/depend:
	cd /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/KsPositivityCheck /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck /home/pangdongwen/文档/程序文件/有限元/KS/KSpositivity/build/KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : KsPositivityCheck/CMakeFiles/Kserrorcheck.dir/depend

