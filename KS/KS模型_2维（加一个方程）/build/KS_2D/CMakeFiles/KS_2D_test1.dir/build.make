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
CMAKE_SOURCE_DIR = /home/pangdongwen/文档/程序文件/有限元/KS模型_2维

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build

# Include any dependencies generated for this target.
include KS_2D/CMakeFiles/KS_2D_test1.dir/depend.make

# Include the progress variables for this target.
include KS_2D/CMakeFiles/KS_2D_test1.dir/progress.make

# Include the compile flags for this target's objects.
include KS_2D/CMakeFiles/KS_2D_test1.dir/flags.make

KS_2D/CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.o: KS_2D/CMakeFiles/KS_2D_test1.dir/flags.make
KS_2D/CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.o: ../KS_2D/KS_2D_test1.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object KS_2D/CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.o"
	cd /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.o -c /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/KS_2D/KS_2D_test1.cpp

KS_2D/CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.i"
	cd /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/KS_2D/KS_2D_test1.cpp > CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.i

KS_2D/CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.s"
	cd /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/KS_2D/KS_2D_test1.cpp -o CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.s

# Object files for target KS_2D_test1
KS_2D_test1_OBJECTS = \
"CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.o"

# External object files for target KS_2D_test1
KS_2D_test1_EXTERNAL_OBJECTS =

KS_2D/KS_2D_test1: KS_2D/CMakeFiles/KS_2D_test1.dir/KS_2D_test1.cpp.o
KS_2D/KS_2D_test1: KS_2D/CMakeFiles/KS_2D_test1.dir/build.make
KS_2D/KS_2D_test1: /usr/lib/x86_64-linux-gnu/libpython3.8.so
KS_2D/KS_2D_test1: KS_2D/CMakeFiles/KS_2D_test1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable KS_2D_test1"
	cd /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KS_2D_test1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
KS_2D/CMakeFiles/KS_2D_test1.dir/build: KS_2D/KS_2D_test1

.PHONY : KS_2D/CMakeFiles/KS_2D_test1.dir/build

KS_2D/CMakeFiles/KS_2D_test1.dir/clean:
	cd /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D && $(CMAKE_COMMAND) -P CMakeFiles/KS_2D_test1.dir/cmake_clean.cmake
.PHONY : KS_2D/CMakeFiles/KS_2D_test1.dir/clean

KS_2D/CMakeFiles/KS_2D_test1.dir/depend:
	cd /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pangdongwen/文档/程序文件/有限元/KS模型_2维 /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/KS_2D /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D /home/pangdongwen/文档/程序文件/有限元/KS模型_2维/build/KS_2D/CMakeFiles/KS_2D_test1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : KS_2D/CMakeFiles/KS_2D_test1.dir/depend

