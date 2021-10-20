# Set the minimum required version of cmake for a project.
cmake_minimum_required(VERSION 3.10)

# Set Variables
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)
set(PROJECT_NAME "calculator")

set(LIBIGL_INCLUDE_DIR /home/sungho/data/libigl/v2.2.0/include)
set(LIBIGL_DIR ${LIBIGL_INCLUDE_DIR}/..)
set(MARIADB_INCLUDE_DIR /usr/local/include/mariadb/)
set(MARIADB_LIB_DIR /usr/local/lib/mariadb/)
set(LIBSSH_INCLUDE_DIR /usr/local/include/libssh/)


# Project
PROJECT(${PROJECT_NAME} CXX)



#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
# libigl
option(LIBIGL_WITH_TETGEN            "Use TetGen"         ON)
find_package(LIBIGL REQUIRED QUIET)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
		    #${MARIADB_INCLUDE_DIR}
		    #${LIBSSH_INCLUDE_DIR}
		    ${LIBIGL_DIR}
                    ${LIBIGL_DIR}/external
                    ${LIBIGL_DIR}/external/glad/include
		    )

# File
FILE(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
FILE(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
	               
# openMP
FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


# MariaDB
FIND_LIBRARY(MARIADB_LIB
        NAMES mariadb
        PATHS ${MARIADB_LIB_DIR}
        )

# libssh
find_package(libssh REQUIRED)
if (NOT TARGET ssh)
     add_library(ssh SHARED IMPORTED)
     set_target_properties(ssh
                           PROPERTIES
                           IMPORTED_LOCATION "${LIBSSH_LIBRARIES}"
                           INTERFACE_INCLUDE_DIRECTORIES "${LIBSSH_INCLUDE_DIR}")
endif()

# The main program
add_executable(${PROJECT_NAME} calculator.cpp ${sources} ${headers})

# Link
target_link_libraries( ${PROJECT_NAME} 
	${MARIADB_LIB} 
	ssh 
	${LIBSSH_LIB} 
	${Geant4_LIBRARIES}
	igl::core
	igl::tetgen
	)
