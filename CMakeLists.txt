# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Optimal_transport_2)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
# find_package(UMFPACK REQUIRED)
# include_directories(${GSL_INCLUDE_DIR})
# include_directories(${UMF_INCLUDE_DIR})

 link_libraries(${GSL_LIBRARY})

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

IF(UNIX)
  IF(CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math")
    ADD_DEFINITIONS(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
  ENDIF(CMAKE_COMPILER_IS_GNUCC)
ENDIF(UNIX)

find_package(CGAL COMPONENTS Qt4)
find_package(Boost COMPONENTS filesystem system serialization program_options FATAL_ERROR)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(Qt4)
include_directories (BEFORE
  ./include
  ./third-party
  ./third-party/liblbfgs/include 
  #  ./third-party/huebner
  # ./third-party/ANN/include
  # ./third-party/quadprog
)
include(${QT_USE_FILE})

#--------------------------------
# The "Regular" demo: Regular_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES Optimal_transportation.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES ./resources/Optimal_transportation.qrc
  ./resources/Optimal_transportation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( "Optimal_transportation.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Optimal_transportation.moc" )

set(OT_SRCS 
  third-party/liblbfgs/lib/lbfgs.c
  third-party/huebner/bundle_method.c
  third-party/huebner/cqpminimizer.c
  third-party/huebner/initial_point.c
  third-party/huebner/convergence_fsdf.c
  third-party/huebner/fsdfminimizer.c
  third-party/huebner/qp_pdip.c
  third-party/huebner/mg_pdip.c
  third-party/quadprog/uQuadProg++.cc)

# SET(MUMPS_LIBRARY ${MUMPS_LIBRARY} -ldmumps)
# link_libraries(${MUMPS_LIBRARY})
# link_libraries(-lumfpack)
add_library(OT SHARED ${OT_SRCS})
link_libraries(OT)
#link_libraries(ANN)
link_libraries(${QT_LIBRARIES}  ${Boost_SERIALIZATION_LIBRARY})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
link_libraries(${Boost_PROGRAM_OPTIONS_LIBRARY})
link_libraries(${GSL_LIBRARY})
link_libraries(${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

# ENABLE_TESTING()
add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(third-party)

# The executable itself.
add_executable  (Optimal_transportation
  Optimal_transportation.cpp
  Optimal_transportation.moc
  ${DT_UI_FILES} ${DT_RESOURCE_FILES} )
