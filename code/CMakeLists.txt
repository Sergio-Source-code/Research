# CMake entry point
cmake_minimum_required (VERSION 2.6)
project (LocalMeshOpt)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

#find_package(OpenGL REQUIRED)
#find_package(Boost QUIET COMPONENTS program_options)
#if(Boost_PROGRAM_OPTIONS_FOUND)
#    if( CGAL_AUTO_LINK_ENABLED )
#        message( STATUS "Boost.ProgramOptions library: found" )
#    else()
#        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
#    endif()
#    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
#    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_LIBRARIES})
#endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(
    /usr/include/
    /usr/include/eigen3
    /usr/local/include/
    /usr/local/include/vtk-9.2/
	./src
    ./verdict/vtkverdict
    .
)

set(ALL_LIBS
#	${OPENGL_LIBRARY}
#	GL
#    boost_program_options
    pthread
    ${VTK_LIBRARIES}
)

add_definitions(
#	-D_CRT_SECURE_NO_WARNINGS
)

set(SRC_OPTIMIZE
    optimize.cpp
    src/Solver.cpp
    src/Solver.h
	src/Mesh.cpp
	src/Mesh.h
	src/MeshFileReader.cpp
	src/MeshFileReader.h
	src/MeshFileWriter.cpp
	src/MeshFileWriter.h
)
add_executable(optimize ${SRC_OPTIMIZE})
target_link_libraries(optimize ${ALL_LIBS})