cmake_minimum_required(VERSION 3.16)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)


# # use linalg from NGSolve, otherwise comment it out
# find_package(NGSolve CONFIG REQUIRED
#   HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
# )

find_package(LAPACK REQUIRED)

include_directories(src "lib/PEP-BLA/src" "lib/PEP-BLA/lib/PEP-HPC/src" "lib/PEP-BLA/lib/PEP-HPC/concurrentqueue")
add_subdirectory (src)

add_executable (test_ode demos/test_ode.cc)
target_link_libraries (test_ode PUBLIC ${LAPACK_LIBRARIES})
# target_link_libraries (test_ode PUBLIC ngbla)

add_executable (test_newmark demos/test_newmark.cc)
target_link_libraries (test_newmark PUBLIC ${LAPACK_LIBRARIES})
# target_link_libraries (test_newmark PUBLIC ngbla)

add_executable (test_alpha demos/test_alpha.cc)
target_link_libraries (test_alpha PUBLIC ${LAPACK_LIBRARIES})
# target_link_libraries (test_alpha PUBLIC ngbla)


# add_subdirectory (mass_spring)
