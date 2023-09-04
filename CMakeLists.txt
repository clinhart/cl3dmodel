cmake_minimum_required(VERSION 3.9)
project(cl3dmodel VERSION 0.0.1 DESCRIPTION "3D model using a flavor of DCEL with quarter-edge-references")
include(GNUInstallDirs)
add_library(cl3dmodel src/cl3dmodel.cpp)
set_target_properties(cl3dmodel PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER api/cl3dmodel.h)
#configure_file(cl3dmodel.pc.in cl3dmodel.pc @ONLY)
target_include_directories(cl3dmodel PRIVATE .)
install(TARGETS cl3dmodel
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(FILES ${CMAKE_BINARY_DIR}/cl3dmodel.pc
#    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


