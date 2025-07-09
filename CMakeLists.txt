cmake_minimum_required(VERSION 3.10)

project(FlexAID)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(MSVC)
    message(STATUS "Configuring for Windows with MSVC")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} /W4 /O2 /MT")
else()
    message(STATUS "Configuring for Unix-like system (Linux, macOS, etc.)")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -O3 -static")
endif()

# Boost path
set(BOOST_INCLUDE_PATH "C:/boost_1_86_0")
include_directories(${BOOST_INCLUDE_PATH})

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/LIB)

# Add source files
set(SOURCES
        LIB/top.cpp
        LIB/boinc.cpp
        LIB/assign_shift.cpp
        LIB/buildic.cpp
        LIB/Vcontacts.cpp
        LIB/ic2cf.cpp
        LIB/vcfunction.cpp
        LIB/create_rebuild_list.cpp
        LIB/read_coor.cpp
        LIB/read_constraints.cpp
        LIB/assign_constraint.cpp
        LIB/modify_pdb.cpp
        LIB/update_constraint.cpp
        LIB/update_optres.cpp
        LIB/update_bonded.cpp
        LIB/read_grid.cpp
        LIB/read_input.cpp
        LIB/read_lig.cpp
        LIB/read_pdb.cpp
        LIB/add2_optimiz_vec.cpp
        LIB/assign_radius.cpp
        LIB/assign_radii.cpp
        LIB/assign_radii_types.cpp
        LIB/assign_types.cpp
        LIB/bondedlist.cpp
        LIB/buildcc.cpp
        LIB/buildlist.cpp
        LIB/geometry.cpp
        LIB/read_conect.cpp
        LIB/residue_conect.cpp
        LIB/wif083.cpp
        LIB/calc_cleftic.cpp
        LIB/write_pdb.cpp
        LIB/gaboom.cpp
        LIB/calc_center.cpp
        LIB/ic_bounds.cpp
        LIB/buildic_point.cpp
        LIB/buildcc_point.cpp
        LIB/read_rmsdst.cpp
        LIB/calc_rmsd_chrom.cpp
        LIB/calc_rmsd.cpp
        LIB/write_rrd.cpp
        LIB/write_sphere.cpp
        LIB/write_grid.cpp
        LIB/write_rrg.cpp
        LIB/read_rotlib.cpp
        LIB/number_of_dihedrals.cpp
        LIB/read_flexscfile.cpp
        LIB/read_rotobs.cpp
        LIB/set_intprob.cpp
        LIB/build_rotamers.cpp
        LIB/check_clash.cpp
        LIB/read_normalgrid.cpp
        LIB/read_eigen.cpp
        LIB/assign_eigen.cpp
        LIB/print_surfmat.cpp
        LIB/alter_mode.cpp
        LIB/dee_pivot.cpp
        LIB/read_emat.cpp
        LIB/read_spheres.cpp
        LIB/cluster.cpp
        LIB/rna_structure.cpp
        LIB/shortest_path.cpp
        LIB/assign_shortflex.cpp
        LIB/slice_grid.cpp
        LIB/partition_grid.cpp
        LIB/generate_grid.cpp
        LIB/maps.cpp
)

# Add executable
add_executable(FlexAID ${SOURCES})

target_link_libraries(FlexAID ${Boost_LIBRARIES})
