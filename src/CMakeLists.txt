# =========================================================================
# FLEXI Options
# =========================================================================

# Build directory list
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
SUBDIRLIST(EQNSYSTEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/equations/)
SUBDIRLIST(LIFTINGS   ${CMAKE_CURRENT_SOURCE_DIR}/src/dg/)
SUBDIRLIST(TESTCASES  ${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/)

# =========================================================================
# Node type
# =========================================================================
SET(FLEXI_NODETYPE "GAUSS" CACHE STRING "Set of nodes used for representing the solution (GAUSS/GAUSS-LOBATTO).")
SET_PROPERTY(CACHE FLEXI_NODETYPE PROPERTY STRINGS "GAUSS" "GAUSS-LOBATTO")
IF(FLEXI_NODETYPE STREQUAL "GAUSS")
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_DEFINITIONS(-DPP_NodeType=2)
ELSE()
  MESSAGE(SEND_ERROR "No such node type known.")
ENDIF()

# =========================================================================
# Equation system
# =========================================================================
SET(FLEXI_EQNSYSNAME "navierstokes" CACHE STRING
                     "Used equation system (navierstokes/linearscalaradvection/rans_sa)")
SET_PROPERTY(CACHE FLEXI_EQNSYSNAME  PROPERTY STRINGS ${EQNSYSTEMS})
IF(   FLEXI_EQNSYSNAME STREQUAL "navierstokes")
  ADD_DEFINITIONS(-DPP_nVar=5 -DPP_nVarPrim=6 -DEQNSYSNR=2)
  SET(WITHEOS 1)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "linearscalaradvection")
  ADD_DEFINITIONS(-DPP_nVar=1 -DPP_nVarPrim=1 -DEQNSYSNR=1)
  SET(WITHEOS 0)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "rans_sa")
  ADD_DEFINITIONS(-DPP_nVar=6 -DPP_nVarPrim=7 -DEQNSYSNR=3)
  SET(WITHEOS 1)
ELSE()
  MESSAGE(SEND_ERROR "Unknown equation system.")
ENDIF()

IF(WITHEOS GREATER 0)
  SET(FLEXI_EOS "idealgas" CACHE STRING "Used equation of state (idealgas/table)")
  SET_PROPERTY(CACHE FLEXI_EOS PROPERTY STRINGS idealgas)
  MARK_AS_ADVANCED(FORCE FLEXI_EOS)
ENDIF()

# =========================================================================
# SplitDG
# =========================================================================
OPTION(FLEXI_SPLIT_DG "Use split DG formulation" OFF)
IF(FLEXI_SPLIT_DG)
  ADD_DEFINITIONS(-DSPLIT_DG)
  IF(FLEXI_NODETYPE STREQUAL "GAUSS")
    MESSAGE(SEND_ERROR "FLEXI_SPLIT_DG only works on FLEXI_NODETYPE==GAUSS-LOBATTO points.")
  ENDIF()
  IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
    MESSAGE(SEND_ERROR "FLEXI_SPLIT_DG only works for FLEXI_EQNSYSNAME==navierstokes.")
  ENDIF()
ENDIF()

# =========================================================================
# PPLimiter
# =========================================================================
OPTION(FLEXI_PP_LIMITER "Enable Positivity Preserving Limiter" OFF)
IF(FLEXI_PP_LIMITER)
  ADD_DEFINITIONS(-DPP_LIMITER)
  IF(FLEXI_NODETYPE STREQUAL "GAUSS")
    MESSAGE(WARNING "Positivity Preserving Limiter only works on volume interpolation points. Positivity on faces is not guaranteed with FLEXI_NODETYPE==GAUSS." )
  ENDIF()
  IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
    MESSAGE(SEND_ERROR "FLEXI_PP_LIMITER only works for FLEXI_EQNSYSNAME==navierstokes.")
  ENDIF()
ENDIF()

# =========================================================================
# Polynomial degree
# =========================================================================
SET(FLEXI_POLYNOMIAL_DEGREE "N" CACHE STRING "N: compile with variable polynomial degree. 1,2,3,...: compile with fixed polynomial degree for performance (untested).")
ADD_DEFINITIONS(-DPP_N=${FLEXI_POLYNOMIAL_DEGREE})

# =========================================================================
# Parabolic terms
# =========================================================================
OPTION(FLEXI_PARABOLIC "Enable parabolic terms of equation-system" ON)
IF (FLEXI_PARABOLIC)
  ADD_DEFINITIONS(-DPARABOLIC=1)

  # Set the lifting
  LIST(APPEND FLEXI_LIFTING_LIST  "br1"
                                  "br2")
  SET(FLEXI_LIFTING "br1" CACHE STRING "Lifting method used for PARABOLIC")
  SET_PROPERTY(CACHE FLEXI_LIFTING PROPERTY STRINGS ${FLEXI_LIFTING_LIST})

  IF(    FLEXI_LIFTING STREQUAL "br1")
    ADD_DEFINITIONS(-DPP_Lifting=1)
  ELSEIF(FLEXI_LIFTING STREQUAL "br2")
    ADD_DEFINITIONS(-DPP_Lifting=2)
  ELSE()
    MESSAGE(SEND_ERROR "No such LIFTING known.")
  ENDIF()

  # Set the law for viscosity
  IF((FLEXI_EQNSYSNAME STREQUAL "navierstokes") OR (FLEXI_EQNSYSNAME STREQUAL "rans_sa"))
    SET(FLEXI_VISCOSITY "constant" CACHE STRING "Method to calculate viscosity (constant/sutherland/powerlaw)")
    SET_PROPERTY(CACHE FLEXI_VISCOSITY PROPERTY STRINGS constant sutherland powerlaw)
    IF(    FLEXI_VISCOSITY STREQUAL "constant")
      ADD_DEFINITIONS(-DPP_VISC=0)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "sutherland")
      ADD_DEFINITIONS(-DPP_VISC=1)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "powerlaw")
      ADD_DEFINITIONS(-DPP_VISC=2)
    ELSE()
      MESSAGE(SEND_ERROR "No such VISCOSITY known.")
    ENDIF()
    IF (FLEXI_EDDYVISCOSITY)
      IF(FLEXI_2D)
        MESSAGE(SEND_ERROR " FLEXI_EDDYVISCOSITY is enabled, please disable FLEXI_2D")
      ENDIF()
    ENDIF()
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DPARABOLIC=0)
ENDIF()

# =========================================================================
# Eddy viscosity models
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_EDDYVISCOSITY "Enable eddy viscosity modells" OFF
                                           "FLEXI_PARABOLIC"               OFF)
IF (FLEXI_EDDYVISCOSITY)
  ADD_DEFINITIONS(-DEDDYVISCOSITY=1)
ELSE()
  ADD_DEFINITIONS(-DEDDYVISCOSITY=0)
ENDIF()

# =========================================================================
# Dimension related preproc variables
# =========================================================================
OPTION(FLEXI_2D "Enable two dimensional computation" OFF)
IF(FLEXI_2D)
  ADD_DEFINITIONS(-DPP_dim=2)
ELSE()
  ADD_DEFINITIONS(-DPP_dim=3)
ENDIF()

# =========================================================================
# FV subcell shock capturing
# =========================================================================
SET(FLEXI_FV "OFF" CACHE STRING "Enable Finite Volume subcell shock capturing")
SET_PROPERTY(CACHE FLEXI_FV PROPERTY STRINGS OFF SWITCH BLEND)

CMAKE_DEPENDENT_OPTION(FLEXI_FV_RECONSTRUCTION "Enable 2nd order reconstruction for FV" ON
	                                             "FLEXI_FV" ON) # Everyhting except "OFF" is logically "ON"

IF(NOT (FLEXI_FV STREQUAL "OFF"))
  IF(    FLEXI_FV STREQUAL "SWITCH")
    ADD_DEFINITIONS(-DFV_ENABLED=1)
  ELSEIF(FLEXI_FV STREQUAL "BLEND")
    IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
      MESSAGE(SEND_ERROR "FLEXI_FV=BLEND only works for FLEXI_EQNSYSNAME==navierstokes.")
    ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS")
      MESSAGE(SEND_ERROR "FLEXI_FV=BLEND only works on FLEXI_NODETYPE==GAUSS-LOBATTO points.")
    ENDIF()
    ADD_DEFINITIONS(-DFV_ENABLED=2)
  ELSE()
    MESSAGE(SEND_ERROR "No such FV shock capturing known.")
  ENDIF()

  IF (FLEXI_FV_RECONSTRUCTION)
    ADD_DEFINITIONS(-DFV_RECONSTRUCT=1)
  ELSE()
    ADD_DEFINITIONS(-DFV_RECONSTRUCT=0)
    IF (FLEXI_PARABOLIC)
      MESSAGE(SEND_ERROR " FLEXI_PARABOLIC is enabled, please also enabled FLEXI_FV_RECONSTRUCTION or disable FLEXI_PARABOLIC")
    ENDIF()
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DFV_ENABLED=0)
  ADD_DEFINITIONS(-DFV_RECONSTRUCT=0)
ENDIF()

# =========================================================================
# Particle tracking
# =========================================================================
OPTION(PARTICLES "Enable particle tracking" OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_CODE_ANALYZE "Output code analyze messages"         OFF
                                              "PARTICLES"                            OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_LOADBALANCE  "Use LoadBalance timemeasurement"      OFF
                                              "PARTICLES;LIBS_USE_MPI"               OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_RW           "Use random walk model for particles"  OFF
                                              "PARTICLES;LIBS_USE_MPI"               OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_INTKIND8     "Enable particle numbers larger standard integer size" OFF
                                              "PARTICLES"                            OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_EXTEND_RHS   "Use extended RHS for particles"       OFF
                                              "PARTICLES;FLEXI_PARABOLIC"            OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_FAXEN_CORR   "Enable Faxen correctio for particles" OFF
                                              "PARTICLES;FLEXI_PARABOLIC"            OFF)
# Extend RHS options
CMAKE_DEPENDENT_OPTION(PARTICLES_UNDISTFLOW   "Use viscous stresses and pressure gradient in RHS for particles" OFF
                                              "PARTICLES_EXTEND_RHS"                 OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_VIRTUALMASS  "Use virtual mass force in RHS for particles" OFF
                                              "PARTICLES_EXTEND_RHS"                 OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_BASSETFORCE  "Use Basset force in RHS for particles" OFF
                                              "PARTICLES_EXTEND_RHS"                 OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_ROTATION     "Consider rotation of particles"       OFF
                                              "PARTICLES_EXTEND_RHS"                 OFF)
CMAKE_DEPENDENT_OPTION(PARTICLES_ANALYZE_RHS  "Write rhs to file for particles"      OFF
                                              "PARTICLES_EXTEND_RHS"                 OFF)
IF(PARTICLES)
  IF(FLEXI_2D)
    MESSAGE(SEND_ERROR "PARTICLES only works on FLEXI_2D=OFF")
  ENDIF()
  IF( NOT FLEXI_PARABOLIC)
    MESSAGE(SEND_ERROR "PARTICLES only works on FLEXI_PARABOLIC=ON")
  ENDIF()
  IF(FLEXI_FV AND PARTICLES_RW)
    MESSAGE(SEND_ERROR "PARTICLES_RW does not work with FLEXI_FV=ON")
  ENDIF()
  ADD_DEFINITIONS(-DUSE_PARTICLES=1)

  # Shared memory split type. Only available with OpenMPI
  SET(FLEXI_SPLIT_TYPE "$ENV{FLEXI_SPLIT_TYPE}")
  IF("${FLEXI_SPLIT_TYPE}" STREQUAL "")
    SET(FLEXI_SPLIT_TYPE "MPI_COMM_TYPE_SHARED")
  ELSE()
    MESSAGE(STATUS "Using user-defined environment variable [FLEXI_SPLIT_TYPE = ${FLEXI_SPLIT_TYPE}] for shared memory communicator splitting. Setting [FLEXI_SHARED_MEMORY = ${FLEXI_SPLIT_TYPE}]")
  ENDIF()
  SET(FLEXI_SHARED_MEMORY "${FLEXI_SPLIT_TYPE}" CACHE STRING "Split type for creating new communicators based on colors and keys (requires MPI 3 or higher). Options with the prefix OMPI_ are specific to Open MPI.")
  SET_PROPERTY(CACHE FLEXI_SHARED_MEMORY PROPERTY STRINGS MPI_COMM_TYPE_SHARED
                                                          OMPI_COMM_TYPE_CORE
                                                          FLEXI_COMM_TYPE_NODE)

  ADD_DEFINITIONS(-DSHARED_MEMORY_METHOD=${FLEXI_SHARED_MEMORY})
  IF(FLEXI_SHARED_MEMORY STREQUAL "MPI_COMM_TYPE_SHARED")
    UNSET(FLEXI_SHARED_MEMORY_CORES CACHE)
    ADD_DEFINITIONS(-DCORE_SPLIT=0)
    MESSAGE(STATUS "Shared memory split type for subcommunicators set to node-level")
  ELSEIF(FLEXI_SHARED_MEMORY STREQUAL "FLEXI_COMM_TYPE_NODE")
    SET(FLEXI_SHARED_MEMORY_CORES "2" CACHE STRING "Number of cores per node when setting FLEXI_SHARED_MEMORY=FLEXI_COMM_TYPE_NODE. All cores must be on the same physical node!")
    ADD_DEFINITIONS(-DCORE_SPLIT=${FLEXI_SHARED_MEMORY_CORES})
    MESSAGE(STATUS "Shared memory split type for subcommunicators set to sub-node-level with user-specific value [FLEXI_SHARED_MEMORY_CORES = ${FLEXI_SHARED_MEMORY_CORES}] cores per node")
  ELSEIF(FLEXI_SHARED_MEMORY STREQUAL "OMPI_COMM_TYPE_CORE")
    UNSET(FLEXI_SHARED_MEMORY_CORES CACHE)
    ADD_DEFINITIONS(-DCORE_SPLIT=1)
    MESSAGE(STATUS "Shared memory split type for subcommunicators set to core-level")
  ENDIF()

  # Code Analyze
  IF(PARTICLES_CODE_ANALYZE)
    ADD_DEFINITIONS(-DCODE_ANALYZE=1)
  ELSE()
    ADD_DEFINITIONS(-DCODE_ANALYZE=0)
  ENDIF()

  # Load Balance
  IF(PARTICLES_LOADBALANCE)
    MESSAGE(STATUS "Building Flexi with LoadBalance")
    ADD_DEFINITIONS(-DUSE_LOADBALANCE=1)
  ELSE()
    ADD_DEFINITIONS(-DUSE_LOADBALANCE=0)
  ENDIF()

  # Random Walk
  IF(PARTICLES_RW)
    ADD_DEFINITIONS(-DUSE_RW=1)
  ELSE()
    ADD_DEFINITIONS(-DUSE_RW=0)
  ENDIF()

  # INTEGER(KIND=8=
  IF (PARTICLES_INTKIND8)
    ADD_DEFINITIONS(-DUSE_INTKIND8=1)
  ELSE()
    ADD_DEFINITIONS(-DUSE_INTKIND8=0)
  ENDIF()

  IF (PARTICLES_FAXEN_CORR)
    ADD_DEFINITIONS(-DUSE_FAXEN_CORR=1)
  ELSE()
    ADD_DEFINITIONS(-DUSE_FAXEN_CORR=0)
  ENDIF()

  IF (PARTICLES_EXTEND_RHS)
    ADD_DEFINITIONS(-DUSE_EXTEND_RHS=1)
    # Undisturbed flow
    IF (PARTICLES_UNDISTFLOW)
      ADD_DEFINITIONS(-DUSE_UNDISTFLOW=1)
    ELSE()
      ADD_DEFINITIONS(-DUSE_UNDISTFLOW=0)
    ENDIF()
    # Virtual mass
    IF (PARTICLES_VIRTUALMASS)
      ADD_DEFINITIONS(-DUSE_VIRTUALMASS=1)
    ELSE()
      ADD_DEFINITIONS(-DUSE_VIRTUALMASS=0)
    ENDIF()
    # Basset force
    IF (PARTICLES_BASSETFORCE)
      ADD_DEFINITIONS(-DUSE_BASSETFORCE=1)
    ELSE()
      ADD_DEFINITIONS(-DUSE_BASSETFORCE=0)
    ENDIF()
    IF (PARTICLES_ROTATION)
      ADD_DEFINITIONS(-DPP_nVarPart=10)
      ADD_DEFINITIONS(-DPP_nVarPartRHS=6)
    ELSE()
      ADD_DEFINITIONS(-DPP_nVarPart=7)
      ADD_DEFINITIONS(-DPP_nVarPartRHS=3)
    ENDIF()
    # Write RHS to file for analyze
    IF (PARTICLES_ANALYZE_RHS)
      ADD_DEFINITIONS(-DANALYZE_RHS=1)
    ELSE()
      ADD_DEFINITIONS(-DANALYZE_RHS=0)
    ENDIF()
  ELSE()
    ADD_DEFINITIONS(-DPP_nVarPart=7)
    ADD_DEFINITIONS(-DPP_nVarPartRHS=3)
    ADD_DEFINITIONS(-DUSE_EXTEND_RHS=0)
    ADD_DEFINITIONS(-DUSE_UNDISTFLOW=0)
    ADD_DEFINITIONS(-DUSE_VIRTUALMASS=0)
    ADD_DEFINITIONS(-DUSE_BASSETFORCE=0)
    ADD_DEFINITIONS(-DANALYZE_RHS=0)
  ENDIF()

ELSE()
  ADD_DEFINITIONS(-DUSE_PARTICLES=0)
  ADD_DEFINITIONS(-DCODE_ANALYZE=0)
  ADD_DEFINITIONS(-DUSE_INTKIND8=0)
  ADD_DEFINITIONS(-DUSE_LOADBALANCE=0)
  ADD_DEFINITIONS(-DUSE_RW=0)
  ADD_DEFINITIONS(-DUSE_EXTEND_RHS=0)
  ADD_DEFINITIONS(-DUSE_FAXEN_CORR=0)
  ADD_DEFINITIONS(-DCORE_SPLIT=1)
ENDIF()

# =========================================================================
# Reduced lifting in performance mode
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_PERFORMANCE_OPTLIFT "Lift only needed variables" OFF
                                                 "FLEXI_PERFORMANCE;NOT POSTI" OFF)
IF(FLEXI_PERFORMANCE_OPTLIFT)
  ADD_DEFINITIONS(-DPP_OPTLIFT=1)
ELSE()
  ADD_DEFINITIONS(-DPP_OPTLIFT=0)
ENDIF()

# =========================================================================
# Test case
# =========================================================================
SET(   FLEXI_TESTCASE   "default"      CACHE STRING
                        "Name of testcase which is used. Must be identical to directory name in testcase folder.")
SET_PROPERTY(CACHE FLEXI_TESTCASE PROPERTY STRINGS ${TESTCASES})

# =========================================================================
# Unit test at end of compilation
# =========================================================================
IF (CMAKE_BUILD_TYPE MATCHES "Sanitize")
  MESSAGE(STATUS "Unit tests automatically disabled with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
  OPTION(FLEXI_UNITTESTS "Enable unit tests after build" OFF)
ELSE()
  OPTION(FLEXI_UNITTESTS "Enable unit tests after build" ON)
ENDIF()
MARK_AS_ADVANCED(FORCE FLEXI_UNITTESTS)
IF(FLEXI_UNITTESTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/unitTests/CMakeLists.txt)
ENDIF()

# =========================================================================
# Flexi sources & targets
# =========================================================================
FUNCTION(add_exec target)
  ADD_EXECUTABLE(       ${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} libflexishared)
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

FUNCTION(add_lib_shared target)
  ADD_LIBRARY(          ${target} SHARED ${ARGN})
  TARGET_LINK_LIBRARIES(${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.o libflexistatic userblocklib stacksizelib ${linkedlibs} ${particlelibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES PREFIX ""
                                             COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran)
  INSTALL(TARGETS ${target} LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                                                GROUP_READ GROUP_EXECUTE
                                                WORLD_READ WORLD_EXECUTE
                            DESTINATION lib)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

# include flexi.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

# collect all F90 files needed for flexilib, except the extra files below
FILE(GLOB_RECURSE flexiF90 ./src/interpolation/*.f90
                           ./src/readintools/*.f90
                           ./src/mortar/*.f90
                           ./src/mesh/*.f90
                           ./src/mpi/*.f90
                           ./src/globals/*.f90
                           ./src/analyze/*.f90
                           ./src/recordpoints/*.f90
                           ./src/restart/*.f90
                           ./src/output/*.f90
                           ./src/sponge/*.f90
                           ./src/filter/*.f90
                           ./src/io_hdf5/*.f90
                           ./src/timedisc/*.f90
                           ./src/precond/*.f90
                           ./src/implicit/*.f90
                           ./src/testcase/${FLEXI_TESTCASE}/*.f90
                           # ./src/flexilib.f90
                           ./unitTests/unittest.f90
                           ./unitTests/unittest_vars.f90)

# The following files/directories where excluded from the above global recursive search,
# since some of the subdirectories should only be included if the respective feature is enabled.

# >> collect F90 files of the equation system (non-recursive) in a new list 'eqnsysF90'
FILE(GLOB eqnsysF90        ./src/equations/${FLEXI_EQNSYSNAME}/*.f90)
# append F90 files of the equation of state to 'eqnsysF90'-list
IF (WITHEOS GREATER 0)
  IF(POSTI)
    FILE(GLOB eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/*.f90
              eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/posti/*.f90)
  ELSE()
    FILE(GLOB eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/*.f90)
  ENDIF()
  LIST(APPEND eqnsysF90 ${eosF90})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/)
ELSE()
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/)
ENDIF()

# >> append F90 files of the eddy viscosity to 'eqnsysF90'-list
IF (FLEXI_EDDYVISCOSITY)
  FILE(GLOB eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/default/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/smagorinsky/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/vreman/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/sigmaModel/*.f90)
  LIST(APPEND eqnsysF90 ${eddyviscF90})
ENDIF()

# >> append 'eqnsysF90'-list to 'flexiF90'-list
LIST(APPEND flexiF90 ${eqnsysF90})

# >> collect F90 files of 'dg' directory (non-recursive) in a new list 'dgF90'
FILE(GLOB dgF90            ./src/dg/*.f90)

# >> append F90 files of lifting to 'dgF90'-list
IF (FLEXI_PARABOLIC)
  #IF (PARTICLES_FAXEN_CORR)
    FILE(GLOB liftingF90     ./src/dg/lifting/*.f90
                             ./src/lifting_gen/*.f90)
  #ELSE()
   # FILE(GLOB liftingF90     ./src/dg/lifting/*.f90)
  #ENDIF()
  LIST(APPEND dgF90 ${liftingF90})
ENDIF()

# >> append 'dgF90'-list to 'flexiF90'-list
LIST(APPEND flexiF90 ${dgF90})

# >> collect F90 files of 'fv' directory (non-recursive) and append to 'flexiF90' list
IF (FLEXI_FV)
  FILE(GLOB fvF90 ./src/indicator/*.f90
                  ./src/fv/*.f90)
  LIST(APPEND flexiF90 ${fvF90})
ENDIF()

# >> collect F90 files of 'hash' directory (non-recursive) and append to 'flexiF90' list
FILE(GLOB   hashF90 ./src/hash/*.f90)
LIST(APPEND flexiF90 ${hashF90})
# GCC<10 has incorrect handling of typeless constants as per https://gcc.gnu.org/gcc-10/changes.html
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  IF(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "10.0.0")
    MESSAGE (STATUS "GFortran v${CMAKE_Fortran_COMPILER_VERSION} with invalid handling of typeless constants detected, setting -fno-range-check")
    SET_SOURCE_FILES_PROPERTIES(./src/hash/sha256.f90 PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays -fno-range-check -Wno-conversion")
  ENDIF()
ENDIF()

# collect F90 files of 'particle' directory (recursive) in new list 'particlesF90'
IF (PARTICLES)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/particles)
  FILE(GLOB_RECURSE particlesF90 ./src/particles/*.f90)
  LIST(APPEND flexiF90 ${particlesF90})

  ADD_LIBRARY(memusagelib STATIC ./src/particles/globals/processmemusage.cxx)
  SET_TARGET_PROPERTIES(memusagelib PROPERTIES LANGUAGE CXX LINKER_LANGUAGE CXX
                                               COMPILE_FLAGS "-fPIC")
  LIST(APPEND particlelibs memusagelib)
ENDIF()


# define libs
ADD_LIBRARY(          libflexistaticF90 OBJECT ${flexiF90})
SET_TARGET_PROPERTIES(libflexistaticF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(          libflexistatic STATIC $<TARGET_OBJECTS:libflexistaticF90>)
SET_TARGET_PROPERTIES(libflexistatic PROPERTIES PREFIX ""
                                     OUTPUT_NAME "libflexif90")
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(libflexistaticF90 ${SELFBUILTEXTERNALS})
ENDIF()

ADD_LIBRARY(userblocklib STATIC ./src/output/read_userblock.c)
ADD_LIBRARY(stacksizelib STATIC ./src/globals/stacksize.c)
SET_TARGET_PROPERTIES(stacksizelib userblocklib PROPERTIES LANGUAGE C LINKER_LANGUAGE C
                                                           COMPILE_FLAGS "-fPIC")

# we link the statically built libs
ADD_LIB_SHARED(       libflexishared ./src/flexilib.f90)
SET_TARGET_PROPERTIES(libflexishared PROPERTIES OUTPUT_NAME "libflexi")
ADD_DEPENDENCIES(     libflexishared libflexistatic userblocklib stacksizelib ${INTERNALLIBS})
ADD_EXEC(        flexi ./src/flexi.f90)
ADD_DEPENDENCIES(flexi libflexistatic userblocklib stacksizelib ${INTERNALLIBS})

# special compiler flags to avoid warnings for empty functions of testcase.f90 and overwrite of intrinsic abort in globals.f90
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLEXI_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-Wno-unused-dummy-argument")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals.f90 PROPERTIES COMPILE_FLAGS "-Wno-intrinsic-shadow")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLEXI_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-warn nounused")
ENDIF()

# =========================================================================
# Ctags
# =========================================================================
MARK_AS_ADVANCED(FORCE CTAGS_PATH)
FIND_PROGRAM(CTAGS_PATH ctags-exuberant)
IF(CTAGS_PATH)
  SET(             EQNS_WO_ACTUAL ${EQNSYSTEMS})
  LIST(REMOVE_ITEM EQNS_WO_ACTUAL ${FLEXI_EQNSYSNAME})
  SET(EXLUDES "")
  FOREACH(item ${EQNS_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/equations/${item}/*" ${EXLUDES})
  ENDFOREACH()

  IF(FLEXI_PARABOLIC)
    SET(             LIFTINGS_WO_ACTUAL ${LIFTINGS})
    LIST(REMOVE_ITEM LIFTINGS_WO_ACTUAL ${FLEXI_LIFTING})
    FOREACH(item ${LIFTINGS_WO_ACTUAL})
       SET(EXLUDES "--exclude=src/dg/${item}/*" ${EXLUDES})
    ENDFOREACH()
  ENDIF()

  SET(             TESTCASES_WO_ACTUAL ${TESTCASES})
  LIST(REMOVE_ITEM TESTCASES_WO_ACTUAL ${FLEXI_TESTCASE})
  FOREACH(item ${TESTCASES_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/testcase/${item}/*" ${EXLUDES})
  ENDFOREACH()

  ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} src && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} --append posti)
  ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
  MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
  ADD_DEPENDENCIES(flexi tags)
ENDIF(CTAGS_PATH)
