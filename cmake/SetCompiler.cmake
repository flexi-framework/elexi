# =========================================================================
# After settings specific compilers, enable named languages for cmake
# =========================================================================
ENABLE_LANGUAGE(Fortran C CXX)
MARK_AS_ADVANCED(FORCE C_PATH CXX_PATH Fortran_PATH)

# =========================================================================
# Set machine-specific definitions and settings
# =========================================================================
# HLRS HAWK
IF (CMAKE_FQDN_HOST MATCHES "hawk\.hww\.hlrs\.de$")
  MESSAGE(STATUS "Compiling on Hawk")
  # Overwrite compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    SET(FLEXI_INSTRUCTION "-march=znver2 -mtune=znver2")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xCORE-AVX2")
  ENDIF()
  # Use AMD Optimized Lapack/BLAS
  # SET(BLA_VENDOR "FLAME")
  # Set LUSTRE definition to account for filesystem and MPI implementation
  ADD_DEFINITIONS(-DLUSTRE)

# SuperMUC
ELSEIF (CMAKE_FQDN_HOST MATCHES "sng\.lrz\.de$")
  MESSAGE(STATUS "Compiling on SuperMUC")
  # Overwrite compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    SET(FLEXI_INSTRUCTION "-march=skylake-avx512 -mtune=skylake-avx512")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xSKYLAKE-AVX512")
    # Explicitely enable usage of AVX512 registers
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-zmm-usage=high")
  ENDIF()
  # Set LUSTRE definition to account for filesystem and MPI implementation
  ADD_DEFINITIONS(-DLUSTRE)

# LUMI
ELSEIF(CMAKE_FQDN_HOST MATCHES "\.can$")
  MESSAGE(STATUS "Compiling on LUMI")
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(FLEXI_INSTRUCTION "-march=znver3 -mtune=znver3")
  ELSE()
    MESSAGE(FATAL_ERROR "LUMI currently only supported using the GNU Compiler Collection (GCC). Please load/swap the following modules: LUMI PrgEnv-gnu cray-hdf5-parallel")
  ENDIF()

# IAG Prandtl
ELSEIF(CMAKE_FQDN_HOST MATCHES "^(prandtl|grafik.*)\.iag\.uni\-stuttgart\.de")
  MESSAGE(STATUS "Compiling on ${CMAKE_HOSTNAME}")
  SET(FLEXI_INSTRUCTION "-march=native -mtune=native")
  # Set LUSTRE definition to account for filesystem
  ADD_DEFINITIONS(-DLUSTRE)

ELSEIF (CMAKE_FQDN_HOST MATCHES "^ila(head.*|cfd.*)\.ila.uni\-stuttgart\.de")
  MESSAGE(STATUS "Compiling on ILA cluster")
  # Overwrite compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    SET(FLEXI_INSTRUCTION "-march=core-avx2 -mtune=core-avx2")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xCORE-AVX2")
  ENDIF()
  # Work around MPI-IO issue 4446 on machines mounting storage via NFS
  ADD_DEFINITIONS(-DNFS)

ELSEIF (CMAKE_FQDN_HOST MATCHES "^(xenon.*|argon.*)\.ila.uni\-stuttgart\.de")
  MESSAGE(STATUS "Compiling on ILA student cluster")
  SET(FLEXI_INSTRUCTION "-march=native -mtune=native")
  # Work around MPI-IO issue 4446 on machines mountng storage via NFS
  ADD_DEFINITIONS(-DNFS)

ELSEIF ("${CMAKE_FQDN_HOST}" MATCHES "gitlab\.ila\.uni\-stuttgart\.de")
  MESSAGE(STATUS "Compiling on ILA Gitlab")
  ADD_DEFINITIONS(-DVDM_ANALYTICAL)
  SET(FLEXI_INSTRUCTION "-march=native -mtune=native")

ELSE()
  MESSAGE(STATUS "Compiling on a generic machine [${CMAKE_HOSTNAME}]")
  # Set compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang" OR CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    SET(FLEXI_INSTRUCTION "-march=native -mtune=native")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xHost")
  ENDIF()
ENDIF()

# =========================================================================
# CHECK SUPPORT FOR VARIOUS FORTRAN (2003,2008) FEATURES
# =========================================================================
INCLUDE(CheckFortranSourceCompiles)
CHECK_FORTRAN_SOURCE_COMPILES(
"MODULE F03MOD
IMPLICIT NONE
PRIVATE
CONTAINS
! Check for MOVE_ALLOC feature
SUBROUTINE F03_MOVE_ALLOC()
REAL,ALLOCATABLE::a(:,:),b(:,:)
ALLOCATE(a(3,3))
CALL MOVE_ALLOC(a,b)
END SUBROUTINE F03_MOVE_ALLOC
END MODULE F03MOD

PROGRAM F03PROG
END"
Fortran2003Check
SRC_EXT F90)

IF(NOT Fortran2003Check)
  MESSAGE(FATAL_ERROR "Failed to compile basic Fortran2003 programm! Please ensure your compiler is up-to-date!")
ENDIF()

MESSAGE(STATUS "Compiling with [${CMAKE_Fortran_COMPILER_ID}] (v${CMAKE_Fortran_COMPILER_VERSION}) fortran compiler using [${FLEXI_INSTRUCTION}] instruction")

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# CMake can always request position independent code
# SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# GNU Compiler Collection (GCC)
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set Flags (disable lto type warnings due to false positives with MATMUL, which is a known bug)
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -finit-real=snan -finit-integer=snan -Wno-lto-type-mismatch -lstdc++ -DGNU")
    # LUMI has an issue with argument types in MPI(CH) calls
    IF(CMAKE_FQDN_HOST MATCHES "\.can$")
      SET (CMAKE_Fortran_FLAGS            "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    ENDIF()
  ENDIF()
  # initialize all variables as signalling NaNs to force the user to correctly initialize these data types
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS} -g  -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays -ffpe-trap=invalid,zero,overflow -fbacktrace")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g  -Og -ggdb3 -ffpe-trap=invalid,zero,overflow -fbounds-check -fbacktrace -Wall")
  SET (CMAKE_Fortran_FLAGS_SANITIZE       "${CMAKE_Fortran_FLAGS} -g  -Og -ggdb3 -ffpe-trap=invalid,zero,overflow,denorm -fbounds-check -fbacktrace -Wall -fsanitize=address,undefined,leak -fno-omit-frame-pointer -Wc-binding-type -Wuninitialized -pedantic")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (FLEXI_COMPILE_FLAGS "-xf95-cpp-input")
  ELSE()
    # Trailing white space required in case variable is unset!
    SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} ")
  ENDIF()

# AMD Optimized LLVM/CLANG
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -std=f2008 -lstdc++ -DFLANG")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions ")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions -ffpe-trap=invalid,zero,overflow -fbacktrace")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLEXI_INSTRUCTION} -finline-functions ")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g  -O0 -ggdb3 -ffpe-trap=invalid,zero,overflow -fbounds-check -finit-real=snan -fbacktrace -Wall")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (FLEXI_COMPILE_FLAGS "-xf95-cpp-input")
  ELSE()
    # Trailing white space required in case variable is unset!
    SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} ")
  ENDIF()

# Intel Compiler
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -lstdc++ -DINTEL")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}    -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS}    -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div -fpe0 -traceback")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -p -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (FLEXI_COMPILE_FLAGS "-fpp -allow nofpp_comments -assume bscc")
  ELSE()
    SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} -allow nofpp_comments -assume bscc")
  ENDIF()

# Cray Compiler
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -ffree -s real64 -s integer64 -em -lstdc++ -hfp0 -DCRAY")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}    -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS}    -O2 -hfp3 -p . -rm -eD")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS}    -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} -F")
ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# Profile-Guided Optimization (PGO)
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_PERFORMANCE_PGO "Enable profile-guided optimization (Only GNU Compiler supported)" OFF
                                             "FLEXI_PERFORMANCE" OFF)
IF (FLEXI_PERFORMANCE_PGO)
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}       -fprofile-use")
    SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_REWITHDEBINFO} -fprofile-use")
    SET(CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS_PROFILE}       -fprofile-generate")
  ELSE()
    MESSAGE(SEND_ERROR "Profile-guided optimization (PGO) currently only supported for GNU compiler. Either set FLEXI_PERFORMANCE_PGO=OFF or use the GNU compiler." )
  ENDIF()
ENDIF()

# Save the current compiler flags to the cache every time cmake configures the project.
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_RELEASE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_PROFILE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_DEBUG)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_SANITIZE)
SET(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRING "Default compiler flags"        FORCE)
SET(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRING "Release compiler flags"        FORCE)
SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "RelWithDebInfo compiler flags" FORCE)
SET(CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS_PROFILE}"        CACHE STRING "Profile compiler flags"        FORCE)
SET(CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRING "Debug compiler flags"          FORCE)
SET(CMAKE_Fortran_FLAGS_SANITIZE       "${CMAKE_Fortran_FLAGS_SANITIZE}"       CACHE STRING "Sanitize compiler flags"       FORCE)
