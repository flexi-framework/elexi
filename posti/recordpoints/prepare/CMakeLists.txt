# =========================================================================
# PREPARE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "preparerecordpoints")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/prepare)

MESSAGE(STATUS "Posti: using prepare record points tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/parameters.f90
                          ${TOOL_DIR}/parametriccoords.f90
                          ${TOOL_DIR}/rpset.f90
                          ${TOOL_DIR}/rpset_vars.f90
                          ${TOOL_DIR}/blprojection.f90
                          ${TOOL_DIR}/hdf5_outputrp.f90
                          ${TOOL_DIR}/visurp.f90 )

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                                    COMPILE_OPTIONS "-fPIC"
                                                    CMAKE_POSITION_INDEPENDENT_CODE ON)
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib libposti)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared libposti)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(${TOOL_NAME}target)
ADD_DEPENDENCIES(${TOOL_NAME}target libflexishared)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 ${TOOL_NAME}target)

# Standalone Posti executable
ADD_EXEC(             posti_${TOOL_NAME} ${TOOL_DIR}/preparerecordpoints.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
