# =========================================================================
# Turbulence calculation
# =========================================================================
SET(TOOL_NAME "calcturb")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/calcturb")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/posti_calcturb.f90
                          ${TOOL_DIR}/calcturb.f90
                          ${TOOL_DIR}/calcturb_exactfunc.f90
                          ${TOOL_DIR}/calcturb_io.f90
                          ${TOOL_DIR}/calcturb_restart.f90
                          ${TOOL_DIR}/calcturb_vars.f90
                          ${TOOL_DIR}/calcturb_visu.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared)

# Standalone Posti executable
ADD_EXEC(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_calcturb.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${SELFBUILTEXTERNALS})
