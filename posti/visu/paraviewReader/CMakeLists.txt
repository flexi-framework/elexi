paraview_add_plugin(visuReader
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES visuReaderModule
  MODULE_FILES "${TOOL_DIR_PV}/Reader/vtk.module")
