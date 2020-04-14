paraview_add_plugin(visuReader${POSTI_PV_SUFFIX}
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES visuReaderModule
  MODULE_FILES "${TOOL_DIR_PV}/Reader/vtk.module")
