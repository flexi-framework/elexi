set(classes
  visuReader)
  
set(headers
  ../../plugin_visu.h
  ../../pluginTypes_visu.h)

vtk_module_add_module(visuReaderModule
  CLASSES ${classes}
  HEADERS ${headers})
  
vtk_module_link(
  visuReaderModule PUBLIC ${TOOL_NAME}lib libflexishared)
  
paraview_add_server_manager_xmls(
  XMLS  visuReader.xml)
