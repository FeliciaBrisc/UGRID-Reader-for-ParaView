cmake_minimum_required(VERSION 2.8.8)

IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

#add_executable( UGRIDReader MACOSX_BUNDLE UGRIDReader.cxx )

ADD_PARAVIEW_PLUGIN(UGRIDReader "0.8.0"
                    SERVER_MANAGER_XML UGRIDServerManagerXml.xml
                    SERVER_MANAGER_SOURCES UGRIDReader.cxx)

target_link_libraries( UGRIDReader PRIVATE ${VTK_LIBRARIES} ${vtkNetCDF_LIBRARIES} ${vtkNetCDF_cxx_LIBRARIES} )


