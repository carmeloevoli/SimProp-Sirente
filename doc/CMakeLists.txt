find_package(Doxygen)
if(DOXYGEN_FOUND)
	MESSAGE(STATUS "Found Doxygen to build documentation")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml COPYONLY)
    add_custom_target(doxy ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)

    set_source_files_properties(${CMAKE_BINARY_DIR}/docstrings_from_doxy.i PROPERTIES GENERATED true )
    add_custom_target(docstrings_from_doxy
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/libs/doxy2swig/doxy2swig.py ${CMAKE_CURRENT_BINARY_DIR}/xml/index.xml -o ${CMAKE_CURRENT_BINARY_DIR}/docstrings_from_doxy.i
      DEPENDS doxy
    )
else()
	MESSAGE(STATUS "Doxygen not found: cannot build documentation")
endif(DOXYGEN_FOUND)