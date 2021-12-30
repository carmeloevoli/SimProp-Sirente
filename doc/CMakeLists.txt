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

# find_package(Doxygen REQUIRED)
# find_package(Sphinx REQUIRED)

# # Find all the public headers
# get_target_property(SIMPROP_PUBLIC_HEADER_DIR simprop INTERFACE_INCLUDE_DIRECTORIES)
# file(GLOB_RECURSE SIMPROP_PUBLIC_HEADERS 
# 	${SIMPROP_PUBLIC_HEADER_DIR}/*.h)

# set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/SimProp)
# set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
# set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
# set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
# set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# # Set the Doxygen input and output directories in the Doxyfile
# configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# # Doxygen won't create this for us
# file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) 

# # Only regenerate Doxygen when the Doxyfile or public headers change
# add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
# 		DEPENDS ${SIMPROP_PUBLIC_HEADERS}
#         	COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
# 		MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
#         	COMMENT "Generating docs"
# 		VERBATIM)

# # Nice named target so we can run the job easily
# add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

# set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
# set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
# set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

# # Only regenerate Sphinx when:
# #  - Doxygen has rerun 
# #  - Our doc files have been updated
# #  - The Sphinx config has been updated
# add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
# 	COMMAND 
# 		${SPHINX_EXECUTABLE} -b html
# 		# Tell Breathe where to find the Doxygen output
# 		-Dbreathe_projects.SimProp=${DOXYGEN_OUTPUT_DIR}/xml
# 		${SPHINX_SOURCE} ${SPHINX_BUILD}
# 	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
# 	DEPENDS 
# 		# Other docs files you want to track should go here (or in some variable)
# 		${CMAKE_CURRENT_SOURCE_DIR}/index.rst
# 		${DOXYGEN_INDEX_FILE}
# 	MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
# 	COMMENT "Generating documentation with Sphinx")

# # Nice named target so we can run the job easily
# add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# include(GNUInstallDirs)
# install(DIRECTORY ${SPHINX_BUILD}
# 	DESTINATION ${CMAKE_INSTALL_DOCDIR})
