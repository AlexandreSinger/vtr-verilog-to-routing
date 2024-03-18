
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_PLACE_SOURCES)

add_library(libvprplace STATIC
    ${VPR_PLACE_SOURCES}
)

set_target_properties(libvprplace PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprplace PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprplace
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprplace PRIVATE ${VPR_COMPILE_DEFINITIONS})

