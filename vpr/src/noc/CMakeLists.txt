
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_NOC_SOURCES)

add_library(libvprnoc STATIC
    ${VPR_NOC_SOURCES}
)

set_target_properties(libvprnoc PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprnoc PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprnoc
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprnoc PRIVATE ${VPR_COMPILE_DEFINITIONS})

