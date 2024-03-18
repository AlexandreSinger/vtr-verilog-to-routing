
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_TIMING_SOURCES)

add_library(libvprtiming STATIC
    ${VPR_TIMING_SOURCES}
)

set_target_properties(libvprtiming PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprtiming PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprtiming
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprtiming PRIVATE ${VPR_COMPILE_DEFINITIONS})

