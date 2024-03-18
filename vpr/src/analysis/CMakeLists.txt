
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_ANALYSIS_SOURCES)

add_library(libvpranalysis STATIC
    ${VPR_ANALYSIS_SOURCES}
)

set_target_properties(libvpranalysis PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvpranalysis PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvpranalysis
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvpranalysis PRIVATE ${VPR_COMPILE_DEFINITIONS})

