
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_BASE_SOURCES)

add_library(libvprbase STATIC
    ${VPR_BASE_SOURCES}
)

set_target_properties(libvprbase PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprbase PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprbase
    ${VPR_DEPENDENCIES}
    libvpranalysis
    libvprdraw
    libvprnoc
    libvprpack
    libvprplace
    libvprpower
    libvprroute
    libvprtiming
    libvprutil
)
target_compile_definitions(libvprbase PRIVATE ${VPR_COMPILE_DEFINITIONS})

