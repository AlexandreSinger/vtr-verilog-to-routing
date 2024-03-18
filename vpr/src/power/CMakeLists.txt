
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_POWER_SOURCES)

add_library(libvprpower STATIC
    ${VPR_POWER_SOURCES}
)

set_target_properties(libvprpower PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprpower PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprpower
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprpower PRIVATE ${VPR_COMPILE_DEFINITIONS})

