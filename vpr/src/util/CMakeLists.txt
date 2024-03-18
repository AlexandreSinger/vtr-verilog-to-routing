
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_UTIL_SOURCES)

add_library(libvprutil STATIC
    ${VPR_UTIL_SOURCES}
)

set_target_properties(libvprutil PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprutil PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprutil
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprutil PRIVATE ${VPR_COMPILE_DEFINITIONS})

