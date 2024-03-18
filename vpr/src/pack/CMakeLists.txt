
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_PACK_SOURCES)

add_library(libvprpack STATIC
    ${VPR_PACK_SOURCES}
)

set_target_properties(libvprpack PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprpack PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprpack
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprpack PRIVATE ${VPR_COMPILE_DEFINITIONS})

