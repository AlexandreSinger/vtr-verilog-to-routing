
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_DRAW_SOURCES)

add_library(libvprdraw STATIC
    ${VPR_DRAW_SOURCES}
)

set_target_properties(libvprdraw PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprdraw PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprdraw
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprdraw PRIVATE ${VPR_COMPILE_DEFINITIONS})

