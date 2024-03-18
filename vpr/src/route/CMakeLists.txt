
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} VPR_ROUTE_SOURCES)

add_library(libvprroute STATIC
    ${VPR_ROUTE_SOURCES}
)

set_target_properties(libvprroute PROPERTIES PREFIX "") #Avoid extra 'lib' prefix
target_include_directories(libvprroute PUBLIC ${LIB_INCLUDE_DIRS})
target_link_libraries(libvprroute
    ${VPR_DEPENDENCIES}
    libvprbase
)
target_compile_definitions(libvprroute PRIVATE ${VPR_COMPILE_DEFINITIONS})

