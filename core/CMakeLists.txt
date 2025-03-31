add_library(glist STATIC
    glist.cc
    glist.h
)
target_link_libraries(glist PRIVATE heap)
target_link_libraries(glist PRIVATE treap)

add_library(Graph INTERFACE)
target_include_directories(Graph INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)