include_directories()

add_library(Index STATIC
    Index.cc
    Index.h
)

add_library(ICPIndex STATIC
    ICPIndex.cc
    ICPIndex.h
)

add_library(OrdIndex STATIC
    OrdIndex.cc
    OrdIndex.h
)

target_link_libraries(ICPIndex PRIVATE Index)
target_link_libraries(OrdIndex PRIVATE Index)
