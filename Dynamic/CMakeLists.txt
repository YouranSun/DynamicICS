add_library(Dynamic STATIC
    Dynamic.cc
    Dynamic.h
)

add_library(ICPDynamic STATIC
    ICPDynamic.cc
    ICPDynamic.h
)

add_library(OrdDynamic STATIC
    OrdDynamic.cc
    OrdDynamic.h
)

target_link_libraries(Dynamic PRIVATE glist)
target_link_libraries(Dynamic PRIVATE Index)
target_link_libraries(ICPDynamic PRIVATE Dynamic)
target_link_libraries(OrdDynamic PRIVATE Dynamic)
target_link_libraries(OrdDynamic PRIVATE treap2)
target_link_libraries(OrdDynamic PRIVATE OrdIndex)