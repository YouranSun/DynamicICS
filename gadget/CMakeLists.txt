add_library(gadget STATIC
    gadget.cc
    gadget.h
)

add_library(heap STATIC
    heap.cc
    heap.h
)

add_library(Tools STATIC
    Tools.cc
    Tools.h
)

add_library(treap STATIC
    treap.cc
    treap.h
)

add_library(treap2 STATIC
    treap2.cc
    treap2.h
)

add_library(defs INTERFACE)
add_library(dsu INTERFACE)