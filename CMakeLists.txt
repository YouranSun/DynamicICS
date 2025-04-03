cmake_minimum_required(VERSION 3.10)
project(ICSDynamic)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-fsanitize=undefined")
set(CMAKE_CXX_FLAGS "-fsanitize=address")
set(EXECUTABLE_OUTPUT_PATH "bin/")

add_subdirectory(gadget)
add_subdirectory(core)
add_subdirectory(Index)
add_subdirectory(Dynamic)

add_executable(${PROJECT_NAME} run/main.cc)
target_link_libraries(${PROJECT_NAME} PRIVATE defs gadget heap
 Tools treap treap2 Graph glist Index ICPIndex OrdIndex Dynamic ICPDynamic OrdDynamic)

add_executable(checker run/checker.cc)
target_link_libraries(checker PRIVATE Graph)

add_executable(case_build casestudy/build.cc)
target_link_libraries(case_build PRIVATE defs gadget heap
    Tools treap treap2 Graph glist Index ICPIndex OrdIndex Dynamic ICPDynamic OrdDynamic)

add_executable(case_query casestudy/query.cc)
target_link_libraries(case_query PRIVATE defs gadget heap
    Tools treap treap2 Graph glist Index ICPIndex OrdIndex Dynamic ICPDynamic OrdDynamic)