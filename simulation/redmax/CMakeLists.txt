set(SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE HEADER ${SOURCE_DIR}/*.h)
file(GLOB_RECURSE SOURCE ${SOURCE_DIR}/*.cpp)
list(REMOVE_ITEM SOURCE ${SOURCE_DIR}/python_interface.cpp)

add_executable(redmax ${HEADER} ${SOURCE})

pybind11_add_module(redmax_py python_interface.cpp ${HEADER} ${SOURCE})

# Includes.
target_include_directories(redmax PRIVATE ${EXTERNAL_HEADER})
target_include_directories(redmax PRIVATE ${SOURCE_DIR}/)
# Add libs.
target_link_libraries(redmax PRIVATE opengl_viewer)
target_link_libraries(redmax PRIVATE pugixml)
target_link_libraries(redmax PRIVATE tiny_obj_loader)
target_link_libraries(redmax PRIVATE sdfgen)

target_include_directories(redmax_py PRIVATE ${EXTERNAL_HEADER})
target_include_directories(redmax_py PRIVATE ${SOURCE_DIR}/)

# Add libs.
target_link_libraries(redmax_py PRIVATE opengl_viewer)
target_link_libraries(redmax_py PRIVATE pugixml)
target_link_libraries(redmax_py PRIVATE tiny_obj_loader)
target_link_libraries(redmax_py PRIVATE sdfgen)

# Add libs.
target_link_libraries(redmax PRIVATE imgui)
target_link_libraries(redmax_py PRIVATE imgui)

set(CMAKE_CXX_STANDARD_LIBRARIES -ldl)
set(CMAKE_VERBOSE_MAKEFILE ON)

