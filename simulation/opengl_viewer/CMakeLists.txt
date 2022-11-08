set(OPENGL_VIEWER_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE OPENGL_VIEWER_HEADER ${OPENGL_VIEWER_DIR}/*.h)
file(GLOB_RECURSE OPENGL_VIEWER_SOURCE ${OPENGL_VIEWER_DIR}/*.cpp ${OPENGL_VIEWER_DIR}/*.glsl)

add_library(opengl_viewer ${OPENGL_VIEWER_HEADER} ${OPENGL_VIEWER_SOURCE})

# Add include directories.
target_include_directories(opengl_viewer PRIVATE ${OPENGL_VIEWER_DIR}/include)

# Add libs.
target_link_libraries(opengl_viewer
  imgui
)
