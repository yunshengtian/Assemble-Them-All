# CMake entry point.
cmake_minimum_required(VERSION 3.1)

# Do not want to generate zeor_check because they cause messy folders in
# visual studio 2017.
set(CMAKE_SUPPRESS_REGENERATION true)

set(CMAKE_CXX_STANDARD 11)

set(EXTERNAL_HEADER
  "${CMAKE_CURRENT_LIST_DIR}/eigen/"
  "${CMAKE_CURRENT_LIST_DIR}/glew/include"
  "${CMAKE_CURRENT_LIST_DIR}/glfw/include/GLFW"
  "${CMAKE_CURRENT_LIST_DIR}/imgui"
  "${CMAKE_CURRENT_LIST_DIR}/libigl/include"
  "${CMAKE_CURRENT_LIST_DIR}/stb"
  "${CMAKE_CURRENT_LIST_DIR}/pugixml/src"
  "${CMAKE_CURRENT_LIST_DIR}/tiny_obj_loader"
  "${CMAKE_CURRENT_LIST_DIR}/sdfgen"
)

include_directories(${EXTERNAL_HEADER})

# libigl
add_subdirectory(libigl)

# tiny obj loader
add_subdirectory(tiny_obj_loader)

# pugixml
add_subdirectory(pugixml)

# pybind11
add_subdirectory(pybind11)

# sdfgen
add_subdirectory(sdfgen)

# viewer related
if(MSVC AND NOT "${MSVC_VERSION}" LESS 1400)
  add_definitions( "/MP" )
endif()

add_definitions(
  -DTW_STATIC
  -DTW_NO_LIB_PRAGMA
  -DTW_NO_DIRECT3D
  -DGLEW_STATIC
  -D_CRT_SECURE_NO_WARNINGS
)

### GLFW ###
add_subdirectory(glfw)

include_directories(
  glfw/include/GLFW/
  glew/include/
)

### GLEW ###
set(GLEW_SOURCE
  glew/src/glew.c
)
set(GLEW_HEADERS
)
add_library(glew STATIC
  ${GLEW_SOURCE}
  ${GLEW_INCLUDE}
)
target_link_libraries(glew
  ${OPENGL_LIBRARY}
)

### IMGUI ###
set(IMGUI_HEADER
  imgui/imgui.h
  imgui/imgui_internal.h
  imgui/imconfig.h
  imgui/stb_rect_pack.h
  imgui/stb_textedit.h
  imgui/stb_truetype.h
)
set(IMGUI_SOURCE
  imgui/imgui.cpp
  imgui/imgui_demo.cpp
  imgui/imgui_draw.cpp
)
add_library(imgui STATIC
  ${IMGUI_HEADER}
  ${IMGUI_SOURCE}
)
# Treat windows and ubuntu differently.
if(MSVC)
  target_link_libraries(imgui
    glew
    glfw
  )
endif(MSVC)
if(UNIX AND NOT APPLE)
  target_link_libraries(imgui
    glew
    glfw
    X11
    Xi
    Xrandr
    Xxf86vm
    Xinerama
    Xcursor
    rt
    m
    pthread
    dl
  )
endif(UNIX AND NOT APPLE)
if(APPLE)
  include_directories(/System/Library/Frameworks)
  find_library(COCOA_LIBRARY Cocoa)
  find_library(IOKIT_LIBRARY IOKit)
  find_library(COREVIDEO_LIBRARY CoreVideo)
  target_link_libraries(imgui
    glew
    glfw
    ${COCOA_LIBRARY}
    ${IOKIT_LIBRARY}
    ${COREVIDEO_LIBRARY}
  )
endif(APPLE)

