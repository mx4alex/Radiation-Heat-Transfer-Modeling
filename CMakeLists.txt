cmake_minimum_required(VERSION 3.27)
project(ThermalSimulation VERSION 1.0)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "11.0")
        message(FATAL_ERROR "GCC 11 or newer is required for C++20 modules")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.30")
        message(FATAL_ERROR "MSVC 19.30 or newer is required for C++20 modules")
    endif()
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPERIMENTAL_CXX_MODULE_CMAKE_API "2182bf5c-ef0d-489a-91da-49dbc3090d2a")
set(CMAKE_EXPERIMENTAL_CXX_MODULE_DYNDEP 1)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenMP REQUIRED)

set(IMGUI_DIR "C:/imgui" CACHE PATH "Path to ImGui")
set(IMGUI_BACKENDS_DIR "${IMGUI_DIR}/backends" CACHE PATH "Path to ImGui backends")

if(NOT EXISTS ${IMGUI_DIR})
    message(FATAL_ERROR "ImGui directory not found at ${IMGUI_DIR}")
endif()

set(MODULE_IFACES
    modules/surface_properties.ixx
    modules/geometry_loader.ixx
    modules/spectral_solver.ixx
    modules/radiation_solver.ixx
    modules/conduction_solver.ixx
    modules/visualizer.ixx
    modules/menu.ixx
)

set(PROJECT_SRCS
    src/geometry_loader.cpp
    src/surface_properties.cpp
    src/spectral_solver.cpp
    src/radiation_solver.cpp
    src/conduction_solver.cpp
    src/visualizer.cpp
    src/menu.cpp
)

foreach(file ${MODULE_IFACES} ${PROJECT_SRCS})
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
        message(FATAL_ERROR "Source file not found: ${file}")
    endif()
endforeach()

set(IMGUI_SRCS
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_BACKENDS_DIR}/imgui_impl_glut.cpp
    ${IMGUI_BACKENDS_DIR}/imgui_impl_opengl2.cpp
)

foreach(file ${IMGUI_SRCS})
    if(NOT EXISTS "${file}")
        message(FATAL_ERROR "ImGui source file not found: ${file}")
    endif()
endforeach()

add_executable(${PROJECT_NAME} 
    ${MODULE_IFACES}
    ${PROJECT_SRCS}
    ${IMGUI_SRCS}
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${IMGUI_DIR}
    ${IMGUI_BACKENDS_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    OpenMP::OpenMP_CXX
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE 
        /experimental:module
        /W4
    )
else()
    target_compile_options(${PROJECT_NAME} PRIVATE 
        -fmodules-ts
        -Wall 
        -Wextra
    )
endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/resources")
    message(FATAL_ERROR "Resources directory not found")
endif()

file(COPY ${CMAKE_SOURCE_DIR}/resources DESTINATION ${CMAKE_BINARY_DIR})

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)

install(DIRECTORY resources/
    DESTINATION bin/resources
)

enable_testing() 