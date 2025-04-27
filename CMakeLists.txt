cmake_minimum_required(VERSION 3.0)

project(a4)
set(CMAKE_CXX_STANDARD 11)


find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)
find_package(nlohmann_json 3.11.2 REQUIRED)  

add_library(col781 src/hw.cpp src/camera.cpp deps/src/gl.c src/skeletal.cpp  src/spring.cpp src/obstacle.cpp)
target_include_directories(col781 PUBLIC deps/include)
target_link_libraries(col781 glm::glm OpenGL::GL SDL2::SDL2 nlohmann_json::nlohmann_json)
# target_link_libraries(col781 glm::glm OpenGL::GL SDL2::SDL2)

add_executable(e2 executables/e2.cpp)
target_include_directories(e2 PUBLIC src)
target_link_libraries(e2 col781)

add_executable(example executables/example.cpp)
target_include_directories(example PUBLIC src)
target_link_libraries(example col781)

add_executable(e1 executables/e1.cpp)
target_include_directories(e1 PUBLIC src)
target_link_libraries(e1 col781)


add_executable(e3 executables/e3.cpp)
target_include_directories(e3 PUBLIC src)
target_link_libraries(e3 col781)

add_executable(e4 executables/e4.cpp)
target_include_directories(e4 PUBLIC src)
target_link_libraries(e4 col781)


# AddressSanitizer Flags
# set(SANITIZE_FLAGS "-fsanitize=address")

# target_compile_options(col781 PRIVATE ${SANITIZE_FLAGS})
# target_link_options(col781 PRIVATE ${SANITIZE_FLAGS})

# target_compile_options(e2 PRIVATE ${SANITIZE_FLAGS})
# target_link_options(e2 PRIVATE ${SANITIZE_FLAGS})

# target_compile_options(e1 PRIVATE ${SANITIZE_FLAGS})
# target_link_options(e1 PRIVATE ${SANITIZE_FLAGS})