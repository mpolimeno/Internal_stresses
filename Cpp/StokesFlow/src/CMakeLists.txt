cmake_minimum_required(VERSION 3.22.2)

project(SingleLayer)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR} ../input_files/FaceCenters.txt
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

 file(COPY ${CMAKE_CURRENT_SOURCE_DIR} ../input_files/NormalDirection.txt
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

 add_executable(${PROJECT_NAME} main_driver.cc funcs.cc)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_14)
