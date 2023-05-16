cmake_minimum_required(VERSION 3.22.2)

project(DLA)

add_executable(${PROJECT_NAME} main_driver.cc funcs.cc)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_14)
