cmake_minimum_required(VERSION 3.19)

project(RTModel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(Reader.exe RTModel/lib/Reader.cpp)
target_include_directories(Reader.exe PRIVATE RTModel/include)

add_executable(InitCond.exe RTModel/lib/InitCond.cpp)
target_include_directories(InitCond.exe PRIVATE RTModel/include)

add_executable(LevMar.exe RTModel/lib/bumper.cpp RTModel/lib/LevMar.cpp RTModel/lib/LevMarFit.cpp
        RTModel/lib/VBBinaryLensingLibrary.cpp)
target_include_directories(LevMar.exe PRIVATE RTModel/include)


add_executable(ModelSelector.exe RTModel/lib/bumper.cpp RTModel/lib/ModelSelector.cpp)
target_include_directories(ModelSelector.exe PRIVATE RTModel/include)

add_executable(Finalizer.exe RTModel/lib/bumper.cpp RTModel/lib/Finalizer.cpp)
target_include_directories(Finalizer.exe PRIVATE RTModel/include)


set(binary_directory "${CMAKE_CURRENT_LIST_DIR}/RTModel/bin")

install(TARGETS Reader.exe DESTINATION ${binary_directory})
install(TARGETS InitCond.exe DESTINATION ${binary_directory})
install(TARGETS LevMar.exe DESTINATION ${binary_directory})
install(TARGETS ModelSelector.exe DESTINATION ${binary_directory})
install(TARGETS Finalizer.exe DESTINATION ${binary_directory})

