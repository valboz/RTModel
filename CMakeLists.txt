cmake_minimum_required(VERSION 3.19)

project(RTModel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(RTModel SHARED
        RTModel/lib/bumper.cpp
        RTModel/lib/Finalizer.cpp
        RTModel/lib/InitCond.cpp
        RTModel/lib/LevMar.cpp
        RTModel/lib/LevMarFit.cpp
        RTModel/lib/ModelSelector.cpp
        RTModel/lib/Reader.cpp
        RTModel/lib/VBBinaryLensingLibrary.cpp
        RTModel/include/LevMarFit.h
        RTModel/include/bumper.h
        RTModel/include/VBBinaryLensingLibrary.h
)
target_include_directories(RTModel PRIVATE RTModel/include)