cmake_minimum_required(VERSION 3.19)

project(RTModel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Build executables
# =================

add_executable(Reader RTModel/lib/Reader.cpp)
target_include_directories(Reader PRIVATE RTModel/include)

add_executable(InitCond RTModel/lib/InitCond.cpp)
target_include_directories(InitCond PRIVATE RTModel/include)

add_executable(LevMar RTModel/lib/bumper.cpp RTModel/lib/LevMar.cpp RTModel/lib/LevMarFit.cpp
        RTModel/lib/VBBinaryLensingLibrary.cpp)
target_include_directories(LevMar PRIVATE RTModel/include)


add_executable(ModelSelector RTModel/lib/bumper.cpp RTModel/lib/ModelSelector.cpp)
target_include_directories(ModelSelector PRIVATE RTModel/include)

add_executable(Finalizer RTModel/lib/bumper.cpp RTModel/lib/Finalizer.cpp)
target_include_directories(Finalizer PRIVATE RTModel/include)

# Install executables
# ===================

set(site_packages_binary_directory "${SKBUILD_PLATLIB_DIR}/RTModel/bin")
message("$SKBUILD_PLATLIB_DIR")

install(TARGETS Reader DESTINATION ${site_packages_binary_directory})
install(TARGETS InitCond DESTINATION ${site_packages_binary_directory})
install(TARGETS LevMar DESTINATION ${site_packages_binary_directory})
install(TARGETS ModelSelector DESTINATION ${site_packages_binary_directory})
install(TARGETS Finalizer DESTINATION ${site_packages_binary_directory})

# Current RTModel pathing results in editable installs searching the source directories.
set(source_binary_directory "${CMAKE_CURRENT_LIST_DIR}/RTModel/bin")

install(TARGETS Reader DESTINATION ${source_binary_directory})
install(TARGETS InitCond DESTINATION ${source_binary_directory})
install(TARGETS LevMar DESTINATION ${source_binary_directory})
install(TARGETS ModelSelector DESTINATION ${source_binary_directory})
install(TARGETS Finalizer DESTINATION ${source_binary_directory})