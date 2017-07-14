cmake_minimum_required(VERSION 3.2.2)
project (jspec)
set (CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 11)
#set (CMAKE_CXX_FLAGS "-Wall")


#For the shared library:
set (PROJECT_LINK_LIBS muparser)
link_directories(lib)

include_directories(include)

file(GLOB SOURCES "src/*.cc")

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_LINK_LIBS})

