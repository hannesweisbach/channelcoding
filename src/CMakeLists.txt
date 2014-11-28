project("Channel coding exercises")
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#add_definitions(-DNDEBUG)

find_package(Threads REQUIRED)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  find_library(LIBCXXABI NAMES c++abi PATHS ${LIBCXXPATH})
  if(${LIBCXXABI} MATCHES "NOTFOUND")
    message(FATAL_ERROR "libc++abi not found. Set LIBCXXPATH for non-standard paths.")
  endif()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lc++abi")
  if(LIBCXXPATH)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBCXXPATH}")
  endif()
  add_definitions(-D_GNU_SOURCE)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(codes codes/codes.c++)

add_executable(exercises exercises.c++)
target_link_libraries(exercises codes)
add_executable(md md.c++)
add_executable(rm rm.c++)

add_library(simulation simulation/simulation.c++)
target_link_libraries(simulation codes ${CMAKE_THREAD_LIBS_INIT})

add_executable(benchmark simulation/benchmark.c++)
target_link_libraries(benchmark simulation)

add_executable(bitflips simulation/bitflips.c++)
target_link_libraries(bitflips simulation)

add_executable(uncoded simulation/uncoded.c++)
target_link_libraries(uncoded simulation)

