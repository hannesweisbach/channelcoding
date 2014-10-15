project("Channel coding exercises")
cmake_minimum_required(VERSION 3.0)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_definitions(-DNDEBUG)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(gf gf.c++)
add_library(bch bch.c++)
add_library(rs rs.c++)

add_executable(tasks tasks.c++)
target_link_libraries(tasks gf bch rs)

add_executable(md md.c++)
add_executable(rm rm.c++)

add_executable(ms_fk minsum_fk.c++)
target_link_libraries(ms_fk bch gf)

add_executable(optimization optimization.c++)
target_link_libraries(optimization bch gf)

add_executable(performance performance.c++)
target_link_libraries(performance bch gf)

add_executable(ln ln.c++)
target_link_libraries(ln bch gf)

