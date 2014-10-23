project("Channel coding exercises")
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_definitions(-DNDEBUG)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(gf gf.c++)
add_library(bch bch.c++)
target_link_libraries(bch gf)
add_library(rs rs.c++)

add_library(eval eval.c++)
target_link_libraries(eval bch)

add_executable(tasks tasks.c++)
target_link_libraries(tasks gf bch rs)

add_executable(md md.c++)
add_executable(rm rm.c++)

add_executable(ms_fk minsum_fk.c++)
target_link_libraries(ms_fk eval gf)

add_executable(optimization optimization.c++)
target_link_libraries(optimization eval gf)

add_executable(optimize_beta optimize_beta.c++)
target_link_libraries(optimize_beta eval gf)

add_executable(optimize_2d optimize_2d.c++)
target_link_libraries(optimize_2d eval gf)

add_executable(performance performance.c++)
target_link_libraries(performance eval gf)

add_executable(uncoded uncoded.c++)
target_link_libraries(uncoded eval gf)

add_executable(ln ln.c++)
target_link_libraries(ln bch gf)

