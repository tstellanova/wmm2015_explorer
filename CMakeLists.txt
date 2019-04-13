cmake_minimum_required(VERSION 3.13)
project(wmm2015)

set(CMAKE_C_STANDARD 11)

include_directories(src)

add_executable(wmm_file
        src/EGM9615.h
        src/GeomagnetismHeader.h
        src/GeomagnetismLibrary.c
        src/wmm_file.c)

add_executable(wmm_grid
        src/EGM9615.h
        src/GeomagnetismHeader.h
        src/GeomagnetismLibrary.c
        src/wmm_grid.c)

add_executable(wmm_point
        src/EGM9615.h
        src/GeomagnetismHeader.h
        src/GeomagnetismLibrary.c
        src/wmm_point.c)

add_executable(wmm_dump
        src/EGM9615.h
        src/GeomagnetismHeader.h
        src/GeomagnetismLibrary.c
        src/wmm_dump.c)