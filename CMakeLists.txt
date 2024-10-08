cmake_minimum_required(VERSION 3.16)
project(bisHash)

set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
#set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}  -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")

add_executable(bisHash
        main.cpp
        configs.h
        smith.h smith.cpp
        params.h params.cpp

        hashing/hash_utils.h hashing/hash_utils.cpp
        hashing/family_min_hash.cpp hashing/family_min_hash.h

        indexing/object_writer.h indexing/object_writer.cpp

        bio_utils/sequence_reader.h bio_utils/sequence_reader.cpp
        bio_utils/sequence.h bio_utils/sequence.cpp
        bio_utils/sam_writer.h bio_utils/sam_writer.cpp

        utils/logger.h utils/logger.cpp
        utils/time_profile.h utils/time_profile.cpp
        utils/multiproc.h utils/multiproc.cpp
)

target_link_libraries(bisHash -lz -lm -lpthread)
