cmake_minimum_required(VERSION 3.13)

project(SEALTinyAlabels VERSION 4.1.1 LANGUAGES CXX)

# If not called from root CMakeLists.txt
if(NOT DEFINED SEAL_BUILD_TINYLABELS)
    set(SEAL_BUILD_TINYLABELS ON)

    # Import Microsoft SEAL
    find_package(SEAL 4.1.1 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

if(SEAL_BUILD_TINYLABELS)
    add_executable(setup)
    target_sources(setup
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/setup.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )

    add_executable(enc1)
    target_sources(enc1
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/enc1.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )

    add_executable(enc2)
    target_sources(enc2
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/enc2.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )

    add_executable(keygen)
    target_sources(keygen
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/keygen.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )

    add_executable(dec)
    target_sources(dec
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/dec.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )

    add_executable(gen_samples)
    target_sources(gen_samples
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/gen_samples.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )

    add_executable(benchmark)
    target_sources(benchmark
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/benchmark.cpp
            ${CMAKE_CURRENT_LIST_DIR}/batchselect.cpp
    )


    if(TARGET SEAL::seal)
        target_link_libraries(setup PRIVATE SEAL::seal)
        target_link_libraries(enc1 PRIVATE SEAL::seal)
        target_link_libraries(enc2 PRIVATE SEAL::seal)
        target_link_libraries(keygen PRIVATE SEAL::seal)
        target_link_libraries(dec PRIVATE SEAL::seal)
        target_link_libraries(gen_samples PRIVATE SEAL::seal)
        target_link_libraries(benchmark PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(setup PRIVATE SEAL::seal_shared)
        target_link_libraries(enc1 PRIVATE SEAL::seal_shared)
        target_link_libraries(enc2 PRIVATE SEAL::seal_shared)
        target_link_libraries(keygen PRIVATE SEAL::seal_shared)
        target_link_libraries(dec PRIVATE SEAL::seal_shared)
        target_link_libraries(gen_samples PRIVATE SEAL::seal_shared)
        target_link_libraries(benchmark PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()
endif()
