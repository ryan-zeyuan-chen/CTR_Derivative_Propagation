cmake_minimum_required(VERSION 3.23.2)
project(ctrOpt VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -DNDEBUG") # -fcatch-undefined-behavior

# Generated binary files will be directed to bin folder 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# addind the library and main executable directories to the project tree
add_subdirectory(ctr_library) # static CTR library
add_subdirectory(executable)  # main executable

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    ON)

find_package(Boost REQUIRED)
find_package(LAPACK REQUIRED)

set(BLA_VENDOR OpenBLAS)
set(BLA_STATIC                   ON)  # use static linkage
find_package(BLAS REQUIRED)

find_package(Blaze REQUIRED)

# Check for the Boost libray, if found print message, include dirs and link libraries.
if (Boost_FOUND)
    message(STATUS "Boost has been found!")
    include_directories(${Boost_INCLUDE_DIRS})
    link_libraries(${Boost_LIBRARIES})

    if (WIN32)
        # disable autolinking in boost
        add_definitions( -DBOOST_ALL_NO_LIB )  
        # asks Boost to report information about automatic linking requests on Visual Studio
        add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
        # force all boost libraries to dynamic link
        add_definitions( -DBOOST_ALL_DYN_LINK )
    endif()
elseif(NOT Boost_FOUND)
    error(STATUS "Boost has not been found!")
endif()

# Check for the LAPACK libray, if found print message, include dirs and link libraries.
if (LAPACK_FOUND)
    message("LAPACK has been found!")
    include_directories(${LAPACK_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES})
elseif(NOT LAPACK_FOUND)
    error("LAPACK API has not been found")
endif()

if (BLAS_FOUND)
    message(STATUS "BLAS has been found!")
    # add_definitions(-DBUILD_WITHOUT_LAPACK=OFF -DBUILD_WITHOUT_CBLAS=OFF -DNOFORTRAN=OFF -DDYNAMIC_ARCH=ON -DDYNAMIC_OLDER=ON)
    include_directories(${BLAS_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${BLAS_LIBRARIES})
elseif(NOT BLAS_FOUND)
    error(STATUS "BLAS library has not been found")
endif()

# Check for the Blaze libray, if found print message, include dirs and link libraries.
find_package( blaze )
if( blaze_FOUND )
   message(STATUS "Blaze has been found!")
   target_link_libraries(${PROJECT_NAME} INTERFACE blaze::blaze)
elseif(NOT Blaze_FOUND)
    error(STATUS "Blaze has not been found!")
endif()

# # Check for the OpenIGTLink libray. If found, print message, include dirs, and link libraries.
# set(OpenIGTLink_DIR "/usr/local/lib/igtl/cmake/igtl-3.1/")
# find_package(OpenIGTLink REQUIRED PATHS OpenIGTLink_DIR)
# if( OpenIGTLink_FOUND )
#    message(STATUS "OpenIGTLink has been found!")
#    target_include_directories(${PROJECT_NAME} PUBLIC ${OpenIGTLink_USE_FILE})
#    target_link_libraries(${PROJECT_NAME} PUBLIC OpenIGTLink)
# elseif(NOT OpenIGTLink_FOUND)
#     error(STATUS "OpenIGTLink has not been found!")
# endif()

# Check for the TBB libray. If found, print message, include dirs, and link libraries.
find_package(TBB REQUIRED)
if( TBB_FOUND )
   message(STATUS "TBB has been found!")
   target_link_libraries(${PROJECT_NAME} PRIVATE TBB::tbb)
elseif(NOT TBB_FOUND)
    error(STATUS "TBB has not been found!")
endif()

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)