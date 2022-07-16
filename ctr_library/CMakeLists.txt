# listing the CTR library source files into a variable
set(CTR_LIB_SOURCES
    src/CTR.cpp
    src/Observer.cpp
    src/ODESystem.cpp
    src/Segment.cpp
    src/Tube.cpp
)

# creating a CTR static library
add_library(CTRlib STATIC ${CTR_LIB_SOURCES})
target_include_directories(CTRlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# linking the library dependency on the NLOPT optimization and math libraries
# target_link_libraries(CTRlib PRIVATE nlopt m)