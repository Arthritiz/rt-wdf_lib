set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(modules/armadillo-code)

set(LIB_NAME "RtWdfLib")

set(SOURCES
        ./Libs/rt-wdf/rt-wdf.cpp
        ./Libs/rt-wdf/rt-wdf_nlModels.cpp
        ./Libs/rt-wdf/rt-wdf_nlSolvers.cpp
)

add_library(${LIB_NAME} STATIC ${SOURCES})

target_link_libraries(${LIB_NAME}
    PRIVATE
        armadillo)

