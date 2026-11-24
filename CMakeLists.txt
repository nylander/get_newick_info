cmake_minimum_required(VERSION 3.10)
project(get_newick_info VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_executable(get_newick_info src/get_newick_info.c)
target_link_libraries(get_newick_info PRIVATE m)
install(TARGETS get_newick_info DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
    add_test(NAME support_values
        COMMAND get_newick_info --support "${CMAKE_CURRENT_SOURCE_DIR}/data/brlens_bp.tree")
    set_tests_properties(support_values PROPERTIES
        PASS_REGULAR_EXPRESSION "60[\r\n]+23[\r\n]+52[\r\n]+53[\r\n]+100")
endif()
