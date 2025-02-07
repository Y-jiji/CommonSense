

# -- Scan Source and Configure Build

function(AutoBuild)
    cmake_parse_arguments(AUTO "STATIC;SHARED" "LIB_DIR;BIN_DIR" "" "${ARGV}")
    # List Lib Source
    file(GLOB_RECURSE SRC ${AUTO_LIB_DIR}/*.cpp ${AUTO_LIB_DIR}/*.c ${AUTO_LIB_DIR}/*.cc)
    file(GLOB_RECURSE INC ${AUTO_LIB_DIR}/*.hpp ${AUTO_LIB_DIR}/*.h)
    # List Unit Tests
    file(GLOB_RECURSE SRC_TEST ${AUTO_LIB_DIR}/*.test.cpp ${AUTO_LIB_DIR}/*.test.cc)
    file(GLOB_RECURSE INC_TEST ${AUTO_LIB_DIR}/*.test.hpp ${AUTO_LIB_DIR}/*.test.h)
    # Exclude Unit Tests
    list(REMOVE_ITEM SRC ${SRC_TEST})
    list(REMOVE_ITEM INC ${INC_TEST})
    # Create Include Directory
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/include/${CMAKE_PROJECT_NAME})
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/${CMAKE_PROJECT_NAME})
    # Configure Header
    # + The header files are symbolicly linked to the build directory. 
    # + So any changes that you did in the build directory is synchronized to your source file. 
    foreach(F ${INC})
        file(RELATIVE_PATH R ${CMAKE_CURRENT_SOURCE_DIR}/${AUTO_LIB_DIR} ${F})
        get_filename_component(D ${R} DIRECTORY)
        file(CREATE_LINK ${AUTO_LIB_DIR}/${R} ${CMAKE_BINARY_DIR}/include/${CMAKE_PROJECT_NAME}/${R})
        install(FILES ${F} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}/${D})
    endforeach(F R D)
    # Configure Build Target
    if(${AUTO_STATIC})
        message(STATUS "STATIC LIBRARY ${CMAKE_PROJECT_NAME}")
        add_library(${CMAKE_PROJECT_NAME} STATIC ${SRC})
    elseif(${AUTO_SHARED})
        message(STATUS "SHARED LIBRARY ${CMAKE_PROJECT_NAME}")
        add_library(${CMAKE_PROJECT_NAME} SHARED ${SRC})
    else()
        message(FATAL_ERROR "PLEASE ADD 'STATIC', 'INTERFACE' OR 'SHARED' FOR LIB TYPE")
    endif()

    # Configure Build Target Directory
    target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_BINARY_DIR}/include)

    # Configure Installation
    include(GNUInstallDirs)
    install(TARGETS ${CMAKE_PROJECT_NAME})

    # Make gtests available in`ctest` test framework
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    message(STATUS "SOURCE OF ${CMAKE_PROJECT_NAME}")
    foreach(F ${SRC})
        message(STATUS "+ " ${F})
    endforeach(F R)

    message(STATUS "UNIT TESTS")
    include(GoogleTest)
    foreach(F ${SRC_TEST}) # unit tests
        file(RELATIVE_PATH R ${CMAKE_CURRENT_SOURCE_DIR} ${F})
        string(REPLACE "/" "-" R ${R})
        string(REPLACE ".cpp" "" R ${R})
        message(STATUS "+ " ${F})
        add_executable(${R} ${F})
        target_link_libraries(${R} GTest::gtest_main GTest::gtest ${CMAKE_PROJECT_NAME})
        set_target_properties(${R} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
        gtest_add_tests(TARGET "${R}")
    endforeach(F N)

    # Remove source variables
    unset(SRC_TEST)
    unset(INC_TEST)
    unset(INC)
    unset(SRC)

    # Process executables
    file(GLOB_RECURSE SRC ${AUTO_BIN_DIR}/*.cpp)
    message(STATUS "EXECUTABLES")
    foreach(F ${SRC})
        message(STATUS "+ " ${F})
        file(RELATIVE_PATH R ${CMAKE_CURRENT_SOURCE_DIR}/${AUTO_BIN_DIR} ${F})
        string(REPLACE ".cpp" "" R ${R})
        add_executable(${R} ${F})
        target_link_libraries(${R} PRIVATE ${CMAKE_PROJECT_NAME})
    endforeach(F R)
    unset(SRC)

endfunction(AutoBuild)

# -- Manage External Dependencies Using Git

function(Execute)
    cmake_parse_arguments(EX "" "WORKING_DIRECTORY" "COMMAND" "${ARGV}")
    list(JOIN EX_COMMAND " " EX_COMMAND)
    message(STATUS "+ EXECUTE " ${EX_COMMAND})
    execute_process(
        ${ARGN}
        WORKING_DIRECTORY ${EX_WORKING_DIRECTORY}
        COMMAND_ERROR_IS_FATAL ANY
        ERROR_QUIET
        OUTPUT_VARIABLE EXECUTION_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(EXECUTION_OUTPUT "${EXECUTION_OUTPUT}" PARENT_SCOPE)
endfunction()

function(CMakeInstall)
    cmake_parse_arguments(CMI "" "PUBLIC;PRIVATE;USER;REPO;PACK;TARGET;NEED_BUILD" "FLAGS" "${ARGV}")
    if("${CMI_PACK}" STREQUAL "")
        set(CMI_PACK ${CMI_REPO})
    endif()
    if("${CMI_TARGET}" STREQUAL "")
        set(CMI_TARGET ${CMI_PACK})
    else()
        set(CMI_TARGET ${CMI_PACK}::${CMI_TARGET})
    endif()
    if (${CMI_NEED_BUILD})
        Execute(
            COMMAND                 cmake -S . -B build ${CMI_FLAGS} 
                -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/3rd_party_install/${CMI_USER}/${CMI_REPO}
                -DCMAKE_INSTALL_BINDIR=bin
                -DCMAKE_INSTALL_LIBDIR=lib
                -DCMAKE_INSTALL_INCLUDEDIR=include
            WORKING_DIRECTORY       ${CMAKE_BINARY_DIR}/3rd_party_download/${CMI_USER}/${CMI_REPO}
        )
        Execute(
            COMMAND                 cmake --build build --target install -j
            WORKING_DIRECTORY       ${CMAKE_BINARY_DIR}/3rd_party_download/${CMI_USER}/${CMI_REPO}
        )
    endif()
    message(STATUS "+ FIND PACKAGE ${CMI_PACK}")
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/3rd_party_install/${CMI_USER}/${CMI_REPO})
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE INTERNAL "CMAKE_PREFIX_PATH")
    find_package(${CMI_PACK} REQUIRED QUIET)
    if(${CMI_PUBLIC})
        message(STATUS "+ ADD ${CMI_TARGET} AS PUBLIC DEPENDENCY")
        target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_BINARY_DIR}/3rd_party_install/${CMI_USER}/${CMI_REPO}/include)
        target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${CMI_TARGET})
        target_link_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_BINARY_DIR}/3rd_party_install/${CMI_USER}/${CMI_REPO}/lib)
    else()
        message(STATUS "+ ADD ${CMI_TARGET} AS PRVIATE DEPENDENCY")
        target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/3rd_party_install/${CMI_USER}/${CMI_REPO}/include)
        target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${CMI_TARGET})
        target_link_directories(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/3rd_party_install/${CMI_USER}/${CMI_REPO}/lib)
    endif()
endfunction()

function(Git)
    cmake_parse_arguments(GH "PUBLIC;PRIVATE" "SITE;USER;REPO;BRANCH;PIPELINE;PACK;DIR;TARGET" "FLAGS" "${ARGV}")
    message(STATUS "GIT " ${GH_USER}/${GH_REPO} @ ${GH_BRANCH})
    # support for ssh
    if ("${GH_SITE}" MATCHES "^git.*")
        set(URL ${GH_SITE}:${GH_USER}/${GH_REPO})
    else()
        set(URL ${GH_SITE}/${GH_USER}/${GH_REPO})
    endif()
    # if the cloned package needs build
    set(NEED_BUILD OFF)
    # download package from using git clone
    if (NOT EXISTS "${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}/.git")
        Execute(
            COMMAND             git clone
                --recurse-submodules -j8 --quiet --depth 1 
                --branch "${GH_BRANCH}" "${URL}" 
                "${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}"
            WORKING_DIRECTORY  "${CMAKE_BINARY_DIR}"
        )
        set(NEED_BUILD ON)
    # pull if needed
    else()
        Execute(
            COMMAND             git rev-parse HEAD
            WORKING_DIRECTORY  "${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}"
        )
        set(GH_HEAD_BEFORE "${EXECUTION_OUTPUT}")
        message(STATUS "+ OLD HEAD ${GH_HEAD_BEFORE}")
        Execute(
            COMMAND             git pull
                --recurse-submodules -j8 --quiet
            WORKING_DIRECTORY  "${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}"
        )
        Execute(
            COMMAND             git rev-parse HEAD
            WORKING_DIRECTORY  "${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}"
        )
        set(GH_HEAD_AFTER "${EXECUTION_OUTPUT}")
        message(STATUS "+ NEW HEAD ${GH_HEAD_AFTER}")
        if (${GH_HEAD_AFTER} STREQUAL ${GH_HEAD_BEFORE})
            set(NEED_BUILD OFF)
        else()
            set(NEED_BUILD ON)
        endif()
    endif()
    # configure how this package is added to current project
    if(${GH_PIPELINE} STREQUAL "INCLUDE")
        message(STATUS "+ DIRECTLY INCLUDE ${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}")
        include_directories(${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO})
    elseif(${GH_PIPELINE} STREQUAL "CMAKE SUBDIR")
        set(CMAKE_MESSAGE_LOG_LEVEL__ ${CMAKE_MESSAGE_LOG_LEVEL})
        set(CMAKE_MESSAGE_LOG_LEVEL ERROR)
        add_subdirectory(
            ${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}
            ${CMAKE_BINARY_DIR}/3rd_party_subdirs/${GH_USER}/${GH_REPO}
            EXCLUDE_FROM_ALL
        )
        set(CMAKE_MESSAGE_LOG_LEVEL ${CMAKE_MESSAGE_LOG_LEVEL__})
    elseif(${GH_PIPELINE} STREQUAL "CMAKE INSTALL")
        # install headers to given directory
        CMakeInstall(
            REPO    ${GH_REPO}
            USER    ${GH_USER}
            PACK    ${GH_PACK}
            TARGET  ${GH_TARGET}
            FLAGS   ${GH_FLAGS}
            PUBLIC  ${GH_PUBLIC}
            PRIVATE ${GH_PRIVATE}
            NEED_BUILD ${NEED_BUILD}
        )
    elseif(${GH_PIPELINE} STREQUAL "AUTOMAKE INSTALL")
        # run the build command if it is a new thing
        if(${NEED_BUILD})
            Execute(
                COMMAND                 ./autogen.sh --prefix ${CMAKE_BINARY_DIR}/3rd_party_install/${GH_USER}/${GH_REPO}
                WORKING_DIRECTORY       ${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}
            )
            Execute(
                COMMAND                 make install -j4
                WORKING_DIRECTORY       ${CMAKE_BINARY_DIR}/3rd_party_download/${GH_USER}/${GH_REPO}
            )
        endif()
        include_directories(${CMAKE_BINARY_DIR}/3rd_party_install/${GH_USER}/${GH_REPO}/include)
        link_directories(${CMAKE_BINARY_DIR}/3rd_party_install/${GH_USER}/${GH_REPO}/lib)
    else()
        message(FATAL_ERROR "UNKNOWN THIRD PARTY PIPELINE ${GH_PIPELINE}")
    endif()
endfunction()