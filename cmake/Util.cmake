include(FetchContent)

function(find_or_fetch DEPENDENCY)
    string(TOUPPER ${DEPENDENCY} _DEP_UPPER)
    # -D<DEP>_PATH=<dir> overrides the local source directory.
    if(DEFINED ${_DEP_UPPER}_PATH AND EXISTS "${${_DEP_UPPER}_PATH}")
        set(${DEPENDENCY}_SOURCE_DIR "${${_DEP_UPPER}_PATH}")
        message(STATUS "Using ${DEPENDENCY} source from ${${_DEP_UPPER}_PATH} (cache override)")
    else()
        set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
        # -D<DEP>_REF=<branch|tag|sha> overrides the fetched ref.
        if(DEFINED ${_DEP_UPPER}_REF AND NOT "${${_DEP_UPPER}_REF}" STREQUAL "")
            set(REMOTE_BRANCH "${${_DEP_UPPER}_REF}")
        else()
            get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
        endif()
        message(STATUS "Using ${DEPENDENCY} ref ${REMOTE_BRANCH} from ${REPO_URL}")

        FetchContent_Declare(
            ${DEPENDENCY}
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            GIT_REPOSITORY ${REPO_URL}
            GIT_TAG ${REMOTE_BRANCH}
        )
        FetchContent_Populate(${DEPENDENCY})
    endif()

    set(_DEP_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY})
    # Parent CTest enablement must not register unbuilt dependency tests.
    set(BUILD_TESTING OFF)
    add_subdirectory(${${DEPENDENCY}_SOURCE_DIR}
        ${_DEP_BINARY_DIR}
        EXCLUDE_FROM_ALL
    )
    # libneo sets CMAKE_RUNTIME_OUTPUT_DIRECTORY to its own binary dir, so its
    # executables (e.g. vacfield.x) land directly in _DEP_BINARY_DIR. Export it
    # for runtime helpers and documentation; this replaces the former
    # $CODE/libneo/build location.
    set(${_DEP_UPPER}_DIR ${_DEP_BINARY_DIR} CACHE PATH
        "${DEPENDENCY} build directory (holds built executables)" FORCE)
endfunction()


function(get_branch_or_main REPO_URL REMOTE_BRANCH)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${BRANCH}
        OUTPUT_VARIABLE BRANCH_EXISTS
        ERROR_VARIABLE GIT_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(BRANCH_EXISTS)
        set(${REMOTE_BRANCH} ${BRANCH} PARENT_SCOPE)
    else()
        set(${REMOTE_BRANCH} "main" PARENT_SCOPE)
    endif()
endfunction()
