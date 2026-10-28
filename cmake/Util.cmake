include(FetchContent)
find_package(Git REQUIRED)

# Populate DEST with REPO_URL at REF for refs a plain clone cannot reach.
#
# FetchContent clones the whole repository and then runs `git checkout <ref>`.
# A clone only carries objects reachable from a branch or tag, so a bare commit
# SHA whose branch has been deleted -- exactly what the downstream libneo gate
# passes via LIBNEO_REF for a merged or fork pull request -- is absent from the
# clone and the checkout dies with "unable to read tree (<sha>)". Fetching the
# commit by name instead works, because GitHub serves fetch-by-SHA requests.
function(fetch_git_ref DEPENDENCY REPO_URL REF DEST)
    if(EXISTS "${DEST}/.git")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${DEST}
            OUTPUT_VARIABLE _HEAD_SHA
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _RESULT
            ERROR_QUIET
        )
        if(_RESULT EQUAL 0 AND _HEAD_SHA STREQUAL "${REF}")
            message(STATUS "Reusing ${DEPENDENCY} checkout at ${REF} in ${DEST}")
            return()
        endif()
    endif()

    file(REMOVE_RECURSE ${DEST})
    file(MAKE_DIRECTORY ${DEST})
    _run_git_step(${DEST} "initialise ${DEPENDENCY} repository" init --quiet)
    _run_git_step(${DEST} "add ${DEPENDENCY} remote" remote add origin ${REPO_URL})
    _run_git_step(${DEST} "fetch ${DEPENDENCY} ref ${REF} from ${REPO_URL}"
        fetch --quiet --depth 1 origin ${REF})
    _run_git_step(${DEST} "check out ${DEPENDENCY} ref ${REF}"
        checkout --quiet --detach FETCH_HEAD)
endfunction()


function(_run_git_step WORK_DIR WHAT)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} ${ARGN}
        WORKING_DIRECTORY ${WORK_DIR}
        RESULT_VARIABLE _RESULT
        ERROR_VARIABLE _STDERR
    )
    if(NOT _RESULT EQUAL 0)
        message(FATAL_ERROR
            "Failed to ${WHAT}.\n"
            "git ${ARGN} exited with ${_RESULT}:\n${_STDERR}\n"
            "The ref must exist in the remote repository; a commit that was "
            "never pushed, or a pull request head from a private fork, cannot "
            "be fetched."
        )
    endif()
endfunction()


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

        execute_process(
            COMMAND ${GIT_EXECUTABLE} ls-remote --heads --tags ${REPO_URL} ${REMOTE_BRANCH}
            OUTPUT_VARIABLE REF_ON_REMOTE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(REF_ON_REMOTE)
            FetchContent_Declare(
                ${DEPENDENCY}
                DOWNLOAD_EXTRACT_TIMESTAMP TRUE
                GIT_REPOSITORY ${REPO_URL}
                GIT_TAG ${REMOTE_BRANCH}
            )
            FetchContent_Populate(${DEPENDENCY})
        else()
            # No branch or tag points at this ref, so a clone would not contain
            # it. Fetch the commit directly instead. See fetch_git_ref above.
            message(STATUS
                "${DEPENDENCY} ref ${REMOTE_BRANCH} is not a branch or tag; "
                "fetching the commit directly")
            set(${DEPENDENCY}_SOURCE_DIR ${FETCHCONTENT_BASE_DIR}/${DEPENDENCY}-src)
            fetch_git_ref(${DEPENDENCY} ${REPO_URL} ${REMOTE_BRANCH}
                ${${DEPENDENCY}_SOURCE_DIR})
        endif()
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
