cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(hermes-2 LANGUAGES CXX)

# Extra CMake scripts in cmake/ subdirectory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Update submodules
# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
#find_package(Git QUIET)
#if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
#  message(STATUS "Submodule update")
#  execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    RESULT_VARIABLE GIT_SUBMOD_RESULT)
#  if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#  endif()
#endif()

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC HERMES_REVISION)
if(HERMES_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(HERMES_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${HERMES_REVISION}")

# BOUT++ is a dependency
#add_subdirectory(external/BOUT-dev)
find_package(bout++ REQUIRED)

set(HERMES_SOURCES
    src/hermes-2.cxx
    include/hermes-2.hxx
    src/diffusion2d.cxx
    include/diffusion2d.hxx
    src/div_ops.cxx
    include/div_ops.hxx
    src/full-velocity.cxx
    include/full-velocity.hxx
    src/loadmetric.cxx
    include/loadmetric.hxx
    src/mixed.cxx
    include/mixed.hxx
    src/neutral-model.cxx
    include/neutral-model.hxx
    include/none.hxx
    src/radiation.cxx
    include/radiation.hxx
    src/recycling.cxx
    include/recycling.hxx
    atomicpp/ImpuritySpecies.cxx
    atomicpp/ImpuritySpecies.hxx
    atomicpp/json.hxx
    atomicpp/Prad.cxx
    atomicpp/Prad.hxx
    atomicpp/RateCoefficient.cxx
    atomicpp/RateCoefficient.hxx
    atomicpp/sharedFunctions.cxx
    atomicpp/sharedFunctions.hxx
    )

# The main executable target
add_executable(hermes-2
               ${HERMES_SOURCES}
               ${CMAKE_CURRENT_BINARY_DIR}/include/revision.hxx)

target_link_libraries(hermes-2 PRIVATE bout++::bout++)

target_include_directories(hermes-2 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/atomicpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/atomicpp>
  $<INSTALL_INTERFACE:include>
  )

# Build the file containing just the commit hash
# This will be rebuilt on every commit!
configure_file(
  "${PROJECT_SOURCE_DIR}/include/revision.hxx.in"
  "${PROJECT_BINARY_DIR}/include/revision.hxx")

# Once built, copy the data and test directories
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/linear-device $<TARGET_FILE_DIR:${PROJECT_NAME}>/linear-device)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/json_database $<TARGET_FILE_DIR:${PROJECT_NAME}>/json_database)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                   ${CMAKE_SOURCE_DIR}/impurity_user_input.json $<TARGET_FILE_DIR:${PROJECT_NAME}>/impurity_user_input.json)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/tokamak $<TARGET_FILE_DIR:${PROJECT_NAME}>/tokamak)
