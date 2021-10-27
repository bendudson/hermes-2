cmake_minimum_required(VERSION 3.9...3.12)

project(hermes-2 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
find_package(bout++ REQUIRED)
add_executable(hermes-2 hermes-2.cxx div_ops.cxx loadmetric.cxx radiation.cxx neutral-model.cxx 		  diffusion2d.cxx recycling.cxx full-velocity.cxx mixed.cxx atomicpp/ImpuritySpecies.cxx atomicpp/Prad.cxx atomicpp/RateCoefficient.cxx atomicpp/sharedFunctions.cxx)
target_link_libraries(hermes-2 PRIVATE bout++::bout++)

