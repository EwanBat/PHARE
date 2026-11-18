cmake_minimum_required (VERSION 3.20.1)

project(test-updater)

set(SOURCES test_updater.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  phare_simulator
  ${GTEST_LIBS})

add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


set(target test-updater-physical-seam)
add_executable(${target} test_updater_physical_seam.cpp)
target_include_directories(${target} PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE phare_core phare_simulator ${GTEST_LIBS})
add_no_mpi_phare_test(${target} ${CMAKE_CURRENT_BINARY_DIR})




set(target test-updater-inflow-flux)
add_executable(${target} test_updater_inflow_flux.cpp)
target_include_directories(${target} PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE phare_core phare_simulator ${GTEST_LIBS})
add_no_mpi_phare_test(${target} ${CMAKE_CURRENT_BINARY_DIR})
