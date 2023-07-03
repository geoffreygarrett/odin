#  # Enable testing functionality
#  enable_testing()
#  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Enable testing of the benchmark library.")
#  include(ExternalProject)
#
#  # Set the Google Test properties
#  set(GTEST_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/googletest")
#  set(GTEST_LIB_DIR "${CMAKE_BINARY_DIR}/lib")
#  set(GTEST_LIBRARIES gtest gtest_main)
#  set(GTEST_INSTALL_DIR "${GTEST_PREFIX}")
#
#  # Include Google Test as an external project
#  ExternalProject_Add(
#          external_gtest
#          GIT_REPOSITORY https://github.com/google/googletest.git
#          GIT_TAG v1.13.0
#          PREFIX ${GTEST_PREFIX}
#          CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${GTEST_LIB_DIR} -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
#          -DCMAKE_INSTALL_PREFIX:PATH=${GTEST_INSTALL_DIR}
#          BUILD_COMMAND ${CMAKE_COMMAND}
#          INSTALL_COMMAND "${CMAKE_MAKE_PROGRAM}" install
#          DEPENDEES run_tests
#          LOG_DOWNLOAD ON
#          LOG_CONFIGURE ON
#          LOG_BUILD ON
#  )
#
#  # Define the include directory and library path for Google Test
#  ExternalProject_Get_Property(external_gtest source_dir binary_dir)
#  set(GTEST_INCLUDE_DIR ${source_dir}/googletest/include)
#  log_message(1 "GTest" "Include directory:" ${_LOG_W} ${GTEST_INCLUDE_DIR})
#
#  foreach (_LIB IN ITEMS ${GTEST_LIBRARIES})
#      set(GTEST_LIBRARY_PATH "${binary_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${_LIB}.a")
#      log_message(1 "GTest" "Library ${_LIB} path:" ${_LOG_W} ${GTEST_LIBRARY_PATH})
#
#      # Create an imported library target for Google Test
#      add_library(${_LIB} UNKNOWN IMPORTED)
#      set_property(TARGET ${_LIB} PROPERTY IMPORTED_LOCATION ${GTEST_LIBRARY_PATH})
#      add_dependencies(${_LIB} external_gtest)
#  endforeach ()
#
#  # Define the runTests executable and the source files associated with it
#  #    add_executable(run_tests tests/tests_state_handler.cpp)
#  file(GLOB TEST_SOURCES "./*.cpp")
#  add_executable(run_tests ${TEST_SOURCES})
#
#  # Add the dependencies to the test executable
#  add_dependencies(run_tests gtest)
#
#  # Link runTests executable with gtest, gtest_main, and external_includes libraries
#  target_link_libraries(run_tests PRIVATE ${GTEST_LIBRARIES} external_includes GSL::gsl)
#
#  # Include the Google Test header files
#  target_include_directories(run_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${GTEST_INCLUDE_DIR} ${AUTODIFF_INCLUDE_DIR} ${CEREAL_INCLUDE_DIR})
#
#  # Use Google Test's module to automatically discover tests
#  include(GoogleTest)
#  gtest_discover_tests(run_tests)
#
#  # Set the GTEST_COLOR environment variable for the runTests executable
#  target_compile_definitions(run_tests PRIVATE GTEST_COLOR=1)



