# Python bindings for nested_cages (built when NESTED_CAGES_BUILD_PYTHON=ON,
# e.g. via `pip install .` which drives this through scikit-build-core).

find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)

# nanobind (fetched, matching the version used by the libigl Python bindings)
FetchContent_Declare(nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG v2.7.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nanobind)

nanobind_add_module(_nested_cages STABLE_ABI
  "${CMAKE_CURRENT_SOURCE_DIR}/nested_cages_bind.cpp")
target_link_libraries(_nested_cages PRIVATE nested_cages_core)

# Generate a .pyi stub next to the compiled module.
nanobind_add_stub(_nested_cages_stub
  MODULE _nested_cages
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_nested_cages.pyi"
  PYTHON_PATH $<TARGET_FILE_DIR:_nested_cages>
  DEPENDS _nested_cages)

# Install into the `nested_cages` package (scikit-build-core wheel layout).
install(TARGETS _nested_cages LIBRARY DESTINATION nested_cages)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/_nested_cages.pyi" DESTINATION nested_cages)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/nested_cages/__init__.py" DESTINATION nested_cages)
