// nanobind bindings for the Nested Cages pipeline.
//
// Exposes a single high-level function that mirrors the command-line tool but
// works directly on numpy arrays, reusing the exact same C++ code path
// (nested_cages() from nested_cages_lib.{h,cpp}).

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/optional.h>

#include <Eigen/Core>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "nested_cages_lib.h"

namespace nb = nanobind;

using CageList = std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>;

static CageList nested_cages_py(
    Eigen::MatrixXd V,
    Eigen::MatrixXi F,
    std::vector<int> num_faces,
    int quad_order,
    std::string energy_expansion,
    std::string energy_final,
    std::optional<std::vector<bool>> regular,
    std::optional<CageList> input_decimations)
{
  const int k = static_cast<int>(num_faces.size());
  if (k == 0)
  {
    throw std::invalid_argument("num_faces must contain at least one level");
  }

  // regular flags (default: adaptive decimation everywhere)
  std::vector<bool> regular_vec(k, false);
  if (regular)
  {
    if (static_cast<int>(regular->size()) != k)
    {
      throw std::invalid_argument("len(regular) must equal len(num_faces)");
    }
    regular_vec = *regular;
  }

  // optional per-level input decimations
  std::vector<Eigen::MatrixXd> V_decim(k);
  std::vector<Eigen::MatrixXi> F_decim(k);
  if (input_decimations)
  {
    if (static_cast<int>(input_decimations->size()) != k)
    {
      throw std::invalid_argument(
          "len(input_decimations) must equal len(num_faces)");
    }
    for (int i = 0; i < k; ++i)
    {
      V_decim[i] = (*input_decimations)[i].first;
      F_decim[i] = (*input_decimations)[i].second;
    }
  }

  CageList out_pairs;
  std::vector<Eigen::MatrixXd> Vout;
  std::vector<Eigen::MatrixXi> Fout;

  bool ok;
  {
    // Release the GIL: the pipeline is a long-running, pure-C++ computation.
    nb::gil_scoped_release release;
    ok = nested_cages(
        V, F, quad_order, num_faces, regular_vec, V_decim, F_decim,
        energy_expansion, energy_final, Vout, Fout);
  }

  if (!ok)
  {
    throw std::runtime_error(
        "nested_cages failed (flow diverged, or a cage self-intersected or "
        "intersected the previous layer). See stdout for details.");
  }

  out_pairs.reserve(Vout.size());
  for (std::size_t i = 0; i < Vout.size(); ++i)
  {
    out_pairs.emplace_back(std::move(Vout[i]), std::move(Fout[i]));
  }
  return out_pairs;
}

NB_MODULE(_nested_cages, m)
{
  m.doc() = "Nested Cages: nested, non-self-intersecting cages for a triangle mesh";

  m.def("nested_cages", &nested_cages_py,
        nb::arg("V"),
        nb::arg("F"),
        nb::arg("num_faces"),
        nb::arg("quad_order") = 1,
        nb::arg("energy_expansion") = "None",
        nb::arg("energy_final") = "None",
        nb::arg("regular") = nb::none(),
        nb::arg("input_decimations") = nb::none(),
        R"(Compute a sequence of nested cages for a triangle mesh.

Parameters
----------
V : (#V, 3) float64 array
    Input (fine) mesh vertex positions.
F : (#F, 3) int32 array
    Input mesh triangle indices into ``V``.
num_faces : list[int]
    Target number of faces for each cage, one entry per level (decreasing).
quad_order : int, optional
    Quadrature order for the shrinking flow (1, 2 or 3). Default 1.
energy_expansion : str, optional
    Energy minimized during re-inflation. One of ``"None"``, ``"DispStep"``,
    ``"DispInitial"``, ``"Volume"``, ``"SurfARAP"``, ``"VolARAP"``. Default "None".
energy_final : str, optional
    Energy minimized after re-inflation. Same choices. Default "None".
regular : list[bool], optional
    Per-level flag: ``True`` for regular decimation, ``False`` (default) for
    adaptive. If given, must have the same length as ``num_faces``.
input_decimations : list[tuple[ndarray, ndarray]], optional
    Per-level explicit input decimations ``(V_i, F_i)`` used instead of
    decimating the previous cage. If given, must have the same length as
    ``num_faces``.

Returns
-------
list[tuple[(#Vc, 3) ndarray, (#Fc, 3) ndarray]]
    One ``(V, F)`` cage per level.

Raises
------
RuntimeError
    If the pipeline fails for any level.
)");
}
