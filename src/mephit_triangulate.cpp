#include "mephit_fem.h"

#include <delaunay32/delaunay.hpp>
#include <delaunay32/quantization.hpp>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

[[noreturn]] void fail(const char *message)
{
  std::fprintf(stderr, "FEM_triangulate_external: %s\n", message);
  std::abort();
}

} // namespace

extern "C" void FEM_triangulate_external(const int npt_inner,
                                           const int npt_outer,
                                           const double *bdry_R,
                                           const double *bdry_Z,
                                           const double R_mid,
                                           const double Z_mid,
                                           const char *fname)
{
  if (npt_inner < 3 || npt_outer < 3 || bdry_R == nullptr ||
      bdry_Z == nullptr || fname == nullptr) {
    fail("invalid boundary input");
  }

  const int npoint = npt_inner + npt_outer;
  std::vector<delaunay32::FloatPoint> source(static_cast<std::size_t>(npoint));
  for (int point = 0; point < npoint; ++point) {
    source[static_cast<std::size_t>(point)] = {
      bdry_R[point], bdry_Z[point]
    };
  }

  try {
    delaunay32::QuantizationOptions quantization_options;
    quantization_options.collision_policy =
      delaunay32::QuantizationCollisionPolicy::Reject;
    const delaunay32::QuantizationResult quantized =
      delaunay32::quantize(source, quantization_options);

    delaunay32::PolygonDomain domain;
    domain.outer_ring.reserve(static_cast<std::size_t>(npt_outer));
    for (int point = 0; point < npt_outer; ++point) {
      domain.outer_ring.push_back(
        static_cast<std::uint32_t>(npt_inner + point));
    }
    std::vector<std::uint32_t> hole;
    hole.reserve(static_cast<std::size_t>(npt_inner));
    for (int point = 0; point < npt_inner; ++point) {
      hole.push_back(static_cast<std::uint32_t>(point));
    }
    domain.holes.push_back(std::move(hole));

    delaunay32::Triangulator triangulator;
    triangulator.set_points(quantized.points);
    triangulator.set_polygons({std::move(domain)});
    const delaunay32::TriangulationResult result = triangulator.triangulate();

    // R_mid and Z_mid are the old Triangle hole marker. The explicit inner
    // ring above carries the same topology; retain the arguments in the C
    // ABI and reject a non-finite marker so callers still get input checks.
    if (!std::isfinite(R_mid) || !std::isfinite(Z_mid)) {
      fail("non-finite hole marker");
    }

    std::ofstream mesh_file(fname);
    if (!mesh_file) {
      fail("could not open output mesh");
    }
    mesh_file.precision(17);
    mesh_file << npoint << ' ' << result.triangles.size() << ' '
              << npoint << '\n';
    for (int point = 0; point < npoint; ++point) {
      mesh_file << bdry_R[point] << ' ' << bdry_Z[point] << " 0\n";
    }
    for (const delaunay32::Triangle &triangle : result.triangles) {
      mesh_file << triangle.i0 + 1U << ' ' << triangle.i1 + 1U << ' '
                << triangle.i2 + 1U << " 1\n";
    }
    for (int point = 0; point < npt_inner; ++point) {
      const int next = (point + 1) % npt_inner;
      mesh_file << point + 1 << ' ' << next + 1 << " 2\n";
    }
    for (int point = 0; point < npt_outer; ++point) {
      const int current = npt_inner + point;
      const int next = npt_inner + (point + 1) % npt_outer;
      mesh_file << current + 1 << ' ' << next + 1 << " 2\n";
    }
    if (!mesh_file) {
      fail("could not write output mesh");
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "FEM_triangulate_external: %s\n", error.what());
    std::abort();
  }
}
