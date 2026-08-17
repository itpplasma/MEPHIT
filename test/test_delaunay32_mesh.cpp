#include "mephit_fem.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Mesh {
  int points = 0;
  int triangles = 0;
  int segments = 0;
  std::vector<std::pair<double, double>> coordinates;
  std::vector<std::array<int, 3>> elements;
  std::vector<std::pair<int, int>> boundary;
};

double area(const std::pair<double, double> &a,
            const std::pair<double, double> &b,
            const std::pair<double, double> &c)
{
  return 0.5 * ((b.first - a.first) * (c.second - a.second) -
                (b.second - a.second) * (c.first - a.first));
}

double polygon_area(const std::vector<double> &x,
                    const std::vector<double> &y,
                    int offset, int count)
{
  double result = 0.0;
  for (int point = 0; point < count; ++point) {
    const int current = offset + point;
    const int next = offset + (point + 1) % count;
    result += x[static_cast<std::size_t>(current)] *
              y[static_cast<std::size_t>(next)] -
              y[static_cast<std::size_t>(current)] *
              x[static_cast<std::size_t>(next)];
  }
  return std::abs(result) * 0.5;
}

bool point_in_polygon(const std::vector<std::pair<double, double>> &points,
                      int offset, int count, double x, double y)
{
  bool inside = false;
  for (int point = 0, previous = count - 1; point < count;
       previous = point++) {
    const auto &current_point =
      points[static_cast<std::size_t>(offset + point)];
    const auto &previous_point =
      points[static_cast<std::size_t>(offset + previous)];
    const bool crosses =
      ((current_point.second > y) != (previous_point.second > y)) &&
      (x < (previous_point.first - current_point.first) *
                 (y - current_point.second) /
                 (previous_point.second - current_point.second) +
               current_point.first);
    if (crosses) inside = !inside;
  }
  return inside;
}

double angle(double opposite_squared, double first_squared,
             double second_squared)
{
  const double cosine =
    (first_squared + second_squared - opposite_squared) /
    (2.0 * std::sqrt(first_squared * second_squared));
  return std::acos(std::max(-1.0, std::min(1.0, cosine))) *
         180.0 / std::acos(-1.0);
}

bool read_mesh(const std::string &path, Mesh &mesh)
{
  std::ifstream input(path);
  if (!(input >> mesh.points >> mesh.triangles >> mesh.segments)) {
    return false;
  }
  mesh.coordinates.resize(static_cast<std::size_t>(mesh.points));
  for (auto &coordinate : mesh.coordinates) {
    int marker;
    if (!(input >> coordinate.first >> coordinate.second >> marker) ||
        marker != 0) {
      return false;
    }
  }
  mesh.elements.resize(static_cast<std::size_t>(mesh.triangles));
  for (auto &element : mesh.elements) {
    int region;
    if (!(input >> element[0] >> element[1] >> element[2] >> region) ||
        region != 1) {
      return false;
    }
  }
  mesh.boundary.resize(static_cast<std::size_t>(mesh.segments));
  for (auto &segment : mesh.boundary) {
    int marker;
    if (!(input >> segment.first >> segment.second >> marker) || marker != 2) {
      return false;
    }
  }
  return true;
}

int canonical_edge(int first, int second)
{
  if (first > second) std::swap(first, second);
  return first * 10000 + second;
}

} // namespace

int main(int argc, char **argv)
{
  if (argc != 2) {
    std::cerr << "usage: test_delaunay32_mesh output.msh\n";
    return 2;
  }

  constexpr int inner_points = 128;
  constexpr int outer_points = 16;
  constexpr double center_R = 170.0;
  constexpr double center_Z = 0.0;
  constexpr double inner_R = 50.0;
  constexpr double inner_Z = 50.0;
  constexpr double outer_R = 100.0;
  constexpr double outer_Z = 100.0;
  std::vector<double> boundary_R(inner_points + outer_points);
  std::vector<double> boundary_Z(inner_points + outer_points);
  for (int point = 0; point < inner_points; ++point) {
    const double theta = 2.0 * std::acos(-1.0) * point / inner_points;
    boundary_R[point] = center_R + inner_R * std::cos(theta);
    boundary_Z[point] = center_Z + inner_Z * std::sin(theta);
  }
  for (int point = 0; point < outer_points; ++point) {
    const double theta = 2.0 * std::acos(-1.0) * point / outer_points;
    boundary_R[inner_points + point] = center_R + outer_R * std::cos(theta);
    boundary_Z[inner_points + point] = center_Z + outer_Z * std::sin(theta);
  }

  FEM_triangulate_external(
    inner_points, outer_points, boundary_R.data(), boundary_Z.data(),
    center_R, center_Z, argv[1]);

  Mesh mesh;
  if (!read_mesh(argv[1], mesh)) {
    std::cerr << "could not parse generated FreeFem mesh\n";
    return 3;
  }
  if (mesh.points != inner_points + outer_points ||
      mesh.segments != inner_points + outer_points || mesh.triangles <= 0) {
    std::cerr << "unexpected mesh counts\n";
    return 4;
  }

  for (int point = 0; point < mesh.points; ++point) {
    if (mesh.coordinates[static_cast<std::size_t>(point)].first !=
          boundary_R[static_cast<std::size_t>(point)] ||
        mesh.coordinates[static_cast<std::size_t>(point)].second !=
          boundary_Z[static_cast<std::size_t>(point)]) {
      std::cerr << "boundary point order or coordinates changed\n";
      return 5;
    }
  }

  std::set<int> actual_edges;
  for (const auto &segment : mesh.boundary) {
    actual_edges.insert(canonical_edge(segment.first, segment.second));
  }
  std::set<int> triangle_edges;
  for (const auto &element : mesh.elements) {
    triangle_edges.insert(canonical_edge(element[0], element[1]));
    triangle_edges.insert(canonical_edge(element[1], element[2]));
    triangle_edges.insert(canonical_edge(element[2], element[0]));
  }
  for (int point = 0; point < inner_points; ++point) {
    const int next = (point + 1) % inner_points + 1;
    if (!actual_edges.count(canonical_edge(point + 1, next)) ||
        !triangle_edges.count(canonical_edge(point + 1, next))) {
      std::cerr << "inner boundary segment was lost\n";
      return 6;
    }
  }
  for (int point = 0; point < outer_points; ++point) {
    const int current = inner_points + point + 1;
    const int next = inner_points + (point + 1) % outer_points + 1;
    if (!actual_edges.count(canonical_edge(current, next)) ||
        !triangle_edges.count(canonical_edge(current, next))) {
      std::cerr << "outer boundary segment was lost\n";
      return 7;
    }
  }

  double mesh_area = 0.0;
  double minimum_angle = std::numeric_limits<double>::max();
  for (const auto &element : mesh.elements) {
    const auto &a = mesh.coordinates[static_cast<std::size_t>(element[0] - 1)];
    const auto &b = mesh.coordinates[static_cast<std::size_t>(element[1] - 1)];
    const auto &c = mesh.coordinates[static_cast<std::size_t>(element[2] - 1)];
    const double signed_area = area(a, b, c);
    if (!(signed_area > 0.0)) {
      std::cerr << "triangle is not counterclockwise\n";
      return 8;
    }
    mesh_area += signed_area;
    const double ab2 = (a.first - b.first) * (a.first - b.first) +
                       (a.second - b.second) * (a.second - b.second);
    const double bc2 = (b.first - c.first) * (b.first - c.first) +
                       (b.second - c.second) * (b.second - c.second);
    const double ca2 = (c.first - a.first) * (c.first - a.first) +
                       (c.second - a.second) * (c.second - a.second);
    minimum_angle = std::min(minimum_angle, angle(ab2, bc2, ca2));
    minimum_angle = std::min(minimum_angle, angle(bc2, ca2, ab2));
    minimum_angle = std::min(minimum_angle, angle(ca2, ab2, bc2));

    const double centroid_R = (a.first + b.first + c.first) / 3.0;
    const double centroid_Z = (a.second + b.second + c.second) / 3.0;
    if (!point_in_polygon(mesh.coordinates, inner_points, outer_points,
                          centroid_R, centroid_Z) ||
        point_in_polygon(mesh.coordinates, 0, inner_points,
                         centroid_R, centroid_Z)) {
      std::cerr << "triangle centroid lies outside the polygonal annulus\n";
      return 9;
    }
  }

  const double expected_area =
    polygon_area(boundary_R, boundary_Z, inner_points, outer_points) -
    polygon_area(boundary_R, boundary_Z, 0, inner_points);
  if (std::abs(mesh_area - expected_area) / expected_area > 1.0e-10) {
    std::cerr << "annulus area mismatch: " << mesh_area << " vs "
              << expected_area << '\n';
    return 10;
  }
  std::cout << "Delaunay32 annulus mesh: " << mesh.points << " points, "
            << mesh.triangles << " triangles, minimum angle "
            << minimum_angle << " degrees\n";
  return 0;
}
