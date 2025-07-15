# MEPHIT Fortran Delaunay Triangulation Implementation Plan

## Overview

This plan outlines the implementation of a custom Fortran Delaunay triangulation library to replace the external TRIANGLE dependency in MEPHIT. The implementation will be based on standard algorithms from computational geometry literature and will be validated against mathematical Delaunay criteria rather than exact output matching.

## Theoretical Foundation

### Delaunay Triangulation Properties
1. **Empty circumcircle property**: No point lies inside the circumcircle of any triangle
2. **Maximizes minimum angle**: Among all triangulations, maximizes the minimum angle
3. **Unique for non-degenerate point sets**: Given a set of points, the Delaunay triangulation is unique
4. **Dual to Voronoi diagram**: Each triangle corresponds to a Voronoi vertex

### Reference Algorithms
- **Bowyer-Watson algorithm** (Bowyer 1981, Watson 1981): Incremental point insertion
- **Lawson's edge flipping** (Lawson 1977): Local optimization
- **Constrained Delaunay** (Chew 1989): Handling boundary constraints
- **Ruppert's algorithm** (Ruppert 1995): Quality mesh generation

## Current TRIANGLE Usage Analysis

### Interface Requirements
- **Input**: Boundary points, segments, hole specifications
- **Output**: Triangulated mesh with points, triangles, segments
- **Constraints**: Minimum angle quality, boundary preservation
- **File format**: Simple text format for mesh data

### Key Function: `FEM_triangulate_external`
- **Purpose**: Creates triangulation for plasma mesh generation
- **Input geometry**: Inner/outer boundaries, hole centers
- **Quality requirements**: Angle constraints, boundary preservation
- **Output**: Mesh file for finite element calculations

## Implementation Architecture

### Module Structure
```
src/
├── triangulation_fortran.f90          # Main API module
├── delaunay_types.f90                 # Data structures
├── geometric_predicates.f90           # Robust geometric tests
├── bowyer_watson.f90                  # Core triangulation algorithm
├── constrained_delaunay.f90           # Boundary constraints
├── mesh_quality.f90                   # Quality improvement
└── triangulation_io.f90               # Input/output handling
```

## Detailed Implementation Plan

### Phase 1: Core Data Structures and Predicates

#### 1.1 Data Structures (`src/delaunay_types.f90`)
```fortran
module delaunay_types
    type :: point_t
        real(dp) :: x, y
        integer :: id
    end type
    
    type :: triangle_t
        integer :: vertices(3)  ! Point indices
        integer :: neighbors(3) ! Neighbor triangle indices
        logical :: valid
    end type
    
    type :: edge_t
        integer :: endpoints(2)
        logical :: constrained
    end type
    
    type :: mesh_t
        type(point_t), allocatable :: points(:)
        type(triangle_t), allocatable :: triangles(:)
        type(edge_t), allocatable :: edges(:)
        integer :: npoints, ntriangles, nedges
    end type
end module
```

#### 1.2 Geometric Predicates (`src/geometric_predicates.f90`)
Based on Shewchuk's robust predicates (1997):
```fortran
module geometric_predicates
    ! Orientation test: CCW, CW, or collinear
    function orientation(p1, p2, p3) result(orient)
    
    ! In-circle test: point d inside circumcircle of abc
    function in_circle(pa, pb, pc, pd) result(inside)
    
    ! Point-in-triangle test
    function point_in_triangle(p, t) result(inside)
    
    ! Circumcenter calculation
    function circumcenter(pa, pb, pc) result(center)
end module
```

### Phase 2: Core Delaunay Algorithm

#### 2.1 Bowyer-Watson Algorithm (`src/bowyer_watson.f90`)
Implementation based on Bowyer (1981) and Watson (1981):
```fortran
module bowyer_watson
    ! Main triangulation routine
    subroutine delaunay_triangulate(points, mesh)
    
    ! Insert single point into existing triangulation
    subroutine insert_point(mesh, point_idx)
    
    ! Find triangles whose circumcircles contain the point
    subroutine find_cavity(mesh, point_idx, cavity_triangles)
    
    ! Create new triangles connecting point to cavity boundary
    subroutine fill_cavity(mesh, point_idx, cavity_boundary)
end module
```

**Algorithm Steps:**
1. Create super-triangle containing all points
2. For each point:
   - Find all triangles whose circumcircles contain the point
   - Remove these triangles (creating a cavity)
   - Connect the point to the cavity boundary
3. Remove super-triangle and its associated triangles

### Phase 3: Constrained Delaunay Implementation

#### 3.1 Constraint Handling (`src/constrained_delaunay.f90`)
Based on Chew (1989) and Sloan (1993):
```fortran
module constrained_delaunay
    ! Insert constrained edge into triangulation
    subroutine insert_constraint(mesh, edge)
    
    ! Recover missing constraint edges
    subroutine recover_edges(mesh, constraints)
    
    ! Split intersecting triangles
    subroutine split_intersecting_triangles(mesh, edge)
end module
```

**Algorithm Steps:**
1. Start with unconstrained Delaunay triangulation
2. For each constraint edge:
   - Find triangles intersecting the edge
   - Split or remove intersecting triangles
   - Insert the constraint edge
3. Locally optimize around constraints

### Phase 4: Quality Improvement

#### 4.1 Mesh Quality (`src/mesh_quality.f90`)
Based on Ruppert (1995) and Chew (1993):
```fortran
module mesh_quality
    ! Refine triangulation to meet quality constraints
    subroutine refine_mesh(mesh, min_angle, max_area)
    
    ! Identify poor-quality triangles
    subroutine find_bad_triangles(mesh, min_angle, bad_triangles)
    
    ! Insert Steiner points to improve quality
    subroutine insert_steiner_points(mesh, bad_triangles)
end module
```

**Quality Criteria:**
- **Minimum angle constraint**: No triangle has angle < threshold
- **Maximum area constraint**: No triangle has area > threshold
- **Aspect ratio**: Limit triangle elongation

### Phase 5: Integration Layer

#### 5.1 Main API (`src/triangulation_fortran.f90`)
```fortran
module triangulation_fortran
    ! Main triangulation interface
    subroutine triangulate_domain(boundary_points, constraints, holes, &
                                 options, result_mesh)
    
    ! Quality triangulation with angle constraints
    subroutine triangulate_quality(boundary_points, constraints, &
                                  min_angle, result_mesh)
    
    ! Triangulation with holes
    subroutine triangulate_with_holes(boundary_points, constraints, &
                                     hole_points, result_mesh)
end module
```

#### 5.2 C Interface (`src/triangulation_c_interface.c`)
Replacement for existing TRIANGLE calls:
```c
// Drop-in replacement for FEM_triangulate_external
void FEM_triangulate_external_fortran(
    const int npt_inner,
    const int npt_outer,
    const double *bdry_R,
    const double *bdry_Z,
    const double R_mid,
    const double Z_mid,
    const char *fname);
```

## Validation Strategy

### Mathematical Validation
1. **Delaunay property**: Verify empty circumcircle property
2. **Convex hull preservation**: Boundary triangles form convex hull
3. **Topology validation**: Euler characteristic V - E + F = 2
4. **Angle constraints**: Verify minimum angle requirements
5. **Area constraints**: Verify maximum area requirements

### Test Cases
```fortran
! Test Delaunay property for all triangles
subroutine test_delaunay_property(mesh)
    do i = 1, mesh%ntriangles
        call verify_empty_circumcircle(mesh, i)
    end do
end subroutine

! Test constraint preservation
subroutine test_constraints(mesh, constraints)
    do i = 1, size(constraints)
        call verify_constraint_present(mesh, constraints(i))
    end do
end subroutine
```

### Performance Validation
- **Complexity**: O(n log n) expected for n points
- **Memory usage**: O(n) space complexity
- **Robustness**: Handle degenerate cases gracefully

## Implementation Schedule

### Week 1: Foundation
- [ ] Implement data structures and types
- [ ] Implement geometric predicates with robust arithmetic
- [ ] Write comprehensive tests for predicates
- [ ] Validate against known geometric cases

### Week 2: Core Algorithm
- [ ] Implement Bowyer-Watson algorithm
- [ ] Handle super-triangle initialization and removal
- [ ] Test with simple point sets
- [ ] Validate Delaunay property

### Week 3: Constraints
- [ ] Implement constraint edge insertion
- [ ] Handle edge recovery algorithms
- [ ] Test with boundary constraints
- [ ] Validate constraint preservation

### Week 4: Quality Improvement
- [ ] Implement bad triangle detection
- [ ] Add Steiner point insertion
- [ ] Test angle and area constraints
- [ ] Validate mesh quality metrics

### Week 5: Advanced Features
- [ ] Implement hole handling
- [ ] Add support for multiple boundaries
- [ ] Test complex geometries
- [ ] Validate topology preservation

### Week 6: Integration
- [ ] Implement C interface layer
- [ ] Update build system (CMake)
- [ ] Test integration with existing MEPHIT code
- [ ] Performance benchmarking

### Week 7: Validation and Documentation
- [ ] Comprehensive test suite
- [ ] Physics validation with MEPHIT
- [ ] Performance optimization
- [ ] Documentation and examples

## Quality Assurance

### Code Quality
- **Fortran standards**: Modern Fortran (2008+) with explicit interfaces
- **Error handling**: Graceful handling of invalid inputs
- **Memory management**: Proper allocation/deallocation
- **Documentation**: Comprehensive inline documentation

### Numerical Robustness
- **Exact arithmetic**: Use rational arithmetic for critical predicates
- **Adaptive precision**: Handle near-degenerate cases
- **Stability**: Consistent results across platforms
- **Tolerance handling**: Appropriate geometric tolerances

### Testing Coverage
- **Unit tests**: Test each module independently
- **Integration tests**: Test complete workflows
- **Edge cases**: Degenerate geometries, boundary cases
- **Performance tests**: Large-scale problem validation

## Risk Mitigation

### Technical Risks
- **Numerical issues**: Use proven robust predicates
- **Performance**: Profile and optimize critical paths
- **Complexity**: Incremental implementation with validation
- **Integration**: Test with existing MEPHIT code early

### Timeline Risks
- **Scope management**: Focus on core functionality first
- **Testing overhead**: Automated testing from day one
- **Documentation**: Document during development
- **Validation**: Continuous validation against mathematical criteria

## Success Criteria

1. **Correctness**: All triangulations satisfy Delaunay property
2. **Robustness**: Handle all reasonable input geometries
3. **Performance**: Acceptable performance for MEPHIT use cases
4. **Integration**: Seamless replacement of TRIANGLE dependency
5. **Maintainability**: Clean, documented, testable code

## References

1. Bowyer, A. (1981). "Computing Dirichlet tessellations"
2. Watson, D.F. (1981). "Computing the n-dimensional Delaunay tessellation"
3. Lawson, C.L. (1977). "Software for C1 surface interpolation"
4. Chew, L.P. (1989). "Constrained Delaunay triangulations"
5. Ruppert, J. (1995). "A Delaunay refinement algorithm for quality 2-dimensional mesh generation"
6. Shewchuk, J.R. (1997). "Adaptive precision floating-point arithmetic and fast robust geometric predicates"

This implementation will provide MEPHIT with a robust, maintainable triangulation capability based on well-established algorithms from computational geometry literature.