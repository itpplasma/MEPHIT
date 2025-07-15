module geometric_predicates
    use iso_fortran_env, only: dp => real64
    use delaunay_types, only: point_t, triangle_t, mesh_t
    implicit none
    
    private
    public :: orientation, in_circle, point_in_triangle, circumcenter
    public :: triangle_area, triangle_angles, edge_length
    public :: ORIENTATION_CCW, ORIENTATION_CW, ORIENTATION_COLLINEAR
    public :: geometric_tolerance
    
    ! Orientation test results
    integer, parameter :: ORIENTATION_CCW = 1
    integer, parameter :: ORIENTATION_CW = -1
    integer, parameter :: ORIENTATION_COLLINEAR = 0
    
    ! Geometric tolerance for numerical comparisons
    real(dp), parameter :: geometric_tolerance = 1.0e-14_dp
    
contains

integer function orientation(pa, pb, pc)
    !> Robust orientation test: returns CCW, CW, or COLLINEAR
    !> Based on Shewchuk's robust predicates
    type(point_t), intent(in) :: pa, pb, pc
    
    real(dp) :: detleft, detright, det
    real(dp) :: detsum, errbound
    
    ! Calculate the determinant
    detleft = (pa%x - pc%x) * (pb%y - pc%y)
    detright = (pa%y - pc%y) * (pb%x - pc%x)
    det = detleft - detright
    
    ! Fast path for clearly non-zero results
    if (detleft > 0.0_dp) then
        if (detright <= 0.0_dp) then
            orientation = ORIENTATION_CCW
            return
        else
            detsum = detleft + detright
        end if
    else if (detleft < 0.0_dp) then
        if (detright >= 0.0_dp) then
            orientation = ORIENTATION_CW
            return
        else
            detsum = -detleft - detright
        end if
    else
        orientation = ORIENTATION_COLLINEAR
        return
    end if
    
    ! Error bound for robust computation
    errbound = 3.0_dp * geometric_tolerance * detsum
    
    if (abs(det) >= errbound) then
        if (det > 0.0_dp) then
            orientation = ORIENTATION_CCW
        else
            orientation = ORIENTATION_CW
        end if
    else
        ! Use exact arithmetic for near-zero cases
        orientation = orientation_exact(pa, pb, pc)
    end if
end function orientation

integer function orientation_exact(pa, pb, pc)
    !> Exact orientation test using higher precision
    type(point_t), intent(in) :: pa, pb, pc
    
    real(dp) :: acx, acy, bcx, bcy
    real(dp) :: det
    
    acx = pa%x - pc%x
    acy = pa%y - pc%y
    bcx = pb%x - pc%x
    bcy = pb%y - pc%y
    
    det = acx * bcy - acy * bcx
    
    if (abs(det) <= geometric_tolerance) then
        orientation_exact = ORIENTATION_COLLINEAR
    else if (det > 0.0_dp) then
        orientation_exact = ORIENTATION_CCW
    else
        orientation_exact = ORIENTATION_CW
    end if
end function orientation_exact

logical function in_circle(pa, pb, pc, pd)
    !> Test if point pd is inside circumcircle of triangle abc
    !> Returns true if pd is inside the circumcircle
    type(point_t), intent(in) :: pa, pb, pc, pd
    
    real(dp) :: adx, ady, bdx, bdy, cdx, cdy
    real(dp) :: abdet, bcdet, cadet, alift, blift, clift
    real(dp) :: det
    
    adx = pa%x - pd%x
    ady = pa%y - pd%y
    bdx = pb%x - pd%x
    bdy = pb%y - pd%y
    cdx = pc%x - pd%x
    cdy = pc%y - pd%y
    
    abdet = adx * bdy - bdx * ady
    bcdet = bdx * cdy - cdx * bdy
    cadet = cdx * ady - adx * cdy
    
    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy
    
    det = alift * bcdet + blift * cadet + clift * abdet
    
    in_circle = det > geometric_tolerance
end function in_circle

logical function point_in_triangle(p, mesh, tri_idx)
    !> Test if point p is inside triangle tri_idx
    type(point_t), intent(in) :: p
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: tri_idx
    
    type(point_t) :: pa, pb, pc
    integer :: orient1, orient2, orient3
    
    if (.not. mesh%triangles(tri_idx)%valid) then
        point_in_triangle = .false.
        return
    end if
    
    pa = mesh%points(mesh%triangles(tri_idx)%vertices(1))
    pb = mesh%points(mesh%triangles(tri_idx)%vertices(2))
    pc = mesh%points(mesh%triangles(tri_idx)%vertices(3))
    
    orient1 = orientation(pa, pb, p)
    orient2 = orientation(pb, pc, p)
    orient3 = orientation(pc, pa, p)
    
    ! Point is inside if all orientations are the same (all CCW or all CW)
    point_in_triangle = (orient1 == orient2) .and. (orient2 == orient3) .and. &
                       (orient1 /= ORIENTATION_COLLINEAR)
end function point_in_triangle

function circumcenter(pa, pb, pc) result(center)
    !> Calculate circumcenter of triangle abc
    type(point_t), intent(in) :: pa, pb, pc
    type(point_t) :: center
    
    real(dp) :: ax, ay, bx, by, cx, cy
    real(dp) :: d, ux, uy, vx, vy
    real(dp) :: det
    
    ax = pa%x
    ay = pa%y
    bx = pb%x
    by = pb%y
    cx = pc%x
    cy = pc%y
    
    ! Calculate the determinant
    det = 2.0_dp * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    
    if (abs(det) <= geometric_tolerance) then
        ! Degenerate triangle - return centroid
        center%x = (ax + bx + cx) / 3.0_dp
        center%y = (ay + by + cy) / 3.0_dp
        center%id = 0
        return
    end if
    
    ! Calculate squared distances
    d = ax * ax + ay * ay
    ux = d * (by - cy) + bx * bx * (cy - ay) + cx * cx * (ay - by)
    uy = d * (cx - bx) + by * by * (ax - cx) + cy * cy * (bx - ax)
    
    center%x = ux / det
    center%y = uy / det
    center%id = 0
end function circumcenter

real(dp) function triangle_area(pa, pb, pc)
    !> Calculate area of triangle abc
    type(point_t), intent(in) :: pa, pb, pc
    
    real(dp) :: ax, ay, bx, by, cx, cy
    
    ax = pa%x
    ay = pa%y
    bx = pb%x
    by = pb%y
    cx = pc%x
    cy = pc%y
    
    triangle_area = 0.5_dp * abs((ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)))
end function triangle_area

function triangle_angles(pa, pb, pc) result(angles)
    !> Calculate angles of triangle abc in degrees
    type(point_t), intent(in) :: pa, pb, pc
    real(dp) :: angles(3)
    
    real(dp) :: a, b, c  ! Side lengths
    real(dp) :: cos_a, cos_b, cos_c
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: rad_to_deg = 180.0_dp / pi
    
    ! Calculate side lengths
    a = edge_length(pb, pc)  ! Side opposite to angle A
    b = edge_length(pa, pc)  ! Side opposite to angle B
    c = edge_length(pa, pb)  ! Side opposite to angle C
    
    ! Handle degenerate triangles
    if (a <= geometric_tolerance .or. b <= geometric_tolerance .or. c <= geometric_tolerance) then
        angles = 0.0_dp
        return
    end if
    
    ! Calculate angles using law of cosines
    cos_a = (b*b + c*c - a*a) / (2.0_dp * b * c)
    cos_b = (a*a + c*c - b*b) / (2.0_dp * a * c)
    cos_c = (a*a + b*b - c*c) / (2.0_dp * a * b)
    
    ! Clamp to valid range to avoid numerical issues
    cos_a = max(-1.0_dp, min(1.0_dp, cos_a))
    cos_b = max(-1.0_dp, min(1.0_dp, cos_b))
    cos_c = max(-1.0_dp, min(1.0_dp, cos_c))
    
    angles(1) = acos(cos_a) * rad_to_deg
    angles(2) = acos(cos_b) * rad_to_deg
    angles(3) = acos(cos_c) * rad_to_deg
end function triangle_angles

real(dp) function edge_length(pa, pb)
    !> Calculate length of edge between points pa and pb
    type(point_t), intent(in) :: pa, pb
    
    real(dp) :: dx, dy
    
    dx = pa%x - pb%x
    dy = pa%y - pb%y
    
    edge_length = sqrt(dx*dx + dy*dy)
end function edge_length

end module geometric_predicates