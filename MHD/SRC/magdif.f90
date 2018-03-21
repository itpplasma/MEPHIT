module magdif
  use from_nrtype, only: dp                                    ! PRELOAD/SRC/from_nrtype.f90
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element,& ! PRELOAD/SRC/mesh_mod.f90 
       mesh_element_rmp, bphicovar
  
  implicit none

  integer log_level
  character (len=1024) :: point_file, tri_file  ! points.dat and triangles.dat from PRELOAD
  logical :: log_err, log_warn, log_info, log_debug

  namelist / settings / log_level, point_file, tri_file
  
  integer, parameter :: log = 6                 ! log to stdout, TODO: make this configurable 

  complex(dp), allocatable :: presn(:)          ! pressure perturbation p_n in each mesh point
  complex(dp), allocatable :: curr(:,:)         ! edge currents R j_n \cdot n weighted by R
  complex(dp), allocatable :: Bnflux(:,:)       ! edge fluxes R B_n \cdot n weighted by R
  complex(dp), allocatable :: Bnphi(:)          ! physical toroidal component of Bn   
  
contains
  
  subroutine magdif_init
    call read_config
    call read_mesh
  end subroutine magdif_init
  
  subroutine read_config
    open(1, file='magdif.in')
    read(1, nml=settings)
    close(1)

    log_err = .false.
    log_warn = .false.
    log_info = .false.
    log_debug = .false.
    if (log_level>0) log_err = .true.
    if (log_level>1) log_warn = .true.
    if (log_level>2) log_info = .true.
    if (log_level>3) log_debug = .true.
  end subroutine read_config
  
  subroutine read_mesh
    open(1, file=point_file, form='unformatted')
    read(1) npoint
    if (log_info) write(log,*) 'npoint = ', npoint
    allocate(mesh_point(npoint))
    allocate(presn(npoint))
    read(1) mesh_point
    close(1)

    open(1, file=tri_file, form='unformatted')
    read(1) ntri
    if (log_info) write(log,*) 'ntri   = ', ntri
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))
    allocate(Bnflux(ntri,3))
    allocate(Bnphi(ntri))
    allocate(curr(ntri,3))
    read(1) mesh_element
    read(1) bphicovar
    close(1)
  end subroutine read_mesh
  
end module magdif
