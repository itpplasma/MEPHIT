from os.path import expanduser
from fffi import FortranLibrary, FortranModule

libmagdif = FortranLibrary('magdif', path=expanduser('~/src/NEO-EQ/build/lib'))
field_eq = FortranModule(libmagdif, name='field_eq_mod')
magdif = FortranModule(libmagdif, name='magdif')
magdif_config = FortranModule(libmagdif, name='magdif_config')

# TODO: add more definitions useful for testing
# Attention: private variables and routines are inaccessible

field_eq.fdef("""
  real(8) :: psif, psib
""")

magdif.fdef("""
  real(8), allocatable    :: B0flux(:,:)
  complex(8), allocatable :: Bnflux(:,:)
  complex(8), allocatable :: Bnphi(:)

  subroutine magdif_init
  end

  !subroutine read_mesh
  !end
""")

magdif_config.fdef("""
  subroutine read_config
  end

  !subroutine log_open
  !end
""")

libmagdif.fdef("""
  subroutine field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  
  real(8), intent(in) :: r,p,z
  real(8), intent(out) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  
  end subroutine field
  
""")

libmagdif.compile(verbose=1)
field_eq.load()
magdif.load()
magdif_config.load()
