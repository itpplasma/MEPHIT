from os.path import expanduser
from fffi import FortranLibrary, FortranModule

libmagdif = FortranLibrary('magdif', path=expanduser('~/src/NEO-EQ/build/lib'))
magdif = FortranModule(libmagdif, name='magdif')
magdif_config = FortranModule(libmagdif, name='magdif_config')

magdif.fdef("""
    subroutine magdif_init
    end

!    subroutine read_mesh
!    end
""")

magdif_config.fdef("""
    subroutine read_config
    end

!    subroutine log_open
!    end
""")

libmagdif.compile(verbose=1)
magdif.load()
magdif_config.load()
