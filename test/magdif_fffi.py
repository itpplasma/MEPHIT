from os.path import expanduser
from fffi import FortranLibrary, FortranModule

libmagdif = FortranLibrary('magdif', path=expanduser('~/src/NEO-EQ/build/lib'))
magdif = FortranModule(libmagdif, name='magdif')

magdif.fdef("""
    subroutine read_mesh
    end subroutine
""")

libmagdif.compile(verbose=1)
libmagdif.load()
