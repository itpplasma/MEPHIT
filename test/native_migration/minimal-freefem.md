# FreeFem runtime without external sparse factors

The tested source is official FreeFem 4.17, commit
`754dae8ae8f858784ad46b2f46dfd7862829cb5a` from
<https://github.com/FreeFem/FreeFem-sources>. No source patches are required.
The optional FortSparse module remains a separate SDK-linked component;
FreeFem's license and SDK terms apply to that component. This recipe installs
into a new local prefix and leaves the installed interpreter untouched.

Prerequisites include GNU C/C++, Autotools, OpenBLAS/LAPACK, and the ordinary
FreeFem parser, graphics, ARPACK and HDF5 build dependencies. The tested host
used GNU 16.2 and OpenBLAS 0.3.34. This is a serial interpreter; HDF5 may still
bring transitive MPI DSOs. The guarantee concerns external sparse factors.

```sh
git clone https://github.com/FreeFem/FreeFem-sources.git freefem-source
git -C freefem-source checkout --detach 754dae8ae8f858784ad46b2f46dfd7862829cb5a
mkdir freefem-build-source
git -C freefem-source archive HEAD | tar -x -C freefem-build-source
freefem_prefix="$PWD/freefem-prefix"
cd freefem-build-source
autoreconf -i
./configure --prefix="$freefem_prefix" --disable-download \
  --disable-system-umfpack --disable-umfpack \
  --without-petsc --without-petsc_complex --without-mpi \
  --disable-mumps --disable-mumps_seq --disable-superlu \
  --disable-scalapack --disable-metis --disable-parmetis --disable-scotch \
  --with-blas=-lopenblas --with-lapack=-lopenblas \
  CXXFLAGS='-O2 -g -DNDEBUG -std=gnu++14' CFLAGS='-O2 -g'
make -j2 -C src/libMesh
make -j2 -C src/lglib
make -j2 -C src/fflib
make -j2 -C src/nw FreeFem++
make -j2 -C plugin/seq ff-c++ include.done
make -j2 -C plugin/seq BinaryIO.so
make -C src/nw install-binPROGRAMS bin_PROGRAMS=FreeFem++
make -j2 -C plugin/seq install-exec-local \
  LIST_COMPILE=BinaryIO.so LIST_COMPILE_PKG= MPIPROG=
mkdir -p "$freefem_prefix/lib/ff++/4.17/idp"
cp idp/*.idp "$freefem_prefix/lib/ff++/4.17/idp/"
```

`--disable-umfpack` alone disables an optional plugin; the separate
`--disable-system-umfpack` is necessary to disable core UMFPACK/CHOLMOD
autodetection. The complex PETSc option uses an underscore. Limiting the
installation lists avoids building unrelated plugins or downloading solvers.

Configure MEPHIT with this exact executable and matching SDK:

```sh
cmake -S . -B build-native \
  -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=OpenBLAS \
  -DFREEFEM_EXECUTABLE="$freefem_prefix/bin/FreeFem++" \
  -DFREEFEM_CXX_EXECUTABLE="$freefem_prefix/bin/ff-c++"
cmake --build build-native --target mephit_run.x fortsparse_freefem -j2
ctest --test-dir build-native -R '^(native_freefem|freefem_protocol)' \
  --output-on-failure -j2
```

Verify actual `ldd` closure of the interpreter, `BinaryIO.so`, the generated
`scripts/fortsparse_native.so`, and `lib/libmephit.so`: none may load UMFPACK,
CHOLMOD, SuperLU or MUMPS. BSD AMD provides ordering only. Dense BLAS/LAPACK and
ARPACK remain permitted. A valid stock interpreter can also execute native
factors, but its unused linked vendor DSOs do not satisfy this stronger closure
claim. Keep the two configurations explicit in test and resource results.
