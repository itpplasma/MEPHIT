  ! Saving non-dummy SuperLU routines for future use
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_superlu_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER :: nrhs,ldb,n,info,iopt
    INTEGER :: talk

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    n = nrow
    nrhs = 1
    ldb = n

    IF (sparse_talk) THEN
       talk = 1
    ELSE
       talk = 0
    END IF

    info = 0

    ! First, factorize the matrix. The factors are stored in *factors* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       iopt = 1
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )


       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
          ELSE
             PRINT *, 'INFO from factorization = ', info
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       iopt = 2
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )

       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
             ! WRITE(*,*) (b(i), i=1, n)
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuperLU
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       iopt = 3
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )
       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Free succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    RETURN
  END SUBROUTINE sparse_solve_superlu_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for complex sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_superluComplex_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER :: nrhs,ldb,n,info,iopt
    INTEGER :: talk

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    n = nrow
    nrhs = 1
    ldb = n

    IF (sparse_talk) THEN
       talk = 1
    ELSE
       talk = 0
    END IF

    info = 0

    ! First, factorize the matrix. The factors are stored in *factors* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       iopt = 1
       CALL c_fortran_zgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )


       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
          ELSE
             PRINT *, 'INFO from factorization = ', info
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       iopt = 2
       CALL c_fortran_zgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )

       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
             ! WRITE(*,*) (b(i), i=1, n)
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuperLU
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       iopt = 3
       CALL c_fortran_zgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )
       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Free succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    RETURN
  END SUBROUTINE sparse_solve_superluComplex_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_superlu_b2(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER :: nrhs,ldb,n,info,iopt
    INTEGER :: talk

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    n = nrow
    nrhs = SIZE(b,2)
    ldb = n

    IF (sparse_talk) THEN
       talk = 1
    ELSE
       talk = 0
    END IF

    info = 0

    ! First, factorize the matrix. The factors are stored in *factors* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       iopt = 1
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )

       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
          ELSE
             PRINT *, 'INFO from factorization = ', info
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       iopt = 2
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )

       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
             ! WRITE(*,*) (b(i), i=1, n)
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuperLU
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       iopt = 3
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            b, ldb, factors, info )
       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Free succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    RETURN
  END SUBROUTINE sparse_solve_superlu_b2
  !-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_superlu_b2_loop(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER :: nrhs,ldb,n,info,iopt
    INTEGER :: i
    INTEGER :: talk
    INTEGER :: info_store

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: bloc

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    n = nrow
    nrhs = 1
    ldb = n

    IF (ALLOCATED(bloc)) DEALLOCATE(bloc)
    ALLOCATE(bloc(nrow))
    bloc = 0.0_dp

    IF (sparse_talk) THEN
       talk = 1
    ELSE
       talk = 0
    END IF

    info_store = 0
    info = 0

    ! First, factorize the matrix. The factors are stored in *factors* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       iopt = 1
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            bloc, ldb, factors, info )

       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
          ELSE
             PRINT *, 'INFO from factorization = ', info
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       iopt = 2
       DO i = 1,SIZE(b,2)
          bloc = b(:,i)
          CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
               bloc, ldb, factors, info )

          IF (sparse_talk) THEN
             IF (info .EQ. 0) THEN
                !PRINT *, 'Solve succeeded',i
                ! WRITE(*,*) (b(i), i=1, n)
             ELSE
                info_store = info_store + info
                PRINT *, 'INFO from triangular solve = ', info
             ENDIF
          END IF

          b(:,i) = bloc
       END DO

       IF (sparse_talk) THEN
          IF (info_store .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve: failed ', info_store, ' times'
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuperLU
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       iopt = 3
       CALL c_fortran_dgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            bloc, ldb, factors, info )
       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Free succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    IF (ALLOCATED(bloc)) DEALLOCATE(bloc)
    RETURN
  END SUBROUTINE sparse_solve_superlu_b2_loop
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_superluComplex_b2_loop(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER :: nrhs,ldb,n,info,iopt
    INTEGER :: i
    INTEGER :: talk
    INTEGER :: info_store

    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: bloc

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    n = nrow
    nrhs = 1
    ldb = n

    IF (ALLOCATED(bloc)) DEALLOCATE(bloc)
    ALLOCATE(bloc(nrow))
    bloc = 0.0_dp

    IF (sparse_talk) THEN
       talk = 1
    ELSE
       talk = 0
    END IF

    info_store = 0
    info = 0

    ! First, factorize the matrix. The factors are stored in *factors* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       iopt = 1
       CALL c_fortran_zgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            bloc, ldb, factors, info )

       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
          ELSE
             PRINT *, 'INFO from factorization = ', info
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       iopt = 2
       DO i = 1,SIZE(b,2)
          bloc = b(:,i)
          CALL c_fortran_zgssv( iopt, n, nz, nrhs, val, irow, pcol, &
               bloc, ldb, factors, info )

          IF (sparse_talk) THEN
             IF (info .EQ. 0) THEN
                !PRINT *, 'Solve succeeded',i
                ! WRITE(*,*) (b(i), i=1, n)
             ELSE
                info_store = info_store + info
                PRINT *, 'INFO from triangular solve = ', info
             ENDIF
          END IF

          b(:,i) = bloc
       END DO

       IF (sparse_talk) THEN
          IF (info_store .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve: failed ', info_store, ' times'
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuperLU
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       iopt = 3
       CALL c_fortran_zgssv( iopt, n, nz, nrhs, val, irow, pcol, &
            bloc, ldb, factors, info )
       IF (sparse_talk) THEN
          IF (info .EQ. 0) THEN
             PRINT *, 'Free succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info
          ENDIF
       END IF
    END IF

    IF (ALLOCATED(bloc)) DEALLOCATE(bloc)
    RETURN
  END SUBROUTINE sparse_solve_superluComplex_b2_loop
  !-------------------------------------------------------------------------------
