MODULE MOD_io_bin

    IMPLICIT NONE

    ! BINARY IO
    INTERFACE prt_bin
        MODULE PROCEDURE store_binary_vector_real4
        MODULE PROCEDURE store_binary_vector_real8
        MODULE PROCEDURE store_binary_vector_integer
        MODULE PROCEDURE store_binary_matrix_real4
        MODULE PROCEDURE store_binary_matrix_real8
        MODULE PROCEDURE store_binary_matrix_integer
        MODULE PROCEDURE store_binary_tensor3_real4
        MODULE PROCEDURE store_binary_tensor3_real8
        MODULE PROCEDURE store_binary_tensor3_integer
    END INTERFACE prt_bin

    INTERFACE load_bin
        MODULE PROCEDURE load_binary_vector_real4
        MODULE PROCEDURE load_binary_vector_real8
        MODULE PROCEDURE load_binary_vector_integer
        MODULE PROCEDURE load_binary_matrix_real4
        MODULE PROCEDURE load_binary_matrix_real8
        MODULE PROCEDURE load_binary_matrix_integer
        MODULE PROCEDURE load_binary_tensor3_real4
        MODULE PROCEDURE load_binary_tensor3_real8
        MODULE PROCEDURE load_binary_tensor3_integer
    END INTERFACE load_bin

    PRIVATE
    ! BINARY IO
    PUBLIC :: prt_bin
    PUBLIC :: load_bin

CONTAINS

! BINARY INPUT/OUTPUT

   SUBROUTINE store_binary_vector_real4(V, filename)
        REAL(4), DIMENSION(:), INTENT(IN) :: V
        CHARACTER(len=*),      INTENT(IN) :: filename
        INTEGER, DIMENSION(1)             :: V_shape
        INTEGER, PARAMETER                :: type_label = 4
        V_shape = shape(V)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) V_shape
        WRITE(6550) V
        CLOSE(6550)
    END SUBROUTINE store_binary_vector_real4

    SUBROUTINE load_binary_vector_real4(filename, V)
        CHARACTER(len=*),                   INTENT(IN)  :: filename
        REAL(4), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V
        INTEGER, DIMENSION(1)                           :: V_shape
        INTEGER                                         :: type_label
        IF (ALLOCATED(V)) DEALLOCATE(V)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) V_shape
        ALLOCATE(V(V_shape(1)))
        READ (6551) V
        CLOSE(6551)
    END SUBROUTINE load_binary_vector_real4

    SUBROUTINE store_binary_matrix_real4(M, filename)
        REAL(4), DIMENSION(:,:), INTENT(IN) :: M
        CHARACTER(len=*),        INTENT(IN) :: filename
        INTEGER, DIMENSION(2)               :: M_shape
        INTEGER, PARAMETER                  :: type_label = 4
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix_real4

    SUBROUTINE load_binary_matrix_real4(filename, M)
        CHARACTER(len=*), INTENT(IN)                      :: filename
        REAL(4), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(2)                             :: M_shape
        INTEGER                                           :: type_label
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix_real4

    SUBROUTINE store_binary_tensor3_real4(M, filename)
        REAL(4), DIMENSION(:,:,:), INTENT(IN) :: M
        CHARACTER(len=*),          INTENT(IN) :: filename
        INTEGER, DIMENSION(3)                 :: M_shape
        INTEGER, PARAMETER                    :: type_label = 4
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_tensor3_real4

    SUBROUTINE load_binary_tensor3_real4(filename, M)
        CHARACTER(len=*),                       INTENT(IN)  :: filename
        REAL(4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(3)                               :: M_shape
        INTEGER                                             :: type_label
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2), M_shape(3)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_tensor3_real4

    SUBROUTINE store_binary_vector_real8(V, filename)
        REAL(8), DIMENSION(:), INTENT(IN) :: V
        CHARACTER(len=*),      INTENT(IN) :: filename
        INTEGER, DIMENSION(1)             :: V_shape
        INTEGER, PARAMETER                :: type_label = 8
        V_shape = shape(V)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) V_shape
        WRITE(6550) V
        CLOSE(6550)
    END SUBROUTINE store_binary_vector_real8

    SUBROUTINE load_binary_vector_real8(filename, V)
        CHARACTER(len=*),                   INTENT(IN)  :: filename
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V
        INTEGER, DIMENSION(1)                           :: V_shape
        INTEGER                                         :: type_label
        IF (ALLOCATED(V)) DEALLOCATE(V)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) V_shape
        ALLOCATE(V(V_shape(1)))
        READ (6551) V
        CLOSE(6551)
    END SUBROUTINE load_binary_vector_real8

    SUBROUTINE store_binary_matrix_real8(M, filename)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: M
        CHARACTER(len=*),        INTENT(IN) :: filename
        INTEGER, DIMENSION(2)               :: M_shape
        INTEGER, PARAMETER                  :: type_label = 8
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix_real8

    SUBROUTINE load_binary_matrix_real8(filename, M)
        CHARACTER(len=*),                     INTENT(IN)  :: filename
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(2)                             :: M_shape
        INTEGER                                           :: type_label
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix_real8

    SUBROUTINE store_binary_tensor3_real8(M, filename)
        REAL(8), DIMENSION(:,:,:), INTENT(IN) :: M
        CHARACTER(len=*),          INTENT(IN) :: filename
        INTEGER, DIMENSION(3)                 :: M_shape
        INTEGER, PARAMETER                    :: type_label = 8
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_tensor3_real8

    SUBROUTINE load_binary_tensor3_real8(filename, M)
        CHARACTER(len=*),                       INTENT(IN)  :: filename
        REAL(8), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(3)                               :: M_shape
        INTEGER                                             :: type_label
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2), M_shape(3)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_tensor3_real8

    SUBROUTINE store_binary_vector_integer(V, filename)
        INTEGER, DIMENSION(:), INTENT(IN) :: V
        CHARACTER(len=*),      INTENT(IN) :: filename
        INTEGER, DIMENSION(1)             :: V_shape
        INTEGER, PARAMETER                :: type_label = 1
        V_shape = shape(V)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) V_shape
        WRITE(6550) V
        CLOSE(6550)
    END SUBROUTINE store_binary_vector_integer

    SUBROUTINE load_binary_vector_integer(filename, V)
        CHARACTER(len=*),                   INTENT(IN)  :: filename
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V
        INTEGER, DIMENSION(1)                           :: V_shape
        INTEGER                                         :: type_label
        IF (ALLOCATED(V)) DEALLOCATE(V)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) V_shape
        ALLOCATE(V(V_shape(1)))
        READ (6551) V
        CLOSE(6551)
    END SUBROUTINE load_binary_vector_integer

    SUBROUTINE store_binary_matrix_integer(M, filename)
        INTEGER, DIMENSION(:,:), INTENT(IN) :: M
        CHARACTER(len=*),        INTENT(IN) :: filename
        INTEGER, DIMENSION(2)               :: M_shape
        INTEGER, PARAMETER                  :: type_label = 1
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix_integer

    SUBROUTINE load_binary_matrix_integer(filename, M)
        CHARACTER(len=*),                     INTENT(IN)  :: filename
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(2)                             :: M_shape
        INTEGER                                           :: type_label
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix_integer

    SUBROUTINE store_binary_tensor3_integer(M, filename)
        INTEGER, DIMENSION(:,:,:), INTENT(IN) :: M
        CHARACTER(len=*),          INTENT(IN) :: filename
        INTEGER, DIMENSION(3)                 :: M_shape
        INTEGER, PARAMETER                    :: type_label = 1
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) type_label
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_tensor3_integer

    SUBROUTINE load_binary_tensor3_integer(filename, M)
        CHARACTER(len=*),                       INTENT(IN)  :: filename
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(3)                               :: M_shape
        INTEGER                                             :: type_label
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) type_label
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2), M_shape(3)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_tensor3_integer

END MODULE MOD_io_bin

! def load_fortran_binary(filename, transpose=False):
!     # Fortran data structure:
!     #     4 byte head --> label
!     #     N byte zeroth record (int32, a label for the type)
!     #     4 byte foot == head
!     #     4 byte head --> N
!     #     N byte first record (int32, the shape of the object contained in the second record)
!     #     4 byte foot == head
!     #     4 byte head --> M
!     #     M byte second record (int32 or float64, the matrix or vector)
!     #     4 byte foot == head
!     #
!     # On the use of transpose=True/False:
!     # Fortran will serialize data in Fortran order (in M(i,j,k), is the fast running index)
!     # if we load this data, and want to keep it as it was in Fortran, then we have to use
!     # transpose=False as default (which counterintuitively has to call transpose())
!     #
!     # additional note: if you have a 2d matrix that you want to load efficiently in python,
!     # then transpose it in Fortran before storing it, and load it with transpose, which will
!     # yield back the matrix as it was before its transpose operation in Fortran.
!     # this allows to skip a transpose call in python.
!     with open(filename) as f:
!         # read shape
!         record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
!         label         = np.fromfile(f, dtype=np.int32, count=record_length//4, sep="")
!         dummy         = np.fromfile(f, dtype=np.int32, count=1, sep="")
!         record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
!         shape         = np.fromfile(f, dtype=np.int32, count=record_length//4, sep="")
!         dummy         = np.fromfile(f, dtype=np.int32, count=1, sep="")
!         # compute matrix size
!         nelem = np.prod(shape)
!         # here we use the record length to determine the type of data (8bit/element -> real, 4bit/element -> integer)
!         record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
!         # select data type
!         if label == 1:
!             dtype = np.int32
!         elif label == 4:
!             dtype = np.float32
!         elif label == 8:
!             dtype = np.float64
!         else:
!             print("Error! wrong type label")
!             return
!         # read data
!         if shape.size > 1:
!             if transpose:
!                 return np.fromfile(f, dtype=dtype, count=nelem, sep="").reshape(shape[::-1])
!             else:
!                 return np.fromfile(f, dtype=dtype, count=nelem, sep="").reshape(shape[::-1]).transpose()
!         else:
!             return np.fromfile(f, dtype=dtype, count=nelem, sep="").flatten()
