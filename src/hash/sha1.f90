! The MIT License (MIT)

! Copyright (c) 2014 Jake Timothy

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!==================================================================================================================================
! This code was provided by Jake Timothy at https://github.com/jaketimothy/fortran-cryptography under the MIT License. A short and
! simple permissive license with conditions only requiring preservation of copyright and license notices. Licensed works,
! modifications, and larger works may be distributed under different terms and without source code.
!>> The version provided here was adjusted to FLEXI style and uses some intrinsic FUNCTIONs of FLEXI.
#include "flexi.h"

! Notes:
! No, fortran does not support unsigned integers.

MODULE MOD_SHA1
! MODULES
USE,INTRINSIC :: ISO_FORTRAN_ENV  ,ONLY: INT8,INT16,INT32,INT64
USE,INTRINSIC :: ISO_FORTRAN_ENV  ,ONLY: REAL32,REAL64,REAL128

IMPLICIT NONE
PRIVATE

INTERFACE SHA1
  MODULE PROCEDURE SHA1scalar, SHA1array
END INTERFACE

PUBLIC :: SHA1
!----------------------------------------------------------------------------------------------------------------------------------

CONTAINS

!==================================================================================================================================
!> SHA-1 interface function.
!==================================================================================================================================
FUNCTION SHA1Scalar(value) RESULT(hash)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(*),         INTENT(IN)  :: value
CHARACTER(LEN=40)             :: hash
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(INT8),DIMENSION(:),ALLOCATABLE :: bytes
INTEGER(INT64)                         :: length
!==================================================================================================================================

SELECT TYPE(value)
  TYPE IS (CHARACTER(LEN=*))
    length = LEN(value)
    ALLOCATE(bytes(((length+8)/64 + 1)*64))
    bytes(:length) = TRANSFER(value, bytes(:length))
  TYPE IS (INTEGER(INT8))
    length = 1
    ALLOCATE(bytes(64))
    bytes(1) = value
  TYPE IS (INTEGER(INT16))
    length = 2
    ALLOCATE(bytes(64))
    bytes(1:length) = TRANSFER(value, bytes(1:length))
  TYPE IS (INTEGER(INT32))
    length = 4
    ALLOCATE(bytes(64))
    bytes(1:length) = TRANSFER(value, bytes(1:length))
  TYPE IS (INTEGER(INT64))
    length = 8
    ALLOCATE(bytes(64))
    bytes(1:length) = TRANSFER(value, bytes(1:length))
  TYPE IS (REAL(REAL32))
    length = 4
    ALLOCATE(bytes(64))
    bytes(1:length) = TRANSFER(value, bytes(1:length))
  TYPE IS (REAL(REAL64))
    length = 8
    ALLOCATE(bytes(64))
    bytes(1:length) = TRANSFER(value, bytes(1:length))
  TYPE IS (REAL(REAL128))
    length = 16
    ALLOCATE(bytes(64))
    bytes(1:length) = TRANSFER(value, bytes(1:length))
  CLASS DEFAULT
    CALL Abort(__STAMP__,'Error: Unsupported type in SHA1.')
END SELECT

hash = SHA1Hash(bytes, length)
DEALLOCATE(bytes)

END FUNCTION SHA1Scalar


!==================================================================================================================================
!> SHA-1 function for arrays
!==================================================================================================================================
FUNCTION SHA1Array(value) RESULT(hash)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(*),DIMENSION(:),INTENT(IN)  :: value
CHARACTER(LEN=40)                 :: hash
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(INT8),DIMENSION(:),ALLOCATABLE :: bytes
INTEGER(INT64)                         :: length
INTEGER                                :: i, width
!==================================================================================================================================

select type (value)
  TYPE IS (CHARACTER(LEN=*))
    length = SIZE(value) * LEN(value(1))
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*LEN(value(i)) + 1:i*LEN(value(i))) = &
        TRANSFER(value(i), bytes((i - 1)*LEN(value(i)) + 1:i*LEN(value(i))))
    END DO
  TYPE IS (INTEGER(INT8))
    length = SIZE(value)
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    bytes(1:length) = value
  TYPE IS (INTEGER(INT16))
    width = 2
    length = SIZE(value) * width
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*width + 1:i*width) = TRANSFER(value(i), bytes((i - 1)*width + 1:i*width))
    END DO
  TYPE IS (INTEGER(INT32))
    width = 4
    length = SIZE(value) * width
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*width + 1:i*width) = TRANSFER(value(i), bytes((i - 1)*width + 1:i*width))
    END DO
  TYPE IS (INTEGER(INT64))
    width = 8
    length = SIZE(value) * width
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*width + 1:i*width) = TRANSFER(value(i), bytes((i - 1)*width + 1:i*width))
    END DO
  TYPE IS (REAL(REAL32))
    width = 4
    length = SIZE(value) * width
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*width + 1:i*width) = TRANSFER(value(i), bytes((i - 1)*width + 1:i*width))
    END DO
  TYPE IS (REAL(REAL64))
    width = 8
    length = SIZE(value) * width
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*width + 1:i*width) = TRANSFER(value(i), bytes((i - 1)*width + 1:i*width))
    END DO
  TYPE IS (REAL(REAL128))
    width = 16
    length = SIZE(value) * width
    ALLOCATE(bytes(((length + 8)/64 + 1)*64))
    DO i = 1,SIZE(value)
        bytes((i - 1)*width + 1:i*width) = TRANSFER(value(i), bytes((i - 1)*width + 1:i*width))
    END DO
  CLASS DEFAULT
    CALL Abort(__STAMP__,'Error: Unsupported type in SHA1.')
END SELECT

hash = SHA1Hash(bytes, length)
DEALLOCATE(bytes)

END FUNCTION SHA1Array


!==================================================================================================================================
!> SHA-1 function for strings
!==================================================================================================================================
FUNCTION SHA1Hash(bytes,length)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(INT8),DIMENSION(:),INTENT(INOUT) :: bytes
INTEGER(INT64),            INTENT(IN)    :: length
CHARACTER(LEN=40)                        :: SHA1Hash
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(INT32)               :: h0,h1,h2,h3,h4,a,b,c,d,e,f,k,temp
INTEGER(INT32),DIMENSION(80) :: w
INTEGER                      :: i,j
!==================================================================================================================================

bytes(length + 1) = IBSET(0_int8, 7)
bytes(length + 2:) = 0_int8
bytes(SIZE(bytes) - 7:) = TRANSFER(length*8_int64, bytes(SIZE(bytes) - 7:))
bytes(SIZE(bytes) - 7:) = bytes(SIZE(bytes):SIZE(bytes) - 7:-1)

h0 =  1732584193_int32 ! z'67452301'
h1 = - 271733879_int32 ! z'EFCDAB89'
h2 = -1732584194_int32 ! z'98BADCFE'
h3 =   271733878_int32 ! z'10325476'
h4 = -1009589776_int32 ! z'C3D2E1F0'

DO i = 1,SIZE(bytes)/64
  DO j = 1,16 ! take 512 bit chunk of string
      w(j) = TRANSFER(bytes((i-1)*64 + j*4:(i-1)*64 + (j-1)*4 + 1:-1), w(j)) ! is the source size less than the result size?
  END DO
  DO j = 17,80 ! Extend the sixteen 32-bit words into eighty 32-bit words
      w(j) = ISHFTC(IEOR(IEOR(IEOR(w(j-3), w(j-8)), w(j-14)), w(j-16)), 1)
  END DO

  a = h0; b = h1; c = h2; d = h3; e = h4
  DO j = 1,80
    SELECT CASE (J)
      CASE(1:20)
        f = IOR(IAND(b, c), IAND(NOT(b), d))
        k = 1518500249_int32  ! k = z'5A827999'
      CASE(21:40)
        f = IEOR(IEOR(b, c), d)
        k = 1859775393_int32  ! k = z'6ED9EBA1'
      CASE(41:60)
        f = IOR(IOR(IAND(b, c), IAND(b, d)), IAND(c, d))
        k = -1894007588_int32 ! k = z'8F1BBCDC'
      CASE(61:80)
        f = IEOR(IEOR(b, c), d)
        k = -899497514_int32  ! k = z'CA62C1D6'
      CASE DEFAULT
        ! Elimate GCC warnings
        f = -1
        k = -1
    END SELECT

    temp = ISHFTC(a, 5) + f + e + w(j) + k
    e = d
    d = c
    c = ISHFTC(b, 30)
    b = a
    a = temp
  END DO

  h0 = h0 + a
  h1 = h1 + b
  h2 = h2 + c
  h3 = h3 + d
  h4 = h4 + e
END DO

WRITE(SHA1Hash(1:8),   "(Z8.8)") h0
WRITE(SHA1Hash(9:16),  "(Z8.8)") h1
WRITE(SHA1Hash(17:24), "(Z8.8)") h2
WRITE(SHA1Hash(25:32), "(Z8.8)") h3
WRITE(SHA1Hash(33:40), "(Z8.8)") h4

END FUNCTION SHA1Hash

END MODULE MOD_SHA1
