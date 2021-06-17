! Copyright (c) 2014 Mikael Leetmaa

! This software is provided 'as-is', without any express or implied
! warranty. In no event will the authors be held liable for any damages
! arising from the use of this software.

! Permission is granted to anyone to use this software for any purpose,
! including commercial applications, and to alter it and redistribute it
! freely, subject to the following restrictions:

!    1. The origin of this software must not be misrepresented; you must not
!    claim that you wrote the original software. If you use this software
!    in a product, an acknowledgment in the product documentation would be
!    appreciated but is not required.

!    2. Altered source versions must be plainly marked as such, and must not be
!    misrepresented as being the original software.

!    3. This notice may not be removed or altered from any source
!    distribution.
!==================================================================================================================================
! This code was provided by Mikael Leetmaa at https://github.com/leetmaa/bitsy under the zlib License. A short permissive license,
! compatible with GPL. Requires altered source versions to be documented as such.
!>> The version provided here was adjusted to FLEXI style and contains minor fixes for conformal variable casting.

MODULE MOD_SHA256
! MODULES
USE ISO_C_BINDING

IMPLICIT NONE
PRIVATE

INTERFACE SHA256
  MODULE PROCEDURE SHA256
END INTERFACE

PUBLIC :: SHA256
! PUBLIC :: DIRTY_SHA256
! Public for the sake of unit-testing.
! PUBLIC :: SHA256B
! PUBLIC :: MS0
! PUBLIC :: MS1
! PUBLIC :: CS0
! PUBLIC :: CS1
! PUBLIC :: MAJ
! PUBLIC :: CH
! PUBLIC :: SWAP32
! PUBLIC :: SWAP64
! PUBLIC :: SWAP64A
! PUBLIC :: CONSUME_CHUNK
!----------------------------------------------------------------------------------------------------------------------------------

CONTAINS

!==================================================================================================================================
!> SHA-256 interface function.
!! @param str : (in) The message to digest.
!! @return    : The SHA-256 digest as a string of length 64.
!==================================================================================================================================
FUNCTION SHA256(STR)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=64) :: SHA256
CHARACTER(LEN=*), INTENT(IN) :: STR
!==================================================================================================================================
! Call the work horse with proper bit swapping.
SHA256 = SHA256B(STR, 1)

END FUNCTION SHA256


!==================================================================================================================================
!> Quick and dirty SHA-256 interface function (no bit-swapping).
!! @param str : (in) The message to digest.
!! @return    : The SHA-256 digest as a string of length 64.
!==================================================================================================================================
FUNCTION DIRTY_SHA256(STR)
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=64) :: DIRTY_SHA256
CHARACTER(LEN=*), INTENT(IN):: STR
!==================================================================================================================================
! Call the work horse - no bit swapping.
DIRTY_SHA256 = SHA256B(STR, 0)

END FUNCTION DIRTY_SHA256


!==================================================================================================================================
!> Calculate the SHA-256 hash of the incomming string.
!! @param str    : (in) The message to take digest.
!! @param swap   : (in) Flag to indicate if swapping to big-endian
!!                 input (swap=1) should be used. swap=1 is needed
!!                 for the routine to pass the standard tests, but
!!                 decreases speed with a factor 2.
!! @return       : The SHA-256 digest as a string of length 64.
!==================================================================================================================================
FUNCTION SHA256B(STR, SWAP)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=64) :: SHA256B
CHARACTER(LEN=*), INTENT(IN) :: str
INTEGER,          INTENT(IN) :: swap
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Helper variables.
INTEGER(KIND=C_INT64_T) :: length
INTEGER(KIND=C_INT32_T) :: i,temp1,temp2
INTEGER                 :: break,pos0
! Parameters for the cruncher.
INTEGER(KIND=C_INT32_T) :: h0_ref(8)
INTEGER(KIND=C_INT32_T) :: k0_ref(64)
! Work areas.
INTEGER(KIND=C_INT32_T) :: h0(8)
INTEGER(KIND=C_INT32_T) :: k0(64)
INTEGER(KIND=C_INT32_T) :: a0(8)
INTEGER(KIND=C_INT32_T) :: w0(64)
!==================================================================================================================================

! Set the initial data.
DATA (h0_ref(i),i=1,8)/&
     & z'6a09e667', z'bb67ae85', z'3c6ef372', z'a54ff53a', z'510e527f', z'9b05688c', z'1f83d9ab', z'5be0cd19'/

DATA (k0_ref(i), i=1,64)/&
     & z'428a2f98', z'71374491', z'b5c0fbcf', z'e9b5dba5', z'3956c25b', z'59f111f1', z'923f82a4', z'ab1c5ed5',&
     & z'd807aa98', z'12835b01', z'243185be', z'550c7dc3', z'72be5d74', z'80deb1fe', z'9bdc06a7', z'c19bf174',&
     & z'e49b69c1', z'efbe4786', z'0fc19dc6', z'240ca1cc', z'2de92c6f', z'4a7484aa', z'5cb0a9dc', z'76f988da',&
     & z'983e5152', z'a831c66d', z'b00327c8', z'bf597fc7', z'c6e00bf3', z'd5a79147', z'06ca6351', z'14292967',&
     & z'27b70a85', z'2e1b2138', z'4d2c6dfc', z'53380d13', z'650a7354', z'766a0abb', z'81c2c92e', z'92722c85',&
     & z'a2bfe8a1', z'a81a664b', z'c24b8b70', z'c76c51a3', z'd192e819', z'd6990624', z'f40e3585', z'106aa070',&
     & z'19a4c116', z'1e376c08', z'2748774c', z'34b0bcb5', z'391c0cb3', z'4ed8aa4a', z'5b9cca4f', z'682e6ff3',&
     & z'748f82ee', z'78a5636f', z'84c87814', z'8cc70208', z'90befffa', z'a4506ceb', z'bef9a3f7', z'c67178f2'/

h0 = h0_ref
k0 = k0_ref

! -----------------------------------
! Function body implementation.
break  = 0
pos0   = 1
length = LEN(TRIM(str))

DO WHILE (break.NE.1)
  ! Get the next 16 32bit words to consume.
  CALL consume_chunk(str, length, w0(1:16), pos0, break, swap)

  ! Extend the first 16 words to fill the work schedule array.
  DO i = 17,64
    w0(i) = ms1(w0(i-2)) + w0(i-16) + ms0(w0(i-15)) + w0(i-7)
  END DO

  ! Initialize the workin variables with the current version of the hash.
  a0 = h0

  ! Run the compression loop.
  DO i = 1,64
    temp1 = a0(8) + cs1(a0(5)) + ch(a0(5),a0(6),a0(7)) + k0(i) + w0(i)
    temp2 = cs0(a0(1)) + maj(a0(1),a0(2),a0(3))

    a0(8) = a0(7)
    a0(7) = a0(6)
    a0(6) = a0(5)
    a0(5) = a0(4) + temp1
    a0(4) = a0(3)
    a0(3) = a0(2)
    a0(2) = a0(1)
    a0(1) = temp1 + temp2
  END DO

  ! Update the state.
  h0 = h0 + a0
END DO

! Write the result to the output variable.
WRITE(SHA256B,'(8z8.8)') h0(1), h0(2), h0(3), h0(4), h0(5), h0(6), h0(7), h0(8)

END FUNCTION SHA256B


!==================================================================================================================================
!> Swap the byte order on a 32bit integer.
!! @param inp : (in) The integer to byte swap.
!! @return    : The byte swapped integer.
!==================================================================================================================================
PURE FUNCTION SWAP32(inp)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: swap32
INTEGER(KIND=C_INT32_T),INTENT(IN) :: inp
!==================================================================================================================================
! Initialize to zero to eliminate GCC warning
swap32 = 0

CALL MVBITS(inp, 24, 8, swap32,  0)
CALL MVBITS(inp, 16, 8, swap32,  8)
CALL MVBITS(inp,  8, 8, swap32, 16)
CALL MVBITS(inp,  0, 8, swap32, 24)

END FUNCTION SWAP32


!==================================================================================================================================
!> Swap the byte order on a 64 bit integer.
!! @param inp : (in) The integer to byte swap.
!! @return    : The byte swapped integer.
!==================================================================================================================================
PURE FUNCTION SWAP64(inp)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT64_T)            :: swap64
INTEGER(KIND=C_INT64_T),INTENT(IN) :: inp
!==================================================================================================================================
CALL MVBITS(inp, 56, 8, swap64,  0)
CALL MVBITS(inp, 48, 8, swap64,  8)
CALL MVBITS(inp, 40, 8, swap64, 16)
CALL MVBITS(inp, 32, 8, swap64, 24)
CALL MVBITS(inp, 24, 8, swap64, 32)
CALL MVBITS(inp, 16, 8, swap64, 40)
CALL MVBITS(inp,  8, 8, swap64, 48)
CALL MVBITS(inp,  0, 8, swap64, 56)

END FUNCTION SWAP64

!==================================================================================================================================
!> Swap the byte order on a 64bit integer as if
!! each half was a 32bit integer to swap.
!! @param inp : (in) The integer to byte swap.
!! @return    : The byte swapped integer.
!==================================================================================================================================
PURE FUNCTION SWAP64A(inp)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT64_T)            :: swap64a
INTEGER(KIND=C_INT64_T),INTENT(IN) :: inp
!==================================================================================================================================
! Initialize to zero to eliminate GCC warning
swap64a = 0

CALL MVBITS(inp,  0, 8, swap64a, 32)
CALL MVBITS(inp,  8, 8, swap64a, 40)
CALL MVBITS(inp, 16, 8, swap64a, 48)
CALL MVBITS(inp, 24, 8, swap64a, 56)
CALL MVBITS(inp, 32, 8, swap64a,  0)
CALL MVBITS(inp, 40, 8, swap64a,  8)
CALL MVBITS(inp, 48, 8, swap64a, 16)
CALL MVBITS(inp, 56, 8, swap64a, 24)

END FUNCTION SWAP64A


!==================================================================================================================================
!> The 'ch' function in SHA-2.
!! @param a : (in) The a input integer.
!! @param b : (in) The b input integer.
!! @param c : (in) The c input integer.
!! @return  : ch(a,b,c), see the code.
!==================================================================================================================================
PURE FUNCTION CH(a, b, c)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: ch
INTEGER(KIND=C_INT32_T),INTENT(IN) :: a
INTEGER(KIND=C_INT32_T),INTENT(IN) :: b
INTEGER(KIND=C_INT32_T),INTENT(IN) :: c
!==================================================================================================================================
ch = IEOR(IAND(a, b), (IAND(NOT(a), c)))

END FUNCTION CH


!==================================================================================================================================
!> The 'maj' function in SHA-2.
!! @param a : (in) The a input integer.
!! @param b : (in) The b input integer.
!! @param c : (in) The c input integer.
!! @return  : maj(a,b,c), see the code.
!==================================================================================================================================
PURE FUNCTION MAJ(a, b, c)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: maj
INTEGER(KIND=C_INT32_T),INTENT(IN) :: a
INTEGER(KIND=C_INT32_T),INTENT(IN) :: b
INTEGER(KIND=C_INT32_T),INTENT(IN) :: c
!==================================================================================================================================
maj = IEOR(IAND(a, b), IEOR(IAND(a, c), IAND(b, c)))

END FUNCTION MAJ

!==================================================================================================================================
!> The '\Sigma_0' function in SHA-2.
!! @param a : (in) The a input integer.
!! @return  : cs0(a), see the code.
!==================================================================================================================================
PURE FUNCTION CS0(a)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: cs0
INTEGER(KIND=C_INT32_T),INTENT(IN) :: a
!==================================================================================================================================
cs0 = IEOR(ISHFTC(a, -2), IEOR(ISHFTC(a, -13), ISHFTC(a, -22)))

END FUNCTION CS0


!==================================================================================================================================
!> The '\Sigma_1' function in SHA-2.
!! @param a : (in) The a input integer.
!! @return  : cs1(a), see the code.
!==================================================================================================================================
PURE FUNCTION CS1(a)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: cs1
INTEGER(KIND=C_INT32_T),INTENT(IN) :: a
!==================================================================================================================================
cs1 = IEOR(ISHFTC(a, -6), IEOR(ISHFTC(a, -11), ISHFTC(a, -25)))

END FUNCTION CS1


!==================================================================================================================================
!> The '\sigma_0' function in SHA-2.
!! @param a : (in) The a input integer.
!! @return  : ms0(a), see the code.
!==================================================================================================================================
PURE FUNCTION MS0(A)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: ms0
INTEGER(KIND=C_INT32_T),INTENT(IN) :: a
!==================================================================================================================================
ms0 = IEOR(ISHFTC(a, -7), IEOR(ISHFTC(a, -18), ISHFT(a, -3)))

END FUNCTION MS0


!==================================================================================================================================
!> The '\sigma_1' function in SHA-2.
!! @param a : (in) The a input integer.
!! @return  : ms1(a), see the code.
!==================================================================================================================================
PURE FUNCTION MS1(a)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=C_INT32_T)            :: ms1
INTEGER(KIND=C_INT32_T),INTENT(IN) :: a
!==================================================================================================================================
ms1 = IEOR(ISHFTC(a, -17), IEOR(ISHFTC(a, -19), ISHFT(a, -10)))

END FUNCTION MS1


!==================================================================================================================================
!> Copy 16 32bit words of data from str(pos0) to inp(1:16). The
!! data is padded a requiered by the SHA-256 algorithm.
!! @param str    : (in) The message to take a chunk from.
!! @param length : (in) The length of the message in 8bit words.
!! @param inp    : (inout) The work area to copy the data to.
!! @param pos0   : (inout) Variable to store the start of the next chunk.
!! @param break  : (inout) Indicates the position in the work flow.
!!                 break=0 on entry -> continue to consume a chunk, pad if needed.
!!                 break=2 on entry -> continue to consume, padding was allready done.
!!                 break=1 one exit -> the last chunk was consumed.
!! @param swap   : (in) Flag to indicate if swapping to big-endian
!!                 input (swap=1) should be used. swap=1 is needed
!!                 for the routine to pass the standard tests, but
!!                 decreases speed with a factor 2.
!==================================================================================================================================
SUBROUTINE consume_chunk(str,length,inp,pos0,break,swap)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),       INTENT(IN)    :: str
INTEGER(KIND=C_INT64_T),INTENT(IN)    :: length
INTEGER(KIND=C_INT32_T),INTENT(INOUT) :: inp(*)
INTEGER,                INTENT(INOUT) :: pos0
INTEGER,                INTENT(INOUT) :: break
INTEGER,                INTENT(IN)    :: swap
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=4)        :: last_word
INTEGER(KIND=C_INT64_T) :: rest
INTEGER(KIND=C_INT32_T) :: to_pad
INTEGER(KIND=C_INT32_T) :: leftover
INTEGER(KIND=C_INT32_T) :: space_left
INTEGER(KIND=C_INT32_T) :: zero
INTEGER(KIND=C_INT8_T)  :: ipad0
INTEGER(KIND=C_INT8_T)  :: ipad1
INTEGER(KIND=C_INT8_T)  :: i
! Initialize data
DATA zero  / b'00000000000000000000000000000000'/
DATA ipad0 / b'00000000' /
DATA ipad1 / b'10000000' /
!==================================================================================================================================

! Calculate the rest.
rest = length - pos0 + 1

! If we are far from the end.
IF (rest.GE.64) THEN
  ! Copy the data over.
  inp(1:16) = TRANSFER(str(pos0:pos0+64-1), inp(1:16))

  ! Big-endian.
  IF (swap.EQ.1) THEN
    DO i = 1,16
      inp(i) = swap32(inp(i))
    END DO
  END IF

  ! Increment the starting position for the next roundx.
  pos0 = pos0 + 64

ELSE
  ! Space left in the input chunk.
  space_left = 16

  ! number of leftover full 32bit words.
  leftover   = INT(rest/4,C_INT32_T)

  ! Copy any leftovers.
  IF (leftover.GT.0) THEN
    inp(1:leftover) = transfer(str(pos0:pos0+leftover*4-1), inp(1:16))

    ! Big-endian.
    IF (swap.EQ.1) THEN
      DO i = 1,INT(leftover,C_INT8_T)
        inp(i) = swap32(inp(i))
      END DO
    END IF

    ! Increment the starting position.
    pos0 = pos0 + leftover*4
    rest = length - pos0 + 1
    space_left = space_left - leftover
  END IF

  IF (space_left.GT.0) THEN

    IF (break.NE.2) THEN
      ! Add any remaining incomplete 32bit word.
      IF (rest.GT.0) THEN
        last_word(1:rest) = str(pos0:pos0+rest-1)
        ! Increment the pos0.
        pos0 = pos0 + INT(rest)
      END IF

      ! Add the '10000000' padding.
      last_word(rest+1:rest+1) = TRANSFER(ipad1, last_word(1:1))

      ! Add zeros for a full 32bit word.
      to_pad = 4 - INT(rest,C_INT32_T) - 1
      DO i = 1,INT(to_pad,C_INT8_T)
        last_word(rest+1+i:rest+1+i) = TRANSFER(ipad0, last_word(1:1))
      END DO

      ! Copy the last full (padded) word over.
      inp(17-space_left) = TRANSFER(last_word(1:4), inp(1))

      IF (swap.EQ.1) THEN
        inp(17-space_left) = swap32(inp(17-space_left))
      END IF

      ! Decrement the space left.
      space_left = space_left - 1

      ! Set the flag to indicate that we have padded.
      break = 2
    END IF

    ! If not enough space to finnish, add zeros.
    IF (space_left.EQ.1) THEN
      inp(16) = zero
      space_left = 0
    END IF

    rest = 0
  END IF

  ! Continue with the last part if there is enough space left.
  IF ((rest.EQ.0) .AND. (space_left.GE.2)) THEN

    ! Add zeros until 64 bits left.
    DO WHILE (space_left.GT.2)
      inp(17-space_left) = zero
      space_left = space_left - 1
    END DO

    ! Add the two last 32bit words.
    inp(15:16) = TRANSFER(swap64a(length*8), inp(15:16))

    ! Set break flag indicating we are done with the whole message.
    break = 1
  END IF
END IF

END SUBROUTINE consume_chunk

END MODULE MOD_SHA256
