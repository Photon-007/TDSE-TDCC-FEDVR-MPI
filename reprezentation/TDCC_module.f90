MODULE TDCC_module
	USE nrtype
	USE gvars

	TYPE TDCC_coeffs
		INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: couple2
		REAL(DP),     DIMENSION(:), ALLOCATABLE :: value
	END TYPE


	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE TDCC_set_tagmax(representation)
	IMPLICIT NONE
	CHARACTER(LEN=4), INTENT(IN) :: representation

	select case (representation)
!	select case (SphericalHarmonicRepresentation)
		case ('Y_lm')
			tagmax = (lmax + 1) * (lmax + 1)

		case ('Y_l0')
			tagmax = lmax + 1
	end select

END SUBROUTINE
!-------------------------------------------------------------------------------
FUNCTION TDCC_get_tag_from_lm(represetation, ll, mm) RESULT(tag)
	IMPLICIT NONE
	CHARACTER(LEN=4), INTENT(IN) :: represetation
	INTEGER(I4B), INTENT(IN)     :: ll, mm
	INTEGER(I4B)                 :: tag

	if((ll<0).or.(ll>lmax).or.(mm<-lmax).or.(mm>lmax))then	! The quantum numers lie outside the studied range of values,
		tag = -1											! so we do not have a tag associated to them
	else
		select case (represetation)
			case ('Y_lm')
				tag = ll**2 + (ll+mm) + 1

			case ('Y_l0')
				tag = ll + 1
		end select
	endif

END FUNCTION
!-------------------------------------------------------------------------------
SUBROUTINE TDCC_get_lm_from_tag(represetation, tag, ll, mm)
	IMPLICIT NONE
	CHARACTER(LEN=4), INTENT(IN) :: represetation
	INTEGER(I4B), INTENT(IN)     :: tag
	INTEGER(I4B), INTENT(OUT)    :: ll, mm

	select case (represetation)
		case ('Y_lm')
			ll = int(dsqrt(1.0d0*(tag - 1)))
			mm = (tag - 1) - ll - ll**2 

		case ('Y_l0')
			ll = tag - 1
			mm = m_init
	end select

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE TDCC_test_mapping(representation, ll, mm)
	IMPLICIT NONE
	CHARACTER(LEN=4), INTENT(IN) :: representation
	INTEGER(I4B), INTENT(IN)     :: ll, mm
	INTEGER(I4B) :: l, m
	INTEGER(I4B) :: tag

	tag = TDCC_get_tag_from_lm(representation, ll, mm)
	write(6,'(I5.1,2x,I5.1,2x,A6,2x,I5.1)') ll, mm, " ---> ", tag
	call TDCC_get_lm_from_tag(representation, tag, l, m)
	write(6,'(I5.1,2x,A6,2x,I5.1,2x,I5.1)') tag, " ---> ", l, m

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE TDCC_set_coupling_coeffs(CC, representation, tag_i, tag_f)
	IMPLICIT NONE
	TYPE(TDCC_coeffs), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: CC
	CHARACTER(LEN=4), INTENT(IN) :: representation
	INTEGER(I4B),     INTENT(IN) :: tag_i, tag_f
	INTEGER(I4B) :: tag, ll, mm

!	if(.not.allocated(CC))allocate(CC(tag_i:tag_f))

	select case (representation)
		case ('Y_lm')
			do tag=tag_i,tag_f
!				if(allocated(CC(tag)%couple2))deallocate(CC(tag)%couple2)
!				if(allocated(CC(tag)%value))  deallocate(CC(tag)%value)
				allocate(CC(tag)%couple2(1:4))	! {l-1, m-1}, {l-1, m+1}, {l+1, m-1}, {l+1, m+1}
				allocate(CC(tag)%value(1:4))

				call TDCC_get_lm_from_tag(representation, tag, ll, mm)

				CC(tag)%couple2(1) = TDCC_get_tag_from_lm(representation, ll-1, mm-1)
				CC(tag)%couple2(2) = TDCC_get_tag_from_lm(representation, ll-1, mm+1)
				CC(tag)%couple2(3) = TDCC_get_tag_from_lm(representation, ll+1, mm-1)
				CC(tag)%couple2(4) = TDCC_get_tag_from_lm(representation, ll+1, mm+1)
				CC(tag)%value(1)   = -0.5 * dsqrt((1.0d0*(ll+mm)*(ll+mm-1))/((2*ll-1)*(2*ll+1)))
				CC(tag)%value(2)   = -0.5 * dsqrt((1.0d0*(ll-mm)*(ll-mm-1))/((2*ll-1)*(2*ll+1)))
				CC(tag)%value(3)   =  0.5 * dsqrt((1.0d0*(ll-mm+1)*(ll-mm+2))/((2*ll+1)*(2*ll+3)))
				CC(tag)%value(4)   =  0.5 * dsqrt((1.0d0*(ll+mm+1)*(ll+mm+2))/((2*ll+1)*(2*ll+3)))
			enddo

		case ('Y_l0')
			do tag=tag_i,tag_f
!				if(allocated(CC(tag)%couple2))deallocate(CC(tag)%couple2)
!				if(allocated(CC(tag)%value))  deallocate(CC(tag)%value)
				allocate(CC(tag)%couple2(1:2))	! l-1, l+1
				allocate(CC(tag)%value(1:2))

				call TDCC_get_lm_from_tag(representation, tag, ll, mm)

				CC(tag)%couple2(1) = TDCC_get_tag_from_lm(representation, ll-1, mm)
				CC(tag)%couple2(2) = TDCC_get_tag_from_lm(representation, ll+1, mm)
				CC(tag)%value(1)   = dsqrt((1.0d0*(ll-mm)*(ll+mm))/((2*ll-1)*(2*ll+1)))
				CC(tag)%value(2)   = dsqrt((1.0d0*(ll-mm+1)*(ll+mm+1))/((2*ll+1)*(2*ll+3)))
			enddo
	end select

END SUBROUTINE
!-------------------------------------------------------------------------------




END MODULE TDCC_module
