MODULE Operator_typedef
	USE nrtype
	USE gvars
	USE FEDVR_module
	USE TDCC_module
	USE HDF5

	INTERFACE Op_create
		MODULE PROCEDURE fedvr_op_full_create
		MODULE PROCEDURE fedvr_op_diag_create
		MODULE PROCEDURE EE_op_create
		MODULE PROCEDURE AA_op_create
	END INTERFACE

	INTERFACE Op_allocate
		MODULE PROCEDURE fedvr_op_full_allocate
		MODULE PROCEDURE fedvr_op_diag_allocate
		MODULE PROCEDURE EE_op_allocate
		MODULE PROCEDURE AA_op_allocate
	END INTERFACE

	INTERFACE Op_mirror
		MODULE PROCEDURE fedvr_op_full_mirror
		MODULE PROCEDURE fedvr_op_diag_mirror
		MODULE PROCEDURE EE_op_mirror
		MODULE PROCEDURE AA_op_mirror
	END INTERFACE

	INTERFACE Op_set_bigmat
		MODULE PROCEDURE set_bigmat_OP_full
		MODULE PROCEDURE set_bigmat_OP_diag
	END INTERFACE

	INTERFACE OP_tofile
		MODULE PROCEDURE fedvr_op_tofile_full
		MODULE PROCEDURE fedvr_op_tofile_diag
	END INTERFACE

	INTERFACE OP_save
		MODULE PROCEDURE fedvr_op_full_save_h5
		MODULE PROCEDURE fedvr_op_diag_save_h5
	END INTERFACE

	TYPE fedvr_matrix
		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matrix
	END TYPE

	TYPE fedvr_op_full
		INTEGER(I4B), DIMENSION(1:2) :: fe_range
		INTEGER(I4B), DIMENSION(1:2) :: gp_range
		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
		INTEGER(I4B) :: Nfun, Nfe, Ngp, ll, mm
		TYPE(fedvr_matrix), DIMENSION(:), ALLOCATABLE :: fe_region
	END TYPE

	TYPE fedvr_op_diag
		INTEGER(I4B), DIMENSION(1:2) :: fe_range
		INTEGER(I4B), DIMENSION(1:2) :: gp_range
		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
		INTEGER(I4B) :: Nfun, Nfe, Ngp, ll, mm
		REAL(DP), DIMENSION(:), ALLOCATABLE :: diag
	END TYPE

	TYPE EE_op
		INTEGER(I4B), DIMENSION(1:2) :: fe_range
		INTEGER(I4B), DIMENSION(1:2) :: gp_range
		INTEGER(I4B), DIMENSION(1:2) :: tag_range
		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
		INTEGER(I4B) :: Nfun, Nfe, Ngp, Ntag, Npulse
		TYPE(pulse), DIMENSION(:), ALLOCATABLE :: pulses		! in case multiple pulses are switched on... have to figure out how to read them from input...
		REAL(DP), DIMENSION(:), ALLOCATABLE :: Ei_t
		REAL(DP), DIMENSION(:), ALLOCATABLE :: Envi_t
		REAL(DP) :: EE_t
		REAL(DP) :: Env_t
		CHARACTER(LEN=4) :: Y_rep
		INTEGER(I4B) :: Nrank
		TYPE(TDCC_coeffs), DIMENSION(:), ALLOCATABLE :: cc
	END TYPE

	TYPE AA_op
		INTEGER(I4B), DIMENSION(1:2) :: fe_range
		INTEGER(I4B), DIMENSION(1:2) :: gp_range
		INTEGER(I4B), DIMENSION(1:2) :: tag_range
		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
		INTEGER(I4B) :: Nfun, Nfe, Ngp, Ntag, Npulse
		TYPE(pulse), DIMENSION(:), ALLOCATABLE :: pulses
		REAL(DP), DIMENSION(:), ALLOCATABLE :: Ai_t
		REAL(DP), DIMENSION(:), ALLOCATABLE :: Envi_t
		REAL(DP) :: AA_t
		REAL(DP) :: Env_t
		CHARACTER(LEN=4) :: Y_rep
		INTEGER(I4B) :: Nrank
		TYPE(TDCC_coeffs), DIMENSION(:), ALLOCATABLE :: cc
		TYPE(fedvr_op_full)	:: deriv1
	END TYPE

!	TYPE OP_BigMatrix
!		INTEGER(I4B) :: Ngp
!		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matrix
!	END TYPE

!	GLOBAL OPERATORS DEFINED IN MODULE SETUP
	TYPE(fedvr_op_full) :: OP_Sym
	TYPE(fedvr_op_full) :: OP_KinVr,  Ham0_tag,  Ham0_init
	TYPE(fedvr_op_diag) :: OP_Vr,     OP_dVr									! used during hhg
	TYPE(fedvr_op_diag) :: OP_Vl_tag, OP_dVl_tag
	TYPE(fedvr_op_diag) :: Abspot_op, Maskfun_op
	TYPE(fedvr_op_full), DIMENSION(:), ALLOCATABLE :: OP_Ham0
	TYPE(fedvr_op_diag), DIMENSION(:), ALLOCATABLE :: OP_Vl, OP_dVl

	TYPE(EE_op) :: OP_EField
	TYPE(AA_OP) :: OP_AField

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE Op_Sym_create(Sym, grid)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(OUT) :: Sym
	TYPE(fedvr_grid), INTENT(IN)  :: grid
	INTEGER(I4B) :: ii, jj, ii_gp, jj_gp
	REAL(DP) :: wi, wj

	call Op_create(Sym, grid, 0, 0)
	call Op_allocate(Sym)

!	do ii_fe=1,OP_Sym%Nfe
!		do jj=1,OP_Sym%Nfun
!			do ii=1,OP_Sym%Nfun
!				OP_Sym%fe_region(ii_fe)%matrix(ii,jj) = grid%ww_fe(ii,ii_fe)/dsqrt(grid%ww_fe(ii,ii_fe)*grid%ww_fe(jj,ii_fe))
!			enddo
!		enddo
!	enddo
!	...DAFUCK!!!!!! What is the difference??????? sqrt(2) factor in the first/last row/col ---> at the FE boundary the grid%ww = sum(grid%ww_fe)
	do ii_fe=1,Sym%Nfe
		do jj=1,Sym%Nfun
			jj_gp = (ii_fe-1)*(Sym%Nfun-1) + jj
			do ii=1,Sym%Nfun
				ii_gp = (ii_fe-1)*(Sym%Nfun-1) + ii
				wi = grid%ww(ii_gp)
				wj = grid%ww(jj_gp)
				Sym%fe_region(ii_fe)%matrix(ii,jj) = wi/dsqrt(wi*wj)
			enddo
		enddo
	enddo

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE fedvr_op_full_create(OP, grid, ll, mm)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(OUT) :: OP
	TYPE(fedvr_grid), INTENT(IN)  :: grid
	INTEGER(I4B), INTENT(IN) :: ll, mm

	OP%fe_range = grid%fe_range
	OP%gp_range = grid%gp_range
	OP%Nfun     = grid%Nfun
	OP%Nfe      = grid%Nfe
	OP%Ngp      = grid%Ngp
	OP%ll       = ll
	OP%mm       = mm

!	if(allocated(OP%fe_region))deallocate(OP%fe_region)		! does this deallocate all structures nested inside???

	if(allocated(OP%fe_gp_range))deallocate(OP%fe_gp_range)
	allocate(OP%fe_gp_range(1:2, 1:grid%Nfe))
	OP%fe_gp_range(:,:) = grid%fe_gp_range(:,:)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_diag_create(OP, grid, ll, mm)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(OUT) :: OP
	TYPE(fedvr_grid), INTENT(IN)  :: grid
	INTEGER(I4B), INTENT(IN) :: ll, mm

	OP%fe_range = grid%fe_range
	OP%gp_range = grid%gp_range
	OP%Nfun     = grid%Nfun
	OP%Nfe      = grid%Nfe
	OP%Ngp      = grid%Ngp
	OP%ll       = ll
	OP%mm       = mm

!	if(allocated(OP%diag))deallocate(OP%diag)

	if(allocated(OP%fe_gp_range))deallocate(OP%fe_gp_range)
	allocate(OP%fe_gp_range(1:2, 1:grid%Nfe))
!	write(6,*)shape(OP%fe_gp_range), shape(grid%fe_gp_range)
	OP%fe_gp_range(:,:) = grid%fe_gp_range(:,:)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE EE_op_create(OP, grid, tag_range, Y_rep, LP)
	IMPLICIT NONE
	TYPE(EE_op), INTENT(OUT) :: OP
	TYPE(fedvr_grid), INTENT(IN)  :: grid
	INTEGER(I4B), DIMENSION(2), INTENT(IN) :: tag_range
	CHARACTER(LEN=4), INTENT(IN) :: Y_rep
	TYPE(pulse), DIMENSION(:), INTENT(IN) :: LP
	INTEGER(I4B) :: ii

	OP%Y_rep    = Y_rep
	OP%fe_range = grid%fe_range
	OP%gp_range = grid%gp_range
	OP%tag_range= tag_range
	OP%Ntag     = tag_range(2) - tag_range(1) + 1
	OP%Nfun     = grid%Nfun
	OP%Nfe      = grid%Nfe
	OP%Ngp      = grid%Ngp
	OP%Npulse   = size(LP)
	OP%EE_t     = 0.0d0
	OP%Env_t    = 0.0d0
	select case (Y_rep)
		case ('Y_lm')
			OP%Nrank = 4
		case ('Y_l0')
			OP%Nrank = 2
	end select

	if(allocated(OP%pulses))deallocate(OP%pulses)
	allocate(OP%pulses(OP%Npulse))

	do ii=1,OP%Npulse
		OP%pulses(ii)%Envelope = LP(ii)%Envelope
		OP%pulses(ii)%EA       = LP(ii)%EA
		OP%pulses(ii)%t0       = LP(ii)%t0
		OP%pulses(ii)%Omega    = LP(ii)%Omega
		OP%pulses(ii)%param    = LP(ii)%param
		OP%pulses(ii)%CEP      = LP(ii)%CEP
		OP%pulses(ii)%t_start  = LP(ii)%t_start
		OP%pulses(ii)%t_end    = LP(ii)%t_end
	enddo

	if(allocated(OP%Ei_t))  deallocate(OP%Ei_t)
	if(allocated(OP%Envi_t))deallocate(OP%Envi_t)
	allocate(OP%Ei_t(OP%Npulse), OP%Envi_t(OP%Npulse))
	OP%Ei_t(:)   = 0.0d0
	OP%Envi_t(:) = 0.0d0

	if(allocated(OP%fe_gp_range))deallocate(OP%fe_gp_range)
	allocate(OP%fe_gp_range(1:2, 1:grid%Nfe))
	OP%fe_gp_range(:,:) = grid%fe_gp_range(:,:)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE AA_op_create(OP, grid, tag_range, Y_rep, LP)
	IMPLICIT NONE
	TYPE(AA_Op), INTENT(OUT) :: OP
	TYPE(fedvr_grid), INTENT(IN)  :: grid
	INTEGER(I4B), DIMENSION(2), INTENT(IN) :: tag_range
	CHARACTER(LEN=4), INTENT(IN) :: Y_rep
	TYPE(pulse), DIMENSION(:), INTENT(IN) :: LP
	INTEGER(I4B) :: ii

	OP%Y_rep    = Y_rep
	OP%fe_range = grid%fe_range
	OP%gp_range = grid%gp_range
	OP%tag_range= tag_range
	OP%Ntag     = tag_range(2) - tag_range(1) + 1
	OP%Nfun     = grid%Nfun
	OP%Nfe      = grid%Nfe
	OP%Ngp      = grid%Ngp
	OP%Npulse   = size(LP)
	OP%AA_t     = 0.0d0
	OP%Env_t    = 0.0d0
	select case (Y_rep)
		case ('Y_lm')
			OP%Nrank = 4
		case ('Y_l0')
			OP%Nrank = 2
	end select
	call fedvr_op_full_create(OP%deriv1, grid, 0, 0)

	if(allocated(OP%pulses))deallocate(OP%pulses)
	allocate(OP%pulses(OP%Npulse))

	do ii = 1,OP%Npulse
		OP%pulses(ii)%Envelope = LP(ii)%Envelope
		OP%pulses(ii)%EA       = LP(ii)%EA
		OP%pulses(ii)%t0       = LP(ii)%t0
		OP%pulses(ii)%Omega    = LP(ii)%Omega
		OP%pulses(ii)%param    = LP(ii)%param
		OP%pulses(ii)%CEP      = LP(ii)%CEP
		OP%pulses(ii)%t_start  = LP(ii)%t_start
		OP%pulses(ii)%t_end    = LP(ii)%t_end
	enddo

	if(allocated(OP%Ai_t))  deallocate(OP%Ai_t)
	if(allocated(OP%Envi_t))deallocate(OP%Envi_t)
	allocate(OP%Ai_t(OP%Npulse), OP%Envi_t(OP%Npulse))
	OP%Ai_t(:)   = 0.0d0
	OP%Envi_t(:) = 0.0d0

	if(allocated(OP%fe_gp_range))deallocate(OP%fe_gp_range)
	allocate(OP%fe_gp_range(1:2, 1:grid%Nfe))
	OP%fe_gp_range(:,:) = grid%fe_gp_range(:,:)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_full_allocate(OP)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(INOUT) :: OP

	if(allocated(OP%fe_region))deallocate(OP%fe_region)
	allocate(OP%fe_region(1:OP%Nfe))
	do ii_fe=1,OP%Nfe
		if(allocated(OP%fe_region(ii_fe)%matrix))deallocate(OP%fe_region(ii_fe)%matrix)
		allocate(OP%fe_region(ii_fe)%matrix(1:OP%Nfun,1:OP%Nfun))
		OP%fe_region(ii_fe)%matrix(:,:) = 0.0d0
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_diag_allocate(OP)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(INOUT) :: OP

	if(allocated(OP%diag))deallocate(OP%diag)
	allocate(OP%diag(1:OP%Ngp))
	OP%diag(:) = 0.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE EE_op_allocate(OP)
	IMPLICIT NONE
	TYPE(EE_op), INTENT(INOUT) :: OP

	if(allocated(OP%cc))deallocate(OP%cc)
	allocate(OP%cc(OP%tag_range(1):OP%tag_range(2)))
	call TDCC_set_coupling_coeffs(OP%cc, OP%Y_rep, OP%tag_range(1), OP%tag_range(2))

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE AA_op_allocate(OP)
	IMPLICIT NONE
	TYPE(AA_Op), INTENT(INOUT) :: OP

	if(allocated(OP%cc))deallocate(OP%cc)
	allocate(OP%cc(OP%tag_range(1):OP%tag_range(2)))
	call TDCC_set_coupling_coeffs(OP%cc, OP%Y_rep, OP%tag_range(1), OP%tag_range(2))
	call fedvr_op_full_allocate(OP%deriv1)
!	call separate routine for first deriv

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_full_mirror(OP, OP_new, behaviour)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN)  :: OP
	TYPE(fedvr_op_full), INTENT(OUT) :: OP_new
	CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: behaviour

	OP_new%fe_range = OP%fe_range
	OP_new%gp_range = OP%gp_range
	OP_new%Nfun     = OP%Nfun
	OP_new%Nfe      = OP%Nfe
	OP_new%Ngp      = OP%Ngp
	OP_new%ll       = OP%ll
	OP_new%mm       = OP%mm

	if(allocated(OP_new%fe_gp_range))deallocate(OP_new%fe_gp_range)
	allocate(OP_new%fe_gp_range(1:2, 1:OP%Nfe))
	OP_new%fe_gp_range(:,:) = OP%fe_gp_range(:,:)

!	call fedvr_op_full_allocate(OP_new)

	if(present(behaviour))then
		if(behaviour == 'alloc')then
			call Op_allocate(OP_new)
		else
			call Op_allocate(OP_new)
			do ii_fe=1,OP_new%Nfe
				OP_new%fe_region(ii_fe)%matrix(:,:) = OP%fe_region(ii_fe)%matrix(:,:)
			enddo
		endif
	else
		call Op_allocate(OP_new)
		do ii_fe=1,OP_new%Nfe
			OP_new%fe_region(ii_fe)%matrix(:,:) = OP%fe_region(ii_fe)%matrix(:,:)
		enddo
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_diag_mirror(OP, OP_new, behaviour)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN)  :: OP
	TYPE(fedvr_op_diag), INTENT(OUT) :: OP_new
	CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: behaviour

	OP_new%fe_range = OP%fe_range
	OP_new%gp_range = OP%gp_range
	OP_new%Nfun     = OP%Nfun
	OP_new%Nfe      = OP%Nfe
	OP_new%Ngp      = OP%Ngp
	OP_new%ll       = OP%ll
	OP_new%mm       = OP%mm

	if(allocated(OP_new%fe_gp_range))deallocate(OP_new%fe_gp_range)
	allocate(OP_new%fe_gp_range(1:2, 1:OP%Nfe))
	OP_new%fe_gp_range(:,:) = OP%fe_gp_range(:,:)

	if(present(behaviour))then
		if(behaviour == 'alloc')then
			call Op_allocate(OP_new)
		else
			call Op_allocate(OP_new)
			OP_new%diag(:) = OP%diag(:)
		endif
	else
		call Op_allocate(OP_new)
		OP_new%diag(:) = OP%diag(:)
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE EE_op_mirror(OP, OP_new)
	IMPLICIT NONE
	TYPE(EE_op), INTENT(IN)  :: OP
	TYPE(EE_op), INTENT(OUT) :: OP_new
	INTEGER(I4B) :: ii !, tag


	OP_new%Y_rep    = OP%Y_rep
	OP_new%fe_range = OP%fe_range
	OP_new%gp_range = OP%gp_range
	OP_new%tag_range= OP%tag_range
	OP_new%Ntag     = OP%Ntag
	OP_new%Nfun     = OP%Nfun
	OP_new%Nfe      = OP%Nfe
	OP_new%Ngp      = OP%Ngp
	OP_new%Npulse   = OP%Npulse
	OP_new%EE_t     = OP%EE_t
	OP_new%Env_t    = OP%Env_t

	if(allocated(OP_new%Ei_t))  deallocate(OP_new%Ei_t)
	if(allocated(OP_new%Envi_t))deallocate(OP_new%Envi_t)
	allocate(OP_new%Ei_t(OP%Npulse), OP_new%Envi_t(OP%Npulse))
	OP_new%Ei_t(:)   = OP%Ei_t(:)
	OP_new%Envi_t(:) = OP%Envi_t(:)

	if(allocated(OP_new%pulses))deallocate(OP_new%pulses)
	allocate(OP_new%pulses(OP%Npulse))
	do ii=1,OP%Npulse
		OP_new%pulses = OP%pulses
	enddo 

	if(allocated(OP_new%fe_gp_range))deallocate(OP_new%fe_gp_range)
	allocate(OP_new%fe_gp_range(1:2, 1:OP_new%Nfe))
	OP_new%fe_gp_range(:,:) = OP%fe_gp_range(:,:)

	call Op_allocate(OP_new)

!	do tag = OP_new%tag_range(1), OP_new%tag_range(2)	! the above OP_allocate takes care
!		OP_new%cc(tag) = OP%cc(tag)
!	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE AA_op_mirror(OP, OP_new)
	IMPLICIT NONE
	TYPE(AA_op), INTENT(IN)  :: OP
	TYPE(AA_op), INTENT(OUT) :: OP_new
	INTEGER(I4B) :: ii !,tag

	OP_new%Y_rep    = OP%Y_rep
	OP_new%fe_range = OP%fe_range
	OP_new%gp_range = OP%gp_range
	OP_new%tag_range= OP%tag_range
	OP_new%Ntag     = OP%Ntag
	OP_new%Nfun     = OP%Nfun
	OP_new%Nfe      = OP%Nfe
	OP_new%Ngp      = OP%Ngp
	OP_New%Npulse   = OP%Npulse
	OP_new%AA_t     = OP%AA_t
	OP_new%Env_t    = OP%Env_t

	if(allocated(OP_new%Ai_t))  deallocate(OP_new%Ai_t)
	if(allocated(OP_new%Envi_t))deallocate(OP_new%Envi_t)
	allocate(OP_new%Ai_t(OP%Npulse), OP_new%Envi_t(OP%Npulse))
	OP_new%Ai_t(:)   = OP%Ai_t(:)
	OP_new%Envi_t(:) = OP%Envi_t(:)

	if(allocated(OP_new%pulses))deallocate(OP_new%pulses)
	allocate(OP_new%pulses(OP%Npulse))
	do ii=1,OP%Npulse
		OP_new%pulses = OP%pulses
	enddo 

	if(allocated(OP_new%fe_gp_range))deallocate(OP_new%fe_gp_range)
	allocate(OP_new%fe_gp_range(1:2, 1:OP_new%Nfe))
	OP_new%fe_gp_range(:,:) = OP%fe_gp_range(:,:)

	call Op_allocate(OP_new)

!	do tag = OP_new%tag_range(1), OP_new%tag_range(2)
!		OP_new%cc(tag) = OP%cc(tag)
!	enddo

	call Op_mirror(OP%deriv1, OP_new%deriv1)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_bigmat_OP_full(OP, bigmat)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN)    :: OP
!	TYPE(OP_BigMatrix),  INTENT(INOUT) :: bigmat
	REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: bigmat

!	bigmat%Ngp = OP%Ngp
!	if(allocated(bigmat%matrix))deallocate(bigmat%matrix)
!	allocate(bigmat%matrix(1:OP%Ngp,1:OP%Ngp))
	if(allocated(bigmat))deallocate(bigmat)
	allocate(bigmat(1:OP%Ngp,1:OP%Ngp))
	bigmat(:,:) = 0.0d0

	do ii_fe = 1,OP%Nfe
!		bigmat%matrix(OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe),OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe)) = bigmat%matrix(OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe),OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe)) + OP%fe_region(ii_fe)%matrix(:,:)
		bigmat(OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe),OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe)) = bigmat(OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe),OP%fe_gp_range(1,ii_fe):OP%fe_gp_range(2,ii_fe)) + OP%fe_region(ii_fe)%matrix(:,:)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_bigmat_OP_diag(OP, bigmat)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN)    :: OP
!	TYPE(OP_BigMatrix),  INTENT(INOUT) :: bigmat
	REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: bigmat
	INTEGER(I4B) :: ii

!	bigmat%Ngp = OP%Ngp
!	if(allocated(bigmat%matrix))deallocate(bigmat%matrix)
!	allocate(bigmat%matrix(1:OP%Ngp,1:OP%Ngp))
	if(allocated(bigmat))deallocate(bigmat)
	allocate(bigmat(1:OP%Ngp,1:OP%Ngp))
	bigmat(:,:) = 0.0d0

	do ii=1,OP%Ngp
!		bigmat%matrix(ii,ii) = OP%diag(ii)
		bigmat(ii,ii) = OP%diag(ii)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_tofile_full(OP, fileName, prec)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN) :: OP
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	CHARACTER(LEN=*), INTENT(IN) :: prec
	CHARACTER(LEN=5) :: n_fmt
	INTEGER(I4B) :: ii, jj, unit_nr
	
	write(n_fmt,'(I5.1)')OP%Nfun

	call get_free_unit_nr(unit_nr)
	open(unit=unit_nr, file=fileName//'.dat', status='replace')
		write(unit_nr,*)'# fe_range = ', OP%fe_range(1), OP%fe_range(2)
		write(unit_nr,*)'# gp_range = ', OP%gp_range(1), OP%gp_range(2)
		write(unit_nr,*)'# Nfe = ', OP%Nfe
		write(unit_nr,*)'# Ngp = ', OP%Ngp
		write(unit_nr,*)'# Nfun= ', OP%Nfun
		write(unit_nr,*)'# ll  = ', OP%ll
		write(unit_nr,*)'# mm  = ', OP%mm
		write(unit_nr,*)''

		do ii_fe=1,OP%Nfe
			jj_fe=OP%fe_range(1) - 1 + ii_fe
			do ii=1,OP%Nfun
				write(unit_nr,'(I5.1,'//trim(n_fmt)//'(2x,'//trim(prec)//'))') jj_fe,(OP%fe_region(ii_fe)%matrix(ii,jj), jj=1,OP%Nfun)
			enddo
			write(unit_nr,*)''
		enddo

	close(unit_nr)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_tofile_diag(OP, fileName, prec)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN) :: OP
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	CHARACTER(LEN=*), INTENT(IN) :: prec
	CHARACTER(LEN=5) :: n_fmt
	INTEGER(I4B) :: ii, jj, unit_nr
	
	write(n_fmt,'(I5.1)')OP%Nfun

	call get_free_unit_nr(unit_nr)
	open(unit=unit_nr, file=fileName//'.dat', status='replace')
		write(unit_nr,*)'# fe_range = ', OP%fe_range(1), OP%fe_range(2)
		write(unit_nr,*)'# gp_range = ', OP%gp_range(1), OP%gp_range(2)
		write(unit_nr,*)'# Nfe = ', OP%Nfe
		write(unit_nr,*)'# Ngp = ', OP%Ngp
		write(unit_nr,*)'# Nfun= ', OP%Nfun
		write(unit_nr,*)'# ll  = ', OP%ll
		write(unit_nr,*)'# mm  = ', OP%mm
		write(unit_nr,*)''

		do ii=1,OP%Ngp
			jj = OP%gp_range(1) - 1 + ii
			write(unit_nr,'(I5.1,2x,'//trim(prec)//')') jj, OP%diag(ii)
		enddo
	close(unit_nr)

END SUBROUTINE
!-------------------------------------------------------------------------------
!	TYPE fedvr_op_full
!		INTEGER(I4B), DIMENSION(1:2) :: fe_range
!		INTEGER(I4B), DIMENSION(1:2) :: gp_range
!		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
!		INTEGER(I4B) :: Nfun, Nfe, Ngp, ll, mm
!		TYPE(fedvr_matrix), DIMENSION(:), ALLOCATABLE :: fe_region
!	END TYPE

!	TYPE fedvr_op_diag
!		INTEGER(I4B), DIMENSION(1:2) :: fe_range
!		INTEGER(I4B), DIMENSION(1:2) :: gp_range
!		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
!		INTEGER(I4B) :: Nfun, Nfe, Ngp, ll, mm
!		REAL(DP), DIMENSION(:), ALLOCATABLE :: diag
!	END TYPE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_full_save_h5(OP, fileName)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN) :: OP
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	CHARACTER(LEN=9), PARAMETER :: dset_name = 'OP_fe_mat'
	INTEGER(HID_T) :: plist_id
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: mspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HSIZE_T), DIMENSION(3) :: dims_f
	INTEGER(HSIZE_T), DIMENSION(2) :: dims_d, dims_m
	INTEGER(HSIZE_T), DIMENSION(1) :: dims_a
	INTEGER(HSIZE_T), DIMENSION(3) :: count_p, offset_p, stride_p, block_p
	INTEGER(I4B) :: error

!	write(6,*)my_id, OP%fe_range
	if(my_id==root_id)then

		call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

		dims_f = (/N_fe, N_fun, N_fun/)
		call h5screate_simple_f(3, dims_f, dspace_id, error)
		call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

		dims_a=(/0/)
		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_fe', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_fe, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_tfe', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_tfe, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_fun', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_fun, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_gp', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_gp, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)
		call h5fclose_f(file_id, error)
	endif


	call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
	call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
	call h5pclose_f(plist_id, error)

	dims_d = (/OP%Nfun, OP%Nfun/)
	dims_m = (/OP%Nfun, OP%Nfun/)
!	write(6,*) my_id, OP%Nfun

	call h5screate_simple_f(2, dims_m, mspace_id, error)
	call h5dopen_f(file_id, dset_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)	! (/N_fe,N_fun,N_fun/)

	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
!	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)			! hangs if the different processes have different number of FE to output!!!
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)			! Use this one!!!

	do ii_fe = OP%fe_range(1), OP%fe_range(2)
		jj_fe = ii_fe - OP%fe_range(1) + 1
!	do jj_fe = 1,OP%Nfe
!		ii_fe = OP%fe_range(1) + jj_fe - 1
!		if(my_id==mpi_size-1)write(6,*)my_id, ii_fe, jj_fe
		count_p  = (/1,1,1/)
		offset_p = (/ii_fe-1,0,0/)
		stride_p = (/1,1,1/)
		block_p  = (/1,OP%Nfun,OP%Nfun/)
!		if(my_id==mpi_size-1)write(6,*) my_id, 'offset = ', offset_p
!		if(my_id==mpi_size-1)write(6,*) my_id, 'block  = ', block_p

		call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
!		if(my_id==mpi_size-1)write(6,*)"chkp1"
		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, OP%fe_region(jj_fe)%matrix, dims_d, error, mspace_id, dspace_id, xfer_prp = plist_id)
!		if(my_id==mpi_size-1)write(6,*)"chkp2"
	enddo

	call h5pclose_f(plist_id, error)
	call h5dclose_f(dset_id, error)
	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)
!	write(6,*)"chkp", my_id

	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_op_diag_save_h5(OP, fileName)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN) :: OP
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	CHARACTER(LEN=7), PARAMETER :: dset_name ='OP_diag'
	INTEGER(HID_T) :: plist_id
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: mspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims_f
	INTEGER(HSIZE_T), DIMENSION(1) :: dims_d, dims_m
	INTEGER(HSIZE_T), DIMENSION(1) :: dims_a
	INTEGER(HSIZE_T), DIMENSION(1) :: count_p, offset_p, stride_p, block_p
	INTEGER(I4B) :: error
	INTEGER(I4B) :: Ngp, il, iu

	if(my_id == root_id)then
		Ngp = OP%Ngp
		il = 1
		iu = OP%Ngp
	else
		Ngp = OP%Ngp - 1
		il = 2
		iu = OP%Ngp
	endif

!	write(6,*)my_id, OP%fe_range
	if(my_id==root_id)then

		call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

		dims_f = (/N_gp/)
		call h5screate_simple_f(1, dims_f, dspace_id, error)
		call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

		dims_a=(/0/)
		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_fe', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_fe, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_tfe', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_tfe, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_fun', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_fun, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5screate_f(H5S_SCALAR_F, aspace_id, error)
		call h5acreate_f(dset_id, 'N_gp', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
		call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_gp, dims_a, error)
		call h5aclose_f(attr_id, error)

		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)
		call h5fclose_f(file_id, error)
	endif


	call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
	call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
	call h5pclose_f(plist_id, error)

	dims_d = (/Ngp/)
	dims_m = (/Ngp/)
!	write(6,*) my_id, OP%Nfun

	call h5screate_simple_f(1, dims_m, mspace_id, error)
	call h5dopen_f(file_id, dset_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)	! (/N_fe,N_fun,N_fun/)

	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)			! in this case it doesn't matter (no loop)
!	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)			! 

	if(my_id == root_id)then
		count_p  = (/1/)
		offset_p = (/OP%gp_range(1)-1/)
		stride_p = (/1/)
		block_p  = (/Ngp/)
	else
		count_p  = (/1/)
		offset_p = (/OP%gp_range(1)/)
		stride_p = (/1/)
		block_p  = (/Ngp/)
	endif

	call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, OP%diag(il:iu), dims_d, error, mspace_id, dspace_id, xfer_prp = plist_id)

	call h5pclose_f(plist_id, error)
	call h5dclose_f(dset_id, error)
	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)
!	write(6,*)"chkp", my_id

	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------


END MODULE Operator_typedef














