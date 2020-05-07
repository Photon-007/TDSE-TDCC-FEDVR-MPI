MODULE Operator_math
	USE nrtype
	USE gvars
	USE Operator_typedef
	USE FEDVR_module
	USE TDCC_module
	USE TDCC_FEDVR

	INTERFACE OPERATOR (+)
		MODULE PROCEDURE ADD_full_full
		MODULE PROCEDURE ADD_full_diag
		MODULE PROCEDURE ADD_diag_full
		MODULE PROCEDURE ADD_diag_diag
		MODULE PROCEDURE ADD_EE_EE				! this would be an elegant way of dealing with multiple pulses: set up the sum at initialization
		MODULE PROCEDURE ADD_AA_AA
	END INTERFACE

	INTERFACE OPERATOR (*)
		MODULE PROCEDURE MULTIPLY_Op1_full_scalar
		MODULE PROCEDURE MULTIPLY_Op2_full_scalar
		MODULE PROCEDURE MULTIPLY_Op1_diag_scalar
		MODULE PROCEDURE MULTIPLY_Op2_diag_scalar
		MODULE PROCEDURE MULTIPLY_Op1_full_psi
		MODULE PROCEDURE MULTIPLY_Op2_full_psi
		MODULE PROCEDURE MULTIPLY_Op1_diag_psi
		MODULE PROCEDURE MULTIPLY_Op2_diag_psi
		MODULE PROCEDURE MULTIPLY_EE_psi
!		MODULE PROCEDURE MULTIPLY_AA_psi
	END INTERFACE


	CONTAINS
!-------------------------------------------------------------------------------
FUNCTION ADD_full_full(OP1, OP2) RESULT(OP_sum)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN)  :: OP1
	TYPE(fedvr_op_full), INTENT(IN)  :: OP2
	TYPE(fedvr_op_full)              :: OP_sum

	call Op_mirror(OP1, OP_sum)
	do ii_fe=1,OP_sum%Nfe
		OP_sum%fe_region(ii_fe)%matrix(:,:) = OP_sum%fe_region(ii_fe)%matrix(:,:) + OP2%fe_region(ii_fe)%matrix(:,:)
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION ADD_full_diag(OP1, OP2) RESULT(OP_sum)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN)  :: OP1
	TYPE(fedvr_op_diag), INTENT(IN)  :: OP2
	TYPE(fedvr_op_full)              :: OP_sum
	INTEGER(I4B) :: ii

	call Op_mirror(OP1, OP_sum)
	do ii_fe=1,OP_sum%Nfe
		jj_fe = OP_sum%fe_range(1) - 1 + ii_fe

		if(jj_fe>1)then
			OP_sum%fe_region(ii_fe)%matrix(1,1) = OP_sum%fe_region(ii_fe)%matrix(1,1) + 0.5d0 * OP2%diag(OP2%fe_gp_range(1, ii_fe))
		else
			OP_sum%fe_region(ii_fe)%matrix(1,1) = OP_sum%fe_region(ii_fe)%matrix(1,1) +         OP2%diag(OP2%fe_gp_range(1, ii_fe))
		endif

		if(jj_fe<N_fe)then
			OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) = OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) + 0.5d0 * OP2%diag(OP2%fe_gp_range(2, ii_fe))
		else
			OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) = OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) +         OP2%diag(OP2%fe_gp_range(2, ii_fe))
		endif

		do ii=2,OP_sum%Nfun-1
			OP_sum%fe_region(ii_fe)%matrix(ii,ii) = OP_sum%fe_region(ii_fe)%matrix(ii,ii) + OP2%diag(OP2%fe_gp_range(1, ii_fe) - 1 + ii)
		enddo
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION ADD_diag_full(OP1, OP2) RESULT(OP_sum)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN)  :: OP1
	TYPE(fedvr_op_full), INTENT(IN)  :: OP2
	TYPE(fedvr_op_full)              :: OP_sum
	INTEGER(I4B) :: ii

	call Op_mirror(OP2, OP_sum)
	do ii_fe=1,OP_sum%Nfe
		jj_fe = OP_sum%fe_range(1) - 1 + ii_fe

		if(jj_fe>1)then
			OP_sum%fe_region(ii_fe)%matrix(1,1) = OP_sum%fe_region(ii_fe)%matrix(1,1) + 0.5d0 * OP1%diag(OP1%fe_gp_range(1, ii_fe))
		else
			OP_sum%fe_region(ii_fe)%matrix(1,1) = OP_sum%fe_region(ii_fe)%matrix(1,1) +         OP1%diag(OP1%fe_gp_range(1, ii_fe))
		endif

		if(jj_fe<N_fe)then
			OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) = OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) + 0.5d0 * OP1%diag(OP1%fe_gp_range(2, ii_fe))
		else
			OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) = OP_sum%fe_region(ii_fe)%matrix(OP_sum%Nfun,OP_sum%Nfun) +         OP1%diag(OP1%fe_gp_range(2, ii_fe))
		endif

		do ii=2,OP_sum%Nfun-1
			OP_sum%fe_region(ii_fe)%matrix(ii,ii) = OP_sum%fe_region(ii_fe)%matrix(ii,ii) + OP1%diag(OP2%fe_gp_range(1, ii_fe) - 1 + ii)
		enddo
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION ADD_diag_diag(OP1, OP2) RESULT(OP_sum)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN)  :: OP1
	TYPE(fedvr_op_diag), INTENT(IN)  :: OP2
	TYPE(fedvr_op_diag)              :: OP_sum


	call Op_mirror(OP1, OP_sum)
	OP_sum%diag(:) = OP_sum%diag(:) + OP2%diag(:)

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION ADD_EE_EE(OP1, OP2) RESULT(OP_sum)
	IMPLICIT NONE
	TYPE(EE_op), INTENT(IN)  :: OP1
	TYPE(EE_op), INTENT(IN)  :: OP2
	TYPE(EE_op)              :: OP_sum
!	INTEGER(I4B) :: ii

	call Op_mirror(OP1, OP_sum)

!	OP_sum%Ei_t   = array_append(OP_sum%Ei_t,  OP2%Ei_t)
	OP_sum%Npulse = OP1%Npulse + OP2%Npulse
	if(allocated(OP_sum%pulses))deallocate(OP_sum%pulses)
	allocate(OP_sum%pulses(OP_sum%Npulse))

	if(allocated(OP_sum%Ei_t))deallocate(OP_sum%Ei_t)
	allocate(OP_sum%Ei_t(OP_sum%Npulse))

	OP_sum%pulses(1:OP1%Npulse)               = OP1%pulses(1:OP1%Npulse)
	OP_sum%pulses(OP1%Npulse+1:OP_sum%Npulse) = OP2%pulses(1:OP2%Npulse)
	OP_sum%Ei_t(1:OP1%Npulse)                 = OP1%Ei_t(1:OP1%Npulse)
	OP_sum%Ei_t(OP1%Npulse+1:OP_sum%Npulse)   = OP2%Ei_t(1:OP2%Npulse)

!	do ii = 1,OP1%Npulse
!		OP_sum%pulses(ii) = OP1%pulses(ii)
!		OP_sum%Ei_t(ii)   = OP1%Ei_t(ii)
!	enddo
!	do ii = 1,OP2%Npulse
!		OP_sum%pulses(OP1%Npulse + ii) = OP2%pulses(ii)
!		OP_sum%Ei_t(OP1%Npulse + ii)   = OP2%Ei_t(ii)
!	enddo

	OP_sum%EE_t = OP1%EE_t + OP2%EE_t

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION ADD_AA_AA(OP1, OP2) RESULT(OP_sum)
	IMPLICIT NONE
	TYPE(AA_op), INTENT(IN)  :: OP1
	TYPE(AA_op), INTENT(IN)  :: OP2
	TYPE(AA_op)              :: OP_sum
!	INTEGER(I4B) :: ii

	call Op_mirror(OP1, OP_sum)

!	OP_sum%Ai_t   = array_append(OP_sum%Ai_t,  OP2%Ai_t)
	OP_sum%Npulse = OP1%Npulse + OP2%Npulse
	if(allocated(OP_sum%pulses))deallocate(OP_sum%pulses)
	allocate(OP_sum%pulses(OP_sum%Npulse))

	if(allocated(OP_sum%Ai_t))deallocate(OP_sum%Ai_t)
	allocate(OP_sum%Ai_t(OP_sum%Npulse))

	OP_sum%pulses(1:OP1%Npulse)               = OP1%pulses(1:OP1%Npulse)
	OP_sum%pulses(OP1%Npulse+1:OP_sum%Npulse) = OP2%pulses(1:OP2%Npulse)
	OP_sum%Ai_t(1:OP1%Npulse)                 = OP1%Ai_t(1:OP1%Npulse)
	OP_sum%Ai_t(OP1%Npulse+1:OP_sum%Npulse)   = OP2%Ai_t(1:OP2%Npulse)
!	do ii = 1,OP1%Npulse
!		OP_sum%pulses(ii) = OP1%pulses(ii)
!		OP_sum%Ai_t(ii)   = OP1%Ai_t(ii)
!	enddo
!	do ii = 1,OP2%Npulse
!		OP_sum%pulses(OP1%Npulse + ii) = OP2%pulses(ii)
!		OP_sum%Ai_t(OP1%Npulse + ii)   = OP2%Ai_t(ii)
!	enddo

	OP_sum%AA_t = OP1%AA_t + OP2%AA_t

END FUNCTION
!-------------------------------------------------------------------------------
!===============================================================================
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op1_full_scalar(scalar, OP) RESULT(OP_new)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN) :: OP
	TYPE(fedvr_op_full)  :: OP_new
	REAL(DP), INTENT(IN) :: scalar

	call Op_mirror(OP, OP_new)
	do ii_fe=1,OP%Nfe
		OP_new%fe_region(ii_fe)%matrix(:,:) = scalar * OP_new%fe_region(ii_fe)%matrix(:,:)
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op2_full_scalar(scalar, OP) RESULT(OP_new)
	IMPLICIT NONE
	TYPE(fedvr_op_full), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: OP
	TYPE(fedvr_op_full), DIMENSION(:), ALLOCATABLE   :: OP_new
	REAL(DP), INTENT(IN) :: scalar
	INTEGER(I4B) :: rank, rank_lb, rank_ub

	rank_lb = lbound(OP,1,I4B)
	rank_ub = ubound(OP,1,I4B)

	if(allocated(OP_new))deallocate(OP_new)
	allocate(OP_new(rank_lb:rank_ub))
	do rank=rank_lb, rank_ub
		call OP_mirror(OP(rank), OP_new(rank))
		do ii_fe=1,OP(rank)%Nfe
			OP_new(rank)%fe_region(ii_fe)%matrix(:,:) = scalar * OP_new(rank)%fe_region(ii_fe)%matrix(:,:)
		enddo
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op1_diag_scalar(scalar, OP) RESULT(OP_new)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN) :: OP
	TYPE(fedvr_op_diag)  :: OP_new
	REAL(DP), INTENT(IN) :: scalar

	call Op_mirror(OP, OP_new)
	OP_new%diag(:) = scalar * OP_new%diag(:)

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op2_diag_scalar(scalar, OP) RESULT(OP_new)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), DIMENSION(:), INTENT(IN) :: OP
	TYPE(fedvr_op_diag), DIMENSION(:), ALLOCATABLE   :: OP_new
	REAL(DP), INTENT(IN) :: scalar
	INTEGER(I4B) :: rank, rank_lb, rank_ub

	rank_lb = lbound(OP,1,I4B)
	rank_ub = ubound(OP,1,I4B)

	if(allocated(OP_new))deallocate(OP_new)
	allocate(OP_new(rank_lb:rank_ub))
	do rank=rank_lb, rank_ub
		call OP_mirror(OP(rank), OP_new(rank))
		OP_new(rank)%diag(:) = scalar * OP_new(rank)%diag(:)
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op1_full_psi(OP, wf) RESULT(wf_dummy)
	IMPLICIT NONE
	TYPE(fedvr_op_full),  INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	TYPE(tdcc_fedvr_wf) :: wf_dummy
	INTEGER(I4B) :: tag


	call wf_mirror(wf, wf_dummy, 'alloc')
	do tag=wf%tag_range(1),wf%tag_range(2)
		do ii_fe=1,OP%Nfe
			call dgemv('N', OP%Nfun, OP%Nfun, 1.0d0, OP%fe_region(ii_fe)%matrix, OP%Nfun, wf%re(wf%fe_gp_range(1,ii_fe), tag), 1, 1.0d0, wf_dummy%re(wf_dummy%fe_gp_range(1,ii_fe), tag), 1)
			call dgemv('N', OP%Nfun, OP%Nfun, 1.0d0, OP%fe_region(ii_fe)%matrix, OP%Nfun, wf%im(wf%fe_gp_range(1,ii_fe), tag), 1, 1.0d0, wf_dummy%im(wf_dummy%fe_gp_range(1,ii_fe), tag), 1)
		enddo
	enddo

	if(wf%Im_a_chunk)then
		call mpi_update_wf_at_bounds(wf_dummy)
	endif

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op2_full_psi(OP, wf) RESULT(wf_dummy)
	IMPLICIT NONE
	TYPE(fedvr_op_full), DIMENSION(:), INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	TYPE(tdcc_fedvr_wf) :: wf_dummy
	INTEGER(I4B) :: tag

	call wf_mirror(wf, wf_dummy, 'alloc')
	do tag=wf%tag_range(1),wf%tag_range(2)
		do ii_fe=1,OP(tag)%Nfe
			call dgemv('N', OP(tag)%Nfun, OP(tag)%Nfun, 1.0d0, OP(tag)%fe_region(ii_fe)%matrix, OP(tag)%Nfun, wf%re(wf%fe_gp_range(1,ii_fe), tag), 1, 1.0d0, wf_dummy%re(wf_dummy%fe_gp_range(1,ii_fe), tag), 1)
			call dgemv('N', OP(tag)%Nfun, OP(tag)%Nfun, 1.0d0, OP(tag)%fe_region(ii_fe)%matrix, OP(tag)%Nfun, wf%im(wf%fe_gp_range(1,ii_fe), tag), 1, 1.0d0, wf_dummy%im(wf_dummy%fe_gp_range(1,ii_fe), tag), 1)
		enddo
	enddo

	if(wf%Im_a_chunk)then
		call mpi_update_wf_at_bounds(wf_dummy)
	endif

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op1_diag_psi(OP, wf) RESULT(wf_dummy)
	IMPLICIT NONE
	TYPE(fedvr_op_diag),  INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	TYPE(tdcc_fedvr_wf) :: wf_dummy
	INTEGER(I4B) :: tag

	call wf_mirror(wf, wf_dummy, 'alloc')
	do tag=wf%tag_range(1),wf%tag_range(2)
		wf_dummy%re(:, tag) = OP%diag(:) * wf%re(:, tag)
		wf_dummy%im(:, tag) = OP%diag(:) * wf%im(:, tag)
	enddo

!	if(wf%Im_a_chunk)then							! NO update for diag ops !!!
!		call mpi_update_wf_at_bounds(wf_dummy)
!	endif

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_Op2_diag_psi(OP, wf) RESULT(wf_dummy)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), DIMENSION(:), INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	TYPE(tdcc_fedvr_wf) :: wf_dummy
	INTEGER(I4B) :: tag

	call wf_mirror(wf, wf_dummy, 'alloc')
	do tag=wf%tag_range(1),wf%tag_range(2)
		wf_dummy%re(:, tag) = OP(tag)%diag(:) * wf%re(:, tag)
		wf_dummy%im(:, tag) = OP(tag)%diag(:) * wf%im(:, tag)
	enddo

!	if(wf%Im_a_chunk)then
!		call mpi_update_wf_at_bounds(wf_dummy)
!	endif

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION MULTIPLY_EE_psi(OP, wf) RESULT(wf_dummy)
	IMPLICIT NONE
	TYPE(EE_op), INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	TYPE(tdcc_fedvr_wf) :: wf_dummy
	INTEGER(I4B) :: tag, rank!, ip
	REAL(DP)     :: EE_t

!	call wf_mirror(wf, wf_dummy)
!	do ip=1,OP%Npulse				! NO need for this loop. It is better to apply the Sum electric field (see below) 
!		do tag=wf%tag_range(1),wf%tag_range(2)
!			do rank=1,OP%Nrank
!				if(OP%cc(tag)%couple2(rank) /= -1)then
!					EE_t = OP%Ei_t(ip) * OP%cc(tag)%value(rank)
!					wf_dummy%re(:,tag) = wf_dummy%re(:,tag) + EE_t * grid%xx(:) * wf%re(:,OP%cc(tag)%couple2(rank))
!					wf_dummy%im(:,tag) = wf_dummy%im(:,tag) + EE_t * grid%xx(:) * wf%im(:,OP%cc(tag)%couple2(rank))
!				endif
!			enddo
!		enddo
!	enddo

	call wf_mirror(wf, wf_dummy)
	do tag=wf%tag_range(1),wf%tag_range(2)
		do rank=1,OP%Nrank
			if(OP%cc(tag)%couple2(rank) /= -1)then
				EE_t = OP%EE_t * OP%cc(tag)%value(rank)
				wf_dummy%re(:,tag) = wf_dummy%re(:,tag) + EE_t * grid%xx(:) * wf%re(:,OP%cc(tag)%couple2(rank))
				wf_dummy%im(:,tag) = wf_dummy%im(:,tag) + EE_t * grid%xx(:) * wf%im(:,OP%cc(tag)%couple2(rank))
			endif
		enddo
	enddo



!	if(wf%Im_a_chunk)then
!		call mpi_update_wf_at_bounds(wf_dummy)
!	endif

END FUNCTION
!-------------------------------------------------------------------------------
!FUNCTION array_append(array, item) RESULT(array)
!	IMPLICIT NONE
!	REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
!	REAL(DP), DIMENSION(:), INTENT(IN)  :: item
!	REAL(DP), DIMENSION(:), ALLOCATABLE :: temp
!	INTEGER(I4B) :: nn, mm

!	mm = size(item)

!	if(allocated(array))then
!		nn = size(array)
!		allocate(temp(1:nn))
!		temp(:) = array(:)
!		deallocate(array)
!		allocate(array(nn+mm))
!		array(1:nn) = temp(1:nn)
!		array(nn+1:nn+mm) = item(1:mm)
!	else
!		allocate(array(mm))
!		array(:) = item(:)
!	endif

!END FUNCTION
!-------------------------------------------------------------------------------










END MODULE Operator_math
