MODULE TDCC_FEDVR
	USE nrtype
	USE gvars
	USE FEDVR_module
	USE TDCC_module
	USE BoundStates
	USE MPI_pars
	USE mpi
	USE HDF5

	INTERFACE OPERATOR (*)
		MODULE PROCEDURE psi_x_psi
	END INTERFACE

	INTERFACE OPERATOR (+)
		MODULE PROCEDURE psi_add_psi
	END INTERFACE

	INTERFACE wf_scale
		MODULE PROCEDURE wf_scale_real
		MODULE PROCEDURE wf_scale_cmplx
	END INTERFACE

	INTERFACE wf_axpy
		MODULE PROCEDURE wf_axpy_real
		MODULE PROCEDURE wf_axpy_cmplx
	END INTERFACE


	TYPE :: tdcc_fedvr_wf
		LOGICAL(LGT) :: Im_a_chunk 
		LOGICAL(LGT) :: w_inc
		INTEGER(I4B), DIMENSION(1:2) :: fe_range
		INTEGER(I4B), DIMENSION(1:2) :: gp_range
		INTEGER(I4B), DIMENSION(1:2) :: tag_range
		INTEGER(I4B) :: Nfun, Nfe, Ngp, Ntag
		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
		REAL(DP),     DIMENSION(:,:), ALLOCATABLE :: re, im
		REAL(DP),     DIMENSION(:),   ALLOCATABLE :: tag_norm, tag_norm_loc		! tag_norm_loc and norm_loc are just dummyes used to calculate the 
		REAL(DP) :: norm, norm_prev, norm_loc, norm_evol									! local values that will be MPI_ALLREDUCE-d to tag_norm and norm (no time consuming allocation)
	END TYPE

	TYPE(tdcc_fedvr_wf) :: psi, psi_G
	TYPE(tdcc_fedvr_wf) :: wf_dummy, wf_temp

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE wf_create(wf, grid, tag_range)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	TYPE(fedvr_grid),    INTENT(IN)    :: grid
	INTEGER(I4B), DIMENSION(2), INTENT(IN) :: tag_range	! this way the same formalism can be used if we are deviding the workload along the tag 

	wf%w_inc     = .TRUE.
	wf%Im_a_chunk= grid%Im_a_chunk
	wf%fe_range  = grid%fe_range
	wf%gp_range  = grid%gp_range
	wf%tag_range = tag_range
	wf%Ntag      = tag_range(2) - tag_range(1) + 1
	wf%Nfun      = grid%Nfun
	wf%Nfe       = grid%Nfe
	wf%Ngp       = grid%Ngp
	wf%norm      = 0.0d0
	wf%norm_prev = 0.0d0
	wf%norm_loc  = 0.0d0
	wf%norm_evol = 1.0d0

	if(allocated(wf%tag_norm_loc))deallocate(wf%tag_norm_loc)
	if(allocated(wf%tag_norm))deallocate(wf%tag_norm)
	if(allocated(wf%re))deallocate(wf%re)
	if(allocated(wf%im))deallocate(wf%im)

	if(allocated(wf%fe_gp_range))deallocate(wf%fe_gp_range)
	allocate(wf%fe_gp_range(1:2, 1:grid%Nfe))
	wf%fe_gp_range(:,:) = grid%fe_gp_range(:,:)


END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_allocate(wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf

	if(allocated(wf%tag_norm_loc))deallocate(wf%tag_norm_loc)
	if(allocated(wf%tag_norm))deallocate(wf%tag_norm)
	if(allocated(wf%re))deallocate(wf%re)
	if(allocated(wf%im))deallocate(wf%im)
	allocate(wf%re(1:wf%Ngp, 1:wf%Ntag))
	allocate(wf%im(1:wf%Ngp, 1:wf%Ntag))
	allocate(wf%tag_norm(1:wf%Ntag))
	allocate(wf%tag_norm_loc(1:wf%Ntag))

	wf%re(:,:)         = 0.0d0
	wf%im(:,:)         = 0.0d0
	wf%tag_norm(:)     = 0.0d0
	wf%tag_norm_loc(:) = 0.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_deallocate(wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf

	if(allocated(wf%tag_norm_loc))deallocate(wf%tag_norm_loc)
	if(allocated(wf%tag_norm))deallocate(wf%tag_norm)
	if(allocated(wf%re))deallocate(wf%re)
	if(allocated(wf%im))deallocate(wf%im)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_clear(wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf

	wf%norm_prev   = wf%norm
	wf%norm        = 0.0d0
	wf%re(:,:)     = 0.0d0
	wf%im(:,:)     = 0.0d0
	wf%tag_norm(:) = 0.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_copy(wf_source, wf_target)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN   ) :: wf_source
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf_target

	wf_target%w_inc       = wf_source%w_inc
	wf_target%re(:,:)     = wf_source%re(:,:)
	wf_target%im(:,:)     = wf_source%im(:,:)
	wf_target%tag_norm(:) = wf_source%tag_norm(:)
	wf_target%norm_evol   = wf_source%norm_evol
	wf_target%norm        = wf_source%norm
	wf_target%norm_prev   = wf_source%norm_prev

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_mirror(wf, wf_new, behaviour)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN   ) :: wf
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf_new
	CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: behaviour
	REAL(DP) :: ti, tf, elapsed, elapsed_max

	N_wf_mirror = N_wf_mirror + 1
	ti = MPI_Wtime()

	wf_new%w_inc     = wf%w_inc
	wf_new%Im_a_chunk= wf%Im_a_chunk
	wf_new%fe_range  = wf%fe_range
	wf_new%gp_range  = wf%gp_range
	wf_new%tag_range = wf%tag_range
	wf_new%Ntag      = wf%Ntag
	wf_new%Nfun      = wf%Nfun
	wf_new%Nfe       = wf%Nfe
	wf_new%Ngp       = wf%Ngp
	wf_new%norm      = wf%norm
	wf_new%norm_prev = wf%norm_prev
	wf_new%norm_loc  = wf%norm_loc
	wf_new%norm_evol = wf%norm_evol

	if(allocated(wf_new%fe_gp_range))deallocate(wf_new%fe_gp_range)
	allocate(wf_new%fe_gp_range(1:2, 1:wf%Nfe))
	wf_new%fe_gp_range(:,:) = wf%fe_gp_range(:,:)

	if(present(behaviour))then
		if(behaviour == 'alloc')then
			call wf_allocate(wf_new)
			wf_new%re(:,:)     = 0.0d0
			wf_new%im(:,:)     = 0.0d0
			wf_new%tag_norm(:) = 0.0d0
		else
			call wf_allocate(wf_new)
			wf_new%re(:,:)     = wf%re(:,:)
			wf_new%im(:,:)     = wf%im(:,:)
			wf_new%tag_norm(:) = wf%tag_norm(:)
		endif
	else
		call wf_allocate(wf_new)
		wf_new%re(:,:)     = wf%re(:,:)
		wf_new%im(:,:)     = wf%im(:,:)
		wf_new%tag_norm(:) = wf%tag_norm(:)
	endif

	tf = MPI_Wtime()
	elapsed = tf - ti
!	call MPI_Allreduce(elapsed, elapsed_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_err)
	Time_wf_mirror = Time_wf_mirror + elapsed

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_norm(wf)	! use this instead of the wf_dotp as it performs half the calcs!!! (no imaginary part)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	INTEGER(I4B) :: tag, lb, ub


	if(wf%fe_range(1)==1)then
		lb = 1
		ub = wf%Ngp
	else
		lb = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		ub = wf%Ngp
	endif

	do tag = 1,wf%Ntag
		wf%tag_norm_loc(tag) = sum(wf%re(lb:ub, tag) * wf%re(lb:ub, tag)) + sum(wf%im(lb:ub, tag) * wf%im(lb:ub, tag))
	enddo
	wf%norm_loc = sum(wf%tag_norm_loc)

	if(wf%Im_a_chunk)then
!		write(6,*)"wf_norm(my_id=", my_id,") do_mpi: ", wf%Im_a_chunk, lb, ub
		call MPI_Allreduce(wf%tag_norm_loc, wf%tag_norm, wf%Ntag, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		call MPI_Allreduce(wf%norm_loc,     wf%norm,     1,       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
!		write(6,*)"norm Allreduce finished", my_id
	else
!		write(6,*)"wf_norm(my_id=", my_id,") do_mpi: ", wf%Im_a_chunk, lb, ub
		wf%tag_norm(:) = wf%tag_norm_loc(:)
		wf%norm        = wf%norm_loc
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_normalize(wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf

	call wf_norm(wf)
	wf%re(:,:) = wf%re(:,:) / sqrt(dabs(wf%norm))
	wf%im(:,:) = wf%im(:,:) / sqrt(dabs(wf%norm))
	wf%norm_evol = wf%norm_evol * wf%norm
	wf%norm_prev = wf%norm
	wf%norm = 1.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_dotp(wf1, wf2, dotp)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf1, wf2
	COMPLEX(DPC), INTENT(INOUT) :: dotp
!	COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: dotp_tag		! will never use this one!!!!
	COMPLEX(DPC) :: cc
	REAL(DP) :: cr, ci
	INTEGER(I4B) :: tag, lb, ub


	dotp = CZERO
	cc   = CZERO

	if(wf1%fe_range(1)==1)then
		lb = 1
		ub = wf1%Ngp
	else
		lb = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		ub = wf1%Ngp
	endif

	do tag = 1,wf1%Ntag
		cr = sum(wf1%re(lb:ub, tag) * wf2%re(lb:ub, tag)) + sum(wf1%im(lb:ub, tag) * wf2%im(lb:ub, tag))
		ci = sum(wf1%re(lb:ub, tag) * wf2%im(lb:ub, tag)) - sum(wf1%im(lb:ub, tag) * wf2%re(lb:ub, tag))
		cc = cc + dcmplx(cr,ci)
	enddo

	if(wf1%Im_a_chunk)then
!		write(6,*)"wf_dotp(my_id=", my_id,") do_mpi: ", wf1%Im_a_chunk, lb, ub
		call MPI_Allreduce(cc, dotp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi_err)
!		write(6,*)"dotp Allreduce finished", my_id
	else
!		write(6,*)"wf_dotp(my_id=", my_id,") do_mpi: ", wf1%Im_a_chunk, lb, ub
		dotp = cc
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
FUNCTION psi_x_psi(wf1, wf2) RESULT(dotp)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf1, wf2
	COMPLEX(DPC) :: dotp
	COMPLEX(DPC) :: cc
	REAL(DP) :: cr, ci
	INTEGER(I4B) :: tag, lb, ub


	dotp = CZERO
	cc   = CZERO

	if(wf1%fe_range(1)==1)then
		lb = 1
		ub = wf1%Ngp
	else
		lb = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		ub = wf1%Ngp
	endif

	do tag = 1,wf1%Ntag
		cr = sum(wf1%re(lb:ub, tag) * wf2%re(lb:ub, tag)) + sum(wf1%im(lb:ub, tag) * wf2%im(lb:ub, tag))
		ci = sum(wf1%re(lb:ub, tag) * wf2%im(lb:ub, tag)) - sum(wf1%im(lb:ub, tag) * wf2%re(lb:ub, tag))
		cc = cc + dcmplx(cr,ci)
	enddo

	if(wf1%Im_a_chunk)then
!		write(6,*)"wf_dotp(my_id=", my_id,") do_mpi: ", wf1%Im_a_chunk, lb, ub
		call MPI_Allreduce(cc, dotp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi_err)
!		write(6,*)"dotp Allreduce finished", my_id
	else
!		write(6,*)"wf_dotp(my_id=", my_id,") do_mpi: ", wf1%Im_a_chunk, lb, ub
		dotp = cc
	endif

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION psi_add_psi(wf1, wf2) RESULT(wf_sum)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf1, wf2
	TYPE(tdcc_fedvr_wf)             :: wf_sum
	INTEGER(I4B) :: tag

	call wf_mirror(wf1, wf_sum, 'alloc')
	do tag = 1,wf1%Ntag
		wf_sum%re(1:wf1%Ngp, tag) = wf1%re(1:wf1%Ngp, tag) + wf2%re(1:wf1%Ngp, tag)
		wf_sum%im(1:wf1%Ngp, tag) = wf1%im(1:wf1%Ngp, tag) + wf2%im(1:wf1%Ngp, tag)
	enddo

END FUNCTION
!-------------------------------------------------------------------------------
SUBROUTINE wf_scale_real(re, wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	REAL(DP), INTENT(IN) :: re

	wf%re(:,:) = re * wf%re(:,:)
	wf%im(:,:) = re * wf%im(:,:)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_scale_cmplx(cc, wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	COMPLEX(DPC), INTENT(IN) :: cc
	REAL(DP) :: re, im

	re = dble(cc)
	im = aimag(cc)

	wf%re(:,:) = re * wf%re(:,:) - im * wf%im(:,:)
	wf%im(:,:) = re * wf%im(:,:) + im * wf%re(:,:)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_axpy_real(re, wf_x, wf_y)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN)    :: wf_x
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf_y
	REAL(DP), INTENT(IN) :: re
	INTEGER(I4B) :: tag
	
	do tag = 1,wf_x%Ntag
		call daxpy(wf_x%Ngp, re, wf_x%re(1,tag), 1, wf_y%re(1,tag), 1)
		call daxpy(wf_x%Ngp, re, wf_x%im(1,tag), 1, wf_y%im(1,tag), 1)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_axpy_cmplx(cc, wf_x, wf_y)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN)    :: wf_x
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf_y
	COMPLEX(DPC), INTENT(IN) :: cc
	REAL(DP) :: re, im
	INTEGER(I4B) :: tag

	re = dble(cc)
	im = aimag(cc)
	
	do tag = 1,wf_x%Ntag
		call daxpy(wf_x%Ngp,  re, wf_x%re(1,tag), 1, wf_y%re(1,tag), 1)
		call daxpy(wf_x%Ngp,  re, wf_x%im(1,tag), 1, wf_y%im(1,tag), 1)
		call daxpy(wf_x%Ngp, -im, wf_x%im(1,tag), 1, wf_y%re(1,tag), 1)
		call daxpy(wf_x%Ngp,  im, wf_x%re(1,tag), 1, wf_y%im(1,tag), 1)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_remove_vec_projection(wf_re, wf_im, vec, prec)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: wf_re, wf_im
	REAL(DP), DIMENSION(:), INTENT(IN)    :: vec
	REAL(DP), INTENT(IN) :: prec
	REAL(DP) :: vec_norm, vec_n, cr, ci, cr_tmp, ci_tmp
	REAL(DP), EXTERNAL :: dot, ddot
	INTEGER(I4B) :: lb, np, k
	LOGICAL(LGT) :: ortho

	if(my_fe(1)==1)then
		lb = 1
		np = NN_gp
	else
		lb = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		np = NN_gp - 1
	endif

	vec_n = ddot(np, vec(lb), 1, vec(lb), 1)
	call MPI_Allreduce(vec_n, vec_norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)

	cr_tmp = ddot(np, wf_re(lb), 1, vec(lb), 1)
	ci_tmp = ddot(np, wf_im(lb), 1, vec(lb), 1)
	call MPI_Allreduce(cr_tmp, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
	call MPI_Allreduce(ci_tmp, ci, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)

	ortho=.FALSE.
	k=0
	do while(.not.ortho)
		k = k+1
		call daxpy(NN_gp, -cr/vec_norm, vec(1), 1, wf_re(1), 1)
		call daxpy(NN_gp, -ci/vec_norm, vec(1), 1, wf_im(1), 1)
		cr_tmp = ddot(np, wf_re(lb), 1, vec(lb), 1)
		ci_tmp = ddot(np, wf_im(lb), 1, vec(lb), 1)
		call MPI_Allreduce(cr_tmp, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		call MPI_Allreduce(ci_tmp, ci, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		if(cdabs(dcmplx(cr,ci)) <= prec) ortho = .TRUE.
		if(k>=10) ortho = .TRUE.		! escape
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE	wf_remove_bState(wf, bState, prec)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	TYPE(BoundState),    INTENT(IN)    :: bState
	REAL(DP), INTENT(IN) :: prec
	REAL(DP), EXTERNAL :: dot, ddot
	REAL(DP) :: cr, ci, cr_tmp, ci_tmp
	INTEGER(I4B) :: Ngp, lb, tag, k
	LOGICAL(LGT) :: ortho

	tag = TDCC_get_tag_from_lm(SphericalHarmonicRepresentation, bState%l, bState%m)

	if(wf%fe_range(1)==1)then
		lb  = 1
		Ngp = wf%Ngp
	else
		lb  = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		Ngp = wf%Ngp - 1
	endif


	cr_tmp = ddot(Ngp, wf%re(lb,tag), 1, bState%vector(lb), 1)
	ci_tmp = ddot(Ngp, wf%im(lb,tag), 1, bState%vector(lb), 1)
	call MPI_Allreduce(cr_tmp, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
	call MPI_Allreduce(ci_tmp, ci, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)

	ortho=.FALSE.
	k=0
	do while(.not.ortho)
		k = k+1
		call daxpy(wf%Ngp, -cr/bState%norm, bState%vector(1), 1, wf%re(1,tag), 1)
		call daxpy(wf%Ngp, -ci/bState%norm, bState%vector(1), 1, wf%im(1,tag), 1)
		cr_tmp = ddot(Ngp, wf%re(lb,tag), 1, bState%vector(lb), 1)
		ci_tmp = ddot(Ngp, wf%im(lb,tag), 1, bState%vector(lb), 1)
		call MPI_Allreduce(cr_tmp, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		call MPI_Allreduce(ci_tmp, ci, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		if(cdabs(dcmplx(cr,ci)) <= prec) ortho = .TRUE.
		if(k>=10) ortho = .TRUE.		! escape
	enddo


END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE	wf_remove_bState_ITP(wf, bState, prec)	!!! NO tag !!!
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	TYPE(BoundState),    INTENT(IN)    :: bState
	REAL(DP), INTENT(IN) :: prec
	REAL(DP), EXTERNAL :: dot, ddot
	REAL(DP) :: cr, ci, cr_tmp, ci_tmp
	INTEGER(I4B) :: Ngp, lb, k
	LOGICAL(LGT) :: ortho


	if(wf%fe_range(1)==1)then
		lb  = 1
		Ngp = wf%Ngp
	else
		lb  = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		Ngp = wf%Ngp - 1
	endif


	cr_tmp = ddot(Ngp, wf%re(lb,1), 1, bState%vector(lb), 1)
	ci_tmp = ddot(Ngp, wf%im(lb,1), 1, bState%vector(lb), 1)
	call MPI_Allreduce(cr_tmp, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
	call MPI_Allreduce(ci_tmp, ci, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)

	ortho=.FALSE.
	k=0
	do while(.not.ortho)
		k = k+1
		call daxpy(wf%Ngp, -cr/bState%norm, bState%vector(1), 1, wf%re(1,1), 1)
		call daxpy(wf%Ngp, -ci/bState%norm, bState%vector(1), 1, wf%im(1,1), 1)
		cr_tmp = ddot(Ngp, wf%re(lb,1), 1, bState%vector(lb), 1)
		ci_tmp = ddot(Ngp, wf%im(lb,1), 1, bState%vector(lb), 1)
		call MPI_Allreduce(cr_tmp, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		call MPI_Allreduce(ci_tmp, ci, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
		if(cdabs(dcmplx(cr,ci)) <= prec) ortho = .TRUE.
		if(k>=10) ortho = .TRUE.		! escape
	enddo


END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_create_output_h5(wf, fileName, Ntime)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf	! this is the whole wf, not the chunk of the mpi processes
	CHARACTER(LEN=*), INTENT(IN) :: fileName
!	TYPE(fedvr_grid), INTENT(IN) :: grid
	INTEGER(I4B), INTENT(IN) :: Ntime
	CHARACTER(LEN=3), PARAMETER :: dset_psi_name  = "Psi"
	CHARACTER(LEN=4), PARAMETER :: dset_grid_name = "Grid"
	CHARACTER(LEN=4), PARAMETER :: dset_norm1_name = "Norm"
	CHARACTER(LEN=8), PARAMETER :: dset_norm2_name = "Norm_tag"
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_psi_id, dset_grid_id, attr_id, dset_norm_id
	INTEGER(HID_T) :: dspace_psi_id, dspace_grid_id, aspace_id, dspace_norm_id
	INTEGER(HSIZE_T), DIMENSION(4) :: psi_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: grid_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: norm_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: attr_dims
	INTEGER(I4B) :: error
	INTEGER(I4B) :: Ngp, Ntag

	Ngp = wf%Ngp
	Ntag= wf%Ntag

	call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

!	#=======================#
!	|	Write grid data		|	!	(..,1) -> xx; (..,2) -> ww
!	#=======================#
	grid_dims = (/Ngp, 2/)
	call h5screate_simple_f(2, grid_dims, dspace_grid_id, error)
	call h5dcreate_f(file_id, dset_grid_name, H5T_NATIVE_DOUBLE, dspace_grid_id, dset_grid_id, error)

	call h5dclose_f(dset_grid_id, error)
	call h5sclose_f(dspace_grid_id, error)

!	#=======================#
!	|	Write WF data		|
!	#=======================#
	psi_dims=(/Ngp,Ntag,2,Ntime/)	!	(..,..,1,..) -> RE(psi); (..,..,2,..) -> IM(psi)
	
	call h5screate_simple_f(4, psi_dims, dspace_psi_id, error)
	call h5dcreate_f(file_id, dset_psi_name, H5T_NATIVE_DOUBLE, dspace_psi_id, dset_psi_id, error)


	norm_dims=(/2,Ntime/)	!	0 -> norm, 1 -> norm_evol
	call h5screate_simple_f(2, norm_dims, dspace_norm_id, error)
	call h5dcreate_f(file_id, dset_norm1_name, H5T_NATIVE_DOUBLE, dspace_norm_id, dset_norm_id, error)
	call h5dclose_f(dset_norm_id, error)
	call h5sclose_f(dspace_norm_id, error)

	norm_dims=(/Ntag,Ntime/)
	call h5screate_simple_f(2, norm_dims, dspace_norm_id, error)
	call h5dcreate_f(file_id, dset_norm2_name, H5T_NATIVE_DOUBLE, dspace_norm_id, dset_norm_id, error)
	call h5dclose_f(dset_norm_id, error)
	call h5sclose_f(dspace_norm_id, error)


!	#=======================#
!	|	Write WF attributes	|
!	#=======================#
	attr_dims=(/0/)
	call h5screate_f(H5S_SCALAR_F, aspace_id, error)
	call h5acreate_f(dset_psi_id, 'Nfun', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Nfun, attr_dims, error)
	call h5aclose_f(attr_id, error)

	call h5acreate_f(dset_psi_id, 'Nfe', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Nfe, attr_dims, error)
	call h5aclose_f(attr_id, error)

	call h5acreate_f(dset_psi_id, 'Ngp', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Ngp, attr_dims, error)
	call h5aclose_f(attr_id, error)

	call h5acreate_f(dset_psi_id, 'Ntag', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Ntag, attr_dims, error)
	call h5aclose_f(attr_id, error)


	call h5sclose_f(aspace_id, error)
	call h5dclose_f(dset_psi_id, error)
	call h5sclose_f(dspace_psi_id, error)

	call h5fclose_f(file_id, error)



END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE wf_dump_output_h5(wf, fileName, grid, idx_time)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	TYPE(fedvr_grid), INTENT(IN) :: grid
	INTEGER(I4B), INTENT(IN) :: idx_time
	CHARACTER(LEN=3), PARAMETER :: dset_psi_name  = "Psi"
	CHARACTER(LEN=4), PARAMETER :: dset_grid_name = "Grid"
	CHARACTER(LEN=4), PARAMETER :: dset_norm1_name = "Norm"
	CHARACTER(LEN=8), PARAMETER :: dset_norm2_name = "Norm_tag"
	INTEGER(HID_T) :: plist_id
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_psi_id, dset_grid_id, dset_norm_id
	INTEGER(HID_T) :: dspace_psi_id, dspace_grid_id, dspace_norm_id
	INTEGER(HID_T) :: mspace_psi_id, mspace_grid_id, mspace_norm_id
	INTEGER(HSIZE_T), DIMENSION(1) :: grid_slab_dims, mspace_grid_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: norm_slab_dims, mspace_norm_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: psi_slab_dims,  mspace_psi_dims
	INTEGER(HSIZE_T), DIMENSION(1:2) :: count_g, offset_g, stride_g, block_g
	INTEGER(HSIZE_T), DIMENSION(1:4) :: count_p, offset_p, stride_p, block_p
	INTEGER(I4B) :: error
	INTEGER(I4B) :: Ngp, Ntag, tag
	INTEGER(I4B) :: il, iu

	call wf_mirror(wf, wf_temp)
	if(wf%w_inc)then
		do tag=wf%tag_range(1),wf%tag_range(2)
			wf_temp%re(:,tag) = wf%re(:,tag) / dsqrt(grid%ww(:))
			wf_temp%im(:,tag) = wf%im(:,tag) / dsqrt(grid%ww(:))
		enddo
	endif

	if(my_id == root_id)then
		Ngp = wf%Ngp
		Ntag= wf%Ntag
		il = 1
		iu = wf%Ngp
	else
		Ngp = wf%Ngp - 1
		Ntag= wf%Ntag
		il = 2
		iu = wf%Ngp
	endif

	call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
	call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
	call h5pclose_f(plist_id, error)

!	#=======================#
!	|	Write grid data		|	!	(..,1) -> xx; (..,2) -> ww
!	#=======================#
	if(idx_time==0)then
		grid_slab_dims = (/Ngp/)
		mspace_grid_dims=(/Ngp/)
		call h5screate_simple_f(1, mspace_grid_dims, mspace_grid_id, error)
		call h5dopen_f(file_id, dset_grid_name, dset_grid_id, error)
		call h5dget_space_f(dset_grid_id, dspace_grid_id, error)

		call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
		call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

		if(my_id == root_id)then
			count_g  = (/1,1/)
			offset_g = (/wf%gp_range(1)-1,0/)
			stride_g = (/1,1/)
			block_g  = (/Ngp,1/)
		else
			count_g  = (/1,1/)
			offset_g = (/wf%gp_range(1),0/)
			stride_g = (/1,1/)
			block_g  = (/Ngp,1/)
		endif
!		write(6,*) my_id, 'offset = ', offset_g
!		write(6,*) my_id, 'block  = ', block_g

		call h5sselect_hyperslab_f(dspace_grid_id, H5S_SELECT_SET_F, offset_g, count_g, error, stride_g, block_g)
		call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid%xx(il:iu), grid_slab_dims, error, mspace_grid_id, dspace_grid_id, xfer_prp = plist_id)

		if(my_id == root_id)then
			count_g  = (/1,1/)
			offset_g = (/wf%gp_range(1)-1,1/)
			stride_g = (/1,1/)
			block_g  = (/Ngp,1/)
		else
			count_g  = (/1,1/)
			offset_g = (/wf%gp_range(1),1/)
			stride_g = (/1,1/)
			block_g  = (/Ngp,1/)
		endif
		call h5sselect_hyperslab_f(dspace_grid_id, H5S_SELECT_SET_F, offset_g, count_g, error, stride_g, block_g)
		call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid%ww(il:iu), grid_slab_dims, error, mspace_grid_id, dspace_grid_id, xfer_prp = plist_id)

		call h5pclose_f(plist_id, error)
		call h5dclose_f(dset_grid_id, error)
		call h5sclose_f(dspace_grid_id, error)
		call h5sclose_f(mspace_grid_id, error)
	endif

!	#=======================#
!	|	Write WF data		|
!	#=======================#
	psi_slab_dims=(/Ngp,Ntag/)
	mspace_psi_dims=(/Ngp,Ntag/)
	
	call h5screate_simple_f(2, mspace_psi_dims, mspace_psi_id, error)
	call h5dopen_f(file_id, dset_psi_name, dset_psi_id, error)
	call h5dget_space_f(dset_psi_id, dspace_psi_id, error)	! (/Ngp,Ntag,2,Ntime/)

	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

	if(my_id == root_id)then
		count_p  = (/1,1,1,1/)
		offset_p = (/wf%gp_range(1)-1,0,0,idx_time/)
		stride_p = (/1,1,1,1/)
		block_p  = (/Ngp,Ntag,1,1/)
	else
		count_p  = (/1,1,1,1/)
		offset_p = (/wf%gp_range(1),0,0,idx_time/)
		stride_p = (/1,1,1,1/)
		block_p  = (/Ngp,Ntag,1,1/)
	endif
!	write(6,*) my_id, 'offset = ', offset_g
!	write(6,*) my_id, 'block  = ', block_g

	call h5sselect_hyperslab_f(dspace_psi_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_psi_id, H5T_NATIVE_DOUBLE, wf_temp%re(il:iu, wf%tag_range(1):wf%tag_range(2)), psi_slab_dims, error, mspace_psi_id, dspace_psi_id, xfer_prp = plist_id)

	if(my_id == root_id)then
		count_p  = (/1,1,1,1/)
		offset_p = (/wf%gp_range(1)-1,0,1,idx_time/)
		stride_p = (/1,1,1,1/)
		block_p  = (/Ngp,Ntag,1,1/)
	else
		count_p  = (/1,1,1,1/)
		offset_p = (/wf%gp_range(1),0,1,idx_time/)
		stride_p = (/1,1,1,1/)
		block_p  = (/Ngp,Ntag,1,1/)
	endif
!	write(6,*) my_id, 'offset = ', offset_g
!	write(6,*) my_id, 'block  = ', block_g

	call h5sselect_hyperslab_f(dspace_psi_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_psi_id, H5T_NATIVE_DOUBLE, wf_temp%im(il:iu, wf%tag_range(1):wf%tag_range(2)), psi_slab_dims, error, mspace_psi_id, dspace_psi_id, xfer_prp = plist_id)

	call h5pclose_f(plist_id, error)
	call h5dclose_f(dset_psi_id, error)
	call h5sclose_f(dspace_psi_id, error)
	call h5sclose_f(mspace_psi_id, error)
	

	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

!	if(my_id == root_id)then
		norm_slab_dims = (/2/)
		mspace_norm_dims=(/2/)
		call h5screate_simple_f(1, mspace_norm_dims, mspace_norm_id, error)
		call h5dopen_f(file_id, dset_norm1_name, dset_norm_id, error)
		call h5dget_space_f(dset_norm_id, dspace_norm_id, error)

		count_g  = (/1,1/)
		offset_g = (/0,idx_time/)
		stride_g = (/1,1/)
		block_g  = (/2,1/)
		call h5sselect_hyperslab_f(dspace_norm_id, H5S_SELECT_SET_F, offset_g, count_g, error, stride_g, block_g)
		call h5dwrite_f(dset_norm_id, H5T_NATIVE_DOUBLE, (/wf%norm, wf%norm_evol/), norm_slab_dims, error, mspace_norm_id, dspace_norm_id, xfer_prp = plist_id)

		call h5dclose_f(dset_norm_id, error)
		call h5sclose_f(dspace_norm_id, error)
		call h5sclose_f(mspace_norm_id, error)


		norm_slab_dims = (/Ntag/)
		mspace_norm_dims=(/Ntag/)
		call h5screate_simple_f(1, mspace_norm_dims, mspace_norm_id, error)
		call h5dopen_f(file_id, dset_norm2_name, dset_norm_id, error)
		call h5dget_space_f(dset_norm_id, dspace_norm_id, error)

		count_g  = (/1,1/)
		offset_g = (/0,idx_time/)
		stride_g = (/1,1/)
		block_g  = (/Ntag,1/)
		call h5sselect_hyperslab_f(dspace_norm_id, H5S_SELECT_SET_F, offset_g, count_g, error, stride_g, block_g)
		call h5dwrite_f(dset_norm_id, H5T_NATIVE_DOUBLE, wf%tag_norm, norm_slab_dims, error, mspace_norm_id, dspace_norm_id, xfer_prp = plist_id)


		call h5dclose_f(dset_norm_id, error)
		call h5sclose_f(dspace_norm_id, error)
		call h5sclose_f(mspace_norm_id, error)
!	endif
	call h5pclose_f(plist_id, error)

	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------

SUBROUTINE wf_save_h5_sp(wf, fileName, grid)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	TYPE(fedvr_grid), INTENT(IN) :: grid
	CHARACTER(LEN=3), PARAMETER :: dset_psi_name  = "Psi"
	CHARACTER(LEN=4), PARAMETER :: dset_grid_name = "Grid"
	CHARACTER(LEN=4), PARAMETER :: dset_norm1_name = "Norm"
	CHARACTER(LEN=8), PARAMETER :: dset_norm2_name = "Norm_tag"
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_psi_id, dset_grid_id, attr_id, dset_norm_id
	INTEGER(HID_T) :: dspace_psi_id, dspace_grid_id, aspace_id, dspace_norm_id
	INTEGER(HID_T) :: mspace_psi_id, mspace_grid_id
	INTEGER(HSIZE_T), DIMENSION(3) :: psi_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: grid_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: grid_slab_dims, mspace_grid_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: psi_slab_dims,  mspace_psi_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: norm_dims
	INTEGER(HSIZE_T), DIMENSION(1:2) :: count_g, offset_g, stride_g, block_g
	INTEGER(HSIZE_T), DIMENSION(1:3) :: count_p, offset_p, stride_p, block_p
	INTEGER(HSIZE_T), DIMENSION(1) :: attr_dims
	INTEGER(I4B) :: error
	INTEGER(I4B) :: Ngp, Ntag, tag

	Ngp = wf%Ngp
	Ntag= wf%Ntag

	call wf_mirror(wf, wf_temp, 'alloc')
	if(wf%w_inc)then
		do tag=wf%tag_range(1),wf%tag_range(2)
			wf_temp%re(:,tag) = wf%re(:,tag) / dsqrt(grid%ww(:))
			wf_temp%im(:,tag) = wf%im(:,tag) / dsqrt(grid%ww(:))
		enddo
	endif

	call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

!	#=======================#
!	|	Write grid data		|	!	(..,1) -> xx; (..,2) -> ww
!	#=======================#
	grid_dims = (/Ngp, 2/)
	grid_slab_dims = (/Ngp/)
	mspace_grid_dims=(/Ngp/)
	call h5screate_simple_f(2, grid_dims, dspace_grid_id, error)
	call h5screate_simple_f(1, mspace_grid_dims, mspace_grid_id, error)
	call h5dcreate_f(file_id, dset_grid_name, H5T_NATIVE_DOUBLE, dspace_grid_id, dset_grid_id, error)

	count_g  = (/1,1/)
	offset_g = (/0,0/)
	stride_g = (/1,1/)
	block_g  = (/Ngp,1/)
	call h5sselect_hyperslab_f(dspace_grid_id, H5S_SELECT_SET_F, offset_g, count_g, error, stride_g, block_g)
	call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid%xx, grid_slab_dims, error, mspace_grid_id, dspace_grid_id)

	count_g  = (/1,1/)
	offset_g = (/0,1/)
	stride_g = (/1,1/)
	block_g  = (/Ngp,1/)
	call h5sselect_hyperslab_f(dspace_grid_id, H5S_SELECT_SET_F, offset_g, count_g, error, stride_g, block_g)
	call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid%ww, grid_slab_dims, error, mspace_grid_id, dspace_grid_id)

	call h5dclose_f(dset_grid_id, error)
	call h5sclose_f(dspace_grid_id, error)
	call h5sclose_f(mspace_grid_id, error)

!	#=======================#
!	|	Write WF data		|
!	#=======================#
	psi_dims=(/Ngp,Ntag,2/)	!	(..,..,1) -> RE(psi); (..,..,2) -> IM(psi)
	psi_slab_dims=(/Ngp,Ntag/)
	mspace_psi_dims=(/Ngp,Ntag/)
	
	call h5screate_simple_f(3, psi_dims, dspace_psi_id, error)
	call h5screate_simple_f(2, mspace_psi_dims, mspace_psi_id, error)
	call h5dcreate_f(file_id, dset_psi_name, H5T_NATIVE_DOUBLE, dspace_psi_id, dset_psi_id, error)

	count_p  = (/1,1,1/)
	offset_p = (/0,0,0/)
	stride_p = (/1,1,1/)
	block_p  = (/Ngp,Ntag,1/)	
	call h5sselect_hyperslab_f(dspace_psi_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_psi_id, H5T_NATIVE_DOUBLE, wf_temp%re, psi_slab_dims, error, mspace_psi_id, dspace_psi_id)

	count_p  = (/1,1,1/)
	offset_p = (/0,0,1/)
	stride_p = (/1,1,1/)
	block_p  = (/Ngp,Ntag,1/)	
	call h5sselect_hyperslab_f(dspace_psi_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_psi_id, H5T_NATIVE_DOUBLE, wf_temp%im, psi_slab_dims, error, mspace_psi_id, dspace_psi_id)


	norm_dims=(/2/)	!	0 -> norm, 1 -> norm_evol
	call h5screate_simple_f(1, norm_dims, dspace_norm_id, error)
	call h5dcreate_f(file_id, dset_norm1_name, H5T_NATIVE_DOUBLE, dspace_norm_id, dset_norm_id, error)
	call h5dwrite_f(dset_norm_id, H5T_NATIVE_DOUBLE, (/wf%norm, wf%norm_evol/), norm_dims, error)
	call h5dclose_f(dset_norm_id, error)
	call h5sclose_f(dspace_norm_id, error)

	norm_dims=(/Ntag/)
	call h5screate_simple_f(1, norm_dims, dspace_norm_id, error)
	call h5dcreate_f(file_id, dset_norm2_name, H5T_NATIVE_DOUBLE, dspace_norm_id, dset_norm_id, error)
	call h5dwrite_f(dset_norm_id, H5T_NATIVE_DOUBLE, wf%tag_norm, norm_dims, error)
	call h5dclose_f(dset_norm_id, error)
	call h5sclose_f(dspace_norm_id, error)

!	#=======================#
!	|	Write WF attributes	|
!	#=======================#
	attr_dims=(/0/)
	call h5screate_f(H5S_SCALAR_F, aspace_id, error)
	call h5acreate_f(dset_psi_id, 'Nfun', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Nfun, attr_dims, error)
	call h5aclose_f(attr_id, error)

	call h5acreate_f(dset_psi_id, 'Nfe', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Nfe, attr_dims, error)
	call h5aclose_f(attr_id, error)

	call h5acreate_f(dset_psi_id, 'Ngp', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Ngp, attr_dims, error)
	call h5aclose_f(attr_id, error)

	call h5acreate_f(dset_psi_id, 'Ntag', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%Ntag, attr_dims, error)
	call h5aclose_f(attr_id, error)

!	fe_range, gp_range, tag_range and fe_gp_range are not saved 
!	attr_dims=(/2/)
!	call h5screate_simple_f(1, attr_dims, aspace_id, error)
!	call h5acreate_f(dset_psi_id, 'fe_range', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
!	call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, wf%fe_range, attr_dims, error)
!	call h5aclose_f(attr_id, error)


	call h5dclose_f(dset_psi_id, error)
	call h5sclose_f(dspace_psi_id, error)

	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------
!===============================================================================
!*******************************************************************************
SUBROUTINE mpi_distribute_local_wf(wf, my_wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(IN)    :: wf
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: my_wf

	if(my_id == root_id)then
		my_wf%norm      = wf%norm
		my_wf%norm_prev = wf%norm_prev
		my_wf%norm_loc  = wf%norm_loc
		my_wf%norm_evol = wf%norm_evol
		my_wf%tag_norm(:) = wf%tag_norm(Idx_tag(1,my_id):Idx_tag(2,my_id))
		my_wf%re(1:my_wf%Ngp, 1:wf%Ntag) = wf%re(Idx_gp(1, my_id):Idx_gp(2, my_id), Idx_tag(1, my_id):Idx_tag(2, my_id))
		my_wf%im(1:my_wf%Ngp, 1:wf%Ntag) = wf%im(Idx_gp(1, my_id):Idx_gp(2, my_id), Idx_tag(1, my_id):Idx_tag(2, my_id))
		do mpi_id=1,mpi_size-1
			call MPI_Send(wf%norm,      1, MPI_DOUBLE_PRECISION, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(wf%norm_prev, 1, MPI_DOUBLE_PRECISION, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(wf%norm_loc,  1, MPI_DOUBLE_PRECISION, mpi_id, 2, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(wf%norm_evol, 1, MPI_DOUBLE_PRECISION, mpi_id, 3, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(wf%tag_norm(Idx_tag(1, my_id):Idx_tag(2, my_id)),                                  N_tag_arr(mpi_id),                  MPI_DOUBLE_PRECISION, mpi_id, 4, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(wf%re(Idx_gp(1, mpi_id):Idx_gp(2, mpi_id), Idx_tag(1, mpi_id):Idx_tag(2, mpi_id)), N_tag_arr(mpi_id)*N_gp_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 5, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(wf%im(Idx_gp(1, mpi_id):Idx_gp(2, mpi_id), Idx_tag(1, mpi_id):Idx_tag(2, mpi_id)), N_tag_arr(mpi_id)*N_gp_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 6, MPI_COMM_WORLD, mpi_err)
		enddo
	else
		call MPI_Recv(my_wf%norm,      1, MPI_DOUBLE_PRECISION, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_wf%norm_prev, 1, MPI_DOUBLE_PRECISION, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_wf%norm_loc,  1, MPI_DOUBLE_PRECISION, root_id, 2, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_wf%norm_evol, 1, MPI_DOUBLE_PRECISION, root_id, 3, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_wf%tag_norm(1:my_wf%Ntag),        my_wf%Ntag,           MPI_DOUBLE_PRECISION, root_id, 4, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_wf%re(1:my_wf%Ngp, 1:my_wf%Ntag), my_wf%Ntag*my_wf%Ngp, MPI_DOUBLE_PRECISION, root_id, 5, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_wf%im(1:my_wf%Ngp, 1:my_wf%Ntag), my_wf%Ntag*my_wf%Ngp, MPI_DOUBLE_PRECISION, root_id, 6, MPI_COMM_WORLD, rstatus, mpi_err)
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE mpi_distribute_local_vec(vec, my_vec)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN)    :: vec
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: my_vec

	if(my_id == root_id)then
		my_vec(my_gp(1):my_gp(2)) = vec(my_gp(1):my_gp(2))
		do mpi_id=1,mpi_size-1
!			write(6,*) my_id, mpi_id, Idx_gp(:,mpi_id)
			call MPI_Send(vec(Idx_gp(1, mpi_id):Idx_gp(2, mpi_id)), N_gp_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
		enddo
	else
!		write(6,*) my_id, my_gp
!		call MPI_Recv(my_vec(my_gp(1):my_gp(2)), NN_gp, MPI_DOUBLE_PRECISION, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_vec(1:NN_gp), NN_gp, MPI_DOUBLE_PRECISION, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE mpi_update_wf_at_bounds(wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	REAL(DP), DIMENSION(:), ALLOCATABLE :: re_buffer, im_buffer

	if(allocated(re_buffer))deallocate(re_buffer)
	if(allocated(im_buffer))deallocate(im_buffer)
	allocate(re_buffer(1:wf%Ntag), im_buffer(1:wf%Ntag))
	re_buffer(:) = 0.0d0
	im_buffer(:) = 0.0d0


	if(my_id /= root_id)then
		call MPI_Send(wf%re(1, 1:wf%Ntag), wf%Ntag, MPI_DOUBLE_PRECISION, my_id-1, 1, MPI_COMM_WORLD, mpi_err)
		call MPI_Send(wf%im(1, 1:wf%Ntag), wf%Ntag, MPI_DOUBLE_PRECISION, my_id-1, 2, MPI_COMM_WORLD, mpi_err)
	endif
	if(my_id /= mpi_size-1)then
		call MPI_Recv(re_buffer,           wf%Ntag, MPI_DOUBLE_PRECISION, my_id+1, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(im_buffer,           wf%Ntag, MPI_DOUBLE_PRECISION, my_id+1, 2, MPI_COMM_WORLD, rstatus, mpi_err)
	endif

	wf%re(wf%Ngp, :) = wf%re(wf%Ngp, :) + re_buffer(:)
	wf%im(wf%Ngp, :) = wf%im(wf%Ngp, :) + im_buffer(:)

	if(my_id /= mpi_size-1)then
		call MPI_Send(wf%re(wf%Ngp, 1:wf%Ntag), wf%Ntag, MPI_DOUBLE_PRECISION, my_id+1, 3, MPI_COMM_WORLD, mpi_err)
		call MPI_Send(wf%im(wf%Ngp, 1:wf%Ntag), wf%Ntag, MPI_DOUBLE_PRECISION, my_id+1, 4, MPI_COMM_WORLD, mpi_err)
	endif
	if(my_id /= root_id)then
		call MPI_Recv(wf%re(1,      1:wf%Ntag), wf%Ntag, MPI_DOUBLE_PRECISION, my_id-1, 3, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(wf%im(1,      1:wf%Ntag), wf%Ntag, MPI_DOUBLE_PRECISION, my_id-1, 4, MPI_COMM_WORLD, rstatus, mpi_err)
	endif

	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE mpi_update_vec_at_bounds(vec)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: vec
	REAL(DP) :: buffer

	buffer = 0.0d0


	if(my_id /= root_id)    call MPI_Send(vec(1), 1, MPI_DOUBLE_PRECISION, my_id-1, 1, MPI_COMM_WORLD, mpi_err)
	if(my_id /= mpi_size-1) call MPI_Recv(buffer, 1, MPI_DOUBLE_PRECISION, my_id+1, 1, MPI_COMM_WORLD, rstatus, mpi_err)

	vec(NN_gp) = vec(NN_gp) + buffer

	if(my_id /= mpi_size-1) call MPI_Send(vec(NN_gp), 1, MPI_DOUBLE_PRECISION, my_id+1, 2, MPI_COMM_WORLD, mpi_err)
	if(my_id /= root_id)    call MPI_Recv(vec(1    ), 1, MPI_DOUBLE_PRECISION, my_id-1, 2, MPI_COMM_WORLD, rstatus, mpi_err)

	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE




































