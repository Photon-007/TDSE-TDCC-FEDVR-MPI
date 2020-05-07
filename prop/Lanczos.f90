MODULE Lanczos
	USE nrtype
	USE gvars
	USE FEDVR_module
	USE TDCC_module
	USE TDCC_FEDVR
	USE Operator_typedef
	USE Operator_math
	USE eigen_calc
	USE MPI_pars
	USE mpi


!	TYPE Krylov_basis
!		TYPE(tdcc_fedvr_wf), DIMENSION(:), ALLOCATABLE :: q_vec
!	END TYPE

	TYPE(tdcc_fedvr_wf), DIMENSION(:),   ALLOCATABLE :: Q_vec					! Krylov basis (orthonormalised)
	REAL(DP),            DIMENSION(:,:), ALLOCATABLE :: Q_mat					! Hamiltonian in the Krylov base
	REAL(DP),            DIMENSION(:),   ALLOCATABLE :: Q_alpha
	REAL(DP),            DIMENSION(:),   ALLOCATABLE :: Q_beta
	COMPLEX(DPC),        DIMENSION(:),   ALLOCATABLE :: coeffs_curr
	COMPLEX(DPC),        DIMENSION(:),   ALLOCATABLE :: coeffs_prev

	REAL(DPC)    :: dt_max
	INTEGER(I4B) :: prop_steps_rt, Krylov_base_size

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE propagator_Re_Time(wf, tt, dt, Energy)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	REAL(DP), INTENT(INOUT) :: tt
	REAL(DP), INTENT(INOUT) :: dt
	REAL(DP), INTENT(INOUT) :: Energy
	REAL(DP) :: Norm_tt
	REAL(DP) :: delta_coeffs, scale_fact
	COMPLEX(DPC) :: dotp
	INTEGER(I4B) :: ii, jj
	INTEGER(I4B) :: GS_starting_qvec, GS_iter

	Q_mat(:,:) = 0.0d0
	Q_alpha(:) = 0.0d0
	Q_beta(:)  = 1.0d0
	Q_mat(:,:) = 0.0d0
	coeffs_curr(:) = CZERO
	coeffs_prev(:) = CZERO

!	write(6,*) 'wf%norm(Lanczos begin) = ', wf%norm, wf%norm_prev
	Norm_tt = wf%norm
	call wf_copy(wf, Q_vec(0))
	call wf_normalize(Q_vec(0))

	do ii = 1, N_Krylov			! Lanczos loop
		Krylov_base_size = ii

		if(memory_save)then
			Q_vec(ii) = ( OP_KinVr * Q_vec(ii-1) ) + ( OP_Vl     * Q_vec(ii-1) ) + ( OP_EField * Q_vec(ii-1) )			
		else
			Q_vec(ii) = ( OP_Ham0  * Q_vec(ii-1) ) + ( OP_EField * Q_vec(ii-1) )
		endif

		if(GS_all_qvecs)then
			GS_starting_qvec = 0
		else
			GS_starting_qvec = MAX(ii-2,0)
		endif

		do jj = ii-1, GS_starting_qvec, -1
			do GS_iter = 1, GS_refinement
				dotp = Q_vec(jj) * Q_vec(ii)
				if(jj==ii-1) Q_alpha(ii-1) = Q_alpha(ii-1) + dble(dotp)
				call wf_axpy(-dble(dotp), Q_vec(jj), Q_vec(ii))
				if((dabs(dble(dotp)) <= Epsilon_ortho).and.(GS_iter >= 2)) EXIT
			enddo
		enddo

		call wf_normalize(Q_vec(ii))
		Q_beta(ii) = Q_beta(ii) * dsqrt(Q_vec(ii)%norm_prev)

!		if(my_id==root_id)write(6,*)"alpha(",ii,")=",Q_alpha(ii-1), " | beta(",ii,")=",Q_beta(ii)
		Q_mat(ii,ii) = Q_alpha(ii-1)
		if(ii/=N_Krylov)then
			Q_mat(ii+1,ii) = Q_beta(ii)
			Q_mat(ii,ii+1) = Q_beta(ii)
		endif


		if(my_id==root_id)then		! in princilpe everyone could do the same and should end up with the same coeffs, but this way we are sure! (just in case the Q_mat diagonalization has some fluctuations on the processes)
			call get_lanczos_coeffs_RealTime(Q_mat, ii, dt, coeffs_curr)
!			write(6,*)coeffs_prev
!			write(6,*)coeffs_curr
			delta_coeffs = sum(cdabs(coeffs_curr(:) - coeffs_prev(:))**2)
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(delta_coeffs,                     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(coeffs_curr(0:N_Krylov-1), N_Krylov, MPI_DOUBLE_COMPLEX,   root_id, MPI_COMM_WORLD, mpi_err)		
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

		coeffs_prev(:) = coeffs_curr(:)
		if(delta_coeffs < Epsilon_conv) EXIT

	enddo	! Lanczos loop
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	Energy = Q_mat(1,1)
!	if(my_id==root_id)write(6,*)Q_mat(1,1)


	if((delta_coeffs/Epsilon_conv) < 0.8d0)then
		tt = tt + dt
		dt = MIN(1.1d0*dt, dt_max)
!		write(6,*)"delta_coeffs:",delta_coeffs," ---> Increase dt:",dt
!		pause
	else
		if(my_id == root_id)then
			do while(delta_coeffs > Epsilon_conv)
				coeffs_curr(:) = CZERO
				coeffs_prev(:) = CZERO
				scale_fact = 0.9d0*(Epsilon_conv/delta_coeffs)**(1.d0/real(N_Krylov))
				dt = scale_fact * dt
				call get_lanczos_coeffs_RealTime(Q_mat, N_Krylov-1, dt, coeffs_prev)
				call get_lanczos_coeffs_RealTime(Q_mat, N_Krylov,   dt, coeffs_curr)
				delta_coeffs = sum(cdabs(coeffs_curr(:) - coeffs_prev(:))**2)
!				write(6,*)"delta_coeffs:",delta_coeffs," ---> Decrease dt:",dt, "(scale_fact=",scale_fact,")"
			enddo
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(dt,                               1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(delta_coeffs,                     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(coeffs_curr(0:N_Krylov-1), N_Krylov, MPI_DOUBLE_COMPLEX,   root_id, MPI_COMM_WORLD, mpi_err)		
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		tt = tt + dt
	endif

!	%=======================================================%
!	|	Build psi(t+dt) as:									|
!	|	psi(t+dt) = coeffs(0)*q(0) + coeffs(1)*q(1) + ...	|
!	%=======================================================%
	call wf_clear(wf)
!	write(6,*)'rescale wf to:', wf%norm_prev
	coeffs_curr(:) = dsqrt(Norm_tt) * coeffs_curr(:)			! no-mather which one I choose, the norm starts to fluctuate (miledly)
!	coeffs_curr(:) = dsqrt(wf%norm_prev) * coeffs_curr(:)

	do ii = 0, N_Krylov-1
		call wf_axpy(coeffs_curr(ii), Q_vec(ii), wf)
	enddo
!	call wf_scale(dsqrt(Norm_tt), wf)
	call wf_norm(wf)
!	write(6,*) 'wf%norm(Lanczos end) = ', wf%norm, wf%norm_prev

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE propagator_Im_Time(wf, tt, dt, Energy)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	REAL(DP), INTENT(INOUT) :: tt
	REAL(DP), INTENT(INOUT) :: dt
	REAL(DP), INTENT(OUT)   :: Energy
!	REAL(DP) :: Norm_tt
	REAL(DP) :: delta_coeffs !, scale_fact
	COMPLEX(DPC) :: dotp
	INTEGER(I4B) :: ii, jj
	INTEGER(I4B) :: GS_starting_qvec, GS_iter

	Q_mat(:,:) = 0.0d0
	Q_alpha(:) = 0.0d0
	Q_beta(:)  = 1.0d0
	Q_mat(:,:) = 0.0d0
	coeffs_curr(:) = CZERO
	coeffs_prev(:) = CZERO

!	Norm_tt = wf%norm
	call wf_copy(wf, Q_vec(0))
	call wf_normalize(Q_vec(0))

	do ii = 1, N_Krylov			! Lanczos loop
		Krylov_base_size = ii

!		if(memory_save)then
!			Q_vec(ii) = ( OP_KinVr * Q_vec(ii-1) ) + ( OP_Vl * Q_vec(ii-1) )
!		else
			Q_vec(ii) =   Ham0_init  * Q_vec(ii-1) 
!		endif

		if(GS_all_qvecs)then
			GS_starting_qvec = 0
		else
			GS_starting_qvec = MAX(ii-2,0)
		endif

		do jj = ii-1, GS_starting_qvec, -1
			do GS_iter = 1, GS_refinement
				dotp = Q_vec(jj) * Q_vec(ii)
				if(jj==ii-1) Q_alpha(ii-1) = Q_alpha(ii-1) + dble(dotp)
				call wf_axpy(-dble(dotp), Q_vec(jj), Q_vec(ii))
				if((dabs(dble(dotp)) <= Epsilon_ortho).and.(GS_iter >= 2)) EXIT
			enddo
		enddo


		call wf_normalize(Q_vec(ii))
		Q_beta(ii) = Q_beta(ii) * dsqrt(Q_vec(ii)%norm_prev)


!		write(6,*)"alpha(",ii,")=",Q_alpha(ii-1), " | beta(",ii,")=",Q_beta(ii)
		Q_mat(ii,ii) = Q_alpha(ii-1)
		if(ii/=N_Krylov)then
			Q_mat(ii+1,ii) = Q_beta(ii)
			Q_mat(ii,ii+1) = Q_beta(ii)
		endif


		if(my_id==root_id)then		! in princilpe everyone could do the same and should end up with the same coeffs, but this way we are sure! (just in case the Q_mat diagonalization has some fluctuations on the processes)
			call get_lanczos_coeffs_ImagTime(Q_mat, ii, dt, coeffs_curr)
!			write(6,*)coeffs_prev
!			write(6,*)coeffs_curr
			delta_coeffs = sum(cdabs(coeffs_curr(:) - coeffs_prev(:))**2)
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(delta_coeffs,                     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_Bcast(coeffs_curr(0:N_Krylov-1), N_Krylov, MPI_DOUBLE_COMPLEX,   root_id, MPI_COMM_WORLD, mpi_err)		
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

!		write(6,*)ii, dt, delta_coeffs
		coeffs_prev(:) = coeffs_curr(:)
		if(delta_coeffs < Epsilon_conv_ITP) EXIT

	enddo	! Lanczos loop
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	Energy = Q_mat(1,1)


	if((delta_coeffs/Epsilon_conv_ITP) < 0.8d0)then
		tt = tt + dt
		dt = MIN(1.1d0*dt, dt_max)
!		write(6,*)"delta_coeffs:",delta_coeffs," ---> Increase dt:",dt
!		pause
	else
!		if(my_id == root_id)then
!			do while(delta_coeffs > Epsilon_conv_ITP)
!				coeffs_curr(:) = CZERO
!				coeffs_prev(:) = CZERO
!				scale_fact = 0.9d0*(Epsilon_conv_ITP/delta_coeffs)**(1.d0/real(N_Krylov))
!				dt = scale_fact * dt
!				call get_lanczos_coeffs_ImagTime(Q_mat, N_Krylov-1, dt, coeffs_prev)
!				call get_lanczos_coeffs_ImagTime(Q_mat, N_Krylov,   dt, coeffs_curr)
!				delta_coeffs = sum(cdabs(coeffs_curr(:) - coeffs_prev(:))**2)
!				write(6,*)"delta_coeffs:",delta_coeffs," ---> Decrease dt:",dt, "(scale_fact=",scale_fact,")"
!				if(ii==5)pause
!			enddo
!		endif
!		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
!		call MPI_Bcast(dt,                               1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
!		call MPI_Bcast(delta_coeffs,                     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
!		call MPI_Bcast(coeffs_curr(0:N_Krylov-1), N_Krylov, MPI_DOUBLE_COMPLEX,   root_id, MPI_COMM_WORLD, mpi_err)		
!		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		tt = tt + dt
	endif

!	%=======================================================%
!	|	Build psi(t+dt) as:									|
!	|	psi(t+dt) = coeffs(0)*q(0) + coeffs(1)*q(1) + ...	|
!	%=======================================================%
	call wf_clear(wf)
!	coeffs_curr(:) = dsqrt(Norm_tt) * coeffs_curr(:)

	do ii = 0, N_Krylov-1
		call wf_axpy(dble(coeffs_curr(ii)), Q_vec(ii), wf)
	enddo

	call wf_normalize(wf)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE get_lanczos_coeffs_RealTime(HH, n, dt, coeff)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: HH
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), INTENT(IN) :: dt
	COMPLEX(DPC), DIMENSION(0:), INTENT(INOUT) :: coeff
	TYPE(eigenpair_obj) :: Eigenp
	REAL(DP) :: re, im
	INTEGER(I4B) :: i, j


	call allocate_eigenpair(Eigenp, n)
	forall (i=1:n, j=1:n) Eigenp%matrix(i,j) = HH(i,j)

	call EigenSolve_ge(Eigenp)

	do i=1,n
		re = sum(Eigenp%l_eigvecs(i,:) * dcos(-dt * Eigenp%eigvals(:)) * Eigenp%r_eigvecs(1,:))
		im = sum(Eigenp%l_eigvecs(i,:) * dsin(-dt * Eigenp%eigvals(:)) * Eigenp%r_eigvecs(1,:))
		coeff(i-1) = dcmplx(re, im)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE get_lanczos_coeffs_ImagTime(HH, n, dt, coeff)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: HH
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), INTENT(IN) :: dt
	COMPLEX(DPC), DIMENSION(0:), INTENT(INOUT) :: coeff
	TYPE(eigenpair_obj) :: Eigenp
!	REAL(DP) :: dtt
	INTEGER(I4B) :: i, j


	call allocate_eigenpair(Eigenp, n)
	forall (i=1:n, j=1:n) Eigenp%matrix(i,j) = HH(i,j)

	call EigenSolve_ge(Eigenp)

	do i=1,n
		coeff(i-1) = dcmplx(sum(Eigenp%l_eigvecs(i,:) * dexp(-dt * Eigenp%eigvals(:)) * Eigenp%r_eigvecs(1,:)), 0.d0)
	enddo

END SUBROUTINE





END MODULE Lanczos

















