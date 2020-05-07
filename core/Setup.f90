MODULE Setup
	USE nrtype
	USE gvars
	USE BoundStates
	USE MPI_pars
	USE FEDVR_module
	USE TDCC_module
	USE TDCC_FEDVR
	USE Operator_typedef
	USE Operator_math
	USE FEDVR_derivatives
	USE FEDVR_potentials
	USE LASER_potentials
	USE eigen_calc
	USE tProp
	USE mpi

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE set_system_operators(grid)
	IMPLICIT NONE
	TYPE(FEDVR_grid), INTENT(IN) :: grid
	TYPE(fedvr_op_full) :: deriv2		! local operators destroyed after subroutine terminates
!	TYPE(fedvr_op_diag) :: Vx, dVx
	INTEGER(I4B) :: tag, ll, mm
	CHARACTER(LEN=5) :: tag_str

	write(proc_tag,'(I5.5)')my_id

	call Op_Sym_create(OP_Sym, grid)

	if(memory_save)then
		if(allocated(OP_Vl)) deallocate(OP_Vl)
		if(allocated(OP_dVl))deallocate(OP_dVl)
		allocate(OP_Vl(my_tag(1):my_tag(2)), OP_dVl(my_tag(1):my_tag(2)))

		call set_2nd_deriv_op(deriv2, grid)
		call set_core_potential_ops(OP_Vr, OP_dVr, grid)

		OP_KinVr = (mHalf * deriv2) + OP_Vr
!		OP_KinVr = (mHalf * deriv2) +    Vr

		do ii_fe=1,OP_Ham0(tag)%Nfe
			OP_KinVr%fe_region(ii_fe)%matrix = OP_Sym%fe_region(ii_fe)%matrix * OP_KinVr%fe_region(ii_fe)%matrix
		enddo

		do tag = my_tag(1),my_tag(2)
			call TDCC_get_lm_from_tag(SphericalHarmonicRepresentation, tag, ll, mm)
			call set_centrifugal_ops(OP_Vl(tag), OP_dVl(tag), grid, ll)
		enddo

		call EE_op_create(OP_EField, grid, my_tag, SphericalHarmonicRepresentation, LaserPulse)
		call EE_op_allocate(OP_EField)

		call set_abspot(Abspot_op, grid)
		call Op_create(Maskfun_op, grid, 0, 0)
		call Op_allocate(Maskfun_op)

	else
		if(allocated(OP_Ham0))deallocate(OP_Ham0)
		if(allocated(OP_Vl))  deallocate(OP_Vl)
		if(allocated(OP_dVl)) deallocate(OP_dVl)
		allocate(OP_Ham0(my_tag(1):my_tag(2)), OP_Vl(my_tag(1):my_tag(2)), OP_dVl(my_tag(1):my_tag(2)))

		call set_2nd_deriv_op(deriv2, grid)
		call set_core_potential_ops(OP_Vr, OP_dVr, grid)

		do tag = my_tag(1),my_tag(2)
			call TDCC_get_lm_from_tag(SphericalHarmonicRepresentation, tag, ll, mm)
			call set_centrifugal_ops(OP_Vl(tag), OP_dVl(tag), grid, ll)

			call Op_create(OP_Ham0(tag), grid, ll, mm)
!			call Op_allocate(OP_Ham0(tag))
			OP_Ham0(tag) = (( mHalf * deriv2 ) + OP_Vr ) + OP_Vl(tag)
!			OP_Ham0(tag) = deriv2
!			OP_Ham0(tag) = OP_Ham0(tag) + OP_Vr

			do ii_fe=1,OP_Ham0(tag)%Nfe
				OP_Ham0(tag)%fe_region(ii_fe)%matrix = OP_Sym%fe_region(ii_fe)%matrix * OP_Ham0(tag)%fe_region(ii_fe)%matrix
			enddo
		enddo

		call EE_op_create(OP_EField, grid, my_tag, SphericalHarmonicRepresentation, LaserPulse)
		call EE_op_allocate(OP_EField)

		call set_abspot(Abspot_op, grid)
		call Op_create(Maskfun_op, grid, 0, 0)
		call Op_allocate(Maskfun_op)
!		call OP_tofile(Abspot_op, 'Abspot_'//proc_tag, 'f12.6')
		do tag=1,tagmax
			write(tag_str,'(I5.5)')tag
			CALL OP_save(OP_Ham0(tag),'OP_H0_'//tag_str)		! have no idea why it hangs on the final worker!!!
		enddo
		CALL OP_save(OP_Vr,'OP_Vr')
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE get_WF_init_Diag(grid, wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(OUT) :: wf
	TYPE(fedvr_grid),    INTENT(IN)  :: grid
	TYPE(fedvr_grid)    :: grid_temp
	TYPE(tdcc_fedvr_wf) :: wf_temp
	TYPE(fedvr_op_full) :: deriv2, Sym
	TYPE(fedvr_op_diag) :: Vx, dVx
	TYPE(fedvr_op_diag) :: Vl, dVl
	TYPE(eigenpair_obj) :: Eigenp
	TYPE(BoundState) :: bState
	INTEGER(I4B) :: n_th, nn
	INTEGER(I4B) :: ll, mm, tag, ii
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: Idx_tag_diag
	INTEGER(I4B), DIMENSION(:)  , ALLOCATABLE :: my_tag_range
	REAL(DP), DIMENSION(:), ALLOCATABLE :: EigenVals, EigenVec
	CHARACTER(LEN=5) :: tag_str


	tag_init = TDCC_get_tag_from_lm(SphericalHarmonicRepresentation, l_init, m_init)
	n_th = n_init - l_init

	call wf_create(wf, grid, my_tag)
	call wf_allocate(wf)


	if(allocated(my_tag_range))deallocate(my_tag_range)
	allocate(my_tag_range(2))
	if(my_id == root_id)then
		if(allocated(Idx_tag_diag))deallocate(Idx_tag_diag)
		allocate(Idx_tag_diag(2,0:mpi_size-1))
		call mpi_distribute(1, tag_init, mpi_size, Idx_tag_diag)
		my_tag_range(:) = Idx_tag_diag(:,my_id)
		do mpi_id=1,mpi_size-1
			call MPI_Send(Idx_tag_diag(:,mpi_id), 2, MPI_INTEGER, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
		enddo
	else
		call MPI_Recv(my_tag_range, 2, MPI_INTEGER, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

!		First of all you need the grid on which the TDSE is discretized
	call fedvr_grid_create(grid_temp, (/1,N_fe/), (/1,N_gp/), N_fun)
	call fedvr_grid_allocate(grid_temp)
	call fedvr_grid_set(grid_temp)

	if(my_tag_range(2)/=0)then	! diagonalize in parallel... if needed (tag_init>1)
		write(6,*) my_id, " - diagonalizes tags: (",my_tag_range(1),",",my_tag_range(2),")"

!			Next, the different operators that make up the Hamiltonian
		call set_2nd_deriv_op(deriv2, grid_temp)
		call set_core_potential_ops(Vx, dVx, grid_temp)
		call Op_Sym_create(Sym, grid_temp)

		do tag = my_tag_range(1),my_tag_range(2)
			write(tag_str,'(I5.5)')tag
			call TDCC_get_lm_from_tag(SphericalHarmonicRepresentation, tag, ll, mm)
			call set_centrifugal_ops(Vl, dVl, grid_temp, ll)

			call OP_create(Ham0_init, grid_temp, l_init, m_init)
			Ham0_init = ((mHalf * deriv2) + Vx) + Vl

			do ii_fe=1,Ham0_init%Nfe
				Ham0_init%fe_region(ii_fe)%matrix = Sym%fe_region(ii_fe)%matrix * Ham0_init%fe_region(ii_fe)%matrix
			enddo

!				Now set up the eigenpair object 
			call allocate_eigenpair(Eigenp,Ham0_init%Ngp)
			call Op_set_bigmat(Ham0_init, Eigenp%matrix)

!			call EigenSolve_sy(Eigenp)
			call EigenSolve_ge(Eigenp)
			call eigenpair_save_h5(Eigenp, 'Ham0_'//tag_str//'_EigenP', grid_temp%xx, grid_temp%ww, .TRUE.)
		enddo
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	if(my_id==root_id)then
		write(tag_str,'(I5.5)')tag_init
		call eigenvals_load_h5(EigenVals, 'Ham0_'//tag_str//'_EigenP')
		call eigenvec_load_h5(EigenVec, n_th, 'Ham0_'//tag_str//'_EigenP')
		EigenVec(:) = EigenVec(:) * dsqrt(grid_temp%ww)
		Energy_sys  = EigenVals(n_th)
		write(6,*)"init energy:", Energy_sys
		call wf_create(wf_temp, grid_temp, (/1, tagmax/))	!just for the creation of 'Psi_init_diag.h5'
		call wf_create_output_h5(wf_temp, 'Psi_init_diag', 1)
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	call mpi_distribute_local_vec(EigenVec, wf%re(:,tag_init))

	call wf_norm(wf)
	call wf_dump_output_h5(wf, 'Psi_init_diag', grid, 0)

	do tag = tag_init, 1, -1
		write(tag_str,'(I5.5)')tag
		if(my_id==root_id)then
			call eigenvals_load_h5(EigenVals, 'Ham0_'//tag_str//'_EigenP')
			nn = 0
			do ii=1,size(EigenVals)							! count the states that need to be frozen
				if(EigenVals(ii)<Energy_sys) nn = nn + 1
			enddo
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(nn, 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)

!		do ii = nn, 1, -1
		do ii = 1, nn		! to be in the same order as for ITP
			call allocate_BoundState(bState, grid%Ngp)
			if(my_id==root_id)then
				call TDCC_get_lm_from_tag(SphericalHarmonicRepresentation, tag, ll, mm)
				call eigenvec_load_h5(EigenVec, ii, 'Ham0_'//tag_str//'_EigenP')
				bState%n = ii
				bState%l = ll
				bState%m = mm
				bState%Energy = EigenVals(ii)
				EigenVec(:) = EigenVec(:) * dsqrt(grid_temp%ww)
			endif
			call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

			call MPI_bcast(bState%n,      1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)
			call MPI_bcast(bState%l,      1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)
			call MPI_bcast(bState%m,      1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)
			call MPI_bcast(bState%Energy, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
			call mpi_distribute_local_vec(EigenVec, bState%vector)
			call normalize_BoundState(bState, grid)
			call append_BoundState(frozen_states, bState)
		enddo
	enddo		! tag
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	if(allocated(frozen_states))call save_BoundState_arr_h5(frozen_states, grid, 'frozen_states_DIAG', .TRUE.)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE get_WF_init_ITP(grid, wf)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(OUT) :: wf
	TYPE(fedvr_grid),    INTENT(IN)  :: grid
	TYPE(tdcc_fedvr_wf) :: wf_temp
	TYPE(fedvr_grid)    :: grid_temp
	TYPE(fedvr_op_full) :: deriv2
	TYPE(fedvr_op_diag) :: Vx, dVx
	TYPE(fedvr_op_diag) :: Vl, dVl
	TYPE(BoundState)    :: bState
	REAL(DP) :: Energy
	INTEGER(I4B) :: n_th, i_th
	INTEGER(I4B) :: tag, ll, mm

	tag_init = TDCC_get_tag_from_lm(SphericalHarmonicRepresentation, l_init, m_init)
	n_th = n_init - l_init

	call wf_create(wf, grid, my_tag)
	call wf_allocate(wf)



	call set_2nd_deriv_op(deriv2, grid)
	call set_core_potential_ops(Vx, dVx, grid)
	call set_centrifugal_ops(Vl, dVl, grid, l_init)

	call OP_create(Ham0_init, grid, l_init, m_init)
	Ham0_init = ((mHalf * deriv2) + Vx) + Vl

	do ii_fe=1,Ham0_init%Nfe
		Ham0_init%fe_region(ii_fe)%matrix = OP_Sym%fe_region(ii_fe)%matrix * Ham0_init%fe_region(ii_fe)%matrix
	enddo

	if(my_id==root_id)then ! tag, ll, mm, i_th, Energy
		write(6,'(A47)')"#============================================#"
		write(6,'(A47)')"|      INVERSE IMAGINARY TIME PROPAGATION    |"
		write(6,'(A47)')"#============================================#"
		write(6,'(A47)')"| tag |  l  |  m  |  n  |       Energy       |"
		write(6,'(A47)')"#--------------------------------------------#"
	endif


	if(allocated(frozen_states_tag))deallocate(frozen_states_tag)
	do i_th=1,n_th
		call wf_create(wf_temp, grid, (/1, 1/))
		call wf_allocate(wf_temp)
		wf_temp%re(:,1) = 1.0d0
		call wf_normalize(wf_temp)
!		write(6,*) wf_temp%norm, dTime_ITP

		call inverse_imag_Time_propag(wf_temp, Energy)
		if(my_id==root_id)write(6,110) tag_init, l_init, m_init, i_th, Energy
		if(i_th /= n_th)then
			call allocate_BoundState(bState, grid%Ngp)
			bState%n = i_th
			bState%l = l_init
			bState%m = m_init
			bState%Energy = Energy
			bState%vector = wf_temp%re(:,1)
			call normalize_BoundState(bState, grid)
			call append_BoundState(frozen_states_tag, bState)
			call append_BoundState(frozen_states,     bState)					! needed for tprop
		else
			Energy_sys = Energy
			wf%re(:,tag_init) = wf_temp%re(:,1)
			wf%im(:,tag_init) = wf_temp%im(:,1)		! should be zero
			call wf_norm(wf)
		endif
	enddo


	if(my_id==root_id)then
		call fedvr_grid_create(grid_temp, (/1,N_fe/), (/1,N_gp/), N_fun)
		call fedvr_grid_allocate(grid_temp)
		call fedvr_grid_set(grid_temp)

		call wf_create(wf_temp, grid_temp, (/1, tagmax/))
		call wf_create_output_h5(wf_temp, 'Psi_init_ITP', 1)
	endif
	call wf_dump_output_h5(wf, 'Psi_init_ITP', grid, 0)


	if(tag_init>1)then		! complete the frozen states
		do tag = tag_init-1, 1, -1
			call TDCC_get_lm_from_tag(SphericalHarmonicRepresentation, tag, ll, mm)
			call set_2nd_deriv_op(deriv2, grid)
			call set_core_potential_ops(Vx, dVx, grid)
			call set_centrifugal_ops(Vl, dVl, grid, ll)

			call OP_create(Ham0_init, grid, ll, mm)
			Ham0_init = ((mHalf * deriv2) + Vx) + Vl

			do ii_fe=1,Ham0_init%Nfe
				Ham0_init%fe_region(ii_fe)%matrix = OP_Sym%fe_region(ii_fe)%matrix * Ham0_init%fe_region(ii_fe)%matrix
			enddo

			if(allocated(frozen_states_tag))deallocate(frozen_states_tag)
!			i_th =0
!			Energy=-10000000000.d0
!			do while(Energy < Energy_sys)
			DO i_th = 1, 10000
				call wf_create(wf_temp, grid, (/1, 1/))
				call wf_allocate(wf_temp)
				wf_temp%re(:,1) = 1.0d0
				call wf_normalize(wf_temp)
		!		write(6,*) wf_temp%norm, dTime_ITP

				call inverse_imag_Time_propag(wf_temp, Energy)
!				i_th = i_th + 1
				if(Energy > Energy_sys) EXIT									! double check: it is possible that in the previous iteration we came close 
				if(my_id==root_id)write(6,110) tag, ll, mm, i_th, Energy
				call allocate_BoundState(bState, grid%Ngp)
				bState%n = i_th
				bState%l = ll
				bState%m = mm
				bState%Energy = Energy
				bState%vector = wf_temp%re(:,1)
				call normalize_BoundState(bState, grid)
				call append_BoundState(frozen_states_tag, bState)
				call append_BoundState(frozen_states,     bState)				! needed for tprop
			enddo

		enddo
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	if(my_id==root_id)write(6,'(A47)')"#================================================#"
	if(allocated(frozen_states))call save_BoundState_arr_h5(frozen_states, grid, 'frozen_states_ITP', .TRUE.)

110	format(" | ",I3.1, " | ", I3.1, " | ", I3.1, " | ", I3.1, " | ", f18.15, " |")
END SUBROUTINE
!-------------------------------------------------------------------------------



END MODULE Setup
