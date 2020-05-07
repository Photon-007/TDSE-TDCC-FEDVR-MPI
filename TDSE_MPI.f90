PROGRAM main
	USE nrtype
	USE gvars
	USE BoundStates
	USE MPI_pars
	USE FEDVR_module
	USE TDCC_module
	USE TDCC_FEDVR
	USE Operator_typedef
	USE Operator_math
	USE Operator_diag
	USE Operator_expect
	USE FEDVR_derivatives
	USE FEDVR_potentials
	USE LASER_potentials
	USE eigen_calc
	USE unit_conversions
	USE Setup
	USE tProp
	USE Lanczos
	USE specfun
	USE HDF5


!	These variables are just temporary for the testing period
	IMPLICIT NONE
	INTEGER(I4B) :: ii !, ll, mm, tag
	REAL(DP), EXTERNAL :: dot, ddot
	INTEGER(I4B) :: h5_error
	INTEGER(I4B) :: iTime
	COMPLEX(DPC) :: matrix_el
	LOGICAL(LGT) :: test_H0_x_psi = .FALSE.

!	#===========================#
!	|	Initialize hdf5 and MPI	|
!	#===========================#
	call h5open_f(h5_error)
	call MPI_init(mpi_err)
	call MPI_COMM_RANK(MPI_COMM_WORLD, my_id,   mpi_err)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_err)
	write(proc_tag,'(I5.5)')my_id
	Main_start = MPI_Wtime()


!	#===============================================#
!	|	Read, set and broadcast global variables	|
!	#===============================================#
	if(my_id==root_id)then
!		call write_Atomic_units(6)
		call read_params('parameters.in')
		call set_dependent()
		call set_propag_params()
		call set_ITP_params()
		call divide_workload()
!		do mpi_id=0,mpi_size-1
!			write(6,*)mpi_id, Idx_fe(1,mpi_id), Idx_fe(2,mpi_id), N_fe_arr(mpi_id), Idx_gp(1,mpi_id), Idx_gp(2,mpi_id), N_gp_arr(mpi_id)
!		enddo
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	call bcast_gvars()
	call bcast_init_state_params()
	call bcast_propag_params()
	call bcast_ITP_params()
	call bcast_local_vars()
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	call allocate_pulse_params(Npulse)
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	if(my_id==root_id)then
		do ii = 1,Npulse
			call read_pulse(Pulse_files(ii))
			call set_pulse_params(ii)
		enddo
		call write_params(6,'all')
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	call bcast_pulse_params()
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

!	write(6,*) my_id, LaserPulse(1)

!	#=======================================#
!	|	If only diagonalization is needed	|
!	#=======================================#
	if(do_diagH0)then
		call diagonalize_H0()
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
		goto 99
	endif

!	#=======================================#
!	|	If only the pulse shape is needed	|
!	#=======================================#

	if(do_pulse)then
		if(my_id==root_id)then
!			write(6,*) E0_au, t0_au, omega_au, swiss_knife, CEP_rad, Envelope_type
			call fedvr_grid_create(grid_temp, (/1, N_fe/), (/1, N_gp/), N_fun)
			call fedvr_grid_allocate(grid_temp)
			call EE_op_create(OP_EField, grid_temp, (/1, tagmax/), SphericalHarmonicRepresentation, LaserPulse)
			call EE_op_allocate(OP_EField)
			call create_LaserField_h5('LaserField', OP_EField, N_step_data+1)
!			open(unit=77,file='Laser_field.dat',status='replace')	! make hdf5 output
			do iTime=0,N_step_data
				Time = iTime * Time_step_data_au
				call get_EE_op(OP_EField, Time)
!				write(77,*)Time, OP_EField%EE_t, OP_EField%Env_t, (OP_EField%Ei_t(ii), OP_EField%Envi_t(ii), ii=1,OP_EField%Npulse)
				call dump_LaserField_h5('LaserField', OP_EField, iTime)
			enddo
!			close(77)
		endif
		goto 99
	endif


!	#=======================================#
!	|	Set up the Schrodinger equation		|
!	|	(grid, wf, OPs on each process) 	|
!	#=======================================#
!	Start with the grid, as it is the 'parent class' for every other object
	call fedvr_grid_create(grid, my_fe, my_gp, N_fun)
	call fedvr_grid_allocate(grid)
!	call mpi_distribute_local_grids(grid_all, grid)								! use either one of  ---> this one only works if the global grid was created and set first
	call fedvr_grid_set(grid)													! these two options
	call fedvr_grid_write_header(grid)
!	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

!	Next, the operators to be used
	call set_system_operators(grid)
!	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	if(test_H0_x_psi)then
		call fedvr_grid_create(grid_temp, (/1,N_fe/), (/1,N_gp/), N_fun)
		call fedvr_grid_allocate(grid_temp)
		call fedvr_grid_set(grid_temp)
		if(my_id==root_id)then
			call wf_create(wf_temp, grid_temp, (/1, 1/))	!just for the creation of 'Psi_test.h5'
			call wf_create_output_h5(wf_temp, 'Psi_test', 4)
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

		call wf_create(psi, grid, (/1,1/))
		call wf_allocate(psi)
		psi%re(:,1) = 1.0d0 * dsqrt(grid%ww)
		call wf_dump_output_h5(psi, 'Psi_test', grid, 0)
		psi = OP_Ham0(1) * psi
		call wf_dump_output_h5(psi, 'Psi_test', grid, 1)

		psi%re(:,1) = grid%xx**2 * dsqrt(grid%ww)
		call wf_dump_output_h5(psi, 'Psi_test', grid, 2)
		psi = OP_Ham0(1) * psi
		call wf_dump_output_h5(psi, 'Psi_test', grid, 3)
		goto 99
	endif

!	Finally, the initial wf
	if(init_by_ITP)then
		call get_WF_init_ITP(grid, psi)
	else
		call get_WF_init_Diag(grid, psi)
	endif
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
!	if(only_initial_wf) goto 99

	call wf_norm(psi)
	write(6,*) my_id, "||psi|| = ", psi%norm, " | <psi|psi> = ", psi*psi

	matrix_el = OP_expect(OP_Ham0, psi)
	do mpi_id=0,mpi_size-1
		if(my_id==mpi_id)then
			write(6,*) my_id, " <psi|H|psi>_mpi = ", matrix_el
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	enddo
	if(only_initial_wf) goto 99

!	%===============================================%
!	|				TIME PROPAGATION				|
!	%===============================================%

	call real_Time_propag(psi, 0.0d0, Time_max_au)



!99	ii=0	!dummy
99	Main_end = MPI_Wtime()
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

	do mpi_id=0,mpi_size-1
		if(my_id==mpi_id)then
!			write(6,*)"Program terminated correctly! (my_id:",my_id,")"
			write(6,*)" my_id:",my_id," - Program terminated correctly in ",Main_end-Main_start," seconds! (",Time_wf_mirror,"/",N_wf_mirror,")"
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	enddo

	call MPI_FINALIZE(mpi_err)

	call h5close_f(h5_error)

END PROGRAM main






















