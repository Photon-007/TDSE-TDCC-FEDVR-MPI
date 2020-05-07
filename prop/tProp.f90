MODULE tProp
	USE nrtype
	USE gvars
	USE Lanczos
	USE FEDVR_potentials
	USE LASER_potentials
	USE MPI_pars
	USE mpi
	USE HDF5

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE inverse_imag_Time_propag(wf, Energy_ground)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	REAL(DP), INTENT(OUT) :: Energy_ground
	REAL(DP) :: Energy_new, Energy_old
!	REAL(DP) :: tt, dt
	REAL(DP) :: rel_diff
	INTEGER(I4B) :: ii, idx_bs
	LOGICAL(LGT) :: stoping_crit



	if(allocated(Q_vec))deallocate(Q_vec)
	if(allocated(Q_mat))deallocate(Q_mat)
	if(allocated(Q_alpha))deallocate(Q_alpha)
	if(allocated(Q_beta)) deallocate(Q_beta)
	if(allocated(coeffs_curr))deallocate(coeffs_curr)
	if(allocated(coeffs_prev))deallocate(coeffs_prev)
	allocate(Q_vec(0:N_Krylov), Q_mat(1:N_Krylov,1:N_Krylov))
	allocate(Q_alpha(0:N_Krylov), Q_beta(0:N_Krylov))
	allocate(coeffs_curr(0:N_Krylov), coeffs_prev(0:N_Krylov))

	do ii = 0,N_Krylov
		call wf_mirror(wf, Q_vec(ii), 'alloc')
	enddo

	Time = 0.0d0
	dt_max = dTime_ITP
!	dt     = dTime_ITP
	stoping_crit = .FALSE.
	Energy_new = 0.0d0

	do while(.not.stoping_crit)
		Energy_old = Energy_new
		call propagator_Im_Time(wf, Time, dTime_ITP, Energy_new)
		if(allocated(frozen_states_tag))then !call remove_States(wf, frozen_states, Epsilon_ortho)
			do idx_bs=1,size(frozen_states_tag)
				call wf_remove_bState_ITP(wf, frozen_states_tag(idx_bs), Epsilon_ortho)
			enddo
			call wf_normalize(wf)
		endif

!		if(my_id==root_id)write(6,*) Time, dTime_ITP, Energy_new

!		PAUSE
		rel_diff = (Energy_new - Energy_old)/Energy_new
		if((Time>Time_max_ITP_au).or.(dabs(rel_diff)<=Epsilon_conv_ITP)) stoping_crit = .TRUE.

	enddo

	Energy_ground = Energy_new
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE real_Time_propag(wf, Time_start, Time_end)
	IMPLICIT NONE
	TYPE(tdcc_fedvr_wf), INTENT(INOUT) :: wf
	REAL(DP), INTENT(IN) :: Time_start
	REAL(DP), INTENT(IN) :: Time_end
	TYPE(fedvr_grid)    :: grid_temp
	TYPE(tdcc_fedvr_wf) :: wf_temp
	REAL(DP) :: Energy
	REAL(DP) :: norm_tt
	REAL(DP) :: dt_trial, rescale
	REAL(DP) :: dt_mean
	INTEGER(I4B) :: ii, idx_bs, Ndt
	INTEGER(I4B) :: idx_wf, idx_data

	write(6,*)
	write(proc_tag,'(I5.5)')my_id

	if(allocated(Q_vec))deallocate(Q_vec)
	if(allocated(Q_mat))deallocate(Q_mat)
	if(allocated(Q_alpha))deallocate(Q_alpha)
	if(allocated(Q_beta)) deallocate(Q_beta)
	if(allocated(coeffs_curr))deallocate(coeffs_curr)
	if(allocated(coeffs_prev))deallocate(coeffs_prev)
	allocate(Q_vec(0:N_Krylov), Q_mat(1:N_Krylov,1:N_Krylov))
	allocate(Q_alpha(0:N_Krylov), Q_beta(0:N_Krylov))
	allocate(coeffs_curr(0:N_Krylov), coeffs_prev(0:N_Krylov))

	do ii = 0,N_Krylov
		call wf_mirror(wf, Q_vec(ii), 'alloc')
	enddo

	dt_max  = dTime
	Time    = Time_start
	rescale = 1.0d0
	Ndt     = 0
	dt_mean = 0.0d0
!	write(6,*)Time_start, Time, Time_end, dTime

	if(restart)then
		ii=1
		do while(ii*Time_step_wf_au < Time_start)
			ii = ii + 1
		enddo
		idx_wf = ii 
		ii=1
		do while(ii*Time_step_data_au < Time_start)
			ii = ii + 1
		enddo
		idx_data = ii 
	else
		idx_wf   = 1
		idx_data = 1
		if(my_id==root_id)then
			call fedvr_grid_create(grid_temp, (/1,N_fe/), (/1,N_gp/), N_fun)
			call fedvr_grid_allocate(grid_temp)
			call fedvr_grid_set(grid_temp)

			call wf_create(wf_temp, grid_temp, (/1, tagmax/))
			call wf_create_output_h5(wf_temp, 'Psi', N_step_wf)
			call fedvr_grid_deallocate(grid_temp)
			call create_LaserField_h5('LaserField', OP_EField, N_step_data)
			call create_TimePropag_timing_h5('Timing_TimePropag', N_step_data)
		endif
	endif

	call wf_norm(wf)
	call wf_dump_output_h5(wf, 'Psi', grid, idx_wf-1)


	if(my_id==root_id)then
		write(6,*)"#=====================================================================#"
		write(6,*)"|    Time    |       ||wf||     |     Sys_energy   |    Laser Field   |"
		write(6,*)"#---------------------------------------------------------------------#"
	endif
	do while(Time<Time_end)

		dt_trial = MIN(idx_wf*Time_step_wf_au-Time, idx_data*Time_step_data_au-Time)
		if(dt_trial < dTime)then
			rescale = dt_trial / dTime
			dTime   = dt_trial
		endif
!		write(6,*)
!		write(6,*)idx_wf, idx_data, dt_trial, dTime, rescale
		if(dTime<0.0d0)then
			write(6,*)"Propagation failed!!!"
			write(6,*)"dt = ", dTime
			Stop
		endif

		tStep_start = MPI_Wtime()
		call get_EE_op(OP_EField, Time)
		call propagator_Re_Time(wf, Time, dTime, Energy)
		norm_tt = wf%norm
		if(allocated(frozen_states))then !call remove_States(wf, frozen_states, Epsilon_ortho)
			do idx_bs=1,size(frozen_states)
				call wf_remove_bState(wf, frozen_states(idx_bs), Epsilon_ortho)
			enddo
			call wf_norm(wf)
			call wf_scale(dsqrt(norm_tt/wf%norm), wf)
		endif

		call get_maskfun(Maskfun_op, Abspot_op, dTime)
!		call OP_tofile(Maskfun_op, 'Maskfun_'//proc_tag, 'f12.6')
		wf = Maskfun_op * wf
		call wf_norm(wf)
		tStep_end = MPI_Wtime()
		tStep_elapsed = tStep_elapsed + (tStep_end - tStep_start)
!		write(6,*)wf%norm
!		pause

		dt_mean = dt_mean + dTime
		Ndt     = Ndt + 1

		if(dabs(Time - idx_wf*Time_step_wf_au)<1e-7)then
			idx_wf = idx_wf + 1
			call wf_dump_output_h5(wf, 'Psi', grid, idx_wf-1)
			dTime   = dTime /rescale
			rescale = 1.0d0
		endif

		if(dabs(Time - idx_data*Time_step_data_au)<1e-7)then
			if(my_id==root_id)then
				write(6,110) Time, dsqrt(wf%norm), Energy, OP_EField%EE_t
!				...
			endif
			call dump_LaserField_h5('LaserField', OP_EField, idx_data)
			call dump_TimePropag_timing_h5('Timing_TimePropag', idx_data, Time, Ndt, dt_mean, tStep_elapsed)

			idx_data = idx_data + 1
			dTime   = dTime /rescale
			rescale = 1.0d0
			dt_mean = 0.0d0
			Ndt     = 0
			tStep_elapsed = 0.0d0
		endif

	enddo

	if(my_id==root_id)write(6,*)"#=====================================================================#"
	call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
110	format(" |",f12.8, "|", f18.15, "|", f18.15, "|", f18.15, "|")

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE create_TimePropag_timing_h5(fileName, Ntime)!(tt, Ndt, dt, elapsed)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B),     INTENT(IN) :: Ntime
	CHARACTER(LEN=4), PARAMETER :: dset_tt_name = "Time"
	CHARACTER(LEN=7), PARAMETER :: dset_dt_name = "<dTime>"
	CHARACTER(LEN=7), PARAMETER :: dset_Ndt_name= "N_dTime"
	CHARACTER(LEN=7), PARAMETER :: dset_el_name = "Elapsed"
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_id, dspace_id
	INTEGER(HSIZE_T), DIMENSION(2) :: dims
	INTEGER(I4B) :: error


	if(my_id==root_id)then
		call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)
		
		dims = (/mpi_size, Ntime/)
!			Time
		call h5screate_simple_f(2, dims, dspace_id, error)
		call h5dcreate_f(file_id, dset_tt_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)

!			N_dTime
		call h5screate_simple_f(2, dims, dspace_id, error)
		call h5dcreate_f(file_id, dset_Ndt_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)

!			<dTime>
		call h5screate_simple_f(2, dims, dspace_id, error)
		call h5dcreate_f(file_id, dset_dt_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)

!			elapsed
		call h5screate_simple_f(2, dims, dspace_id, error)
		call h5dcreate_f(file_id, dset_el_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)

		call h5fclose_f(file_id, error)
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE dump_TimePropag_timing_h5(fileName, idx_time, tt, Ndt, dt, elapsed)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B),     INTENT(IN) :: idx_time, Ndt
	REAL(DP),         INTENT(IN) :: tt, dt, elapsed
	CHARACTER(LEN=4), PARAMETER :: dset_tt_name = "Time"
	CHARACTER(LEN=7), PARAMETER :: dset_dt_name = "<dTime>"
	CHARACTER(LEN=7), PARAMETER :: dset_Ndt_name= "N_dTime"
	CHARACTER(LEN=7), PARAMETER :: dset_el_name = "Elapsed"
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: plist_id
	INTEGER(HID_T) :: dset_id, dspace_id, mspace_id
	INTEGER(HSIZE_T), DIMENSION(1) :: slab_dims, mspace_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: count_p, offset_p, stride_p, block_p
	INTEGER(I4B) :: error

!	if(allocated(tmp_real))deallocate(tmp_real)
!	if(allocated(tmp_int)) deallocate(tmp_int)
!	allocate(tmp_real(mpi_size), tmp_int(mpi_size))

	slab_dims   = (/1/)
	mspace_dims = (/1/)

	count_p  = (/1,1/)
	offset_p = (/my_id,idx_time/)
	stride_p = (/1,1/)
	block_p  = (/1,1/)

	call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
	call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
	call h5pclose_f(plist_id, error)


	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

!			Time
	call h5screate_simple_f(1, mspace_dims, mspace_id, error)
	call h5dopen_f(file_id, dset_tt_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)

	call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tt, slab_dims, error, mspace_id, dspace_id, xfer_prp = plist_id)

	call h5dclose_f(dset_id, error)
	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)

!			N_dTime
	call h5screate_simple_f(1, mspace_dims, mspace_id, error)
	call h5dopen_f(file_id, dset_Ndt_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)

	call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, Ndt, slab_dims, error, mspace_id, dspace_id, xfer_prp = plist_id)

	call h5dclose_f(dset_id, error)
	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)

!			<dTime>
	call h5screate_simple_f(1, mspace_dims, mspace_id, error)
	call h5dopen_f(file_id, dset_dt_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)

	call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dt, slab_dims, error, mspace_id, dspace_id, xfer_prp = plist_id)

	call h5dclose_f(dset_id, error)
	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)

!			elapsed
	call h5screate_simple_f(1, mspace_dims, mspace_id, error)
	call h5dopen_f(file_id, dset_el_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)

	call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, elapsed, slab_dims, error, mspace_id, dspace_id, xfer_prp = plist_id)

	call h5dclose_f(dset_id, error)
	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)



	call h5pclose_f(plist_id, error)

	call h5fclose_f(file_id, error)




END SUBROUTINE
!-------------------------------------------------------------------------------










END MODULE tProp
