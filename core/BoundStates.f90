MODULE BoundStates
	USE nrtype
	USE gvars
	USE FEDVR_module
	USE HDF5

	TYPE BoundState
		INTEGER(I4B) :: n, l, m
		REAL(DP) :: Energy
		REAL(DP) :: norm
		REAL(DP), DIMENSION(:), ALLOCATABLE :: vector
	END TYPE

	TYPE(BoundState), DIMENSION(:), ALLOCATABLE :: frozen_states
	TYPE(BoundState), DIMENSION(:), ALLOCATABLE :: frozen_states_tag

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE allocate_BoundState(bState, Ngp)	! Ngp so you can use it on either the full grid or on a chunk (ex. if you build the frozen_states with ITP)
	IMPLICIT NONE
	TYPE(BoundState), INTENT(INOUT) :: bState
	INTEGER(I4B),     INTENT(IN)    :: Ngp

	bState%n = 0
	bState%l = 0
	bState%m = 0
	bState%Energy	= 0.0d0
	if(allocated(bState%vector))deallocate(bState%vector)
	allocate(bState%vector(1:Ngp))
	bState%vector(:) = 0.0d0
	bState%norm      = 0.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE allocate_BoundState_arr(bState_arr, Nst, Ngp)
	IMPLICIT NONE
	TYPE(BoundState), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: bState_arr
	INTEGER(I4B), INTENT(IN) :: Nst, Ngp
	INTEGER(I4B) :: ii

	if(allocated(bState_arr))deallocate(bState_arr)
	allocate(bState_arr(1:Nst))
	do ii = 1, Nst
		call allocate_BoundState(bState_arr(ii), Ngp)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE append_BoundState(bState_arr, bState)
	IMPLICIT NONE
	TYPE(BoundState), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: bState_arr
	TYPE(BoundState), INTENT(IN) :: bState
	TYPE(BoundState), DIMENSION(:), ALLOCATABLE :: temp
	INTEGER(I4B) :: Nst, Ngp, ii

	Ngp = size(bState%vector)
	if(.not.allocated(bState_arr))then
		call allocate_BoundState_arr(bState_arr, 1, Ngp)
		bState_arr(1) = bState
	else
		Nst = size(bState_arr)
		call allocate_BoundState_arr(temp, Nst, Ngp)
		do ii = 1, Nst
			temp(ii) = bState_arr(ii)
		enddo
		call allocate_BoundState_arr(bState_arr, Nst+1, Ngp)
		do ii = 1, Nst
			bState_arr(ii) = temp(ii)
		enddo
		bState_arr(Nst+1) = bState
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE normalize_BoundState(bState, grid)
	IMPLICIT NONE
	TYPE(BoundState), INTENT(INOUT) :: bState
	TYPE(fedvr_grid), INTENT(IN)    :: grid
	REAL(DP), EXTERNAL :: dot, ddot
	REAL(DP) :: vec_norm, bState_norm
	INTEGER(I4B) :: Ngp, lb

	if(grid%fe_range(1)==1)then
		lb  = 1
		Ngp = grid%Ngp
	else
		lb  = 2				! so you don't count the overlapping gp twice when runing on multiple nodes
		Ngp = grid%Ngp - 1
	endif

	vec_norm = ddot(Ngp, bState%vector(lb), 1, bState%vector(lb), 1)
	call MPI_Allreduce(vec_norm, bState_norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
	bState%vector(:) = bState%vector(:) / dsqrt(bState_norm)
	bState%norm = 1.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE save_BoundState_arr_h5(bState_arr, grid, fileName, remove_weights)
	IMPLICIT NONE
	TYPE(BoundState), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: bState_arr
	TYPE(fedvr_grid), INTENT(IN) :: grid
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	LOGICAL(LGT), INTENT(IN), OPTIONAL :: remove_weights
	CHARACTER(LEN=6), PARAMETER  :: dset_grid_name = "Grid"
	CHARACTER(LEN=12),PARAMETER  :: dset_vec_name  = "State_vector"
	CHARACTER(LEN=12),PARAMETER  :: dset_E_name = "State_Energy"
	CHARACTER(LEN=7),PARAMETER   :: dset_n_name = "State_n"
	CHARACTER(LEN=7),PARAMETER   :: dset_l_name = "State_l"
	CHARACTER(LEN=7),PARAMETER   :: dset_m_name = "State_m"
	INTEGER(HID_T) :: plist_id
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_grid_id,   dset_vec_id,  dset_prop_id
	INTEGER(HID_T) :: dspace_grid_id, dspace_vec_id, dspace_prop_id
	INTEGER(HID_T) :: mspace_grid_id, mspace_vec_id
	INTEGER(HSIZE_T), DIMENSION(2) :: grid_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: vec_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: prop_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: grid_slab_dims, mspace_grid_dims
	INTEGER(HSIZE_T), DIMENSION(1) :: vec_slab_dims, mspace_vec_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: count_p, offset_p, stride_p, block_p
	REAL(DP),     DIMENSION(:), ALLOCATABLE :: tmp_real
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: tmp_int
	INTEGER(I4B) :: error
	INTEGER(I4B) :: Nst, Ngp, ii
	INTEGER(I4B) :: il, iu

	Nst = size(bState_arr)
	if(my_id == root_id)then
		Ngp = grid%Ngp
		il = 1
		iu = grid%Ngp
	else
		Ngp = grid%Ngp - 1
		il = 2
		iu = grid%Ngp
	endif
!	write(6,*) my_id, Nst, Ngp, grid%Ngp

	if(my_id==root_id)then	! create the file and the groups
		if(allocated(tmp_real))deallocate(tmp_real)
		if(allocated(tmp_int)) deallocate(tmp_int)
		allocate(tmp_real(Nst), tmp_int(Nst))
		call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

		grid_dims = (/N_gp, 2/)
		call h5screate_simple_f(2, grid_dims, dspace_grid_id, error)
		call h5dcreate_f(file_id, dset_grid_name, H5T_NATIVE_DOUBLE, dspace_grid_id, dset_grid_id, error)
		call h5dclose_f(dset_grid_id, error)
		call h5sclose_f(dspace_grid_id, error)

!		write(6,*)(bState_arr(ii)%Energy,ii=1,Nst)
		prop_dims = (/Nst/)
		do ii = 1, Nst
			tmp_real(ii) = bState_arr(ii)%Energy
		enddo
		call h5screate_simple_f(1, prop_dims, dspace_prop_id, error)
		call h5dcreate_f(file_id, dset_E_name, H5T_NATIVE_DOUBLE, dspace_prop_id, dset_prop_id, error)
		call h5dwrite_f(dset_prop_id, H5T_NATIVE_DOUBLE, tmp_real, prop_dims, error)
		call h5dclose_f(dset_prop_id, error)
		call h5sclose_f(dspace_prop_id, error)

		do ii = 1, Nst
			tmp_int(ii) = bState_arr(ii)%n
		enddo
		call h5screate_simple_f(1, prop_dims, dspace_prop_id, error)
		call h5dcreate_f(file_id, dset_n_name, H5T_NATIVE_INTEGER, dspace_prop_id, dset_prop_id, error)
		call h5dwrite_f(dset_prop_id, H5T_NATIVE_INTEGER, tmp_int, prop_dims, error)
		call h5dclose_f(dset_prop_id, error)
		call h5sclose_f(dspace_prop_id, error)

		do ii = 1, Nst
			tmp_int(ii) = bState_arr(ii)%l
		enddo
		call h5screate_simple_f(1, prop_dims, dspace_prop_id, error)
		call h5dcreate_f(file_id, dset_l_name, H5T_NATIVE_INTEGER, dspace_prop_id, dset_prop_id, error)
		call h5dwrite_f(dset_prop_id, H5T_NATIVE_INTEGER, tmp_int, prop_dims, error)
		call h5dclose_f(dset_prop_id, error)
		call h5sclose_f(dspace_prop_id, error)

		do ii = 1, Nst
			tmp_int(ii) = bState_arr(ii)%m
		enddo
		call h5screate_simple_f(1, prop_dims, dspace_prop_id, error)
		call h5dcreate_f(file_id, dset_m_name, H5T_NATIVE_INTEGER, dspace_prop_id, dset_prop_id, error)
		call h5dwrite_f(dset_prop_id, H5T_NATIVE_INTEGER, tmp_int, prop_dims, error)
		call h5dclose_f(dset_prop_id, error)
		call h5sclose_f(dspace_prop_id, error)



		vec_dims = (/N_gp, Nst/)
		call h5screate_simple_f(2, vec_dims, dspace_vec_id, error)
		call h5dcreate_f(file_id, dset_vec_name, H5T_NATIVE_DOUBLE, dspace_vec_id, dset_vec_id, error)
		call h5dclose_f(dset_vec_id, error)
		call h5sclose_f(dspace_vec_id, error)

		call h5fclose_f(file_id, error)
	endif



	call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
	call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
	call h5pclose_f(plist_id, error)


!	#=======================#
!	|	Write grid data		|	!	(..,1) -> xx; (..,2) -> ww
!	#=======================#
	grid_slab_dims = (/Ngp/)
	mspace_grid_dims=(/Ngp/)
	call h5screate_simple_f(1, mspace_grid_dims, mspace_grid_id, error)
	call h5dopen_f(file_id, dset_grid_name, dset_grid_id, error)
	call h5dget_space_f(dset_grid_id, dspace_grid_id, error)

	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

	if(my_id == root_id)then
		count_p  = (/1,1/)
		offset_p = (/grid%gp_range(1)-1,0/)
		stride_p = (/1,1/)
		block_p  = (/Ngp,1/)
	else
		count_p  = (/1,1/)
		offset_p = (/grid%gp_range(1),0/)
		stride_p = (/1,1/)
		block_p  = (/Ngp,1/)
	endif
!	write(6,*) my_id, 'offset = ', offset_g
!	write(6,*) my_id, 'block  = ', block_g

	call h5sselect_hyperslab_f(dspace_grid_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid%xx(il:iu), grid_slab_dims, error, mspace_grid_id, dspace_grid_id, xfer_prp = plist_id)

	if(my_id == root_id)then
		count_p  = (/1,1/)
		offset_p = (/grid%gp_range(1)-1,1/)
		stride_p = (/1,1/)
		block_p  = (/Ngp,1/)
	else
		count_p  = (/1,1/)
		offset_p = (/grid%gp_range(1),1/)
		stride_p = (/1,1/)
		block_p  = (/Ngp,1/)
	endif
	call h5sselect_hyperslab_f(dspace_grid_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid%ww(il:iu), grid_slab_dims, error, mspace_grid_id, dspace_grid_id, xfer_prp = plist_id)

	call h5pclose_f(plist_id, error)
	call h5dclose_f(dset_grid_id, error)
	call h5sclose_f(dspace_grid_id, error)
	call h5sclose_f(mspace_grid_id, error)

!	#=======================#
!	|	Write State vector	|
!	#=======================#

	vec_slab_dims   = (/Ngp/)
	mspace_vec_dims = (/Ngp/)

	call h5screate_simple_f(1, mspace_vec_dims, mspace_vec_id, error)
	call h5dopen_f(file_id, dset_vec_name, dset_vec_id, error)
	call h5dget_space_f(dset_vec_id, dspace_vec_id, error)	! (/N_gp,Nst/)

	call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
	call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

	if(allocated(tmp_real))deallocate(tmp_real)
	allocate(tmp_real(Ngp))
	do ii = 1,Nst
		if(my_id == root_id)then
			count_p  = (/1,1/)
			offset_p = (/grid%gp_range(1)-1, ii-1/)
			stride_p = (/1,1/)
			block_p  = (/Ngp,1/)
		else
			count_p  = (/1,1/)
			offset_p = (/grid%gp_range(1), ii-1/)
			stride_p = (/1,1/)
			block_p  = (/Ngp,1/)
		endif
!		write(6,*) my_id, 'offset = ', offset_p
!		write(6,*) my_id, 'block  = ', block_p
		tmp_real = bState_arr(ii)%vector(il:iu)
		if(present(remove_weights))then
			if(remove_weights) tmp_real = tmp_real / dsqrt(grid%ww(il:iu))
		endif
		call h5sselect_hyperslab_f(dspace_vec_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
		call h5dwrite_f(dset_vec_id, H5T_NATIVE_DOUBLE, tmp_real, vec_slab_dims, error, mspace_vec_id, dspace_vec_id, xfer_prp = plist_id)

	enddo

	call h5pclose_f(plist_id, error)
	call h5dclose_f(dset_vec_id, error)
	call h5sclose_f(dspace_vec_id, error)
	call h5sclose_f(mspace_vec_id, error)

	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE BoundStates

































































