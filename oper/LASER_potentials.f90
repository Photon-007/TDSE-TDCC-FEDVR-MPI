MODULE LASER_potentials
	USE nrtype
	USE gvars
	USE Operator_typedef
	USE TDCC_module
	USE HDF5


	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE get_EE_op(OP, tt)
	IMPLICIT NONE
	TYPE(EE_op), INTENT(INOUT) :: OP
	REAL(DP), INTENT(IN) :: tt
	REAL(DP) :: tt_t0
	INTEGER(I4B) :: ii

	OP%EE_t  = 0.0d0
	OP%Env_t = 0.0d0
	do ii=1,OP%Npulse
		tt_t0 = tt - OP%pulses(ii)%t0
		select case (OP%pulses(ii)%Envelope)
			case ('sin2')
				if((tt>=OP%pulses(ii)%t_start).and.(tt<=OP%pulses(ii)%t_end))then
					OP%Envi_t(ii) = OP%pulses(ii)%EA * dsin(OP%pulses(ii)%param * tt_t0 + PIO2)**2
					OP%Ei_t(ii)   = dsin(OP%pulses(ii)%Omega * tt_t0 + OP%pulses(ii)%CEP) * OP%Envi_t(ii)
				else
					OP%Ei_t(ii) = 0.0d0
				endif
			case ('cos2')
				if((tt>=OP%pulses(ii)%t_start).and.(tt<=OP%pulses(ii)%t_end))then
					OP%Envi_t(ii) = OP%pulses(ii)%EA * dcos(OP%pulses(ii)%param * tt_t0)**2
					OP%Ei_t(ii)   = dsin(OP%pulses(ii)%Omega * tt_t0 + OP%pulses(ii)%CEP) * OP%Envi_t(ii)
				else
					OP%Ei_t(ii) = 0.0d0
				endif
			case ('gauss')
					OP%Envi_t(ii) = OP%pulses(ii)%EA * dexp(OP%pulses(ii)%param * tt_t0**2)
					OP%Ei_t(ii)   = dsin(OP%pulses(ii)%Omega * tt_t0 + OP%pulses(ii)%CEP) * OP%Envi_t(ii)
			case ('flat')
					OP%Envi_t(ii) = OP%pulses(ii)%EA * 1.0d0
					OP%Ei_t(ii)   = dsin(OP%pulses(ii)%Omega * tt_t0 + OP%pulses(ii)%CEP) * OP%Envi_t(ii)
		end select
		OP%EE_t  = OP%EE_t  + OP%Ei_t(ii)
		OP%Env_t = OP%Env_t + OP%Envi_t(ii)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE create_LaserField_h5(fileName, LF, Ntime)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	TYPE(EE_op),      INTENT(IN) :: LF
	INTEGER(I4B),     INTENT(IN) :: Ntime
	CHARACTER(LEN=10), PARAMETER :: dset_EEsum_name = "EField_sum"
	CHARACTER(LEN=9)             :: dset_EEith_name
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HID_T) :: dset_id, dspace_id
	INTEGER(HID_T) :: attr_id, aspace_id
	INTEGER(HSIZE_T), DIMENSION(2) :: dims
	INTEGER(HSIZE_T), DIMENSION(1) :: attr_dims
	INTEGER(SIZE_T) :: attr_len
	INTEGER(I4B) :: error, idx_pulse
	CHARACTER(LEN=2) :: pulse_id


	if(my_id==root_id)then
		call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

		dims = (/2, Ntime/)
!			EField_sum
		call h5screate_simple_f(2, dims, dspace_id, error)
		call h5dcreate_f(file_id, dset_EEsum_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)

!			EField_ii
		do idx_pulse = 1,LF%Npulse
			write(pulse_id, '(I2.2)')idx_pulse
			dset_EEith_name = "EField_"//pulse_id

			call h5screate_simple_f(2, dims, dspace_id, error)
			call h5dcreate_f(file_id, dset_EEith_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

!				Attributes (i.e. pulse parameters)
			attr_dims=(/1/)
			attr_len = 10
			call h5screate_simple_f(1, attr_dims, aspace_id, error)
			call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)		! create a 'string' type from character
			call h5tset_size_f(atype_id, attr_len, error)

			call h5acreate_f(dset_id, 'Envelope_type', atype_id, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, atype_id, LF%pulses(idx_pulse)%Envelope, attr_dims, error)
			call h5aclose_f(attr_id, error)
			call h5sclose_f(aspace_id, error)


			attr_dims=(/0/)
			call h5screate_f(H5S_SCALAR_F, aspace_id, error)
			call h5acreate_f(dset_id, 'Amplitude', H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LF%pulses(idx_pulse)%EA, attr_dims, error)
			call h5aclose_f(attr_id, error)

			call h5acreate_f(dset_id, 'Time_center', H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LF%pulses(idx_pulse)%t0, attr_dims, error)
			call h5aclose_f(attr_id, error)

			call h5acreate_f(dset_id, 'Omega', H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LF%pulses(idx_pulse)%Omega, attr_dims, error)
			call h5aclose_f(attr_id, error)

			call h5acreate_f(dset_id, 'CEP', H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LF%pulses(idx_pulse)%CEP, attr_dims, error)
			call h5aclose_f(attr_id, error)

			call h5acreate_f(dset_id, 'Time_start', H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LF%pulses(idx_pulse)%t_start, attr_dims, error)
			call h5aclose_f(attr_id, error)

			call h5acreate_f(dset_id, 'Time_end', H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
			call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LF%pulses(idx_pulse)%t_end, attr_dims, error)
			call h5aclose_f(attr_id, error)

			call h5sclose_f(aspace_id, error)
			call h5dclose_f(dset_id, error)
			call h5sclose_f(dspace_id, error)

		enddo


		call h5fclose_f(file_id, error)
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE dump_LaserField_h5(fileName, LF, idx_time)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	TYPE(EE_op),      INTENT(IN) :: LF
	INTEGER(I4B),     INTENT(IN) :: idx_time
	CHARACTER(LEN=10), PARAMETER :: dset_EEsum_name = "EField_sum"
	CHARACTER(LEN=9)             :: dset_EEith_name
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_id, dspace_id, mspace_id
	INTEGER(HSIZE_T), DIMENSION(1) :: slab_dims, mspace_dims
	INTEGER(HSIZE_T), DIMENSION(2) :: count_p, offset_p, stride_p, block_p
	INTEGER(I4B) :: error, idx_pulse
	CHARACTER(LEN=2) :: pulse_id
	REAL(DP), DIMENSION(2) :: tiny_vec


	if(my_id==root_id)then
		call h5fopen_f(fileName//'.h5', H5F_ACC_RDWR_F, file_id, error)
		slab_dims   = (/2/)
		mspace_dims = (/2/)

		count_p  = (/1,1/)
		offset_p = (/0,idx_time/)
		stride_p = (/1,1/)
		block_p  = (/2,1/)

		tiny_vec = (/LF%EE_t, LF%Env_t/)
!		call h5screate_simple_f(2, slab_dims, dspace_id, error)
!		call h5screate_simple_f(1, mspace_dims, mspace_id, error)
!		call h5dcreate_f(file_id, dset_EEsum_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

		call h5screate_simple_f(1, mspace_dims, mspace_id, error)
		call h5dopen_f(file_id, dset_EEsum_name, dset_id, error)
		call h5dget_space_f(dset_id, dspace_id, error)

		call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tiny_vec, slab_dims, error, mspace_id, dspace_id)

		call h5dclose_f(dset_id, error)
		call h5sclose_f(dspace_id, error)
		call h5sclose_f(mspace_id, error)

!			EField_ii
		do idx_pulse = 1,LF%Npulse
			write(pulse_id, '(I2.2)')idx_pulse
			dset_EEith_name = "EField_"//pulse_id

			tiny_vec = (/LF%Ei_t(idx_pulse), LF%Envi_t(idx_pulse)/)
!			call h5screate_simple_f(2, slab_dims, dspace_id, error)
!			call h5screate_simple_f(1, mspace_dims, mspace_id, error)
!			call h5dcreate_f(file_id, dset_EEith_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

			call h5screate_simple_f(1, mspace_dims, mspace_id, error)
			call h5dopen_f(file_id, dset_EEith_name, dset_id, error)
			call h5dget_space_f(dset_id, dspace_id, error)

			call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
			call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tiny_vec, slab_dims, error, mspace_id, dspace_id)

			call h5dclose_f(dset_id, error)
			call h5sclose_f(dspace_id, error)
			call h5sclose_f(mspace_id, error)

		enddo

		call h5fclose_f(file_id, error)
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE











































