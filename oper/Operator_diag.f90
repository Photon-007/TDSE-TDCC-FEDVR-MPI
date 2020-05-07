MODULE Operator_diag
	USE nrtype
	USE gvars
	USE MPI_pars
	USE FEDVR_module
	USE TDCC_module
	USE TDCC_FEDVR
	USE Operator_typedef
	USE Operator_math
	USE FEDVR_derivatives
	USE FEDVR_potentials
	USE Laser_potentials
	USE eigen_calc


	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE diagonalize_H0(ll_in, mm_in)
	IMPLICIT NONE
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: ll_in
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: mm_in
	TYPE(fedvr_op_full) :: Ham0
	TYPE(fedvr_op_full) :: deriv2
	TYPE(fedvr_op_diag) :: Vx, dVx
	TYPE(fedvr_op_diag) :: Vl, dVl
	TYPE(eigenpair_obj) :: Eigenp
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: Idx_tag_diag
	INTEGER(I4B), DIMENSION(:)  , ALLOCATABLE :: my_tag_range
	INTEGER(I4B) :: tag, ll, mm
	CHARACTER(LEN=5) :: tag_str



	if(present(ll_in))then
		if(my_id==root_id)then
			ll = ll_in
			if(present(mm_in))then
				mm = mm_in
			else
				mm = m_init
			endif
			tag = TDCC_get_tag_from_lm(SphericalHarmonicRepresentation, ll, mm)
			write(tag_str,'(I5.5)')tag
			write(6,*)"l=",ll
			write(6,*)"m=",mm

!			First of all you need the grid on which the TDSE is discretized
			call fedvr_grid_create(grid, (/1,N_fe/), (/1,N_gp/), N_fun)
			call fedvr_grid_allocate(grid)
			call fedvr_grid_set(grid)

			call set_2nd_deriv_op(deriv2, grid)
			call set_core_potential_ops(Vx, dVx, grid)
			call set_centrifugal_ops(Vl, dVl, grid, ll)

			call Op_create(Ham0, grid, ll, mm)
			call Op_allocate(Ham0)

			Ham0 = ((mHalf * deriv2) + Vx) + Vl

			call allocate_eigenpair(Eigenp,Ham0%Ngp)
			call Op_set_bigmat(Ham0, Eigenp%matrix)

!			call EigenSolve_sy(Eigenp)
			call EigenSolve_ge(Eigenp)
			call eigenpair_save_h5(Eigenp, 'Ham0_'//tag_str//'_EigenP', grid%xx, grid%ww)

		endif	!my_id
	else
		if(allocated(my_tag_range))deallocate(my_tag_range)
		allocate(my_tag_range(2))
		if(my_id == root_id)then
			if(allocated(Idx_tag_diag))deallocate(Idx_tag_diag)
			allocate(Idx_tag_diag(2,0:mpi_size-1))
			call mpi_distribute(1, tagmax, mpi_size, Idx_tag_diag)
			my_tag_range(:) = Idx_tag_diag(:,my_id)
			do mpi_id=1,mpi_size-1
				call MPI_Send(Idx_tag_diag(:,mpi_id), 2, MPI_INTEGER, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
			enddo
		else
			call MPI_Recv(my_tag_range, 2, MPI_INTEGER, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)


!		First of all you need the grid on which the TDSE is discretized
		call fedvr_grid_create(grid, (/1,N_fe/), (/1,N_gp/), N_fun)
		call fedvr_grid_allocate(grid)
		call fedvr_grid_set(grid)

		call set_2nd_deriv_op(deriv2, grid)
		call set_core_potential_ops(Vx, dVx, grid)

		do tag = my_tag_range(1), my_tag_range(2)
			write(tag_str,'(I5.5)')tag
			call TDCC_get_lm_from_tag(SphericalHarmonicRepresentation, tag, ll, mm)
			write(6,*)" my_id:",my_id," | l=",ll," | m=",mm," | tag=",tag

			call set_centrifugal_ops(Vl, dVl, grid, ll)

			call Op_create(Ham0, grid, ll, mm)
			call Op_allocate(Ham0)

			Ham0 = ((mHalf * deriv2) + Vx) + Vl

			call allocate_eigenpair(Eigenp,Ham0%Ngp)
			call Op_set_bigmat(Ham0, Eigenp%matrix)


!			call EigenSolve_sy(Eigenp)
			call EigenSolve_ge(Eigenp)
			call eigenpair_save_h5(Eigenp, 'Ham0_'//tag_str, grid%xx, grid%ww)

		enddo
	endif


END SUBROUTINE
!-------------------------------------------------------------------------------



END MODULE Operator_diag








































