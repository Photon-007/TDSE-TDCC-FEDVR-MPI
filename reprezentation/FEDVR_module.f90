MODULE FEDVR_module
	USE nrtype
	USE gvars
	USE gridsub
	USE MPI_pars
	USE mpi



	TYPE fedvr_grid
		LOGICAL(LGT) :: Im_a_chunk
		INTEGER(I4B), DIMENSION(1:2) :: fe_range
		INTEGER(I4B), DIMENSION(1:2) :: gp_range
		INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: fe_gp_range
		INTEGER(I4B) :: Nfun, Nfe, Ngp
		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xx_fe, ww_fe					! (jj,ii_fe) -> fortran is column major (leftmost index runs fastest)!
		REAL(DP), DIMENSION(:),   ALLOCATABLE :: xx, ww
		REAL(DP), DIMENSION(:),   ALLOCATABLE :: ww_lo, ww_up					! weight of the lower/upper neighbour
	END TYPE


	TYPE(fedvr_grid) :: grid													! local grid of each process
	TYPE(fedvr_grid) :: grid_G													! global grid with all points
	TYPE(fedvr_grid) :: grid_temp

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_grid_write_header(grid)
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(INOUT) :: grid

	do mpi_id=0,mpi_size-1
		if(my_id==mpi_id)then
			write(6,*) 'my_id = ', my_id, 'grid : ', grid%Nfun, grid%Nfe, grid%Ngp, grid%fe_range, grid%gp_range
		endif
		call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_grid_create(grid, fe_range, gp_range, Nfun)
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(INOUT) :: grid
	INTEGER(I4B), DIMENSION(2), INTENT(IN) :: fe_range
	INTEGER(I4B), DIMENSION(2), INTENT(IN) :: gp_range
	INTEGER(I4B), INTENT(IN) :: Nfun

	if((fe_range(1)==1).and.(fe_range(2)==N_fe))then
		grid%Im_a_chunk = .FALSE.
	else
		grid%Im_a_chunk = .TRUE.
	endif
	grid%fe_range = fe_range
	grid%gp_range = gp_range
	grid%Nfun     = Nfun
	grid%Nfe      = fe_range(2) - fe_range(1) + 1
	grid%Ngp      = gp_range(2) - gp_range(1) + 1

	if(allocated(grid%xx_fe))deallocate(grid%xx_fe)
	if(allocated(grid%ww_fe))deallocate(grid%ww_fe)
	if(allocated(grid%xx))deallocate(grid%xx)
	if(allocated(grid%ww))deallocate(grid%ww)
	if(allocated(grid%ww_lo))deallocate(grid%ww_lo)
	if(allocated(grid%ww_up))deallocate(grid%ww_up)
	if(allocated(grid%fe_gp_range))deallocate(grid%fe_gp_range)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_grid_allocate(grid)
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(INOUT) :: grid

	if(allocated(grid%xx_fe))deallocate(grid%xx_fe)
	if(allocated(grid%ww_fe))deallocate(grid%ww_fe)
	if(allocated(grid%xx))deallocate(grid%xx)
	if(allocated(grid%ww))deallocate(grid%ww)
	if(allocated(grid%ww_lo))deallocate(grid%ww_lo)
	if(allocated(grid%ww_up))deallocate(grid%ww_up)
	if(allocated(grid%fe_gp_range))deallocate(grid%fe_gp_range)
	allocate(grid%xx_fe(1:grid%Nfun, 1:grid%Nfe))
	allocate(grid%ww_fe(1:grid%Nfun, 1:grid%Nfe))
	allocate(grid%xx(1:grid%Ngp))
	allocate(grid%ww(1:grid%Ngp))
	allocate(grid%ww_lo(1:grid%Nfe))
	allocate(grid%ww_up(1:grid%Nfe))
	allocate(grid%fe_gp_range(1:2, 1:grid%Nfe))

	grid%xx_fe(:,:) = 0.0d0
	grid%ww_fe(:,:) = 0.0d0
	grid%xx(:) = 0.0d0
	grid%ww(:) = 0.0d0
	grid%ww_lo(:) = 0.0d0
	grid%ww_up(:) = 0.0d0

!	do ii_fe=grid%fe_range(1),grid%fe_range(2)
	do ii_fe=1,grid%Nfe
		grid%fe_gp_range(1, ii_fe) = (ii_fe-1)*(grid%Nfun-1) + 1
		grid%fe_gp_range(2, ii_fe) = (ii_fe-1)*(grid%Nfun-1) + grid%Nfun
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_grid_deallocate(grid)
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(INOUT) :: grid

	if(allocated(grid%xx_fe))deallocate(grid%xx_fe)
	if(allocated(grid%ww_fe))deallocate(grid%ww_fe)
	if(allocated(grid%xx))deallocate(grid%xx)
	if(allocated(grid%ww))deallocate(grid%ww)
	if(allocated(grid%ww_lo))deallocate(grid%ww_lo)
	if(allocated(grid%ww_up))deallocate(grid%ww_up)
	if(allocated(grid%fe_gp_range))deallocate(grid%fe_gp_range)
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_grid_set(grid)
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(INOUT) :: grid
	REAL(DP), DIMENSION(1:2) :: endpoints
	REAL(DP), DIMENSION(:), ALLOCATABLE :: scratch 
	REAL(DP), DIMENSION(:), ALLOCATABLE :: gauss_weights
	REAL(DP), DIMENSION(:), ALLOCATABLE :: roots
	REAL(DP), DIMENSION(:), ALLOCATABLE :: xma, xmi
	INTEGER(I4B) :: jj


	if(allocated(xma))deallocate(xma)
	if(allocated(xmi))deallocate(xmi)
	allocate(xma(1:N_fe), xmi(1:N_fe))
	xma(:) = 0.0d0
	xmi(:) = 0.0d0
	do jj_fe=1,N_fe-1
		if(jj_fe <= N_tfe)then
			xma(jj_fe) = xmi(jj_fe) + dx_min + (dx_max-dx_min)*(jj_fe-0.5d0)/N_tfe
		else
			xma(jj_fe) = xmi(jj_fe) + dx_max
		endif
		xmi(jj_fe+1) = xma(jj_fe)
	enddo
	xma(N_fe) = xmi(N_fe) + dx_max


	endpoints(1) = -1.d0
	endpoints(2) =  1.d0

	do ii_fe = 1,grid%Nfe							! local
		jj_fe = grid%fe_range(1) - 1 + ii_fe		! global
		if(jj_fe==1)then
			if(allocated(scratch))deallocate(scratch)
			if(allocated(gauss_weights))deallocate(gauss_weights)
			if(allocated(roots))deallocate(roots)            
			allocate(scratch(1:grid%Nfun+1), roots(1:grid%Nfun+1), gauss_weights(1:grid%Nfun+1))
			call gaussq('legendre', grid%Nfun+1, 0.d0, 1.d0, 2, endpoints, scratch, roots, gauss_weights)
			do jj=1,grid%Nfun
				grid%xx_fe(jj, ii_fe) = xmi(jj_fe) + (roots(jj+1)+1.d0)*(xma(jj_fe)-xmi(jj_fe))*0.5d0
				grid%ww_fe(jj, ii_fe) = gauss_weights(jj+1)*0.5d0*(xma(jj_fe)-xmi(jj_fe))
			enddo
		else
			if(allocated(scratch))deallocate(scratch)
			if(allocated(gauss_weights))deallocate(gauss_weights)
			if(allocated(roots))deallocate(roots)            
			allocate(scratch(1:grid%Nfun),roots(1:grid%Nfun),gauss_weights(1:grid%Nfun))
			call gaussq('legendre', grid%Nfun, 0.d0, 1.d0, 2, endpoints, scratch, roots, gauss_weights)
			do jj=1,grid%Nfun
				grid%xx_fe(jj, ii_fe) = xmi(jj_fe) + (roots(jj)+1.d0)*(xma(jj_fe)-xmi(jj_fe))*0.5d0
				grid%ww_fe(jj, ii_fe) = gauss_weights(jj)*0.5d0*(xma(jj_fe)-xmi(jj_fe))
			enddo
		endif

		do jj=1,grid%Nfun
			grid%xx((ii_fe-1)*(grid%Nfun-1) + jj) = grid%xx_fe(jj, ii_fe)
			grid%ww((ii_fe-1)*(grid%Nfun-1) + jj) = grid%ww((ii_fe-1)*(grid%Nfun-1) + jj) + grid%ww_fe(jj, ii_fe)
		enddo
	enddo

!-----------------------------------------------------------#
!		Get the weights of the neighbouring FE first/last	|
!		point used for the bridge functions!				|
!-----------------------------------------------------------#	
	do ii_fe = 2,grid%Nfe
		grid%ww_lo(ii_fe) = grid%ww_fe(grid%Nfun,ii_fe-1)
	enddo
	do ii_fe = 1,grid%Nfe-1
		grid%ww_up(ii_fe) = grid%ww_fe(1,ii_fe+1)
	enddo

	if(grid%fe_range(1)==2)then			! most certainly never to be used
		if(allocated(scratch))deallocate(scratch)
		if(allocated(gauss_weights))deallocate(gauss_weights)
		if(allocated(roots))deallocate(roots)            
		allocate(scratch(1:grid%Nfun+1), roots(1:grid%Nfun+1), gauss_weights(1:grid%Nfun+1))
		call gaussq('legendre', grid%Nfun+1, 0.d0, 1.d0, 2, endpoints, scratch, roots, gauss_weights)
		grid%ww_lo(1) = gauss_weights(grid%Nfun+1)*0.5d0*(xma(grid%fe_range(1)-1)-xmi(grid%fe_range(1)-1))
	elseif(grid%fe_range(1)>2)then
		if(allocated(scratch))deallocate(scratch)
		if(allocated(gauss_weights))deallocate(gauss_weights)
		if(allocated(roots))deallocate(roots)            
		allocate(scratch(1:grid%Nfun), roots(1:grid%Nfun), gauss_weights(1:grid%Nfun))
		call gaussq('legendre', grid%Nfun, 0.d0, 1.d0, 2, endpoints, scratch, roots, gauss_weights)
		grid%ww_lo(1) = gauss_weights(grid%Nfun)*0.5d0*(xma(grid%fe_range(1)-1)-xmi(grid%fe_range(1)-1))
	endif

	if(grid%fe_range(2)<N_fe)then
		if(allocated(scratch))deallocate(scratch)
		if(allocated(gauss_weights))deallocate(gauss_weights)
		if(allocated(roots))deallocate(roots)            
		allocate(scratch(1:grid%Nfun), roots(1:grid%Nfun), gauss_weights(1:grid%Nfun))
		call gaussq('legendre', grid%Nfun, 0.d0, 1.d0, 2, endpoints, scratch, roots, gauss_weights)
		grid%ww_up(grid%Nfe) = gauss_weights(1)*0.5d0*(xma(grid%fe_range(2)+1)-xmi(grid%fe_range(2)+1))
	endif

!!!!!!!!! Takes care of the mpi boundaries !!!!!!!!
	if(grid%fe_range(1) /= 1)    grid%ww(1)        = grid%ww(1)        + grid%ww_lo(1)
	if(grid%fe_range(2) /= N_fe) grid%ww(grid%Ngp) = grid%ww(grid%Ngp) + grid%ww_up(grid%Nfe)
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE fedvr_grid_tofile(grid,fileName)	! change to h5
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(IN) :: grid
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B) :: ii, jj, unit_nr

	call get_free_unit_nr(unit_nr)
	open(unit=unit_nr, file=fileName//'_fe.dat', status='replace')
		write(unit_nr,*)'# fe_range = ', grid%fe_range(1), grid%fe_range(2)
		write(unit_nr,*)'# gp_range = ', grid%gp_range(1), grid%gp_range(2)
		write(unit_nr,*)'# Nfe = ', grid%Nfe
		write(unit_nr,*)'# Ngp = ', grid%Ngp
		write(unit_nr,*)'# Nfun= ', grid%Nfun
		write(unit_nr,*)''

		write(unit_nr,*)'# ife    ifun    xx_fe    ww_fe'
		do ii_fe=1,grid%Nfe
			jj_fe = grid%fe_range(1) - 1 + ii_fe
			do ii=1,grid%Nfun
				write(unit_nr,'(2(I12.1), 2(4x,g25.17E3))') jj_fe, ii, grid%xx_fe(ii,ii_fe), grid%ww_fe(ii,ii_fe)
			enddo
		enddo
		write(unit_nr,*)''

		write(unit_nr,*)'# ife    lo    up'
		do ii_fe=1,grid%Nfe
			jj_fe = grid%fe_range(1) - 1 + ii_fe
			write(unit_nr,'(I12.1, 2(4x,g25.17E3))') jj_fe, grid%ww_lo(ii_fe), grid%ww_up(ii_fe)
		enddo

	close(unit_nr)

	open(unit=unit_nr, file=fileName//'_gp.dat', status='replace')
		write(unit_nr,*)'# fe_range = ', grid%fe_range(1), grid%fe_range(2)
		write(unit_nr,*)'# gp_range = ', grid%gp_range(1), grid%gp_range(2)
		write(unit_nr,*)'# Nfe = ', grid%Nfe
		write(unit_nr,*)'# Ngp = ', grid%Ngp
		write(unit_nr,*)'# Nfun= ', grid%Nfun
		write(unit_nr,*)''

		write(unit_nr,*)'# igp    xx    ww'
		do ii = 1,grid%Ngp
			jj = grid%gp_range(1) - 1 + ii
			write(unit_nr,'(I12.1, 2(4x,g25.17E3))') jj, grid%xx(ii), grid%ww(ii)
		enddo

	close(unit_nr)

END SUBROUTINE
!-------------------------------------------------------------------------------

!===============================================================================
!*******************************************************************************
SUBROUTINE mpi_distribute_local_grids(grid_global, grid_local)
	IMPLICIT NONE
	TYPE(fedvr_grid), INTENT(IN)    :: grid_global
	TYPE(fedvr_grid), INTENT(INOUT) :: grid_local

	if(my_id == root_id)then
		grid_local%xx_fe(:,:) = grid_global%xx_fe(1:grid_local%Nfun, grid_local%fe_range(1):grid_local%fe_range(2))
		grid_local%ww_fe(:,:) = grid_global%ww_fe(1:grid_local%Nfun, grid_local%fe_range(1):grid_local%fe_range(2))
		grid_local%xx(:) = grid_global%xx(grid_local%gp_range(1):grid_local%gp_range(2))
		grid_local%ww(:) = grid_global%ww(grid_local%gp_range(1):grid_local%gp_range(2))
		grid_local%ww_lo(:) = grid_global%ww_lo(grid_local%fe_range(1):grid_local%fe_range(2))
		grid_local%ww_up(:) = grid_global%ww_up(grid_local%fe_range(1):grid_local%fe_range(2))
		do mpi_id=1,mpi_size-1
			call MPI_Send(grid_global%xx_fe(:,Idx_fe(1, mpi_id):Idx_fe(2, mpi_id)), N_fe_arr(mpi_id)*N_fun, MPI_DOUBLE_PRECISION, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(grid_global%ww_fe(:,Idx_fe(1, mpi_id):Idx_fe(2, mpi_id)), N_fe_arr(mpi_id)*N_fun, MPI_DOUBLE_PRECISION, mpi_id, 2, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(grid_global%xx(Idx_gp(1, mpi_id):Idx_gp(2, mpi_id)), N_gp_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 3, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(grid_global%ww(Idx_gp(1, mpi_id):Idx_gp(2, mpi_id)), N_gp_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 4, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(grid_global%ww_lo(Idx_fe(1, mpi_id):Idx_fe(2, mpi_id)), N_fe_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 5, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(grid_global%ww_up(Idx_fe(1, mpi_id):Idx_fe(2, mpi_id)), N_fe_arr(mpi_id), MPI_DOUBLE_PRECISION, mpi_id, 6, MPI_COMM_WORLD, mpi_err)
		enddo
	else
		call MPI_Recv(grid_local%xx_fe, grid_local%Nfe*grid_local%Nfun, MPI_DOUBLE_PRECISION, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(grid_local%ww_fe, grid_local%Nfe*grid_local%Nfun, MPI_DOUBLE_PRECISION, root_id, 2, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(grid_local%xx, grid_local%Ngp, MPI_DOUBLE_PRECISION, root_id, 3, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(grid_local%ww, grid_local%Ngp, MPI_DOUBLE_PRECISION, root_id, 4, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(grid_local%ww_lo, grid_local%Nfe, MPI_DOUBLE_PRECISION, root_id, 5, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(grid_local%ww_up, grid_local%Nfe, MPI_DOUBLE_PRECISION, root_id, 6, MPI_COMM_WORLD, rstatus, mpi_err)
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------

END MODULE FEDVR_module











































