MODULE MPI_pars
	USE nrtype
	USE mpi

	INTEGER(I4B) :: mpi_err
	INTEGER(I4B) :: mpi_size
	INTEGER(I4B) :: my_id 
!	INTEGER(I4B) :: nr_proc 
	INTEGER(I4B) :: root_id = 0 ! the rank of the "root" process
	INTEGER(I4B) :: mpi_id		! used throughout the code wherever a loop ower id is needed
	INTEGER(I4B), DIMENSION(MPI_STATUS_SIZE) :: rstatus
	CHARACTER(LEN=5) :: proc_tag

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE mpi_para_range(ni, nf, nproc, irank, ista, iend)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ni, nf
	INTEGER(I4B), INTENT(IN) :: nproc, irank
	INTEGER(I4B), INTENT(INOUT) :: ista, iend
	INTEGER(I4B) :: iwork1, iwork2

	iwork1 = (nf - ni + 1) / nproc
	iwork2 = MOD(nf - ni + 1, nproc)
	ista = irank * iwork1 + ni + MIN(irank, iwork2)
	iend = ista + iwork1 - 1
	if(iwork2 > irank) iend = iend + 1
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE mpi_distribute(ni, nf, nproc, bins)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ni, nf
	INTEGER(I4B), INTENT(IN) :: nproc
	INTEGER(I4B), DIMENSION(2,0:nproc-1), INTENT(INOUT) :: bins
	INTEGER(I4B) :: remainder, chunk, ii, nn

	nn = (nf - ni + 1)
	chunk = int(nn / nproc)
	remainder = MOD(nn, nproc)
!	write(6,*)chunk, remainder
	bins(1,0) = ni
	bins(2,0) = ni+chunk-1
	if(nproc>=2)then
		do ii=1,nproc-1
			if(ii <= (nproc-1-remainder))then
				bins(1,ii) = bins(2,ii-1) + 1
				bins(2,ii) = bins(1,ii) + chunk - 1
			else
				bins(1,ii) = bins(2,ii-1)  + 1
				bins(2,ii) = bins(1,ii) + chunk 
			endif
		enddo
	endif
END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE MPI_pars

