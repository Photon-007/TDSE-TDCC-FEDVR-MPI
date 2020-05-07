MODULE eigen_calc
	USE nrtype
	USE HDF5


	TYPE :: eigenpair_obj
		INTEGER(I4B) :: nn
		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matrix
		REAL(DP), DIMENSION(:,:), ALLOCATABLE :: l_eigvecs, r_eigvecs
		REAL(DP), DIMENSION(:),   ALLOCATABLE :: eigvals
	END TYPE eigenpair_obj


	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE allocate_eigenpair(Eigenp,n)
	IMPLICIT NONE
	TYPE(eigenpair_obj), INTENT(INOUT) :: Eigenp
	INTEGER(I4B), INTENT(IN) :: n

	Eigenp%nn = n
	if(allocated(Eigenp%matrix))deallocate(Eigenp%matrix)
	if(allocated(Eigenp%l_eigvecs))deallocate(Eigenp%l_eigvecs)
	if(allocated(Eigenp%r_eigvecs))deallocate(Eigenp%r_eigvecs)
	if(allocated(Eigenp%eigvals))deallocate(Eigenp%eigvals)
	allocate(Eigenp%matrix(1:n,1:n))
	allocate(Eigenp%l_eigvecs(1:n,1:n))
	allocate(Eigenp%r_eigvecs(1:n,1:n))
	allocate(Eigenp%eigvals(1:n))

	Eigenp%matrix(:,:)    = 0.0d0
	Eigenp%l_eigvecs(:,:) = 0.0d0
	Eigenp%r_eigvecs(:,:) = 0.0d0
	Eigenp%eigvals(:)     = 0.0d0

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE deallocate_eigenpair(Eigenp)
	IMPLICIT NONE
	TYPE(eigenpair_obj), INTENT(INOUT) :: Eigenp

	if(allocated(Eigenp%matrix))deallocate(Eigenp%matrix)
	if(allocated(Eigenp%l_eigvecs))deallocate(Eigenp%l_eigvecs)
	if(allocated(Eigenp%r_eigvecs))deallocate(Eigenp%r_eigvecs)
	if(allocated(Eigenp%eigvals))deallocate(Eigenp%eigvals)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE EigenSolve_ge(Eigenp)
	IMPLICIT NONE
	TYPE(eigenpair_obj), INTENT(INOUT) :: Eigenp
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a
	REAL(DP), DIMENSION(:), ALLOCATABLE :: work 
	REAL(DP), DIMENSION(:), ALLOCATABLE :: wr, wi
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vl, vr, vn
	INTEGER(I4B) :: ii, nn
	INTEGER(I4B) :: lda, ldvl, ldvr, lwork, info
	CHARACTER(LEN=1) :: which
	LOGICAL(LGT) :: allReal


	nn = Eigenp%nn

	lda = nn
	ldvl = nn
	ldvr = nn
	lwork = 4*nn

	if(allocated(a))deallocate(a)
	if(allocated(work))deallocate(work)
	if(allocated(wr))deallocate(wr)
	if(allocated(wi))deallocate(wi)
	if(allocated(vl))deallocate(vl)
	if(allocated(vr))deallocate(vr)
	if(allocated(vn))deallocate(vn)
	allocate(a(1:nn,1:nn))
	allocate(work(1:lwork))
	allocate(wr(1:nn))
	allocate(wi(1:nn))
	allocate(vl(1:nn,1:nn))
	allocate(vr(1:nn,1:nn))
	allocate(vn(1:nn,1:nn))

	a(:,:)=Eigenp%matrix(:,:)

!	http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga66e19253344358f5dee1e60502b9e96f.html
!	call dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
!	call dgeev('N', 'V',  n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
	call dgeev('V', 'V', nn, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

	if(info<0)then
		write(6,*)"  ERROR!!! The ",-INFO,"-th argument of \'dgeev\'"
		write(6,*)"       in \'EigenSolve_ge\' had an illega value!"
		STOP
	elseif(info>0)then
		write(6,*)"  ERROR!!! The QR algorithm failed to compute all the"
		write(6,*)"       eigenvals, and no eigenvecs have been computed!"
		STOP
	endif


	allReal=.true.
	do ii=1,nn
		if(abs(wi(ii)).gt. 1.d-10)then
			allReal=.false.
			cycle
		endif
	enddo


	if(allReal)then
		call sort_eigenp(wr, vl, vr)
		which = 'l'
		call normalize_eigenvectors(vl, vr, which, vn)
	else
		write(6,*)"  ERROR!!! Matrix not hermitic in \'EigenSolve\'!"
		STOP
	endif

	Eigenp%eigvals(:) = wr(:)
	if(which=='l')then
		Eigenp%l_eigvecs(:,:) = vn(:,:)
		Eigenp%r_eigvecs(:,:) = vr(:,:)
	elseif(which=='r')then
		Eigenp%l_eigvecs(:,:) = vl(:,:)
		Eigenp%r_eigvecs(:,:) = vn(:,:)
	else
		write(6,*)"Incorrect value for 'which' parameter"
		stop
	endif

	deallocate(a, work, wr, wi, vl, vr, vn)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE EigenSolve_sy(Eigenp)
	IMPLICIT NONE
	TYPE(eigenpair_obj), INTENT(INOUT) :: Eigenp
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a
	REAL(DP), DIMENSION(:), ALLOCATABLE :: work 
	REAL(DP), DIMENSION(:), ALLOCATABLE :: w
	INTEGER(I4B) :: nn
	INTEGER(I4B) :: lda, lwork, info


	nn = Eigenp%nn

	lda = nn
	lwork = 4*nn

	if(allocated(a))deallocate(a)
	if(allocated(work))deallocate(work)
	if(allocated(w))deallocate(w)
	allocate(a(1:nn,1:nn))
	allocate(work(1:lwork))
	allocate(w(1:nn))

	a(:,:)=Eigenp%matrix(:,:)

!	http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html
!	call dsyev(jobz, uplo, n,  a, lda, w, work, lwork, info)
!	call dsyev( 'N',  'U', n,  a, lda, w, work, lwork, info)
	call dsyev( 'V',  'U', nn, a, lda, w, work, lwork, info)

	if(info<0)then
		write(6,*)"  ERROR!!! The ",-INFO,"-th argument of \'dsyev\'"
		write(6,*)"       in \'EigenSolve_sy\' had an illega value!"
		STOP
	elseif(info>0)then
		write(6,*)"  ERROR!!! The \'dsyev\' algorithm failed to converge!"
		write(6,*)"       ",info," off-diagonal elements did not converge to zero!"
		STOP
	endif


	Eigenp%eigvals(:) = w(:)
	Eigenp%r_eigvecs(:,:) = a(:,:)


	deallocate(a, work, w)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE normalize_eigenvectors(vl, vr, which, vn)
!	%---------------------------------------------------------------%
!	| if which = 'l' --> the left eigenvalues will be renormalised  |
!	| if which = 'r' --> the right eigenvalues will be renormalised |
!	%---------------------------------------------------------------%
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: vl, vr
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: vn
	CHARACTER(LEN=1), INTENT(IN) :: which
	REAL(DP) :: sum_v, alpha
	INTEGER(I4B) :: ii, jj, nn
	INTEGER(I4B), DIMENSION(2) :: matrix_shape

	matrix_shape = shape(vl)
	if(matrix_shape(1) /= matrix_shape(2)) then
		write(6,*)"Incorrect matrix shape"
		stop
	else
		nn = matrix_shape(1)
	endif

!	write(6,*)"###           The ", which," eigenvectors are renormalised         ###"
!	write(6,*)"#############################################################"

	do ii=1,nn		!loop over the eigenvectors
		sum_v = 0.0d0
		do jj=1,nn
			sum_v = sum_v + vl(jj,ii) * vr(jj,ii)
		enddo
		alpha = 1.0d0/sum_v
!		alpha = alpha/n
		if(which=='l')then
			vn(:,ii) = vl(:,ii)*alpha
		elseif(which=='r')then
			vn(:,ii) = vr(:,ii)*alpha
		else
			write(6,*)"Incorrect argument for NORMALIZE_EIGENVECTOR(...) subroutine"
		endif
	enddo
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE matrix_inversion(matrix, inv_matrix)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: matrix
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: inv_matrix
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tmp_matrix
	REAL(DP), DIMENSION(:), ALLOCATABLE :: work, ipiv
	INTEGER(I4B) :: nn, lda, lwork, info
	INTEGER(I4B), DIMENSION(2) :: matrix_shape

	matrix_shape = shape(matrix)
	if(matrix_shape(1) /= matrix_shape(2)) then
		write(6,*)"Incorrect matrix shape"
	else
		nn = matrix_shape(1)

		lwork = 2*nn
		lda = nn
		if(allocated(tmp_matrix))deallocate(tmp_matrix)
		if(allocated(work))deallocate(work)
		if(allocated(ipiv))deallocate(ipiv)
		allocate(tmp_matrix(1:nn,1:nn), work(1:lwork), ipiv(1:lwork))
		tmp_matrix = matrix

!		%---------------------------------%
!		|  LU factorization of the matrix |
!		%---------------------------------%

!		call dgetrf(  m,  n, a, lda, ipiv, info )
		call dgetrf( nn, nn, tmp_matrix, lda, ipiv, info )
		if(info/=0)then
!			if(verbose)then
				write(6,*)"LU factorization info = ", info
!			endif
		endif

!		%--------------------------------------%
!		|  calculation of the inverse matrix   |
!		%--------------------------------------%

!		call dgetri( n, a, lda, ipiv, work, lwork, info )
		call dgetri(nn, tmp_matrix, lda, ipiv, work, lwork, info )
		if(info/=0)then
!			if(verbose)then
				write(6,*)"Matrix inversion info = ", info
!			endif
		endif

		inv_matrix(:,:) = tmp_matrix(:,:)

	endif

END SUBROUTINE
!-------------------------------------------------------------------------------

SUBROUTINE sort_eigenp(array, matrix1, matrix2)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:),   INTENT(INOUT) :: array
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: matrix1, matrix2
	REAL(DP), DIMENSION(:), ALLOCATABLE :: v_temp1, v_temp2
	REAL(DP) :: temp
	INTEGER(I4B) :: n, i, swap(1), iswap
!	INTRINSIC MAXLOC
	INTRINSIC MINLOC

	n = size(array)

	if(allocated(v_temp1))deallocate(v_temp1)
	if(allocated(v_temp2))deallocate(v_temp2)
	allocate(v_temp1(1:n),v_temp2(1:n))

	do i = 1, n-1
		swap = MINLOC(array(i:n))
		iswap = swap(1) + i-1
		if(iswap.ne.i)then
			temp = array(i)
			v_temp1(:) = matrix1(:,i)
			v_temp2(:) = matrix2(:,i)
			array(i) = array(iswap)
			matrix1(:,i) = matrix1(:,iswap)
			matrix2(:,i) = matrix2(:,iswap)
			array(iswap) = temp
			matrix1(:,iswap) = v_temp1(:)
			matrix2(:,iswap) = v_temp2(:)
		endif
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE sort_array(array)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:),   INTENT(INOUT) :: array
	REAL(DP) :: temp
	INTEGER(I4B) :: n, i, swap(1), iswap
!	INTRINSIC MAXLOC
	INTRINSIC MINLOC

	n = size(array)
!	write(6,*) n

	do i = 1, n-1
		swap = MINLOC(array(i:n))
		iswap = swap(1) + i-1
		if(iswap.ne.i)then
			temp = array(i)
			array(i) = array(iswap)
			array(iswap) = temp
		endif

	enddo


END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE eigenpair_save(Eigenp,fileName, xx)
	IMPLICIT NONE
	TYPE(eigenpair_obj), INTENT(INOUT) :: Eigenp
	REAL(DP), DIMENSION(:), INTENT(IN) :: xx
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B) :: nn, ii, jj
	CHARACTER(LEN=12) :: str_nn

!	nn=size(Eigenp%eigvals)
	nn=Eigenp%nn

	open(unit=99, file=fileName//'_eigVals.dat', status='replace')
		do ii=1,nn
			write(99,'(I12.1,2x,g25.17E3)')ii,Eigenp%eigvals(ii)
		enddo
	close(99)

	if(size(xx)==nn)then
		write(str_nn,'(I12.1)')nn+1
		open(unit=99, file=fileName//'_eigVecs.dat', status='replace')
			do ii=1,nn
				write(99,'('//trim(str_nn)//'(2x,g25.17E3))') xx(ii), (Eigenp%r_eigvecs(ii,jj), jj=1,nn)
			enddo
		close(99)
	else
		write(str_nn,'(I12.1)')nn
		open(unit=99, file=fileName//'_eigVecs.dat', status='replace')
			do ii=1,nn
				write(99,'('//trim(str_nn)//'(2x,g25.17E3))') (Eigenp%r_eigvecs(ii,jj), jj=1,nn)
			enddo
		close(99)
	endif
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE eigenpair_save_h5(Eigenp, fileName, xx, ww, remove_weights)
	IMPLICIT NONE
	TYPE(eigenpair_obj), INTENT(INOUT) :: Eigenp
	REAL(DP), DIMENSION(:), INTENT(IN) :: xx, ww
	LOGICAL(LGT), INTENT(IN), OPTIONAL :: remove_weights
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	CHARACTER(LEN=7), PARAMETER :: eigval_name = "EigVals"
	CHARACTER(LEN=7), PARAMETER :: eigvec_name = "EigVecs"
	CHARACTER(LEN=4), PARAMETER :: grid_name = "Grid"
	INTEGER(I4B) :: error, rank
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_eigval_id,   dset_eigvec_id,   dset_grid_id
	INTEGER(HID_T) :: dspace_eigval_id, dspace_eigvec_id, dspace_grid_id
	INTEGER(HSIZE_T), DIMENSION(1) :: eigval_dim
	INTEGER(HSIZE_T), DIMENSION(2) :: eigvec_dim
	INTEGER(HSIZE_T), DIMENSION(2) :: grid_dim
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: grid
	REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: vectors
	INTEGER(I4B) :: ii



	call h5fcreate_f(fileName//'.h5', H5F_ACC_TRUNC_F, file_id, error)

!	%===================%
!	|	write grid		|
!	%===================%
	rank=2
	grid_dim=(/Eigenp%nn, 2/)
	if(allocated(grid))deallocate(grid)
	allocate(grid(Eigenp%nn, 2))
	grid(:,1) = xx(:)
	grid(:,2) = ww(:)
	call h5screate_simple_f(rank, grid_dim, dspace_grid_id, error)
	call h5dcreate_f(file_id, grid_name, H5T_NATIVE_DOUBLE, dspace_grid_id, dset_grid_id, error)
	call h5dwrite_f(dset_grid_id, H5T_NATIVE_DOUBLE, grid, grid_dim, error)
	call h5dclose_f(dset_grid_id, error)
	call h5sclose_f(dspace_grid_id, error)
	
!	%=======================%
!	|	write eigenvals		|
!	%=======================%
	rank=1
	eigval_dim(1) = Eigenp%nn
	call h5screate_simple_f(rank, eigval_dim, dspace_eigval_id, error)
	call h5dcreate_f(file_id, eigval_name, H5T_NATIVE_DOUBLE, dspace_eigval_id, dset_eigval_id, error)
	call h5dwrite_f(dset_eigval_id, H5T_NATIVE_DOUBLE, Eigenp%eigvals, eigval_dim, error)
	call h5dclose_f(dset_eigval_id, error)
	call h5sclose_f(dspace_eigval_id, error)

!	%=======================%
!	|	write eigenvecs		|
!	%=======================%

	rank=2
	eigvec_dim(1) = Eigenp%nn
	eigvec_dim(2) = Eigenp%nn
	call h5screate_simple_f(rank, eigvec_dim, dspace_eigvec_id, error)
	call h5dcreate_f(file_id, eigvec_name, H5T_NATIVE_DOUBLE, dspace_eigvec_id, dset_eigvec_id, error)
	if(present(remove_weights))then
		if(remove_weights)then
			if(allocated(vectors))deallocate(vectors)
			allocate(vectors(Eigenp%nn, Eigenp%nn))
			do ii = 1,Eigenp%nn
				vectors(:,ii) = Eigenp%r_eigvecs(:,ii) / dsqrt(ww)
			enddo
			call h5dwrite_f(dset_eigvec_id, H5T_NATIVE_DOUBLE, vectors, eigvec_dim, error)
		else
			call h5dwrite_f(dset_eigvec_id, H5T_NATIVE_DOUBLE, Eigenp%r_eigvecs, eigvec_dim, error)
		endif
	else
		call h5dwrite_f(dset_eigvec_id, H5T_NATIVE_DOUBLE, Eigenp%r_eigvecs, eigvec_dim, error)
	endif
	call h5dclose_f(dset_eigvec_id, error)
	call h5sclose_f(dspace_eigvec_id, error)


	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE eigenpair_load_h5(Eigenp, fileName)
	TYPE(eigenpair_obj), INTENT(OUT) :: Eigenp
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B) :: nn
	CHARACTER(LEN=7), PARAMETER :: eigval_name = "EigVals"
	CHARACTER(LEN=7), PARAMETER :: eigvec_name = "EigVecs"
	INTEGER(I4B) :: error
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_eigval_id,   dset_eigvec_id
	INTEGER(HID_T) :: dspace_eigval_id, dspace_eigvec_id
	INTEGER(HSIZE_T), DIMENSION(1) :: eigval_dim, eigval_dim_max
	INTEGER(HSIZE_T), DIMENSION(2) :: eigvec_dim, eigvec_dim_max


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDONLY_F, file_id, error)

	call h5dopen_f(file_id, eigval_name, dset_eigval_id, error)
	call h5dopen_f(file_id, eigvec_name, dset_eigvec_id, error)
	call h5dget_space_f(dset_eigval_id, dspace_eigval_id, error)
	call h5dget_space_f(dset_eigvec_id, dspace_eigvec_id, error)
	call h5sget_simple_extent_dims_f(dspace_eigval_id, eigval_dim, eigval_dim_max, error)
	call h5sget_simple_extent_dims_f(dspace_eigvec_id, eigvec_dim, eigvec_dim_max, error)
!	write(6,*)eigval_dim, ' / ', eigval_dim_max
!	write(6,*)eigvec_dim, ' / ', eigvec_dim_max
	nn = eigval_dim(1)

	call allocate_eigenpair(Eigenp, nn)
	call h5dread_f(dset_eigval_id, H5T_NATIVE_DOUBLE, Eigenp%eigvals, eigval_dim, error)
	call h5dread_f(dset_eigvec_id, H5T_NATIVE_DOUBLE, Eigenp%r_eigvecs, eigvec_dim, error)


	call h5sclose_f(dspace_eigval_id, error)
	call h5sclose_f(dspace_eigvec_id, error)
	call h5dclose_f(dset_eigval_id, error)
	call h5dclose_f(dset_eigvec_id, error)
	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE eigenvals_load_h5(EigenVals, fileName)
	REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: EigenVals
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B) :: nn
	CHARACTER(LEN=7), PARAMETER :: eigval_name = "EigVals"
	INTEGER(I4B) :: error
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_eigval_id
	INTEGER(HID_T) :: dspace_eigval_id
	INTEGER(HSIZE_T), DIMENSION(1) :: eigval_dim, eigval_dim_max


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDONLY_F, file_id, error)

	call h5dopen_f(file_id, eigval_name, dset_eigval_id, error)
	call h5dget_space_f(dset_eigval_id, dspace_eigval_id, error)
	call h5sget_simple_extent_dims_f(dspace_eigval_id, eigval_dim, eigval_dim_max, error)
!	write(6,*)eigval_dim, ' / ', eigval_dim_max
!	write(6,*)eigvec_dim, ' / ', eigvec_dim_max
	nn = eigval_dim(1)

	if(allocated(EigenVals))deallocate(EigenVals)
	allocate(EigenVals(nn))
	call h5dread_f(dset_eigval_id, H5T_NATIVE_DOUBLE, EigenVals, eigval_dim, error)


	call h5sclose_f(dspace_eigval_id, error)
	call h5dclose_f(dset_eigval_id, error)
	call h5fclose_f(file_id, error)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE eigenvec_load_h5(EigenVec, idx, fileName)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: EigenVec
	CHARACTER(LEN=*), INTENT(IN) :: fileName
	INTEGER(I4B), INTENT(IN) :: idx
	INTEGER(I4B) :: nn
	CHARACTER(LEN=7), PARAMETER :: eigvec_name = "EigVecs"
	INTEGER(I4B) :: error
	INTEGER(HID_T) :: file_id
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: mspace_id
	INTEGER(HSIZE_T), DIMENSION(2) :: eigvec_dim, eigvec_dim_max
	INTEGER(HSIZE_T), DIMENSION(1) :: slab_dim, mspace_dim
	INTEGER(HSIZE_T), DIMENSION(1:2) :: count_p, offset_p, stride_p, block_p


	call h5fopen_f(fileName//'.h5', H5F_ACC_RDONLY_F, file_id, error)

	call h5dopen_f(file_id, eigvec_name, dset_id, error)
	call h5dget_space_f(dset_id, dspace_id, error)
	call h5sget_simple_extent_dims_f(dspace_id, eigvec_dim, eigvec_dim_max, error)
!	write(6,*)eigval_dim, ' / ', eigval_dim_max
!	write(6,*)eigvec_dim, ' / ', eigvec_dim_max
	nn = eigvec_dim(2)

	if(allocated(EigenVec))deallocate(EigenVec)
	allocate(EigenVec(nn))
	slab_dim   = (/nn/)
	mspace_dim = (/nn/)

	count_p  = (/1,1/)
	offset_p = (/0,idx-1/)
	stride_p = (/1,1/)
	block_p  = (/nn,1/)	
	call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_p, count_p, error, stride_p, block_p)
	call h5screate_simple_f(1, mspace_dim, mspace_id, error)
	call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, EigenVec, slab_dim, error, mspace_id, dspace_id)


	call h5sclose_f(dspace_id, error)
	call h5sclose_f(mspace_id, error)
	call h5dclose_f(dset_id, error)
	call h5fclose_f(file_id, error)


END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE eigen_calc






























