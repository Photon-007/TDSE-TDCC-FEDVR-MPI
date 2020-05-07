MODULE FEDVR_derivatives
	USE nrtype
	USE gvars
	USE FEDVR_module
	USE Operator_typedef



	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE set_1st_deriv_op(df_op, grid)	
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(OUT) :: df_op
	TYPE(fedvr_grid),    INTENT(IN)    :: grid
	REAL(DP), DIMENSION(:),   ALLOCATABLE :: x, y
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: f, df, ddf, d_mat
	INTEGER(I4B) :: ii, jj, np

	call Op_create(df_op,  grid, 0, 0)
	call Op_allocate(df_op)
!

	do ii_fe = 1,grid%Nfe						! local
		jj_fe = grid%fe_range(1) - 1 + ii_fe	! global
		if(jj_fe == 1)then
			np = grid%Nfun+1
			if(allocated(x))deallocate(x)
			if(allocated(y))deallocate(y)
			if(allocated(f))deallocate(f)
			if(allocated(df)) deallocate(df)
			if(allocated(ddf))deallocate(ddf)
			if(allocated(d_mat))deallocate(d_mat)
			allocate(x(1:np), y(1:np))
			allocate(f(1:np,1:np), df(1:np,1:np), ddf(1:np,1:np))
			allocate(d_mat(1:np,1:np))

			x(1) = 0.0d0
			x(2:np) = grid%xx_fe(1:np-1, ii_fe)
			y(:) = x(:)
!-----------------------------------------------------------------------------!
!         CALCULATE THE LAGRANGE INTERPOLATION POLYNOMIALS AND                !
!             THEIR FIRST AND SECOND ORDER DERIVATES                          !
!      p(i,j) =           L_j(x_i)                                            !
!     dp(i,j) = d/dx      L_j(x_i)                                            !
!    ddp(i,j) = d^2/dx^2  L_j(x_i)                                            !
!-----------------------------------------------------------------------------!
			call lgngr(f,df,ddf,x,y,np,np,"all")

			d_mat(:,:) = 0.0d0
			do jj=1,np		! fun-index
				d_mat( 1,jj) =  f( 1, 1) * df( 1,jj)
				d_mat(np,jj) = -f(np,np) * df(np,jj)
				do ii=1,np	! point-index
					if(ii==1)then
						d_mat(ii,jj) = d_mat(ii,jj) + ddf(ii,jj) * grid%ww_fe(np-1,1)
					else
						d_mat(ii,jj) = d_mat(ii,jj) + ddf(ii,jj) * grid%ww_fe(ii-1,1)
					endif
				enddo
			enddo
			df_op%fe_region(ii_fe)%matrix(:,:) = d_mat(2:np, 2:np)

		else !jj_fe/=1
			np = grid%Nfun
			if(allocated(x))deallocate(x)
			if(allocated(y))deallocate(y)
			if(allocated(f))deallocate(f)
			if(allocated(df)) deallocate(df)
			if(allocated(ddf))deallocate(ddf)
			if(allocated(d_mat))deallocate(d_mat)
			allocate(x(1:np), y(1:np))
			allocate(f(1:np,1:np), df(1:np,1:np), ddf(1:np,1:np))
			allocate(d_mat(1:np,1:np))

			x(:) = grid%xx_fe(:, ii_fe)
			y(:) = x(:)

			call lgngr(f,df,ddf,x,y,np,np,"all")

			d_mat(:,:) = 0.0d0
			do jj=1,np		! fun-index
				d_mat(1,jj)     =  f( 1, 1) * df( 1,jj)
				d_mat(N_fun,jj) = -f(np,np) * df(np,jj)
				do ii=1,np	! point-index
					d_mat(ii,jj) = d_mat(ii,jj) + ddf(ii,jj) * grid%ww_fe(ii,ii_fe)
				enddo
			enddo
			df_op%fe_region(ii_fe)%matrix(:,:) = d_mat(:,:)

		endif ! jj_fe ? 1
	enddo ! ii_fe


	np = grid%Nfun

	do ii_fe = 1,grid%Nfe
		do jj=1,np
			do ii=1,np
				df_op%fe_region(ii_fe)%matrix(ii,jj) = 0.5 * (df_op%fe_region(ii_fe)%matrix(ii,jj) + df_op%fe_region(ii_fe)%matrix(jj,ii))
				df_op%fe_region(ii_fe)%matrix(jj,ii) = df_op%fe_region(ii_fe)%matrix(ii,jj)
			enddo
		enddo
	enddo


	do ii_fe = 1,grid%Nfe
		jj_fe = grid%fe_range(1) - 1 + ii_fe
		do ii = 1, np
			if(ii==1 .and. jj_fe/=1)then
				df_op%fe_region(ii_fe)%matrix(ii,:) = df_op%fe_region(ii_fe)%matrix(ii,:) / (grid%ww_lo(ii_fe) + grid%ww_fe(1,ii_fe))
			elseif(ii==np .and. jj_fe/=N_fe)then
				df_op%fe_region(ii_fe)%matrix(ii,:) = df_op%fe_region(ii_fe)%matrix(ii,:) / (grid%ww_fe(np,ii_fe) + grid%ww_up(ii_fe))
			else
				df_op%fe_region(ii_fe)%matrix(ii,:) = df_op%fe_region(ii_fe)%matrix(ii,:) / grid%ww_fe(ii,ii_fe)
			endif
		enddo
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------

SUBROUTINE set_2nd_deriv_op(ddf_op, grid)	!OK
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(OUT) :: ddf_op
	TYPE(fedvr_grid),    INTENT(IN)  :: grid
	REAL(DP), DIMENSION(:),   ALLOCATABLE :: x, y
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: f, df, ddf, dd_mat
	INTEGER(I4B) :: ii, jj, np

	call Op_create(ddf_op, grid, 0, 0)
	call Op_allocate(ddf_op)

	do ii_fe = 1,grid%Nfe						! local
		jj_fe = grid%fe_range(1) - 1 + ii_fe	! global
		if(jj_fe == 1)then
			np = grid%Nfun+1
			if(allocated(x))deallocate(x)
			if(allocated(y))deallocate(y)
			if(allocated(f))deallocate(f)
			if(allocated(df)) deallocate(df)
			if(allocated(ddf))deallocate(ddf)
			if(allocated(dd_mat))deallocate(dd_mat)
			allocate(x(1:np), y(1:np))
			allocate(f(1:np,1:np), df(1:np,1:np), ddf(1:np,1:np))
			allocate(dd_mat(1:np,1:np))

			x(1) = 0.0d0
			x(2:np) = grid%xx_fe(1:np-1, ii_fe)
			y(:) = x(:)
!-----------------------------------------------------------------------------!
!         CALCULATE THE LAGRANGE INTERPOLATION POLYNOMIALS AND                !
!             THEIR FIRST AND SECOND ORDER DERIVATES                          !
!      p(i,j) =           L_j(x_i)                                            !
!     dp(i,j) = d/dx      L_j(x_i)                                            !
!    ddp(i,j) = d^2/dx^2  L_j(x_i)                                            !
!-----------------------------------------------------------------------------!
			call lgngr(f,df,ddf,x,y,np,np,"all")

			dd_mat(:,:) = 0.0d0
			do jj=1,np		! fun-index
				dd_mat( 1,jj) =  f( 1, 1) * df( 1,jj)
				dd_mat(np,jj) = -f(np,np) * df(np,jj)
				do ii=1,np	! point-index
					if(ii==1)then
						dd_mat(ii,jj) = dd_mat(ii,jj) + ddf(ii,jj) * grid%ww_fe(np-1,1)
					else
						dd_mat(ii,jj) = dd_mat(ii,jj) + ddf(ii,jj) * grid%ww_fe(ii-1,1)
					endif
				enddo
			enddo
			ddf_op%fe_region(ii_fe)%matrix(:,:) = dd_mat(2:np, 2:np)

		else !jj_fe/=1
			np = grid%Nfun
			if(allocated(x))deallocate(x)
			if(allocated(y))deallocate(y)
			if(allocated(f))deallocate(f)
			if(allocated(df)) deallocate(df)
			if(allocated(ddf))deallocate(ddf)
			if(allocated(dd_mat))deallocate(dd_mat)
			allocate(x(1:np), y(1:np))
			allocate(f(1:np,1:np), df(1:np,1:np), ddf(1:np,1:np))
			allocate(dd_mat(1:np,1:np))

			x(:) = grid%xx_fe(:, ii_fe)
			y(:) = x(:)

			call lgngr(f,df,ddf,x,y,np,np,"all")

			dd_mat(:,:) = 0.0d0
			do jj=1,np		! fun-index
				dd_mat(1,jj)     =  f( 1, 1) * df( 1,jj)
				dd_mat(N_fun,jj) = -f(np,np) * df(np,jj)
				do ii=1,np	! point-index
					dd_mat(ii,jj) = dd_mat(ii,jj) + ddf(ii,jj) * grid%ww_fe(ii,ii_fe)
				enddo
			enddo
			ddf_op%fe_region(ii_fe)%matrix(:,:) = dd_mat(:,:)

		endif ! jj_fe ? 1
	enddo ! ii_fe


	np = grid%Nfun

	do ii_fe = 1,grid%Nfe
		do jj=1,np
			do ii=1,np
				ddf_op%fe_region(ii_fe)%matrix(ii,jj) = 0.5 * (ddf_op%fe_region(ii_fe)%matrix(ii,jj) + ddf_op%fe_region(ii_fe)%matrix(jj,ii))
				ddf_op%fe_region(ii_fe)%matrix(jj,ii) = ddf_op%fe_region(ii_fe)%matrix(ii,jj)
			enddo
		enddo
	enddo


	do ii_fe = 1,grid%Nfe
		jj_fe = grid%fe_range(1) - 1 + ii_fe
		do ii = 1, np
			if(ii==1 .and. jj_fe/=1)then
				ddf_op%fe_region(ii_fe)%matrix(ii,:) = ddf_op%fe_region(ii_fe)%matrix(ii,:) / (grid%ww_lo(ii_fe) + grid%ww_fe(1,ii_fe))
			elseif(ii==np .and. jj_fe/=N_fe)then
				ddf_op%fe_region(ii_fe)%matrix(ii,:) = ddf_op%fe_region(ii_fe)%matrix(ii,:) / (grid%ww_fe(np,ii_fe) + grid%ww_up(ii_fe))
			else
				ddf_op%fe_region(ii_fe)%matrix(ii,:) = ddf_op%fe_region(ii_fe)%matrix(ii,:) / grid%ww_fe(ii,ii_fe)
			endif
		enddo
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------

END MODULE
















































