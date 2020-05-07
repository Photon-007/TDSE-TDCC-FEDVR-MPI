MODULE FEDVR_potentials
	USE nrtype
	USE gvars
	USE FEDVR_module
	USE Operator_typedef


	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE set_centrifugal_ops(VV_op, DV_op, grid, ll)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(INOUT) :: VV_op, DV_op
	TYPE(fedvr_grid),    INTENT(IN)    :: grid
	INTEGER(I4B), INTENT(IN) :: ll
	INTEGER(I4B) :: ii
	
	call Op_create(VV_op, grid, ll, 0)
	call Op_create(DV_op, grid, ll, 0)
	call Op_allocate(VV_op)
	call Op_allocate(DV_op)

	do ii=1,grid%Ngp
		VV_op%diag(ii) =  1.0d0 * ll*(ll+1) / (2.0d0 * grid%xx(ii)**2)
		DV_op%diag(ii) = -1.0d0 * ll*(ll+1) / (1.0d0 * grid%xx(ii)**3)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_core_potential_ops(VV_op, DV_op, grid)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(INOUT) :: VV_op, DV_op
	TYPE(fedvr_grid), INTENT(IN)  :: grid
	INTEGER(I4B) :: ii
	REAL(DP), DIMENSION(6) :: aa_ng
	REAL(DP) :: aa, bb, cc, rr


	call Op_create(VV_op, grid, 0, 0)
	call Op_create(DV_op, grid, 0, 0)
	call Op_allocate(VV_op)
	call Op_allocate(DV_op)

	select case (target_atom)
		case ('H')
			aa_ng(:) = (/ 0.000d0, 0.000d0,   0.000d0, 0.000d0,  0.000d0, 0.000d0 /)
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('H_like')
			aa_ng(:) = (/ 0.000d0, 0.000d0,   0.000d0, 0.000d0,  0.000d0, 0.000d0 /)
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('He')				! J. Phys. B, 38 2593-2600 (2005)
			aa_ng(:) = (/ 1.231d0, 0.662d0,  -1.325d0, 1.236d0, -0.231d0, 0.480d0 /)
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('Ne')				! J. Phys. B, 38 2593-2600 (2005)
			aa_ng(:) = (/ 8.069d0, 2.148d0,  -3.570d0, 1.986d0,  0.931d0, 0.602d0 /)
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('Ar')				! J. Phys. B, 38 2593-2600 (2005)
			aa_ng(:) = (/16.039d0, 2.007d0, -25.543d0, 4.525d0,  0.961d0, 0.443d0 /)	
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('Kr')				! Need to find the coeffs!!!
			aa_ng(:) = (/ 0.000d0, 0.000d0,   0.000d0, 0.000d0,  0.000d0, 0.000d0 /)
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('Xe')				! Phys. Rev. A 84 043420 (2011)
			aa_ng(:) = (/51.356d0, 2.112d0, -99.927d0, 3.737d0,  1.644d0, 0.431d0 /)
			call get_potential_ng(Z_core, grid%xx, aa_ng, VV_op%diag, DV_op%diag)

		case ('Li')				! Klapisch model potential PRA 100 (2019) 023407, Comput. Phys. Commun.2 (1971) 239, J. Phys. B 32 (1999) 5639
			aa = 7.90875d0
			bb = 3.90006d0
			cc = 10.3210d0
			do ii=1,grid%Ngp
				rr = grid%xx(ii)
				VV_op%diag(ii) = -1.0d0 * (1.0d0 + (Z_core-1.0d0)*dexp(-aa*rr) + cc*rr*dexp(-bb*rr)) / rr
				DV_op%diag(ii)= (1.0d0 + (Z_core-1.0d0)*dexp(-aa*rr) + cc*rr*dexp(-bb*rr)) / rr**2 + (aa*(Z_core-1.0d0)*dexp(-aa*rr) + bb*cc*rr*dexp(-bb*rr) - cc*dexp(-bb*rr)) / rr
			enddo

		case ('H2p')			! diatomic molecules modelled by a central (spherical) potential.....
			do ii=1,grid%Ngp
				rr = grid%xx(ii)
				if(rr.gt.RR_core/2.0d0)then
					VV_op%diag(ii) = -1.0d0 * Z_core / rr
				else
					VV_op%diag(ii) = -2.0d0 * Z_core / RR_core
				endif
			enddo

		case ('H2p_SNZ')		! central model potential for H2 and H2p based on Journal of Modern Optics 55 (2008) 2667
			do ii=1,grid%Ngp
				rr = grid%xx(ii)
				VV_op%diag(ii) = -1.0d0 * Z_core * (1.0d0 + SIGN(1.0d0,ascreen)*dexp(-2.0d0*rr/dsqrt(dabs(ascreen)))) / rr
			enddo

		case default
			write(6,*)"ERROR!!! \'target_atom\' not recognized: ", target_atom
			STOP
	end select

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE get_potential_ng(ZZ, rr, aa, VV, DV)
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: ZZ
	REAL(DP), DIMENSION(:), INTENT(IN) :: rr
	REAL(DP), DIMENSION(6), INTENT(IN) :: aa
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: VV, DV
	INTEGER(I4B) :: ii, lb, ub

	lb = lbound(VV,1,I4B)
	ub = ubound(VV,1,I4B)
	do ii =lb,ub
		VV(ii) = -(ZZ + aa(1)*dexp(-aa(2)*rr(ii)) + aa(3)*rr(ii)*dexp(-aa(4)*rr(ii)) + aa(5)*dexp(-aa(6)*rr(ii))) / rr(ii)
		DV(ii) =  (ZZ + aa(1)*(1.0d0+aa(2)*rr(ii))*dexp(-aa(2)*rr(ii)) + aa(3)*aa(4)*dexp(-aa(4)*rr(ii))*rr(ii)**2 + aa(5)*(1.0d0+aa(6)*rr(ii))*dexp(-aa(6)*rr(ii))) / (rr(ii)**2)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_abspot(VA_op, grid)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(INOUT) :: VA_op
	TYPE(fedvr_grid),    INTENT(IN)    :: grid
	INTEGER(I4B) :: ii
	REAL(DP) :: arg

	abspot_r0 = xmax - (xmax * abspot_frac / 100.0d0)
!	write(6,*) my_id, abspot_r0, xmax

	call Op_create(VA_op, grid, 0, 0)
	call Op_allocate(VA_op)

	do ii = 1,VA_op%Ngp
		if(grid%xx(ii) >= abspot_r0)then
			arg = (grid%xx(ii) - abspot_r0) / (xmax - abspot_r0)
			VA_op%diag(ii) = -abspot_alpha * dlog(dcos(arg))
		endif
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE get_maskfun(VM_op, VA_op, dt)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(INOUT) :: VM_op
	TYPE(fedvr_op_diag), INTENT(IN)    :: VA_op
	REAL(DP) :: dt

	VM_op%diag(:) = dexp(-VA_op%diag(:)*dt)

END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE FEDVR_potentials















