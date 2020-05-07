MODULE gvars
	USE nrtype
	USE unit_conversions
	USE MPI_pars
	USE mpi
	IMPLICIT NONE

	INTERFACE write_params
		MODULE PROCEDURE write_params_tounit
		MODULE PROCEDURE write_params_tofile
	END INTERFACE

	TYPE pulse
		CHARACTER(LEN=10) :: Envelope
		REAL(DP) :: EA
		REAL(DP) :: t0
		REAL(DP) :: Omega
		REAL(DP) :: param
		REAL(DP) :: CEP
		REAL(DP) :: t_start
		REAL(DP) :: t_end
	END TYPE

!	%===============================%
!	|		GENERAL SWITCHES		|
!	%===============================%

	LOGICAL(LGT) :: do_chp_calcs       = .FALSE.
	LOGICAL(LGT) :: linearly_pol_laser = .FALSE.
	LOGICAL(LGT) :: memory_save        = .FALSE.
	LOGICAL(LGT) :: restart
	LOGICAL(LGT) :: only_initial_wf
	LOGICAL(LGT) :: init_by_ITP
	LOGICAL(LGT) :: do_diagH0
	LOGICAL(LGT) :: do_pulse
	LOGICAL(LGT) :: do_propag
	LOGICAL(LGT) :: do_hhg
	LOGICAL(LGT) :: do_autocorr
	namelist / general_switches / do_chp_calcs, linearly_pol_laser, memory_save, restart, only_initial_wf, init_by_ITP, do_diagH0, do_pulse, do_propag, do_hhg, do_autocorr

!	%===============================%
!	|	NUMERICAL GRID PARAMETERS	|
!	%===============================%

	INTEGER(I4B) :: N_fun					! the number of quadrature points inside a finite element
	INTEGER(I4B) :: N_fe					! number of the finite elements
	INTEGER(I4B) :: N_tfe					! the number of transient finite elements
	REAL(DP) :: dx_min						! the smalles finite element
	REAL(DP) :: dx_max						! the largest finite element
	REAL(DP) :: xmax						! the size of the R simulation box
	INTEGER(I4B) :: N_gp					! the total number of radial gridpoints
	INTEGER(I4B) :: NN_gp, NN_fe, NN_fun	! local values of each process

	CHARACTER(LEN=4) :: SphericalHarmonicRepresentation
	INTEGER(I4B) :: lmax
	INTEGER(I4B) :: tagmax
	INTEGER(I4B) :: NN_tag
	namelist / grid_parameters / N_fun, N_fe, N_tfe, dx_min, dx_max, lmax

	INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: N_fe_arr
	INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: N_gp_arr
	INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: N_tag_arr
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: Idx_fe
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: Idx_gp
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: Idx_tag
	INTEGER(I4B), DIMENSION(2) :: my_fe, my_gp, my_tag

	INTEGER(I4B) :: ii_fe, jj_fe			! local and global index of the FE-s 

!	#=======================#
!	|	ABSORBING POTENTIAL	|
!	#=======================#
	REAL(DP) :: abspot_r0
	REAL(DP) :: abspot_alpha
	REAL(DP) :: abspot_frac
	namelist / abspot_params / abspot_alpha, abspot_frac


!	%===========================================%
!	|		INITIAL STATE (Default: H 1s)		|
!	%===========================================%
	CHARACTER(LEN=10) :: target_atom = 'H'
	INTEGER(I4B) :: n_init = 1
	INTEGER(I4B) :: l_init = 0
	INTEGER(I4B) :: m_init = 0
	INTEGER(I4B) :: tag_init
	REAL(DP) :: Z_core  = 1.0d0
	REAL(DP) :: RR_core = 2.0d0
	REAL(DP) :: ascreen = -0.01d0
	REAL(DP) :: Energy_sys
	namelist / initialization_params / target_atom, Z_core, n_init, l_init, m_init, RR_core, ascreen

!	%===============================%
!	|		LASER PARAMETERS		|
!	%===============================%
	CHARACTER(LEN=30), DIMENSION(99) :: Pulse_files = ''
	INTEGER(I4B) :: Npulse
	namelist / laser_input / Pulse_files

	TYPE(pulse), DIMENSION(:), ALLOCATABLE :: LaserPulse

	TYPE(phys_variable) :: Intensity
	TYPE(phys_variable) :: Amplitude
	TYPE(phys_variable) :: PhotonEnergy
	TYPE(phys_variable) :: Wavelength
	TYPE(phys_variable) :: FWHM
	TYPE(phys_variable) :: Duration
	TYPE(phys_variable) :: Time_center
	TYPE(phys_variable) :: CEP
	TYPE(phys_variable) :: Oscillations 	! unit: 'FWHM', 'Duration'
	CHARACTER(LEN=10)   :: Envelope_type
	namelist / laser_params / Intensity, Amplitude, PhotonEnergy, WaveLength, CEP, FWHM, Duration, Time_center, Oscillations, Envelope_type


	REAL(DP), DIMENSION(:), ALLOCATABLE :: II0_au, II0_Wcm2				! peak intensity of the pulse
	REAL(DP), DIMENSION(:), ALLOCATABLE :: E0_au, E0_Vcm  				! amplitude of the electric field
	REAL(DP), DIMENSION(:), ALLOCATABLE :: A0_au						! amplitude of the vector potential
	REAL(DP), DIMENSION(:), ALLOCATABLE :: freq_au,   freq_Hz			! central frequency
	REAL(DP), DIMENSION(:), ALLOCATABLE :: omega_au,  omega_eV			! angular frequency, i.e. central energy of the pulse
	REAL(DP), DIMENSION(:), ALLOCATABLE :: WvLen_au,  WvLen_nm			! wavelength 
	REAL(DP), DIMENSION(:), ALLOCATABLE :: period_au, period_fs
	REAL(DP), DIMENSION(:), ALLOCATABLE :: t0_au,     t0_fs				! time moment of the center of the pulse
	REAL(DP), DIMENSION(:), ALLOCATABLE :: CEP_rad
	REAL(DP), DIMENSION(:), ALLOCATABLE :: FWHM_au
	REAL(DP), DIMENSION(:), ALLOCATABLE :: FWHM_fs
	REAL(DP), DIMENSION(:), ALLOCATABLE :: Duration_au					! pulse-duration
	REAL(DP), DIMENSION(:), ALLOCATABLE :: Duration_fs					! pulse-duration
	REAL(DP), DIMENSION(:), ALLOCATABLE :: swiss_knife					! universal parameter defining the pulse used throughout the propagation (saving an if construct each tprop step). Either: duration_au or FWHM_au
	REAL(DP), DIMENSION(:), ALLOCATABLE :: t_start, t_end				! starting and ending time moments of the pulse
	CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: Envelope_fun

!	%===============================%
!	|	PROPAGATION PARAMETERS		|
!	%===============================%

	CHARACTER(LEN=10)   :: Gauge = 'length'		! length/velocity
	TYPE(phys_variable) :: Time_max
	TYPE(phys_variable) :: Time_step_wf
	TYPE(phys_variable) :: Time_step_data
	TYPE(phys_variable) :: Delta_t
	REAL(DP) :: Epsilon_ortho = 1.0d-14
	REAL(DP) :: Epsilon_conv  = 1.0d-24
	INTEGER(I4B) :: GS_refinement = 1
	LOGICAL(LGT) :: GS_all_qvecs  = .FALSE.
	INTEGER(I4B) :: N_Krylov = 8
	namelist / propag_params / Gauge, Time_max, Time_step_wf, Time_step_data, Delta_t, Epsilon_ortho, Epsilon_conv, GS_refinement, GS_all_qvecs, N_Krylov

	REAL(DP) :: Time_max_au,       Time_max_fs
	REAL(DP) :: Time_step_wf_au,   Time_step_wf_fs
	REAL(DP) :: Time_step_data_au, Time_step_data_fs
	REAL(DP) :: dTime
	REAL(DP) :: Time
	INTEGER(I4B) :: N_step_wf, N_step_data

!	%===================================%
!	|	IMAG TIME PROPAG PARAMETERS		|
!	%===================================%

	TYPE(phys_variable) :: Time_max_ITP
	TYPE(phys_variable) :: Delta_t_ITP
	REAL(DP) :: Epsilon_conv_ITP  = 1.0d-24
	namelist / ITP_params / Time_max_ITP, Delta_t_ITP, Epsilon_conv_ITP

	REAL(DP) :: Time_max_ITP_au, Time_max_ITP_fs
	REAL(DP) :: dTime_ITP = 1.0d-2

!	#=======================#
!	|	TIMING VARIABLES	|
!	#=======================#

	INTEGER(I8B) :: N_wf_mirror = 0
	REAL(DP) :: Time_wf_mirror  = 0.0d0
	REAL(DP) :: Main_start,  Main_end
	REAL(DP) :: tStep_start, tStep_end, tStep_elapsed


	CONTAINS
!===============================================================================
!-------------------------------------------------------------------------------
SUBROUTINE set_propag_params()
	IMPLICIT NONE

	if(Time_max%units=='')then
		write(6,*)" ERROR!!! You should specify \'Time_max\' (no defalt value)!"
		write(6,*)" Available units: \'fs\', \'au\'"
		STOP
	else
		select case (Time_max%units)
			case ('au')
				Time_max_au = Time_max%value
				Time_max_fs = Atomic2SI(Time_max%value, 'fs')
			case ('fs')
				Time_max_fs = Time_max%value
				Time_max_au = SI2Atomic(Time_max%value, 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Time_max\': ", Time_max%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

	if(Time_step_wf%units=='')then
		write(6,*)" ERROR!!! You should specify \'Time_step_wf\' (no defalt value)!"
		write(6,*)" Available units: \'fs\', \'au\'"
		STOP
	else
		select case (Time_step_wf%units)
			case ('au')
				Time_step_wf_au = Time_step_wf%value
				Time_step_wf_fs = Atomic2SI(Time_step_wf%value, 'fs')
			case ('fs')
				Time_step_wf_fs = Time_step_wf%value
				Time_step_wf_au = SI2Atomic(Time_step_wf%value, 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Time_step_wf\': ", Time_step_wf%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

	if(Time_step_data%units=='')then
		write(6,*)" ERROR!!! You should specify \'Time_step_data\' (no defalt value)!"
		write(6,*)" Available units: \'fs\', \'au\'"
		STOP
	else
		select case (Time_step_data%units)
			case ('au')
				Time_step_data_au = Time_step_data%value
				Time_step_data_fs = Atomic2SI(Time_step_data%value, 'fs')
			case ('fs')
				Time_step_data_fs = Time_step_data%value
				Time_step_data_au = SI2Atomic(Time_step_data%value, 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Time_step_data\': ", Time_step_data%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

	if(Delta_t%units=='')then
		write(6,*)" ERROR!!! You should specify \'Delta_t\' (no defalt value)!"
		write(6,*)" Available units: \'fs\', \'au\'"
		STOP
	else
		select case (Delta_t%units)
			case ('au')
				dTime = Delta_t%value
			case ('fs')
				dTime = SI2Atomic(Delta_t%value, 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Delta_t\': ", Delta_t%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

	N_step_wf   = int(Time_max_au/Time_step_wf_au)   + 1
	N_step_data = int(Time_max_au/Time_step_data_au) + 1

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_ITP_params()
	IMPLICIT NONE

	if(Time_max_ITP%units=='')then
		write(6,*)" ERROR!!! You should specify \'Time_max_ITP\' (no defalt value)!"
		write(6,*)" Available units: \'fs\', \'au\'"
		STOP
	else
		select case (Time_max_ITP%units)
			case ('au')
				Time_max_ITP_au = Time_max_ITP%value
				Time_max_ITP_fs = Atomic2SI(Time_max_ITP%value, 'fs')
			case ('fs')
				Time_max_ITP_fs = Time_max_ITP%value
				Time_max_ITP_au = SI2Atomic(Time_max_ITP%value, 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Time_max_ITP\': ", Time_max_ITP%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

	if(Delta_t_ITP%units=='')then
		write(6,*)" ERROR!!! You should specify \'Delta_t_ITP\' (no defalt value)!"
		write(6,*)" Available units: \'fs\', \'au\'"
		STOP
	else
		select case (Delta_t%units)
			case ('au')
				dTime_ITP = Delta_t_ITP%value
			case ('fs')
				dTime_ITP = SI2Atomic(Delta_t_ITP%value, 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Delta_t_ITP\': ", Delta_t_ITP%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_pulse_params(ii)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ii

	if((Intensity%units=='').and.(Amplitude%units==''))then
		write(6,*)' ERROR!!! You should specify either the pulse Intensity or Amplitude!'
		STOP
	endif
	if((PhotonEnergy%units=='').and.(Wavelength%units==''))then
		write(6,*)' ERROR!!! You should specify either the PhotonEnergy or the Wavelength!'
		STOP
	endif
	if((Duration%units=='').and.(FWHM%units=='').and.(Oscillations%units==''))then
		write(6,*)' ERROR!!! You should specify either the puled Duration, FWHM or Oscillations!'
		STOP
	endif
	if(Time_center%units=='')then
		write(6,*)' ERROR!!! You should specify the time moment of the center of the pulse!'
		STOP
	endif


	if(Intensity%units /= '')then
		select case (Intensity%units)
			case ('W/cm2')
				II0_Wcm2(ii) = Intensity%value
				II0_au(ii)   = SI2Atomic(Intensity%value * 100**2, 'W/m2')
				E0_Vcm(ii)   = dsqrt(II0_Wcm2(ii) / convert_I2E) !/ 100
				E0_au(ii)    = SI2Atomic(E0_Vcm(ii) * 100, 'V/m') 
			case ('W/m2')
				II0_Wcm2(ii) = Intensity%value / 100**2
				II0_au(ii)   = SI2Atomic(Intensity%value , 'W/m2')
				E0_Vcm(ii)   = dsqrt(II0_Wcm2(ii) / convert_I2E)
				E0_au(ii)    = SI2Atomic(E0_Vcm(ii) * 100, 'V/m') 
			case ('au')	!!!NOT WORKING PROPERLY!!!
				II0_au(ii)   = Intensity%value
				II0_Wcm2(ii) = Atomic2SI(Intensity%value, 'W/m2') / 100**2
				E0_Vcm(ii)   = dsqrt(II0_Wcm2(ii) / convert_I2E)
				E0_au(ii)    = SI2Atomic(E0_Vcm(ii) * 100, 'V/m') 
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Intensity\': ", Intensity%units
				write(6,*)" Available options: \'W/cm2\', \'W/m2\',\'au\'"
				STOP
		end select
	elseif(Amplitude%units /= '')then
		select case (Amplitude%units)
			case ('V/cm')
				E0_Vcm(ii)   = Amplitude%value
				E0_au(ii)    = SI2Atomic(Amplitude%value * 100, 'V/m')
				II0_Wcm2(ii) = convert_I2E * E0_Vcm(ii)**2
				II0_au(ii)   = SI2Atomic(II0_Wcm2(ii) * 100**2, 'W/m2')
			case ('V/m')
				E0_Vcm(ii)   = Amplitude%value / 100
				E0_au(ii)    = SI2Atomic(Amplitude%value , 'V/m')
				II0_Wcm2(ii) = convert_I2E * E0_Vcm(ii)**2
				II0_au(ii)   = SI2Atomic(II0_Wcm2(ii) * 100**2, 'W/m2')
			case ('au')
				E0_au(ii)    = Amplitude%value
				E0_Vcm(ii)   = Atomic2SI(Amplitude%value, 'V/m') / 100
				II0_Wcm2(ii) = convert_I2E * E0_Vcm(ii)**2
				II0_au(ii)   = SI2Atomic(II0_Wcm2(ii) * 100**2, 'W/m2')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Amplitude\': ", Amplitude%units
				write(6,*)" Available options: \'V/cm\', \'V/m\',\'au\'"
				STOP
		end select
	endif

	if(Wavelength%units /= '')then
		select case (Wavelength%units)
			case ('m')
				WvLen_au(ii)  = SI2Atomic(Wavelength%value, 'm')
				WvLen_nm(ii)  = Wavelength%value * 1.0d9
				freq_Hz(ii)   = constSI_lightSpeed / (WvLen_nm(ii) * 1.0d-9)
				freq_au(ii)   = SI2Atomic(freq_Hz(ii), 'Hz')
				omega_au(ii)  = TWOPI * freq_au(ii)
				omega_eV(ii)  = constSI_hbar * freq_Hz(ii) * convert_J2eV * TWOPI
				period_au(ii) = 1.0 / freq_au(ii)
				period_fs(ii) = 1.0 / freq_Hz(ii) * 1.0d15
			case ('nm')
				WvLen_au(ii)  = SI2Atomic(Wavelength%value, 'nm')
				WvLen_nm(ii)  = Wavelength%value 
				freq_Hz(ii)   = constSI_lightSpeed / (WvLen_nm(ii) * 1.0d-9)
				freq_au(ii)   = SI2Atomic(freq_Hz(ii), 'Hz')
				omega_au(ii)  = TWOPI * freq_au(ii)
				omega_eV(ii)  = constSI_hbar * freq_Hz(ii) / convert_J2eV * TWOPI		! why do I need the 2pi???
				period_au(ii) = 1.0 / freq_au(ii)
				period_fs(ii) = 1.0 / freq_Hz(ii) * 1.0d15
			case ('au')
				WvLen_au(ii)  = Wavelength%value
				WvLen_nm(ii)  = Atomic2SI(Wavelength%value, 'nm')
				freq_Hz(ii)   = constSI_lightSpeed / (WvLen_nm(ii) * 1.0d-9)
				freq_au(ii)   = SI2Atomic(freq_Hz(ii), 'Hz')
				omega_au(ii)  = TWOPI * freq_au(ii)
				omega_eV(ii)  = constSI_hbar * freq_Hz(ii) * convert_J2eV * TWOPI
				period_au(ii) = 1.0 / freq_au(ii)
				period_fs(ii) = 1.0 / freq_Hz(ii) * 1.0d15
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Wavelength\': ", Wavelength%units
				write(6,*)" Available options: \'m\', \'nm\',\'au\'"
				STOP
		end select
	elseif(PhotonEnergy%units /= '')then
		select case (PhotonEnergy%units)
			case ('J')
				omega_au(ii)  = SI2Atomic(PhotonEnergy%value, 'J')
				omega_eV(ii)  = PhotonEnergy%value / convert_J2eV
				freq_au(ii)   = omega_au(ii) / TWOPI
				freq_Hz(ii)   = Atomic2SI(freq_au(ii), 'Hz')
				period_au(ii) = 1.0 / freq_au(ii)
				period_fs(ii) = 1.0 / freq_Hz(ii) * 1.0d15
				WvLen_nm(ii)  = constSI_lightSpeed / freq_Hz(ii) * 1.0d9
				WvLen_au(ii)  = SI2Atomic(WvLen_nm(ii), 'nm')
			case ('eV')
				omega_au(ii)  = SI2Atomic(PhotonEnergy%value * convert_J2eV, 'J')
				omega_eV(ii)  = PhotonEnergy%value
				freq_au(ii)   = omega_au(ii) / TWOPI
				freq_Hz(ii)   = Atomic2SI(freq_au(ii), 'Hz')
				period_au(ii) = 1.0 / freq_au(ii)
				period_fs(ii) = 1.0 / freq_Hz(ii) * 1.0d15
				WvLen_nm(ii)  = constSI_lightSpeed / freq_Hz(ii) * 1.0d9
				WvLen_au(ii)  = SI2Atomic(WvLen_nm(ii), 'nm')
			case ('au')
				omega_au(ii)  = PhotonEnergy%value
				omega_eV(ii)  = Atomic2SI(PhotonEnergy%value, 'J') / convert_J2eV
				freq_au(ii)   = omega_au(ii) / TWOPI
				freq_Hz(ii)   = Atomic2SI(freq_au(ii), 'Hz')
				period_au(ii) = 1.0 / freq_au(ii)
				period_fs(ii) = 1.0 / freq_Hz(ii) * 1.0d15
				WvLen_nm(ii)  = constSI_lightSpeed / freq_Hz(ii) * 1.0d9
				WvLen_au(ii)  = SI2Atomic(WvLen_nm(ii), 'nm')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'PhotonEnergy\': ", PhotonEnergy%units
				write(6,*)" Available options: \'J\', \'eV\',\'au\'"
				STOP
		end select
	endif

	if(Time_center%units /= '')then
		select case (Time_center%units)
			case ('fs')
				t0_fs(ii) = Time_center%value
				t0_au(ii) = SI2Atomic(t0_fs(ii) * 1.0d-15, 's')
			case ('au')
				t0_au(ii) = Time_center%value
				t0_fs(ii) = Atomic2SI(t0_au(ii), 'fs') 
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Time_center\': ", Time_center%units
				write(6,*)" Available options: \'fs\',\'au\'"
				STOP
		end select
	endif

	if(CEP%units /= '')then
		select case (CEP%units)
			case ('rad')
				CEP_rad(ii) = CEP%value
			case ('deg')
				CEP_rad(ii) = CEP%value * PI / 180.0d0
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'CEP\': ", CEP%units,"|"
				write(6,*)" Available options: \'rad\', \'deg\'"
				STOP
		end select
	endif

	if(Duration%units /= '')then
		select case (Duration%units)
			case ('fs')
				Duration_fs(ii) = Duration%value
				Duration_au(ii) = SI2Atomic(Duration_fs(ii), 'fs')
			case ('au')
				Duration_au(ii) = Duration%value
				Duration_fs(ii) = Atomic2SI(Duration_au(ii), 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Duration\': ", Duration%units
				write(6,*)" Available options: \'fs\', \'au\'"
				STOP
		end select
	elseif(FWHM%units /= '')then
		select case (FWHM%units)
			case ('fs')
				FWHM_fs(ii) = FWHM%value
				FWHM_au(ii) = SI2Atomic(FWHM_fs(ii), 'fs')
			case ('au')
				FWHM_au(ii) = FWHM%value
				FWHM_fs(ii) = Atomic2SI(FWHM_au(ii), 'fs')
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'FWHM\': ", FWHM%units
				write(6,*)" Available options: \'fs\', \'au\'"
				STOP
		end select
	elseif(Oscillations%units /= '')then
		select case (Oscillations%units)
			case ('Duration')
				Duration_au(ii) = Oscillations%value * period_au(ii)
				Duration_fs(ii) = Oscillations%value * period_fs(ii)
			case ('FWHM')
				FWHM_au(ii) = Oscillations%value * period_au(ii)
				FWHM_fs(ii) = Oscillations%value * period_fs(ii)
			case default
				write(6,*)" ERROR!!! Unknown unit for the \'Oscillations\': ", Oscillations%units
				write(6,*)" Available options: \'Duration\', \'FWHM\'"
				STOP
		end select
	endif

	select case (Envelope_type)
		case ('sin2')
			Envelope_fun(ii) = 'sin2'
!			write(6,*)" I'm here! "
			if(Duration_au(ii)/=0.0d0)then
				swiss_knife(ii) = PI / Duration_au(ii)
				t_start(ii) = t0_au(ii) - 0.5d0 * Duration_au(ii) 
				t_end(ii)   = t0_au(ii) + 0.5d0 * Duration_au(ii) 
			else
				write(6,*)" ERROR!!! \'Duration\' is not specified for a sin^2  envelope!"
!				STOP
			endif
		case ('cos2')
			Envelope_fun(ii) = 'cos2'
			if(FWHM_au(ii)/=0.0d0)then
				swiss_knife(ii) = 1.14372 / FWHM_au(ii)
				t_start(ii) = t0_au(ii) - 0.5 * PI * FWHM_au(ii) / 1.14372d0 
				t_end(ii)   = t0_au(ii) + 0.5 * PI * FWHM_au(ii) / 1.14372d0 
			else
				write(6,*)" ERROR!!! \'FWHM\' is not specified for a cos^2  envelope!"
			endif
		case ('gauss')
			Envelope_fun(ii) = 'gauss'
			if(FWHM_au(ii)/=0.0d0)then
				swiss_knife(ii) = -1.0d0 / (2.0d0 * (FWHM_au(ii) / dsqrt(4.0d0*dlog(2.0d0)))**2)
				t_start(ii) = t0_au(ii) - 1.5d0 * FWHM_au(ii)
				t_end(ii)   = t0_au(ii) + 1.5d0 * FWHM_au(ii)
			else
				write(6,*)" ERROR!!! \'FWHM\' is not specified for a gauss  envelope!"
			endif
		case ('flat')
			Envelope_fun(ii) = 'flat'
			swiss_knife(ii) = 1.0d0
			if(FWHM_au(ii) /= 0.0d0)then
				t_start(ii) = t0_au(ii) - 0.5d0 * FWHM_au(ii)
				t_end(ii)   = t0_au(ii) + 0.5d0 * FWHM_au(ii)
			elseif(Duration_au(ii) /= 0.0d0)then
				t_start(ii) = t0_au(ii) - 0.5d0 * Duration_au(ii) 
				t_end(ii)   = t0_au(ii) + 0.5d0 * Duration_au(ii) 
			endif
		case default
			write(6,*)" ERROR!!! Unknown \'Envelope_type\': ", trim(Envelope_type)
			write(6,*)" Available options: \'sin^2\', \'cos^2\',  \'gauss\' and  \'flat\' "
			STOP
	end select
	
	if (t_start(ii)<0.0d0)then
		write(6,*)" "
		write(6,*)"--------------------------------------------------------------------------------"
		write(6,*)" ERROR! The pulse is switched on befor the propagation starts!"
		write(6,*)" t_lim: [",t_start(ii),":",t_end(ii),"]"
		write(6,*)" Shift t_center (",t0_au(ii),") to a later time moment!"
		write(6,*)" (at least: ",t0_au(ii)+dabs(t_start(ii))," au)"
		write(6,*)"--------------------------------------------------------------------------------"
		STOP
	endif

	LaserPulse(ii)%Envelope = Envelope_fun(ii)
	LaserPulse(ii)%EA       = E0_au(ii)
	LaserPulse(ii)%t0       = t0_au(ii)
	LaserPulse(ii)%Omega    = omega_au(ii)
	LaserPulse(ii)%param    = swiss_knife(ii)
	LaserPulse(ii)%CEP      = CEP_rad(ii)
	LaserPulse(ii)%t_start  = t_start(ii)
	LaserPulse(ii)%t_end    = t_end(ii)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE set_dependent()
	IMPLICIT NONE

	xmax = N_tfe * (dx_max +dx_min) * 0.5d0 + dx_max * (N_fe - N_tfe)
	N_gp = N_fun * N_fe - (N_fe - 1)  ! the no. of total gridpoints = (no. of local gridpoints inside each finet element) X (no. of finite element) - (no. of overlaping gridpoints)
	if(do_chp_calcs)then
		SphericalHarmonicRepresentation = 'Y_lm'
		tagmax = (lmax+1) * (lmax+1)
	elseif(linearly_pol_laser)then
		SphericalHarmonicRepresentation = 'Y_l0'
		tagmax = (lmax + 1)
	else
		SphericalHarmonicRepresentation = 'Y_lm'
		tagmax = (lmax+1) * (lmax+1)
	endif

END SUBROUTINE set_dependent
!-------------------------------------------------------------------------------
SUBROUTINE divide_workload()
	IMPLICIT NONE
	INTEGER(I4B) :: id

	if(allocated(N_fe_arr)) deallocate(N_fe_arr)
	if(allocated(N_gp_arr)) deallocate(N_gp_arr)
	if(allocated(N_tag_arr))deallocate(N_tag_arr)
	if(allocated(Idx_fe)) deallocate(Idx_fe)	
	if(allocated(Idx_gp)) deallocate(Idx_gp)
	if(allocated(Idx_tag))deallocate(Idx_tag)
	allocate(N_fe_arr(0:mpi_size-1), N_gp_arr(0:mpi_size-1), N_tag_arr(0:mpi_size-1))
	allocate(Idx_fe(2,0:mpi_size-1), Idx_gp(2,0:mpi_size-1), Idx_tag(2,0:mpi_size-1))	! mpi_id as a second dimension, hence no runtime warvnings (contiguous memory locations are sent to the processes)

	call mpi_distribute(1, N_fe,   mpi_size, Idx_fe)
!	call mpi_distribute(1, tagmax, mpi_size, Idx_tag)	! everyone is working on all tags
	Idx_tag(1,:) = 1
	Idx_tag(2,:) = tagmax

	do id=0,mpi_size-1
		Idx_gp(1,id) = (Idx_fe(1,id) - 1) * (N_fun - 1) + 1
		Idx_gp(2,id) = Idx_fe(2,id) * (N_fun - 1)  + 1
		N_fe_arr(id) = Idx_fe(2,id) - Idx_fe(1,id) + 1
		N_gp_arr(id) = Idx_gp(2,id) - Idx_gp(1,id) + 1
		N_tag_arr(id)= Idx_tag(2,id)- Idx_tag(1,id)+ 1
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE read_params(param_file)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: param_file
	INTEGER(I4B) :: ii

	write(6,*)'mpi_size=', mpi_size
	write(6,*)'my_id   =', my_id
	write(6,*)"########################################################"
	write(6,*)"###     THE INPUT PARAMETERS ARE READ FROM FILE.     ###"
	open(unit=99, file = param_file, delim='apostrophe')
		read(99,nml=general_switches)
		rewind(99)
		read(99,nml=grid_parameters)
		rewind(99)
		read(99,nml=initialization_params)
		rewind(99)
		read(99,nml=ITP_params)
		rewind(99)
		read(99,nml=propag_params)
		rewind(99)
		read(99,nml=abspot_params)
		rewind(99)
		read(99,nml=laser_input)
		rewind(99)
	close(99)

	do ii = 1, size(Pulse_files)
		if(len(trim(Pulse_files(ii)))/=0)Npulse=ii
	enddo

END SUBROUTINE read_params
!------------------------------------------------------------------------------
SUBROUTINE write_params_tofile(filename, what)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: filename, what

	open(unit=99,file=filename,status="UNKNOWN")
	call write_params_tounit(99, what)
	close(99)

END SUBROUTINE write_params_tofile
!------------------------------------------------------------------------------
SUBROUTINE write_params_tounit(unit_no, what)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: unit_no
	CHARACTER(LEN=*)  :: what
	CHARACTER(LEN=21) :: ich
	INTEGER(I4B) :: ii

	if((what=='all').or.(what=='grid'))then
!		write(unit_no,*)"#######################################################"
		write(unit_no,*)"############# RADIAL GRID PARAMETERS: ##################"
		write(ich,'(I8.1)')N_fun
		write(unit_no,*)"|       N_fun = ",ich,"                  |"
		write(ich,'(I8.1)')N_fe
		write(unit_no,*)"|       N_fe  = ",ich,"                  |"
		write(ich,'(I8.1)')N_tfe
		write(unit_no,*)"|       N_tfe = ",ich,"                  |"
		write(ich,'(I8.1)')N_gp
		write(unit_no,*)"|       N_gp  = ",ich,"                  |"
		write(ich,'(F12.4)')dx_min
		write(unit_no,*)"|      dx_min = ",ich,"                  |"
		write(ich,'(F12.4)')dx_max
		write(unit_no,*)"|      dx_max = ",ich,"                  |"
!		write(ich,'(F12.3)')x_arr(N_gp)
		write(ich,'(F12.4)')xmax
		write(unit_no,*)"|     Box_dim = ",ich,"                  |"
		write(unit_no,*)"############# ANGULAR GRID PARAMETERS: #################"
		write(ich,'(I8.1)')lmax
		write(unit_no,*)"|       lmax  = ",ich,"                  |"
		write(ich,'(I8.1)')tagmax
		write(unit_no,*)"|     tagmax  = ",ich,"                  |"
		write(unit_no,*)"########################################################"
	endif
	if((what=='all').or.(what=='system'))then
!		write(unit_no,*)"########################################################"
		write(unit_no,*)"###############  INITIAL WF PARAMETERS: ################"
		write(unit_no,*)"| target atom =        ",target_atom,"                      |"
		write(ich,'(F12.4)')Z_core
		write(unit_no,*)"|      Z_core = ",ich,"                  |"
		write(ich,'(I8.1)')n_init
		write(unit_no,*)"|      n_init = ",ich,"                  |"
		write(ich,'(I8.1)')l_init
		write(unit_no,*)"|      l_init = ",ich,"                  |"
		write(ich,'(I8.1)')m_init
		write(unit_no,*)"|      m_init = ",ich,"                  |"
		write(unit_no,*)"########################################################"
	endif
	if((what=='all').or.(what=='prop'))then
!		write(unit_no,*)"########################################################"
		write(unit_no,*)"############ TIME PROPAGATION PARAMETERS: ##############"
	!	write(unit_no,*)"|      gauge  = ",gauge,"                         |"
	!	write(ich,'(F5.2)')t_max
!		write(unit_no,*)"|      t_max  = ",t_max,"          |"
	!	write(ich,'(F11.7)')dt
!		write(unit_no,*)"|         dt  = ",dt,"          |"
		write(unit_no,*)"########################################################"
	endif
	if((what=='all').or.(what=='laser'))then
		do ii=1,Npulse
!			write(unit_no,*)"#######################################################"
			write(unit_no,'(A26,I2.1,A29)')" #################  LASER(",ii,") PARAMETERS: ################"
			write(ich,'(1p,d21.14)')II0_Wcm2(ii)
			write(unit_no,*)"|        Intensity = ",ich," W/cm^2      |"
			write(ich,'(1p,d21.14)')II0_au(ii)
			write(unit_no,*)"|        Intensity = ",ich," au          |"
			write(ich,'(1p,d21.14)')E0_Vcm(ii)
			write(unit_no,*)"| El. Field Stren. = ",ich," V/cm        |"
			write(ich,'(1p,d21.14)')E0_au(ii)
			write(unit_no,*)"| El. Field Stren. = ",ich," au          |"
			write(unit_no,*)"|                                                      |"
			write(ich,'(1p,d21.14)')omega_eV(ii)
			write(unit_no,*)"|    Photon Energy = ",ich," eV          |"
			write(ich,'(1p,d21.14)')omega_au(ii)
			write(unit_no,*)"|    Photon Energy = ",ich," au          |"
			write(ich,'(1p,d21.14)')freq_Hz(ii)
			write(unit_no,*)"|        Frequency = ",ich," Hz          |"
			write(ich,'(1p,d21.14)')freq_au(ii)
			write(unit_no,*)"|        Frequency = ",ich," au          |"
			write(ich,'(1p,d21.14)')WvLen_nm(ii)
			write(unit_no,*)"|      Wave Length = ",ich," nm          |"
			write(ich,'(1p,d21.14)')WvLen_au(ii)
			write(unit_no,*)"|      Wave Length = ",ich," au          |"
			write(ich,'(1p,d21.14)')period_fs(ii)
			write(unit_no,*)"|           Period = ",ich," fs          |"
			write(ich,'(1p,d21.14)')period_au(ii)
			write(unit_no,*)"|           Period = ",ich," au          |"
			write(unit_no,*)"|                                                      |"
			write(ich,'(1p,d21.14)')t0_fs(ii)
			write(unit_no,*)"|           Center = ",ich," fs          |"
			write(ich,'(1p,d21.14)')t0_au(ii)
			write(unit_no,*)"|           Center = ",ich," au          |"
			write(ich,'(1p,d21.14)')FWHM_fs(ii)
			write(unit_no,*)"|             FWHM = ",ich," fs          |"
			write(ich,'(1p,d21.14)')FWHM_au(ii)
			write(unit_no,*)"|             FWHM = ",ich," au          |"
			write(ich,'(1p,d21.14)')Duration_fs(ii)
			write(unit_no,*)"|         Duration = ",ich," fs          |"
			write(ich,'(1p,d21.14)')Duration_au(ii)
			write(unit_no,*)"|         Duration = ",ich," au          |"
			write(ich,'(1p,d21.14)')CEP_rad(ii)
			write(unit_no,*)"|              CEP = ",ich," rad         |"
			write(ich,'(1p,d21.14)')CEP_rad(ii) * 180.0d0 / PI
			write(unit_no,*)"|              CEP = ",ich," deg         |"
			write(ich,'(1p,d21.14)')t_start(ii)
			write(unit_no,*)"|          t_start = ",ich," au          |"
			write(ich,'(1p,d21.14)')t_end(ii)
			write(unit_no,*)"|          t_end   = ",ich," au          |"
			write(unit_no,*)"########################################################"
		enddo
	endif

END SUBROUTINE write_params_tounit
!-------------------------------------------------------------------------------
SUBROUTINE read_pulse(pulse_file)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: pulse_file

	Intensity%value = 0.0d0
	Intensity%units = ''
	Amplitude%value = 0.0d0
	Amplitude%units = ''
	PhotonEnergy%value = 0.0d0
	PhotonEnergy%units = ''
	Wavelength%value = 0.0d0
	Wavelength%units = ''
	FWHM%value = 0.0d0
	FWHM%units = ''
	Duration%value = 0.0d0
	Duration%units = ''
	Time_center%value = 0.0d0
	Time_center%units = ''
	CEP%value = 0.0d0
	CEP%units = ''
	Oscillations%value = 0.0d0
	Oscillations%units = ''
	Envelope_type = ''

	open(unit=99, file = pulse_file, delim='apostrophe')
		read(99,nml=laser_params)
	close(99)
END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE allocate_pulse_params(nn)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nn

	if(allocated(LaserPulse))deallocate(LaserPulse)
	if(allocated(II0_au))deallocate(II0_au)
	if(allocated(II0_Wcm2))deallocate(II0_Wcm2)
	if(allocated(E0_au))deallocate(E0_au)
	if(allocated(E0_Vcm))deallocate(E0_Vcm)
	if(allocated(A0_au))deallocate(A0_au)
	if(allocated(freq_au))deallocate(freq_au)
	if(allocated(freq_Hz))deallocate(freq_Hz)
	if(allocated(omega_au))deallocate(omega_au)
	if(allocated(omega_eV))deallocate(omega_eV)
	if(allocated(WvLen_au))deallocate(WvLen_au)
	if(allocated(WvLen_nm))deallocate(WvLen_nm)
	if(allocated(period_au))deallocate(period_au)
	if(allocated(period_fs))deallocate(period_fs)
	if(allocated(t0_au))deallocate(t0_au)
	if(allocated(t0_fs))deallocate(t0_fs)
	if(allocated(CEP_rad))deallocate(CEP_rad)
	if(allocated(FWHM_au))deallocate(FWHM_au)
	if(allocated(FWHM_fs))deallocate(FWHM_fs)
	if(allocated(Duration_au))deallocate(Duration_au)
	if(allocated(Duration_fs))deallocate(Duration_fs)
	if(allocated(swiss_knife))deallocate(swiss_knife)
	if(allocated(t_start))deallocate(t_start)
	if(allocated(t_end))deallocate(t_end)
	if(allocated(Envelope_fun))deallocate(Envelope_fun)

	allocate(LaserPulse(nn))
	allocate(II0_au(nn),   E0_au(nn),  freq_au(nn), omega_au(nn), WvLen_au(nn), period_au(nn), t0_au(nn), FWHM_au(nn), Duration_au(nn))
	allocate(II0_Wcm2(nn), E0_Vcm(nn), freq_Hz(nn), omega_eV(nn), WvLen_nm(nn), period_fs(nn), t0_fs(nn), FWHM_fs(nn), Duration_fs(nn))
	allocate(A0_au(nn), CEP_rad(nn), swiss_knife(nn), t_start(nn), t_end(nn), Envelope_fun(nn))
	CEP_rad(:) = 0.0d0
	FWHM_au(:) = 0.0d0
	FWHM_fs(:) = 0.0d0
	Duration_au(:) = 0.0d0
	Duration_fs(:) = 0.0d0

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_free_unit_nr(nr)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: nr
	LOGICAL(LGT) :: UNITOK, UNITOP
	INTEGER(I4B) :: i

	do i = 1,100000
		INQUIRE(unit=i, exist=UNITOK,opened=UNITOP)
		if(UNITOK .and. .not.UNITOP)then
			nr = i
			EXIT
		endif
	enddo

END SUBROUTINE
!===============================================================================
!*******************************************************************************
SUBROUTINE bcast_gvars()
	IMPLICIT NONE

	call MPI_bcast(only_initial_wf, 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(init_by_ITP    , 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(memory_save    , 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(do_diagH0      , 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(do_pulse       , 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(do_propag      , 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(restart        , 1, MPI_LOGICAL, root_id, MPI_COMM_WORLD, mpi_err)

	call MPI_bcast(N_fun , 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(N_tfe , 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(N_fe  , 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(N_gp  , 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(lmax  , 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(tagmax, 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(dx_min, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(dx_max, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(xmax  , 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(SphericalHarmonicRepresentation, 4, MPI_CHARACTER, root_id, MPI_COMM_WORLD, mpi_err)

	call MPI_bcast(abspot_alpha, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(abspot_frac,  1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)

	call MPI_bcast(Npulse, 1, MPI_INTEGER, root_id, MPI_COMM_WORLD, mpi_err)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE bcast_local_vars()
	IMPLICIT NONE

	if(my_id == root_id)then
		NN_fun= N_fun
		NN_fe = N_fe_arr(my_id)
		NN_gp = N_gp_arr(my_id)
		my_fe(:) = Idx_fe(:,my_id)
		my_gp(:) = Idx_gp(:,my_id)
		my_tag(:)= Idx_tag(:, my_id)
		do mpi_id=1,mpi_size-1
			call MPI_Send(N_fun,             1, MPI_INTEGER, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(N_fe_arr(mpi_id),  1, MPI_INTEGER, mpi_id, 1, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(N_gp_arr(mpi_id),  1, MPI_INTEGER, mpi_id, 2, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(N_tag_arr(mpi_id), 1, MPI_INTEGER, mpi_id, 3, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(Idx_fe(:,mpi_id),  2, MPI_INTEGER, mpi_id, 4, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(Idx_gp(:,mpi_id),  2, MPI_INTEGER, mpi_id, 5, MPI_COMM_WORLD, mpi_err)
			call MPI_Send(Idx_tag(:,mpi_id), 2, MPI_INTEGER, mpi_id, 6, MPI_COMM_WORLD, mpi_err)
		enddo
	else
		call MPI_Recv(NN_fun, 1, MPI_INTEGER, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(NN_fe,  1, MPI_INTEGER, root_id, 1, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(NN_gp,  1, MPI_INTEGER, root_id, 2, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(NN_tag, 1, MPI_INTEGER, root_id, 3, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_fe,  2, MPI_INTEGER, root_id, 4, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_gp,  2, MPI_INTEGER, root_id, 5, MPI_COMM_WORLD, rstatus, mpi_err)
		call MPI_Recv(my_tag, 2, MPI_INTEGER, root_id, 6, MPI_COMM_WORLD, rstatus, mpi_err)
	endif

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE bcast_pulse_params() ! only those which are used during the propagation: E0_au, t0_au, omega_au, swiss_knife, CEP_rad
	IMPLICIT NONE
	INTEGER(I4B) :: ii

	do ii = 1,Npulse
		call MPI_bcast(LaserPulse(ii)%Envelope, 10, MPI_CHARACTER,        root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%EA,        1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%t0,        1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%Omega,     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%param,     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%CEP,       1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%t_start,   1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
		call MPI_bcast(LaserPulse(ii)%t_end,     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	enddo

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE bcast_propag_params()	! maybe remove some if needed only by root_id
	IMPLICIT NONE

	call MPI_bcast(Gauge,            10, MPI_CHARACTER,        root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Epsilon_ortho,     1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Epsilon_conv,      1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Time_max_au,       1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Time_max_fs,       1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Time_step_wf_au,   1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Time_step_wf_fs,   1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Time_step_data_au, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Time_step_data_fs, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(dTime,             1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(N_step_wf,         1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(N_step_data,       1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(GS_all_qvecs,      1, MPI_LOGICAL,          root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(GS_refinement,     1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(N_Krylov,          1, MPI_INTEGER,          root_id, MPI_COMM_WORLD, mpi_err)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE bcast_ITP_params()
	IMPLICIT NONE

	call MPI_bcast(Time_max_ITP_au,  1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(dTime_ITP,        1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(Epsilon_conv_ITP, 1, MPI_DOUBLE_PRECISION, root_id, MPI_COMM_WORLD, mpi_err)

END SUBROUTINE
!-------------------------------------------------------------------------------
SUBROUTINE bcast_init_state_params()
	IMPLICIT NONE

	call MPI_bcast(target_atom, 10, MPI_CHARACTER,          root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(n_init,       1, MPI_INTEGER,            root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(l_init,       1, MPI_INTEGER,            root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(m_init,       1, MPI_INTEGER,            root_id, MPI_COMM_WORLD, mpi_err)
!	call MPI_bcast(tag_init,     1, MPI_INTEGER,            root_id, MPI_COMM_WORLD, mpi_err)	! not assigned when this routine is called 
	call MPI_bcast(Z_core,       1, MPI_DOUBLE_PRECISION,   root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(RR_core,      1, MPI_DOUBLE_PRECISION,   root_id, MPI_COMM_WORLD, mpi_err)
	call MPI_bcast(ascreen,      1, MPI_DOUBLE_PRECISION,   root_id, MPI_COMM_WORLD, mpi_err)

END SUBROUTINE
!-------------------------------------------------------------------------------
END MODULE







































