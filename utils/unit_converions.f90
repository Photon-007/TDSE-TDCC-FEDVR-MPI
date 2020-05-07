MODULE unit_conversions
	USE nrtype
	IMPLICIT NONE


	REAL(DP), PARAMETER :: constSI_mass_e = 9.1093837015d-31	! kg
	REAL(DP), PARAMETER :: constSI_charge = 1.602176634d-19		! C
	REAL(DP), PARAMETER :: constSI_hbar   = 1.054571817d-34		! Js
	REAL(DP), PARAMETER :: constSI_permittivity = 1.11265e-10	! 1/Fm

	REAL(DP), PARAMETER :: constAU_alpha = 1.0d0 / 137.03599911d0
	REAL(DP), PARAMETER :: constAU_lightSpeed = 1.0d0 / constAU_alpha
	REAL(DP), PARAMETER :: constSI_lightSpeed = 299792458		! m/s

	REAL(DP), PARAMETER :: constSI_length    = constSI_permittivity * constSI_hbar**2 / (constSI_mass_e * constSI_charge**2)	! m
	REAL(DP), PARAMETER :: constSI_energy    = constSI_hbar**2 / (constSI_mass_e * constSI_length**2)							! J
	REAL(DP), PARAMETER :: constSI_time      = constSI_hbar / constSI_energy													! s
	REAL(DP), PARAMETER :: constSI_velocity  = constSI_length / constSI_time													! m/s
	REAL(DP), PARAMETER :: constSI_frequency = constSI_velocity / ( constSI_length)												! Hz
	REAL(DP), PARAMETER :: constSI_electricfield = constSI_charge / (constSI_permittivity * constSI_length**2)					! V/m
!	REAL(DP), PARAMETER :: constSI_intensity = constSI_energy**2 / (constSI_hbar * constSI_length**2)							! W/m**2
	REAL(DP), PARAMETER :: constSI_intensity = constSI_energy / constSI_time / constSI_length**2
!	REAL(DP), PARAMETER :: constSI_intensity = constSI_hbar / constSI_time**2 / constSI_length**2


	REAL(DP), PARAMETER :: convert_au2eV = 27.2113878129
	REAL(DP), PARAMETER :: convert_J2eV  = 1.60217656535d-19
	REAL(DP), PARAMETER :: convert_eV2cm1= convert_J2eV / (TWOPI * constSI_hbar * constSI_lightSpeed)
	REAL(DP), PARAMETER :: convert_I2E   = constSI_permittivity / (4.0d0 * PI) * constSI_lightSpeed / 2.0d0


	CHARACTER(LEN=2), PARAMETER :: UNIT_MASS           = 'kg'
	CHARACTER(LEN=1), PARAMETER :: UNIT_CHARGE         = 'C'
	CHARACTER(LEN=1), PARAMETER :: UNIT_LENGTH         = 'm'
	CHARACTER(LEN=2), PARAMETER :: UNIT_ANGSTROM       = 'AA'
	CHARACTER(LEN=2), PARAMETER :: UNIT_NANOMETER      = 'nm'
	CHARACTER(LEN=1), PARAMETER :: UNIT_ENERGY         = 'J'
	CHARACTER(LEN=1), PARAMETER :: UNIT_TIME           = 's'
	CHARACTER(LEN=2), PARAMETER :: UNIT_TIME_FS        = 'fs'
	CHARACTER(LEN=2), PARAMETER :: UNIT_TIME_AS        = 'as'
	CHARACTER(LEN=2), PARAMETER :: UNIT_FREQUENCY      = 'Hz'
	CHARACTER(LEN=3), PARAMETER :: UNIT_FIELD_STRENGTH = 'V/m' 
	CHARACTER(LEN=4), PARAMETER :: UNIT_FIELD_INTENSITY= 'W/m2'


	TYPE phys_variable
		REAL(DP)         :: value = 0.0d0
		CHARACTER(LEN=9) :: units = ''
	END TYPE

	CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE write_Atomic_units(unit_no)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: unit_no

	write(unit_no,*)"########################################################"	! 55
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| mass_e:        ", constSI_mass_e,       " [kg]            |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| charge_e:      ", constSI_charge,       " [C]             |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| hbar:          ", constSI_hbar,         " [Js]            |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| permittivity:  ", constSI_permittivity, " [F/m]           |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| length:        ", constSI_length,       " [m]             |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| energy:        ", constSI_energy,       " [J]             |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| time:          ", constSI_time,         " [s]             |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| velocity:      ", constSI_velocity,     " [m/s]           |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| frequency:     ", constSI_frequency,    " [Hz]            |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| electricfield: ", constSI_electricfield," [V/m]           |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| intensity:     ", constSI_intensity,    " [W/m^2]         |"
	write(unit_no,'(1x,A17,1p,d21.14,A18)') "| convert_I2E:   ", convert_I2E * constSI_electricfield**2, "                 |"
	write(unit_no,*)"########################################################"
	write(unit_no,*)""

	write(unit_no,*) TWOPI / (constSI_hbar * constSI_lightSpeed) * convert_J2eV
	write(unit_no,*) convert_J2eV / (TWOPI * constSI_hbar * constSI_lightSpeed) 
END SUBROUTINE
!-------------------------------------------------------------------------------
FUNCTION Atomic2SI(value, what)
	IMPLICIT NONE
	REAL(DP) :: Atomic2SI
	REAL(DP), INTENT(IN) :: value
	CHARACTER(LEN=*), INTENT(IN) :: what

	select case (what)
		case(UNIT_MASS)
			Atomic2SI = value * constSI_mass_e

		case(UNIT_CHARGE)
			Atomic2SI = value * constSI_charge

		case(UNIT_LENGTH)
			Atomic2SI = value * constSI_length 

		case(UNIT_ANGSTROM)
			Atomic2SI = value * constSI_length / 1d-10

		case(UNIT_NANOMETER)
			Atomic2SI =  value * constSI_length / 1d-9

		case(UNIT_ENERGY)
			Atomic2SI = value * constSI_energy

		case(UNIT_TIME)
			Atomic2SI = value * constSI_time

		case(UNIT_TIME_FS)
			Atomic2SI = value * constSI_time * 1.0d15

		case(UNIT_TIME_AS)
			Atomic2SI = value * constSI_time * 1.0d18

		case(UNIT_FREQUENCY)
			Atomic2SI = value * constSI_frequency

		case(UNIT_FIELD_STRENGTH)
			Atomic2SI = value * constSI_electricfield

		case(UNIT_FIELD_INTENSITY)
			Atomic2SI = value * constSI_intensity

	end select

END FUNCTION Atomic2SI
!-------------------------------------------------------------------------------
FUNCTION SI2Atomic(value, what)
	IMPLICIT NONE
	REAL(DP) :: SI2Atomic
	REAL(DP), INTENT(IN) :: value
	CHARACTER(LEN=*), INTENT(IN) :: what

	select case (what)
		case(UNIT_MASS)
			SI2Atomic = value / constSI_mass_e

		case(UNIT_CHARGE)
			SI2Atomic = value / constSI_charge

		case(UNIT_ANGSTROM)
			SI2Atomic = value / constSI_length * 1d-10

		case(UNIT_NANOMETER)
			SI2Atomic = value / constSI_length * 1d-9

		case(UNIT_ENERGY)
			SI2Atomic = value / constSI_energy

		case(UNIT_TIME)
			SI2Atomic = value / constSI_time

		case(UNIT_TIME_FS)
			SI2Atomic = value / constSI_time * 1.0d-15

		case(UNIT_TIME_AS)
			SI2Atomic = value / constSI_time * 1.0d-18

		case(UNIT_FREQUENCY)
			SI2Atomic = value / constSI_frequency

		case(UNIT_FIELD_STRENGTH)
			SI2Atomic = value / constSI_electricfield

		case(UNIT_FIELD_INTENSITY)
			SI2Atomic = value / constSI_intensity

	end select

END FUNCTION SI2Atomic
!-------------------------------------------------------------------------------
SUBROUTINE PROUST(T)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: T				!The time, in seconds. Positive only, please.
!	INTEGER(I4B),     DIMENSION(5), PARAMETER :: USIZE = (/7*24*60*60, 24*60*60, 60*60,   60,    1/)				!Size of the time unit.
!	CHARACTER(LEN=3), DIMENSION(5), PARAMETER :: UNAME = (/      "wk",      "d",  "hr","min","sec"/)				!Name of the time unit.
!	CHARACTER(LEN=28) :: TEXT	!A scratchpad.
!	INTEGER(I4B) I,L,N,S		!Assistants.

1	FORMAT (I0,1X,A)							!I0 means I only: variable-length, no leading spaces.
!	S = T										!A copy I can mess with.
!	L = 0										!No text has been generated.
!	DO I = 1,NTYPES								!Step through the types to do so.
!		N = S/USIZE(I)							!Largest first.
!		IF (N.GT.0) THEN						!Above the waterline?
!			S = S - N*USIZE(I)					!Yes! Remove its contribution.
!			IF (L.GT.0) THEN					!Is this the first text to be rolled?
!				L = L + 2						!No.
!				TEXT(L - 1:L) = ", "			!Cough forth some punctuation.
!			END IF								!Now ready for this count.
!			WRITE (TEXT(L + 1:),1) N,UNAME(I)	!Place, with the unit name.
!			L = LEN_TRIM(TEXT)					!Find the last non-blank resulting.
!		END IF									!Since I'm not keeping track.
!	END DO										!On to the next unit.
!	Cast forth the result.
!	WRITE (6,*) T,">",TEXT(1:L),"<"				!With annotation.
END SUBROUTINE									!Simple enough with integers.
!-------------------------------------------------------------------------------

END MODULE unit_conversions
