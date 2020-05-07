MODULE Operator_expect
	USE nrtype
	USE gvars
	USE MPI_pars
	USE FEDVR_module
	USE TDCC_module
	USE TDCC_FEDVR
	USE Operator_typedef
	USE Operator_math


	INTERFACE OP_expect
		MODULE PROCEDURE fedvr_op1_full_expect
		MODULE PROCEDURE fedvr_op2_full_expect
		MODULE PROCEDURE fedvr_op1_diag_expect
		MODULE PROCEDURE fedvr_op2_diag_expect
	END INTERFACE

	CONTAINS
!-------------------------------------------------------------------------------
FUNCTION fedvr_op1_full_expect(OP, wf) RESULT(expect)
	IMPLICIT NONE
	TYPE(fedvr_op_full), INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	COMPLEX(DPC) :: expect

	expect = wf * (OP * wf)
	expect = expect / wf%norm

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION fedvr_op2_full_expect(OP, wf) RESULT(expect)
	IMPLICIT NONE
	TYPE(fedvr_op_full), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	COMPLEX(DPC) :: expect

	expect = wf * (OP * wf)
	expect = expect / wf%norm

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION fedvr_op1_diag_expect(OP, wf) RESULT(expect)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	COMPLEX(DPC) :: expect

	expect = wf * (OP * wf)
	expect = expect / wf%norm

END FUNCTION
!-------------------------------------------------------------------------------
FUNCTION fedvr_op2_diag_expect(OP, wf) RESULT(expect)
	IMPLICIT NONE
	TYPE(fedvr_op_diag), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: OP
	TYPE(tdcc_fedvr_wf), INTENT(IN) :: wf
	COMPLEX(DPC) :: expect

	expect = wf * (OP * wf)
	expect = expect / wf%norm

END FUNCTION
!-------------------------------------------------------------------------------

END MODULE Operator_expect










