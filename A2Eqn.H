
if (FiberClosure == "free")
{
	// solve BEqn
	fvSymmTensorMatrix B2Eqn(
		fvm::ddt(B2) + fvm::div(phi, B2, "div(phi,B2)")
	);
	solve(B2Eqn == SourceB2*alpha1);

	// reconstruct A
	#include "A2fromB2.H"
}
else
{
	// solve A2 equation
	fvSymmTensorMatrix A2Eqn(
		fvm::ddt(A2) + fvm::div(phi, A2, "div(phi,A2)")
	);
	solve(A2Eqn == SourceA2*alpha1);

	if (FiberClosure == "FEC") {
		// solve B2 equation
		fvSymmTensorMatrix B2Eqn(
			fvm::ddt(B2) + fvm::div(phi, B2, "div(phi,A2)")
		);
		solve(B2Eqn == SourceB2*alpha1);
	}

	#include "A2Eig.H"
}


/*
if (FiberClosure == "IBOF") {
	Sp_imp1 = 6.0*Ci*(pow(0.5*(D && D.T()), 0.5));
	A2Eqn += fvm::Sp(Sp_imp1, A2) + fvm::Sp(ClosureA4_imp, A2);
}
*/

