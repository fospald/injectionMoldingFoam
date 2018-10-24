{

// see:
// http://www.cfd-online.com/Forums/openfoam-solving/69103-diverging-result-temperature-field-interfoam.html

surfaceScalarField alpha1f = min(max(fvc::interpolate(alpha1_zg), scalar(0)), scalar(1));
surfaceScalarField kappaf = k1*alpha1f + k2*(1.0 - alpha1f);
volScalarField rhoCp = rho1*Cp1*alpha1_zg + rho2*Cp2*alpha2_zg;
//surfaceScalarField rhoCpf = rho1*Cp1*alpha1f + rho2*Cp2*(1.0 - alpha1f);
//volVectorField grad_alpha = fvc::grad(alpha1_zg);
//dimensionedScalar small("small", pow(dimLength, -2), SMALL);
//volVectorField U_proj = U - (1 - 4.0*alpha1_zg*alpha2_zg)*((U & grad_alpha) * grad_alpha)/((grad_alpha & grad_alpha) + small);

    fvScalarMatrix TEqn
    (
	fvm::ddt(rhoCp, T)

//     	+ fvm::div(rhoCpPhi, T, "div(rho*phi,T)")

	// rewriting:
	// div(rhoCp*U*T) = (grad(rhoCp) & U)*T + rhoCp*div(U*T)
	// div(U*T) = div(U)*T + U*grad(T)
	// if there should be no energy transfer accross the interface then (grad(rhoCp) & n) should be zero
	// i.e. we get (grad(rhoCp) & U)*T = (grad(rhoCp) & (U - (U&n)*n))*T
	// since rho1*Cp1 is currently constant and we do not any heat generation at the interface, we set this term zero
//	+ fvm::Sp(fvc::grad(rhoCp) & U_proj, T)

	+ rhoCp*fvm::div(phi, T)

	- fvm::laplacian(kappaf, T)

// TODO: add energy dissipation, currenty too much heat is generated at the boundary?
//	- pow(alpha1_zg, 10.0)*mu*pow(magnitudeD, 2.0)
	- dissipationScale*mu*magnitudeD*magnitudeD
    );

    TEqn.relax();
    TEqn.solve();

    twoPhaseProperties.correct();
}
