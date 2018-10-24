
//tmp<volScalarField> ex = -viscosityA1*(T - viscosityTref)/(viscosityA3 + viscosityD3*p + (T - viscosityTref));
//Info << "ex min=" << gMin(ex()) << " max = " << gMax(ex()) << endl;

// compute the pointwise viscosity depending on the viscosity law
if (viscosityModel == "cross-wlf") {
	// FIXME: the Temperature Jump accross the interface causes instabilities
	tmp<volScalarField> dT = Foam::max(T - viscosityTref, dimensionedScalar("0", viscosityTref.dimensions(), 0));
	tmp<volScalarField> mu1T = mu1*viscosityD1*Foam::exp(-viscosityA1*dT()/(viscosityA3 + viscosityD3*p + dT()));
	mu1c = Foam::min(Foam::max(mu1T()/(1 + Foam::pow(mu1T()*magnitudeD/viscosityTauStar, 1-viscosityN)), mu1*viscosityLowerLimitFactor), mu1*viscosityUpperLimitFactor);
}
else if (viscosityModel == "carreau-wlf") {
	// FIXME: the Temperature Jump accross the interface causes instabilities
	mu1c = mu1
		*Foam::exp(-viscosityA2*(T - viscosityTref))
		*Foam::pow(1 + Foam::pow(viscosityA0*magnitudeD, 2), 0.5*(viscosityA1 - 1.0));
}
else if (viscosityModel == "constant") {
	mu1c = mu1;
}
else {
	FatalErrorIn("UEqn.H")
		<< "unknown viscosity model: '" << viscosityModel << "'"
		<< abort(FatalError);
}

mu = alpha1*mu1c + alpha2*mu2;


// Info << "viscosityA1=" << viscosityA1 << "viscosityA3=" << viscosityA3 << "viscosityD1=" << viscosityD1 << " mu1=" << mu1 << " mu min = " << gMin(mu) << " max = " << gMax(mu) << " viscosityLimitFactor=" << viscosityLimitFactor << endl;


// create the matrix for the momentum predictor
fvVectorMatrix UEqn
(
	fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
);

// add the isotropic part of the viscous stress div(2*mu*D)
// there are several equivalent forms of the equation, with different implicity/explicit splitting
if (UEqnForm == "standard") {
	UEqn += (
	      -fvm::laplacian(mu, U)
	      -fvc::div(mu*Foam::T(fvc::grad(U)))
	);
}
else if (UEqnForm == "stabilized") {
	UEqn += (
	      -symLaplace(mu, U, UEqnStabParam)
	);
}
else if (UEqnForm == "curl") {
	UEqn += (
	      -fvm::laplacian(2.0*mu, U)
	      -fvc::div(mu*Foam::T(fvc::grad(U))-fvc::grad(U))
	);
}
else {
	FatalErrorIn("UEqn.H")
		<< "unknown form for UEqn: '" << UEqnForm << "'"
		<< abort(FatalError);
}

// add the anisotropic part of the viscous stress div(2*mu*Np*(A4:D))
if (FiberSuspension && (FiberNp != 0.0))
{
	if (FiberClosure == "free") {
		throw "ClosureA4 not implemented for free closure";
	}

	dimensionedScalar stab = (2.0*FiberNp)*UEqnNpStabParam1;

	UEqn += // stress term:
		- fvc::div((2.0*FiberNp)*mu*dev(ClosureA4), "div(2*mu*Np*dev(ClosureA4))")
		- fvm::laplacian(stab*mu, U, "laplacian(Np,U)")
		+ fvc::laplacian(stab*mu, U, "div(Np,U)");

#if 0
	dimensionedScalar stab = (2.0*FiberNp*mu1)*UEqnNpStabParam1;

	UEqn += // stress term:
		- fvc::div((2.0*FiberNp*UEqnNpCoupling)*alpha1*mu*dev(ClosureA4), "div(2*mu*Np*ClosureA4)")
		// stabilization terms:
#if 0
		- fvm::laplacian(stab*alpha1, U, "laplacian(Np,U)")
		+ fvc::laplacian(stab*pow(alpha1, UEqnNpStabParam2), U, "div(Np,U)");
#else
		- fvm::laplacian(stab*alpha1_zg, U, "laplacian(Np,U)")
		+ fvc::laplacian(stab*pow(alpha1_zg, UEqnNpStabParam2), U, "div(Np,U)");
#endif
//		+ fvc::div(stab*pow(alpha1, UEqnNpStabParam2)*dev(fvc::grad(U)), "div(Np,U)");
#endif
}

// relaxation is not appropriate in a transient simulation 
// UEqn.relax();

// possibly solve the momentum predictor
if (pimple.momentumPredictor())
{
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
}

