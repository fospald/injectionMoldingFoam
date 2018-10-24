{
// TODO: It would be better to use MULES for the concentration flux...

/*
	dimensionSet:
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles,
    const scalar current,
    const scalar luminousIntensity
*/

// calculate boundary field values

//c.boundaryField() = (-FiberKd) / (2*FiberAe*FiberA2*magnitudeD.boundaryField());

//(-Kd + Ae^2 gamma Kc phi - Ae^2 d gammab Keta sngradlnetar + d Sqrt[(4 Ae^2 gammab Kc Kd phi)/ d^2 + (Kd/d - (Ae^2 gamma Kc phi)/d + Ae^2 gammab Keta sngradlnetar)^2])/(2 Ae^2 gammab Kc)



    dimensionedScalar small("small", dimensionSet(0, 0, -1, 0, 0, 0, 0), SMALL);
    volScalarField etar = 0.5*mu1c*(magnitudeD + small);	// TODO: viscous energy does not consider anisotropy
    dimensionedScalar dim_etar("dim_etar", etar.dimensions(), 1.0); // dimensionSet(1, -1, -3, 0, 0, 0, 0), 1.0);

    fvScalarMatrix phiFiberEqn
    (
	fvm::ddt(phiFiber)
	+ fvm::div(phi, phiFiber)
	- fvm::laplacian(alpha1*FiberKd*magnitudeD, phiFiber)
	+ fvc::div(((FiberKc*FiberAe*FiberAe)*fvc::grad(magnitudeD*phiFiber) + (FiberKeta*FiberAe*FiberAe)*fvc::grad(log(SMALL+etar/dim_etar)))*phiFiber*alpha1)
    );

    phiFiberEqn.relax();
    phiFiberEqn.solve();

// fix negative values (TODO: possible mass conservation issue)
        // integrate the overshoots and limit phiFiber

        scalar phiFiberVolFrac = phiFiber.weightedAverage(mesh.Vsc()).value();
        Foam::tmp< Field<scalar> > phiFiberminusOne = phiFiber.internalField() - 1.0;
        Foam::tmp< Field<scalar> > phiFiberOvers = pos(phiFiberminusOne())*phiFiberminusOne();
        Foam::tmp< Field<scalar> > phiFiberUnders = -neg(phiFiber.internalField())*phiFiber.internalField();
        phiFiberOvershoots += Foam::gSum(phiFiberOvers()*mesh.V());
        phiFiberUndershoots += Foam::gSum(phiFiberUnders()*mesh.V());

        phiFiber.internalField() -= phiFiberOvers();
        phiFiber.internalField() += phiFiberUnders();

        scalar phiFiberOvershootFrac = 100.0 * phiFiberOvershoots / meshV;
        scalar phiFiberUndershootFrac = 100.0 * phiFiberUndershoots / meshV;

    Info<< "Concentration volume fraction = " << (phiFiberVolFrac*100.0) << "%"
 //       << "  Min(phiFiber) = " << phiFiberMin
 //       << "  Max(phiFiber) = " << phiFiberMax
        << "  overshoots = " << (phiFiberOvershootFrac) << "%"
        << "  undershoots = " << (phiFiberUndershootFrac) << "%"
        << endl;

}
