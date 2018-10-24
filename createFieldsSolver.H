
// define some variables

    scalar alpha1Overshoots = 0.0;
    scalar alpha1Undershoots = 0.0;
    scalar meshV = Foam::gSum(mesh.V());

    int lastWriteIndex = 0, nextWriteIndex = 0;

// read solver properties

    Info<< "Reading solverProperties\n" << endl;

    IOdictionary solverProperties
    (
        IOobject
        (
            "solverProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalarList fillPercentWritePoints
    (
        solverProperties.lookup("fillPercentWritePoints")
    );

    bool outputStats(
        solverProperties.lookupOrDefault<bool>("outputStats", true)
    );

    scalar fillPercentBreak(
        solverProperties.lookupOrDefault<scalar>("fillPercentBreak", 100.0)
    );

    scalar fillPercentTol(
	solverProperties.lookupOrDefault<scalar>("fillPercentTol", 1e-4)
    );

    scalar fillPercentWriteInterval(
	solverProperties.lookupOrDefault<scalar>("fillPercentWriteInterval", 10.0)
    );

    scalar fillSingularityBreak(
        solverProperties.lookupOrDefault<scalar>("fillSingularityBreak", 0.95)
    );

    scalar fillPercentOvershootBreak(
        solverProperties.lookupOrDefault<scalar>("fillPercentOvershootBreak", 5.0)
    );

    word UEqnForm(
        solverProperties.lookupOrDefault<word>("UEqnForm", "standard")
    );

    scalar UEqnStabParam(
        solverProperties.lookupOrDefault<scalar>("UEqnStabParam", 1.0e-3)
    );

    scalar UEqnNpCoupling(
        solverProperties.lookupOrDefault<scalar>("UEqnNpCoupling", 1.0)
    );

    scalar UEqnNpStabParam1(
        solverProperties.lookupOrDefault<scalar>("UEqnNpStabParam1", 1.0)
    );

    scalar UEqnNpStabParam2(
        solverProperties.lookupOrDefault<scalar>("UEqnNpStabParam2", 3.0)
    );

    scalar alphaCutoffU(
        solverProperties.lookupOrDefault<scalar>("alphaCutoffU", 0.001)
    );


// read the viscosity model

    word viscosityModel(
        solverProperties.lookupOrDefault<word>("viscosityModel", "constant")
    );

    Info<< "Viscosity Model: " << viscosityModel << endl;

    dimensionedScalar viscosityA0("A0", dimensionSet(0, 0, 1, 0, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("viscosityA0", 0.0)
    );

    scalar viscosityUpperLimitFactor(
	solverProperties.lookupOrDefault<scalar>("viscosityUpperLimitFactor", 10.0)
    );

    scalar viscosityLowerLimitFactor(
	solverProperties.lookupOrDefault<scalar>("viscosityLowerLimitFactor", 0.1)
    );

    scalar viscosityD1(
	solverProperties.lookupOrDefault<scalar>("viscosityD1", 0.0)
    );

    scalar viscosityA1(
	solverProperties.lookupOrDefault<scalar>("viscosityA1", 0.0)
    );

    scalar dissipationScale(
	solverProperties.lookupOrDefault<scalar>("dissipationScale", 1.0)
    );

    scalar viscosityN(
	solverProperties.lookupOrDefault<scalar>("viscosityN", 0.0)
    );

    dimensionedScalar viscosityA2("A2", dimensionSet(0, 0, 0, -1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("viscosityA2", 0.0)
    );

    dimensionedScalar viscosityA3("A3", dimensionSet(0, 0, 0, 1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("viscosityA3", INFINITY)
    );

    dimensionedScalar viscosityD3("D3", dimensionSet(-1, 1, 2, 1, 0, 0, 0), // K/Pa
	solverProperties.lookupOrDefault<scalar>("viscosityD3", 0.0)
    );

    dimensionedScalar viscosityTauStar("TauStar", dimensionSet(1, -1, -2, 0, 0, 0, 0), // Pa
	solverProperties.lookupOrDefault<scalar>("viscosityTauStar", 1.0)
    );

    dimensionedScalar viscosityTref("Tref", dimensionSet(0, 0, 0, 1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("viscosityTref", 0.0)
    );

// read thermal properties

    dimensionedScalar wallTemperature("Twall", dimensionSet(0, 0, 0, 1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("wallTemperature", 0.0)
    );

    dimensionedScalar airTemperature("Tair", dimensionSet(0, 0, 0, 1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("airTemperature", 0.0)
    );

    dimensionedScalar k1("k1", dimensionSet(1, 1, -3, -1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("heatConductivity1", 0.0)
    );

    dimensionedScalar k2("k2", dimensionSet(1, 1, -3, -1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("heatConductivity2", 0.0)
    );

    dimensionedScalar Cp1("Cp1", dimensionSet(0, 2, -2, -1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("specificHeat1", 0.0)
    );

    dimensionedScalar Cp2("Cp2", dimensionSet(0, 2, -2, -1, 0, 0, 0),
	solverProperties.lookupOrDefault<scalar>("specificHeat2", 0.0)
    );

    Cp2 = rho1*Cp1/rho2;

/*
    scalar fillAlphaBoundary(
	solverProperties.lookupOrDefault<scalar>("fillAlphaBoundary", 0.5)
    );

    scalar fillAlphaBoundaryContact(
	solverProperties.lookupOrDefault<scalar>("fillAlphaBoundaryContact", 0.90)
    );

    Switch compressible(
	solverProperties.lookupOrDefault<Switch>("compressible", true)
    );

    Switch stabUEqn(
        solverProperties.lookup("stabUEqn")
    );

    scalar stabUEqnParam(
	solverProperties.lookupOrDefault<scalar>("stabUEqnParam", 1.5)
    );

    Switch useDevClosure(
        solverProperties.lookup("useDevClosure")
    );

    scalar coolTemperatureBreak(
	readScalar(solverProperties.lookup("coolTemperatureBreak"))
    );

    scalar deltaTcool(
        readScalar(solverProperties.lookup("deltaTcool"))
    );

    scalar alphaWallExp(
        readScalar(solverProperties.lookup("alphaWallExp"))
    );

    scalar sourceA2Exp(
        readScalar(solverProperties.lookup("sourceA2Exp"))
    );

    Switch harmonicMu(
        solverProperties.lookup("harmonicMu")
    );
    
    Switch implicitAlpha(
        solverProperties.lookup("implicitAlpha")
    );

    Switch initAlpha(
        solverProperties.lookup("initAlpha")
    );

    
    Switch ficticiousDom(
        solverProperties.lookup("ficticiousDom")
    );

    Switch useFluxBcP(
        solverProperties.lookup("useFluxBcP")
    );

    Switch useCurlStab(
        solverProperties.lookup("useCurlStab")
    );

    Switch useFluxBcU(
        solverProperties.lookup("useFluxBcU")
    );
*/

