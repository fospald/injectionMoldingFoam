/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "SVD.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "processorFvPatchField.H"
#include "codedMixedFvPatchField.H"
#include "eigenCalc.H"
#include "fastExactClosure.H"

#if 1
tmp<fvVectorMatrix> symLaplace(tmp<volScalarField> nu, volVectorField& U, scalar lambdaReg)
{
	tmp<volTensorField> tmp_gradU(fvc::grad(U));
	volTensorField& gradU = tmp_gradU();
	tmp<volTensorField> tmp_gradUT(Foam::T(gradU));
	volTensorField& gradUT = tmp_gradUT();
	dimensionedScalar lambda("lambda", dimensionSet(0, 0, -2, 0, 0, 0, 0), lambdaReg);

//    volTensorField At_inv2_field("At_inv2", Foam::inv((gradU & gradUT) + lambda*Foam::I));
//	At_inv2_field.write();

	volTensorField At_inv_field("At_inv", Foam::inv((gradU & gradUT) + lambda*Foam::I));

/*
	scalarRectangularMatrix Ar(3, 3);
	scalar eps = 0.01;

	forAll (gradU, i)
	{
		tensor At = gradU[i]&gradUT[i];
		tensor& At_inv = At_inv_field[i];

		for (char j = 0; j < 3; j++) {
		for (char k = 0; k < 3; k++) {
			Ar[j][k] = At.v_[j*3+k];
		}}
		
		SVD svdA(Ar, eps);
		const scalarRectangularMatrix& Ar_inv = svdA.VSinvUt();

		for (char j = 0; j < 3; j++) {
		for (char k = 0; k < 3; k++) {
			At_inv.v_[j+3*k] = Ar_inv[j][k];
		}}
	}

		forAll(At_inv_field.boundaryField(), patchI)
		{
			if (!isA< calculatedFvPatchField<tensor> >(At_inv_field.boundaryField()[patchI])) {
				// skip processor boundaries, etc.
				// std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
				continue;
			}

                        Foam::fvPatchField<tensor>& bf = At_inv_field.boundaryField()[patchI];
                        bf == bf.patchInternalField();
                }

*/


	tmp<volTensorField> G = nu()*(Foam::I + (gradUT & gradUT & At_inv_field));

	tmp<volTensorField> diff = nu()*(gradU + gradUT) - (G() & gradU);

	dimensionedScalar small("small", dimensionSet(1, -1, -2, 0, 0, 0, 0), SMALL);
	tmp<volScalarField> diffMag = mag(diff())/(small + mag(G() & gradU));

	volScalarField diffw("saveDiff", diffMag());
	diffw.write();
	

//	Info << "weightedAverage(diffMag) " << gSum(U.mesh().V()) << endl;
	Info << "weightedAverage(diffMag) " << gSum((diffMag()*U.mesh().V())())/gSum(U.mesh().V()) << endl;

/*
	At_inv_field.write();


	dimensionedScalar q("q", dimensionSet(0, -1, 1, 0, 0, 0, 0), 1.0);
	tmp<volScalarField> scale = q*mag(U)*inv(SMALL + gMax(mag(U)()));

//	tmp<volTensorField> G = nu()*(Foam::I + (gradUT & gradUT & Foam::inv((gradU & gradUT) + lambda*Foam::I)));


   Info << "gMax(gradU) " << gMax(((nu()*gradU) && (nu()*gradU))()) << " diff " << gMax((diff() && diff())()) << endl;

	
    volTensorField diffw("saveDiff", diff());
	diffw.write();

    volTensorField diffwa("nugradu", nu()*(gradU + gradUT));
	diffwa.write();

    volTensorField diffwG("G", G());
	diffwG.write();

    volTensorField diffwF("GgradU", G() & gradU);
	diffwF.write();

    volScalarField diffwA("detGradu", det(gradU));
	diffwA.write();
*/
	return (
	      fvm::laplacian(G(), U) +
	      fvc::div(diff)
	);

	return (
              fvm::laplacian(nu, U) +
              fvc::div(nu*dev(T(gradU)))
	);
}
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    (void) maxDeltaT; // Prevent unused variable warnings

//    #include "createPrghCorrTypes.H"
//    #include "correctPhi.H"
//    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

	// Correct value of alpha2
	alpha1.boundaryField().updateCoeffs();
	alpha2 = scalar(1) - alpha1;

	// init temperature field to air temperature
	for(label i=0; i<T.size(); i++)
	{
		T[i] = airTemperature.value();
//		p_rgh[i] = pref.value();
//		p[i] = pref.value();
	}


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    runTime.write();

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info << "\n----------------------------------------------------" << endl;
        Info << "Time Step " << runTime.timeIndex() << endl;
        Info<< "Time = " << runTime.timeName() << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            (void)MULESCorr; (void)alphaApplyPrevCorr; // Prevent unused variable warnings

            if (pimple.firstIter())// || alphaOuterCorrectors)
            {
                twoPhaseProperties.correct();

                #include "alphaEqnSubCycle.H"
                interface.correct();
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

	#include "TEqn.H"

	if (FiberSuspension) {        
		// solve A2 equation
		#include "fiberClosure.H"
		#include "A2Eqn.H"

		if (FiberConcentration) {
			#include "PhiEqn.H"
		}
	}

	// Check if we need to write results

	if ((nextWriteIndex > lastWriteIndex) || runTime.end()) {
		Info << "Writing data for t = " << runTime.value() << endl;
		runTime.writeNow();
		lastWriteIndex = nextWriteIndex;
	}
	else {
		runTime.write();
	}

	if (outputStats) {
		#include "stats.H"
	}

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

