{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phic.boundaryField()[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> tphiAlpha;

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(phic*interface.nHatf());

        tmp<surfaceScalarField> tphiAlphaUn
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );



        surfaceScalarField tphiAlphaUn_NoOut = tphiAlphaUn();

#if 1
                // zero outflow for phase alpha1 (NEW)
                // however this may cause overshooting the value above 1.0,
                // especially in corners
                // we therefore limit alpha1 and report the overshoots as integral value

		forAll(injectionPointMask.boundaryField(), patchI)
		{
			if (!isA< calculatedFvPatchField<scalar> >(injectionPointMask.boundaryField()[patchI])) {
				// skip processor boundaries, etc.
				// std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
				continue;
			}

                        Foam::fvsPatchField<scalar>& bf = tphiAlphaUn_NoOut.boundaryField()[patchI];
                        bf == min(0.0, injectionPointMask.boundaryField()[patchI]*bf);
                }
#endif



        tphiAlpha = tphiAlphaUn;
        MULES::explicitSolve(alpha1, phi, tphiAlphaUn_NoOut, 1, 0);

        alpha2 = 1.0 - alpha1;

        interface.correct();
    }

    rhoPhi = tphiAlpha()*(rho1 - rho2) + phi*rho2;
//    rhoCpPhi = tphiAlpha()*(rho1*Cp1 - rho2*Cp2) + phi*rho2*Cp2;

#if 1
        // integrate the overshoots and limit alpha1
        Foam::tmp< Field<scalar> > alpha1minusOne = alpha1.internalField() - 1.0;
        Foam::tmp< Field<scalar> > alpha1Overs = pos(alpha1minusOne())*alpha1minusOne();
        Foam::tmp< Field<scalar> > alpha1Unders = -neg(alpha1.internalField())*alpha1.internalField();
        alpha1Overshoots += Foam::gSum(alpha1Overs()*mesh.V());
        alpha1Undershoots += Foam::gSum(alpha1Unders()*mesh.V());

        alpha1.internalField() -= alpha1Overs();
        alpha1.internalField() += alpha1Unders();
        alpha2 = 1.0 - alpha1;
#endif

	// compute zero gradient alpha
	alpha1_zg.internalField() = alpha1.internalField();
	forAll(injectionPointMask.boundaryField(), patchI)
	{
		if (!isA< calculatedFvPatchField<scalar> >(injectionPointMask.boundaryField()[patchI])) {
			// skip processor boundaries, etc.
			// std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
			continue;
		}

		alpha1_zg.boundaryField()[patchI] = alpha1_zg.boundaryField()[patchI].patchInternalField();
	}
	alpha2_zg = 1.0 - alpha1_zg;

        scalar alpha1VolFrac = alpha1.weightedAverage(mesh.Vsc()).value();
        scalar alpha1OvershootFrac = 100.0 * alpha1Overshoots / meshV;
        scalar alpha1UndershootFrac = 100.0 * alpha1Undershoots / meshV;
	scalar alpha1Min = gMin(alpha1);
	scalar alpha1MinBoundary = 1.0;
	scalar alpha1Max = gMax(alpha1);

	// calculate singularity indicator alpha1MinBoundary
	forAll(injectionPointMask.boundaryField(), patchI)
	{
		if (!isA< calculatedFvPatchField<scalar> >(injectionPointMask.boundaryField()[patchI])) {
			// skip processor boundaries, etc.
			// std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
			continue;
		}

                Foam::fvPatchField<scalar>& pf = alpha1.boundaryField()[patchI];
		alpha1MinBoundary = min(alpha1MinBoundary, min(pf.patchInternalField()));
	}

	reduce(alpha1MinBoundary, minOp<scalar>());

	// determine next write index
	scalar alpha1VolPercent = 100.0 * alpha1VolFrac;
	nextWriteIndex = static_cast<int>(alpha1VolPercent / fillPercentWriteInterval);

	// check if we need to write now
	scalarList newFillPercentWritePoints;
	for (Foam::label i=0; i<fillPercentWritePoints.size(); i++) {
		if (alpha1VolPercent >= fillPercentWritePoints[i]) {
			lastWriteIndex = -1;
			continue;
		}
		newFillPercentWritePoints.append(fillPercentWritePoints[i]);
	}
	fillPercentWritePoints = newFillPercentWritePoints;

	bool endRun = false;
        if (alpha1VolFrac >= fillPercentBreak || ((1.0 - alpha1VolFrac*100.0/fillPercentBreak) < fillPercentTol)) {
    		Info << "Ending fill because exceeding fill volume" << endl;
		endRun = true;
	}
	else if (alpha1MinBoundary >= fillSingularityBreak) {
    		Info << "Ending fill because exceeding singularity threshold" << endl;
		endRun = true;
        }
	else if (alpha1OvershootFrac >= fillPercentOvershootBreak) {
    		Info << "Ending fill because exceeding overshoot threshold" << endl;
		endRun = true;
        }

	if (endRun) {
		runTime.setEndTime(runTime.value());
		nextWriteIndex = lastWriteIndex + 1;
	}

    Info<< "Liquid phase volume fraction = " << (alpha1VolFrac*100.0) << "%"
        << "  Min(alpha1) = " << alpha1Min
        << "  Max(alpha1) = " << alpha1Max
        << "  singularity = " << (alpha1MinBoundary*100.0) << "%"
        << "  overshoots = " << (alpha1OvershootFrac) << "%"
        << "  undershoots = " << (alpha1UndershootFrac) << "%"
        << endl;
}
