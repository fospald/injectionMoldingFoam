{
	// estimate initial timestep based on injection velocities and
	// initialize alpha1 field next to boundary
	scalar dtMin = HUGE;

	forAll(alpha1.boundaryField(), patchI)
	{
		Foam::fvPatchField<scalar>& pf = alpha1.boundaryField()[patchI];

			if (!isA< calculatedFvPatchField<scalar> >(injectionPointMask.boundaryField()[patchI])) {
				// skip processor boundaries, etc.
				// std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
				continue;
			}

			tmp< Field<scalar> > p = mag(injectionVelocities.boundaryField()[patchI]);
			scalar vMax = SMALL + max(p());
/*
			if (initAlpha) {
				Field<scalar>& intFld = const_cast<Field<scalar>&>(alpha1.internalField());
				UIndirectList<scalar>(intFld, pf.patch().faceCells()) =
					injectionPointMask.boundaryField()[patchI]*p()/vMax;
			}
*/
			dtMin = min(dtMin, min(1.0/(pf.patch().deltaCoeffs()*vMax)));

/*
		std::cout << "patch " << patchI
			<< " " << alpha1.boundaryField()[patchI].patch().name() 
			<< " " << U.boundaryField()[patchI].patch().name() 
			<< " " << p_rgh.boundaryField()[patchI].patch().name() 
			<< " " << T.boundaryField()[patchI].patch().name()
			<< " " << injectionPointMask.boundaryField()[patchI].patch().name()
			<< std::endl;
*/
	}

	reduce(dtMin, minOp<scalar>());
	dtMin *= maxCo;

	Info<< "dtMin = " << dtMin << endl;
	runTime.setDeltaT(dtMin);
}

