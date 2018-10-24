
Info<< "Reading injectionPoints\n" << endl;

IOdictionary injectionPoints
(
    IOobject
    (
        "injectionPoints",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

wordList boundaryTypes;
forAll(alpha1.boundaryField(), f) {
    boundaryTypes.append("calculated");
}

volScalarField injectionPointMask
(
    IOobject
    (
        "injectionPointMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0,
    boundaryTypes
);


volScalarField injectionTemperatures
(
    IOobject
    (
        "injectionTemperatures",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    wallTemperature.value(),
    boundaryTypes
);

volVectorField injectionVelocities
(
    IOobject
    (
        "injectionVelocities",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    pTraits<vector>::zero,
    boundaryTypes
);

const symmTensor A2iso = symmTensor(1.0/3.0, 0, 0, 1.0/3.0, 0, 1.0/3.0);
const symmTensor B2iso = symmTensor(1.0, 0, 0, 1.0, 0, 1.0);

volSymmTensorField injectionA2
(
    IOobject
    (
        "injectionA2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    A2iso,
    boundaryTypes
);

volSymmTensorField injectionB2
(
    IOobject
    (
        "injectionB2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    B2iso,
    boundaryTypes
);


scalar sumFlowRates = 0.0;

forAll(injectionPoints.toc(), key)
{
    word name = injectionPoints.toc()[key];

    Info << "Reading injection point: " << name << endl;

    dictionary& ip = injectionPoints.subDict(name);
    word type = ip.lookup("type");
    scalar flowRate(readScalar(ip.lookup("flowRate")));
    scalar profileExp(readScalar(ip.lookup("profileExp")));
    scalar profileExpA2(readScalar(ip.lookup("profileExpA2")));
    scalar Temp(readScalar(ip.lookup("T")));

    sumFlowRates += flowRate;

    if (type == "geodesic")
    {
        // read dimensions
        vector center(ip.lookup("center"));
        scalar radius(readScalar(ip.lookup("radius")));

	volScalarField d
	(
	    IOobject
	    (
		name,
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

        // mark the injection regions in each boundary field
        forAll(injectionPointMask.boundaryField(), f)
        {
            Foam::fvPatchField<scalar>& bf = injectionPointMask.boundaryField()[f];

            if (!isA< calculatedFvPatchField<scalar> >(bf)) {
		// skip processor boundaries, etc.
		std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
                continue;
            }

            Foam::fvPatchField<scalar>& dbf = d.boundaryField()[f];

            // create mask
            tmp< Field<scalar> > r = dbf;
            tmp< Field<scalar> > mask = pos(radius - r());    // indicator function for circle interior
            bf == Foam::max(bf, mask());
            
            // add velocity field
            Foam::fvPatchField<vector>& vf = injectionVelocities.boundaryField()[f];
            // compute velocity profile
            tmp< Field<scalar> > p = mask()*max(0.0, (1.0 - Foam::pow(r()/radius, profileExp)));
            tmp< Field<scalar> > pSf = p() * bf.patch().magSf();
            tmp< Field<vector> > v = (-flowRate / (SMALL + Foam::gSum(pSf()))) * p() * bf.patch().nf();
            vf == vf + v();

            // add temperature field
            Foam::fvPatchField<scalar>& Tf = injectionTemperatures.boundaryField()[f];
            Tf == (mask() * Temp) + (1.0 - mask()) * Tf;

            // compute initial fiber oritentation at injection points
            tmp < Field<vector> > nf = injectionPointMask.boundaryField()[f].patch().nf();
            tmp < Foam::Field<symmTensor> > dir = sqr(nf());
            Foam::fvPatchField<symmTensor>& A2f = injectionA2.boundaryField()[f];
            tmp< Field<scalar> > pA2 = Foam::pow(p(), profileExpA2);
            tmp< Field<symmTensor> > s = (1.0 - pA2())*dir() + pA2()*A2iso;
            A2f == (1 - mask())*A2f + mask()*s();

            // compute initial tensor B for fast exact closure
            // FIXME: Assuming isotropic state for B2
            Foam::fvPatchField<symmTensor>& B2f = injectionB2.boundaryField()[f];
            tmp< Field<symmTensor> > q = (1.0 - pA2())*B2iso + pA2()*B2iso;
            B2f == (1 - mask())*B2f + mask()*q();
        }
    }
    else if (type == "circle")
    {
        // read cylinder dimensions
        vector axis(ip.lookup("axis"));
        vector center(ip.lookup("center"));
        scalar radius(readScalar(ip.lookup("radius")));
        scalar height(readScalar(ip.lookup("height")));

        // normalize axis vector
        axis /= Foam::mag(axis);
        
        // mark the injection regions in each boundary field
        forAll(injectionPointMask.boundaryField(), f)
        {
            Foam::fvPatchField<scalar>& bf = injectionPointMask.boundaryField()[f];

            if (!isA< calculatedFvPatchField<scalar> >(bf)) {
		// skip processor boundaries, etc.
		std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
                continue;
            }

            // create mask
            tmp< Field<vector> > d = bf.patch().Cf() - center;  // vector to cylinder center
            tmp< Field<scalar> > y = d() & axis;  // distance parallel to axis
            tmp< Field<vector> > r = (center + y()*axis - bf.patch().Cf());  // distance perpendicular to axis
            tmp< Field<scalar> > mask = pos(height - Foam::mag(y()))*pos(radius - Foam::mag(r()));    // indicator function for circle interior
            bf == Foam::max(bf, mask());
            
            // add velocity field
            Foam::fvPatchField<vector>& vf = injectionVelocities.boundaryField()[f];
            // compute velocity profile
            tmp< Field<scalar> > p = mask()*max(0.0, (1.0 - Foam::pow(Foam::mag(r())/radius, profileExp)));
            tmp< Field<scalar> > pSf = p() * bf.patch().magSf();
            tmp< Field<vector> > v = (-flowRate / (SMALL + Foam::gSum(pSf()))) * p() * bf.patch().nf();
            vf == vf + v();

            // add temperature field
            Foam::fvPatchField<scalar>& Tf = injectionTemperatures.boundaryField()[f];
            Tf == (mask() * Temp) + (1.0 - mask()) * Tf;

            // compute initial fiber oritentation at injection points
            tmp < Field<vector> > nf = injectionPointMask.boundaryField()[f].patch().nf();
            tmp < Foam::Field<symmTensor> > dir = sqr(nf());
            Foam::fvPatchField<symmTensor>& A2f = injectionA2.boundaryField()[f];
            tmp< Field<scalar> > pA2 = Foam::pow(p(), profileExpA2);
            tmp< Field<symmTensor> > s = (1.0 - pA2())*dir() + pA2()*A2iso;
            A2f == (1 - mask())*A2f + mask()*s();

            // compute initial tensor B for fast exact closure
            // FIXME: Assuming isotropic state for B2
            Foam::fvPatchField<symmTensor>& B2f = injectionB2.boundaryField()[f];
            tmp< Field<symmTensor> > q = (1.0 - pA2())*B2iso + pA2()*B2iso;
            B2f == (1 - mask())*B2f + mask()*q();
        }
    }
    else if (type == "patch")
    {
        word patchName = ip.lookup("patchName");
        bool patchFound = false;
        
        // mark the selected patch as injection region
        forAll(injectionPointMask.boundaryField(), f)
        {
            Foam::fvPatchField<scalar>& bf = injectionPointMask.boundaryField()[f];

            if (!isA< calculatedFvPatchField<scalar> >(bf)) {
		// skip processor boundaries, etc.
		// std::cout << "skipping " << bf.patch().name() << " of type " << bf.type() << std::endl;
                continue;
            }

            // create mask
            if (bf.patch().name() != patchName) continue;
            patchFound = true;
            bf == 1.0;
            
            // add velocity field
            Foam::fvPatchField<vector>& vf = injectionVelocities.boundaryField()[f];
            const pointField& pf = mesh.boundaryMesh()[f].points();
            scalar xMin = Foam::gMin(pf.component(0)), xMax = Foam::gMax(pf.component(0));
            scalar yMin = Foam::gMin(pf.component(1)), yMax = Foam::gMax(pf.component(1));
            scalar zMin = Foam::gMin(pf.component(2)), zMax = Foam::gMax(pf.component(2));
            // compute velocity profile
            tmp< Field<scalar> > p = max(0.0, (
                Foam::pow(4.0*(bf.patch().Cf().component(0) - xMin + SMALL)*(xMax + SMALL - bf.patch().Cf().component(0))/pow(xMax - xMin + 2*SMALL, 2), 0.5*profileExp)*
                Foam::pow(4.0*(bf.patch().Cf().component(1) - yMin + SMALL)*(yMax + SMALL - bf.patch().Cf().component(1))/pow(yMax - yMin + 2*SMALL, 2), 0.5*profileExp)*
                Foam::pow(4.0*(bf.patch().Cf().component(2) - zMin + SMALL)*(zMax + SMALL - bf.patch().Cf().component(2))/pow(zMax - zMin + 2*SMALL, 2), 0.5*profileExp)));
            tmp< Field<vector> > v = bf.patch().nf() * p();
            tmp< Field<vector> > vs = (-flowRate / Foam::gSum(v() & bf.patch().Sf())) * v();
            vf == vf + vs();

            // add temperature field
            Foam::fvPatchField<scalar>& Tf = injectionTemperatures.boundaryField()[f];
            Tf == Temp;

            // compute initial fiber oritentation at injection points
            tmp < Field<vector> > nf = injectionPointMask.boundaryField()[f].patch().nf();
            tmp < Foam::Field<symmTensor> > dir = sqr(nf());
            Foam::fvPatchField<symmTensor>& A2f = injectionA2.boundaryField()[f];
            tmp< Field<scalar> > pA2 = Foam::pow(p(), profileExpA2);
            tmp< Field<symmTensor> > s = (1.0 - pA2())*dir() + pA2()*A2iso;
            A2f == s();

            // compute initial tensor B for fast exact closure
            // FIXME: Assuming isotropic state for B2
            Foam::fvPatchField<symmTensor>& B2f = injectionB2.boundaryField()[f];
            tmp< Field<symmTensor> > q = (1.0 - pA2())*B2iso + pA2()*B2iso;
            B2f == q();
        }
        
        if (!patchFound) {
            FatalErrorIn("injectionPoints")
                << "patch not found: " << patchName
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("injectionPoints")
            << "invalid injection point type: " << type
            << abort(FatalError);
    }
}

#if 0
injectionPointMask.write();
injectionVelocities.write();
injectionTemperatures.write();
injectionA2.write();
injectionB2.write();
#endif

U.boundaryField() = injectionVelocities.boundaryField();
A2.boundaryField() = injectionA2.boundaryField();
B2.boundaryField() = injectionB2.boundaryField();
T.boundaryField() = injectionTemperatures.boundaryField();
alpha1.boundaryField() = injectionPointMask.boundaryField();

Info << "Expected fill time: " << (meshV/sumFlowRates) << endl;

