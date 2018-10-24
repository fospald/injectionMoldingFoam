
// compute eigenvalues and eigenvectors for A2
// symmetrize A2 and renormalize A2 by its trace

// SMALL = 1e-15

// TODO: make this more efficient using symmetry of A2

#if 0
{
	// Testing eigenvalue decomposition: ok

	scalar* pA2[3];
	scalar* pEvecs[3];

	symmTensor A2s(1,0,0,0,0,0);
	tensor A2t(A2s.xx(), A2s.xy(), A2s.xz(), A2s.xy(), A2s.yy(), A2s.yz(), A2s.xz(), A2s.yz(), A2s.zz());
	tensor evec;
	vector eval;
	scalar* evals = eval.v_;
	scalar* evecs = evec.v_;

	pA2[0] = A2t.v_; 
	pA2[1] = A2t.v_ + 3; 
	pA2[2] = A2t.v_ + 6; 

	pEvecs[0] = evecs;
	pEvecs[1] = evecs + 3;
	pEvecs[2] = evecs + 6;

//	vtkMath::Jacobi(pA2, evals, pEvecs);
	eigen_decomposition(pA2, pEvecs, evals);

	Info << "#############1'" << ((evec & DiagTensor<scalar>(eval[0], eval[1], eval[2]) & evec.T()) - A2s) << endl;
	Info << "#############2'" << (evec.T() & evec) << endl;
	Info << "#############3'" << (evec & evec.T()) << endl;
	Info << "#############5'" << (evec.T() & A2s & evec) << endl;
	Info << "#############5'" << (evec & A2s & evec.T()) << endl;
	exit(1);
}
#endif

try
{
	scalar* pA2[3];
	scalar* pEvecs[3];

	forAll(mesh.C(), i)
	{
		symmTensor& A2s = A2[i];
		tensor A2t(A2s.xx(), A2s.xy(), A2s.xz(), A2s.xy(), A2s.yy(), A2s.yz(), A2s.xz(), A2s.yz(), A2s.zz());
		scalar* evals = eigenValuesA2[i].v_;
		scalar* evecs = eigenVectorsA2[i].v_;

		pA2[0] = A2t.v_; 
		pA2[1] = A2t.v_ + 3; 
		pA2[2] = A2t.v_ + 6; 

		pEvecs[0] = evecs;
		pEvecs[1] = evecs + 3;
		pEvecs[2] = evecs + 6;

	//	vtkMath::Jacobi(pA2, evals, pEvecs);
		eigen_decomposition(pA2, pEvecs, evals);

		evals[0] = max(evals[0], SMALL);
		evals[1] = max(evals[1], SMALL);
		evals[2] = max(evals[2], SMALL);

		scalar escale = 1.0 / (evals[0] + evals[1] + evals[2]);

		evals[0] *= escale;
		evals[1] *= escale;
		evals[2] *= escale;

		A2s = symm(eigenVectorsA2[i] & (DiagTensor<scalar>(evals[0], evals[1], evals[2]) & eigenVectorsA2[i].T()));

#if 0
		continue;


		tensor A2t(A2s.xx(), A2s.xy(), A2s.xz(), A2s.xy(), A2s.yy(), A2s.yz(), A2s.xz(), A2s.yz(), A2s.zz());
		scalar* evals = eigenValuesA2[i].v_;
		scalar* evecs = eigenVectorsA2[i].v_;

		pA2[0] = A2t.v_; 
		pA2[1] = A2t.v_ + 3; 
		pA2[2] = A2t.v_ + 6; 

		pEvecs[0] = evecs;
		pEvecs[1] = evecs + 3;
		pEvecs[2] = evecs + 6;

	//	vtkMath::Jacobi(pA2, evals, pEvecs);
		eigen_decomposition(pA2, pEvecs, evals);

#if 0
		evals[0] = max(evals[0], SMALL);
		evals[1] = max(evals[1], SMALL);
		evals[2] = max(evals[2], SMALL);

#endif

		symmTensor& B2s = B2[i];
		eigenValuesB2[i] = eigenValues2(B2s);
		scalar* evalsB = eigenValuesB2[i].v_;


//Info << "i: " << i << "e: " << eigenValuesA2[i] << "e: " << eigenValuesB2[i] << " " << B2s << endl;

/*
		evals[0] = max(evals[0], 0);
		evals[1] = max(evals[1], 0);
		evals[2] = max(evals[2], 0);
*/

		if ((evals[0] + evals[1] + evals[2]) < 1e-10) Info << "scaling problem i: " << i << " " << eigenValuesA2[i] << endl;
		scalar escale = 1.0 / (evals[0] + evals[1] + evals[2]);


		evals[0] *= escale;
		evals[1] *= escale;
		evals[2] *= escale;

		A2s = symm(eigenVectorsA2[i] & (DiagTensor<scalar>(evals[0], evals[1], evals[2]) & eigenVectorsA2[i].T()));
//		scalar ascale = 1.0 / (A2s.xx() + A2s.yy() + A2s.zz());
//		A2s *= ascale;

/*		evalsB[0] = max(evalsB[0], 0);
		evalsB[1] = max(evalsB[1], 0);
		evalsB[2] = max(evalsB[2], 0);
*/

		if (evalsB[0] * evalsB[1] * evalsB[2] < 1e-10) Info << "scaling problem B i: " << i << " " << eigenValuesB2[i] << endl;
		scalar bscale = 1.0 / std::pow(evalsB[0] * evalsB[1] * evalsB[2], 1.0/3.0);

		evalsB[0] *= bscale;
		evalsB[1] *= bscale;
		evalsB[2] *= bscale;


		B2s = symm(eigenVectorsA2[i] & (DiagTensor<scalar>(evalsB[0], evalsB[1], evalsB[2]) & eigenVectorsA2[i].T()));
		
/*
		scalar A[3][3];
		scalar V[3][3];
		scalar d[3];

		A[0][0] = A2s.xx();
		A[1][1] = A2s.yy();
		A[2][2] = A2s.zz();
		A[0][1] = A[1][0] = A2s.xy();
		A[0][2] = A[2][0] = A2s.xz();
		A[1][2] = A[2][1] = A2s.yz();

//		Info << "I: " << A2s << endl;
		eigen_decomposition(A, V, d);

		tensor evecs(V[0][0], V[0][1], V[0][2], V[1][0], V[1][1], V[1][2], V[2][0], V[2][1], V[2][2]); 
//		symmTensor evecs(V[0][0], V[0][1], V[0][2], V[1][1], V[1][2], V[2][2]); 
		vector evals(d[0], d[1], d[2]);
//		Info << evals << endl;
//		Info << evecs << endl;

		evals.x() = max(evals.x(), SMALL);
		evals.y() = max(evals.y(), SMALL);
		evals.z() = max(evals.z(), SMALL);

		scalar scale = 1.0 / (evals.x() + evals.y() + evals.z());

		evals.x() *= scale;
		evals.y() *= scale;
		evals.z() *= scale;

		eigenValuesA2[i] = evals;
		eigenVectorsA2[i] = evecs;
//		A2[i] = symm(evecs & diagTensor(evals.x(), evals.y(), evals.z()) & evecs.T());
		A2[i] = A2s*scale;

//		Info << "O: " << A2[i] << endl;


*/



/*
//		Info << A2s << endl;
		vector evals = eigenValues2(A2s);
//		Info << evals << endl;

		vector ev1 = eigenVector(A2s, evals.x());
//		Info << ev1 << endl;
		vector ev2 = eigenVector(A2s, evals.y());
//		Info << ev2 << endl;

		if (mag(ev1) < SMALL) {
			ev1 = vector(1.0, 0.0, 0.0);
		}

		if (fabs(ev1 & ev2) > SMALL || mag(ev2) < SMALL) {

			ev2 = ev2 - (ev1 & ev2) * ev1 / sqr(mag(ev1));

			if (mag(ev2) < SMALL) {
				ev2.x() = -ev1.y();
				ev2.y() = ev1.x() + ev1.z();
				ev2.z() = -ev1.y();
			}

			ev2 /= mag(ev2);
		}

		vector ev3 = ev1 ^ ev2;
//		Info << ev3 << endl;

		evals.x() = max(evals.x(), SMALL);
		evals.y() = max(evals.y(), SMALL);
		evals.z() = max(evals.z(), SMALL);

		scalar scale = 1.0 / (evals.x() + evals.y() + evals.z());

		evals.x() *= scale;
		evals.y() *= scale;
		evals.z() *= scale;

		eigenValuesA2[i] = evals;
		eigenVectorsA2[i] = tensor(ev1, ev2, ev3);

		A2[i] = eigenVectorsA2[i] & diagTensor(evals.x(), evals.y(), evals.z()) & eigenVectorsA2[i].T();

//		Info << A2[i] << endl;
*/		

/*
		A2s.xx() = max(A2s.xx(), SMALL);
		A2s.yy() = max(A2s.yy(), SMALL);
		A2s.zz() = max(A2s.zz(), SMALL);

		scalar maxOff = max(max(fabs(A2s.xy()), fabs(A2s.xz())), fabs(A2s.yz()));
		scalar scaleOff = 1.0/max(maxOff, 1.0);

		A2s.xy() *= scaleOff;
		A2s.xz() *= scaleOff;
		A2s.yz() *= scaleOff;

		if (CorrectA2) {
			scalar traceA2 = A2s.xx() + A2s.yy() + A2s.zz();
			A2s *= 1.0 / traceA2;
		}

		eig3sym(A2s, eigenVectorsA2[i], eigenValuesA2[i]);
		A2[i] = A2s;
*/
#endif
	}
}
catch (...)
{
	Info << "Eigenvalue exception !!!!!" << endl;
}


