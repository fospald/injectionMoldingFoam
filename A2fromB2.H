
// compute eigenvalues and eigenvectors for A2
// symmetrize A2 and renormalize A2 by its trace

// SMALL = 1e-15

// TODO: make this more efficient using symmetry of A2

try
{
	scalar* pB2[3];
	scalar* pEvecs[3];

	forAll(mesh.C(), i)
	{
		symmTensor& A2s = A2[i];
		symmTensor& B2s = B2[i];
		tensor B2t(B2s.xx(), B2s.xy(), B2s.xz(), B2s.xy(), B2s.yy(), B2s.yz(), B2s.xz(), B2s.yz(), B2s.zz());
		scalar* evalsA = eigenValuesA2[i].v_;
		scalar* evalsB = eigenValuesB2[i].v_;
		scalar* evecs = eigenVectorsA2[i].v_;

		pB2[0] = B2t.v_; 
		pB2[1] = B2t.v_ + 3; 
		pB2[2] = B2t.v_ + 6; 

		pEvecs[0] = evecs;
		pEvecs[1] = evecs + 3;
		pEvecs[2] = evecs + 6;

	//	vtkMath::Jacobi(pA2, evalsB, pEvecs);
		eigen_decomposition(pB2, pEvecs, evalsB);
	//	eigenVectorsB2[i] = eigenVectorsA2[i];
	
		scalar b1 = evalsB[0], b2 = evalsB[1], b3 = evalsB[2];

		scalar minb = std::min(b1, std::min(b2, b3));

		if (minb < 0) {
			minb *= 1.001;
			b1 -= minb;
			b2 -= minb;
			b3 -= minb;
		}

		scalar detB2inv = std::pow(b1*b2*b3, -1.0/3.0);

		b1 *= detB2inv;
		b2 *= detB2inv;
		b3 *= detB2inv;


		/*
		if (detB2inv < 1e-6 || detB2inv > 1e6 || std::isnan(detB2inv) || std::isinf(detB2inv)) {
			std::cout << "b1=" << b1 << " b2=" << b2 << " b3=" << b3 << std::endl;
			std::cout << "detB2inv=" << detB2inv << std::endl;
		}
		*/

	//	try {
			B2s *= detB2inv;
	//	}
	//	catch (...) {
		//	std::cout << "b1=" << b1 << " b2=" << b2 << " b3=" << b3 << std::endl;
		//	std::cout << "detB2inv=" << detB2inv << std::endl;
	//	}

		evalsA[0] = (1.0/3.0)*boost::math::ellint_rj(b1, b2, b3, b1);
		evalsA[1] = (1.0/3.0)*boost::math::ellint_rj(b1, b2, b3, b2);
		evalsA[2] = (1.0/3.0)*boost::math::ellint_rj(b1, b2, b3, b3);

		scalar tr = evalsA[0]+evalsA[1]+evalsA[2];
		symmTensor S(evalsA[0]/tr, 0.0, 0.0, evalsA[1]/tr, 0.0, evalsA[2]/tr);
		//A2s = symm(eigenVectorsA2[i].T() & (S & eigenVectorsA2[i]));
		A2s = symm(eigenVectorsA2[i] & (S & eigenVectorsA2[i].T()));
	}
}
catch (...)
{
	Info << "Eigenvalue exception !!!!!" << endl;
}


