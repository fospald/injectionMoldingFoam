
#include <cassert>
#include <cmath>

#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>

void fec_dyad(const symmTensor& C, const symmTensor& N, symmTensor& res)
{
	res.xx() = C.xx()*N.xx() + C.xy()*N.yy() + C.xz()*N.zz();
	res.yy() = C.xy()*N.xx() + C.yy()*N.yy() + C.yz()*N.zz();
	res.zz() = C.xz()*N.xx() + C.yz()*N.yy() + C.zz()*N.zz();
}

void fec_dyad(const symmTensor& D, const symmTensor& DT, const symmTensor& N, symmTensor& res)
{
	fec_dyad(D, N, res);

	res.xy() = 2.0*DT.xy()*N.xy();
	res.xz() = 2.0*DT.xz()*N.xz();
	res.yz() = 2.0*DT.yz()*N.yz();
}

// Assumptions
// - eigenValuesA2 are sorted descending
// - eigenValuesB2 are sorted ascending
// - Omega = skew(gradU)
// - Gamma = symm(gradU)
void fec_source(const volSymmTensorField& A2, const volSymmTensorField& B2,
	const volTensorField eigenVectorsA2, const volVectorField& eigenValuesA2, const volVectorField& eigenValuesB2,
	const volTensorField& Omega, const volSymmTensorField& Gamma, scalar lambda, scalar Ci, const volScalarField& magD,
	volSymmTensorField& sourceA2, volSymmTensorField& sourceB2)
{
	scalar eps = 1e-4;

label i0 = 225;
Info << "eps = " << eps << endl;

	for(label i = 0; i < A2.size(); i++)
	{
		// FO tensors
		const symmTensor& A2i = A2[i];
		const symmTensor& B2i = B2[i];

		// eigenvalues of A2i and B2i
		const vector& eigenValuesA2i = eigenValuesA2[i];
		const vector& eigenValuesB2i = eigenValuesB2[i];

		// eigenvector change of basis matrix
		const tensor& V = eigenVectorsA2[i];

		// Eigenvalues of A2: a1 >= a2 >= a3
		scalar a1 = eigenValuesA2i.x();
		scalar a2 = eigenValuesA2i.y();
		scalar a3 = eigenValuesA2i.z();
		//std::assert((a1 >= a2) && (a2 >= a3));

		// Eigenvalues of B2: b1 <= b2 <= b3
		scalar b1 = eigenValuesB2i.x();
		scalar b2 = eigenValuesB2i.y();
		scalar b3 = eigenValuesB2i.z();
		//std::assert((b1 <= b2) && (b2 <= b3));

Info << "-------------------" << endl;
Info << "i: " << i << endl;
Info << "V: " << V << endl;
Info << "A2: " << A2i << endl;
Info << "B2: " << B2i << endl;
Info << "eigenValues(B2): " << (eigenValues2(B2i)) << endl;
Info << "a1 >= a2 >= a3: " << a1 << " " << a2 << " " << a3 << endl;
Info << "b1 <= b2 <= b3: " << b1 << " " << b2 << " " << b3 << endl;
Info << "V.T() & A & V: " << (V.T() & A2i & V) << endl;
Info << "V.T() & B & V: " << (V.T() & B2i & V) << endl;

if (!((a1 >= a2) && (a2 >= a3))) {
	FatalErrorIn("fec_source") << "invalid eigenvalue order " << (a1-a2) << " " << (a2-a3) << abort(FatalError);
}


if (!((b1 <= b2) && (b2 <= b3))) {
	FatalErrorIn("fec_source") << "invalid eigenvalue order" << abort(FatalError);
}

		scalar small = 1e-3;

		scalar b2_b1 = b2 - b1;
		scalar b3_b1 = b3 - b1;
		scalar b3_b2 = b3 - b2;

		scalar c1122, c1133, c2233, c1111, c2222, c3333;

		if (b2_b1 > eps && b3_b2 > eps)
		{
Info << "eigenvalues are not close together" << endl;
			// eigenvalues are not close together
			c1122 = 0.5*(a1 - a2)/b2_b1;
			c1133 = 0.5*(a1 - a3)/b3_b1;
			c2233 = 0.5*(a2 - a3)/b3_b2;
		}
		else if (b2_b1 <= eps && b3_b2 <= eps)
		{
Info << "eigenvalues are close together" << endl;
			// all eigenvalues are close together
			scalar c1 = b1 - 1.0;
			scalar c2 = b2 - 1.0;
			scalar c3 = b3 - 1.0;
			scalar* ci[3] = {&c1122, &c1133, &c2233};

			for (int i = 0; i < 3; i++) {
				if (i == 1) std::swap(c2, c3);
				if (i == 2) std::swap(c1, c3);
				*ci[i] = 1.0/10.0 - 3.0/28.0*c1 - 3.0/28.0*c2 - 1.0/28.0*c3 + 5/48*c1*c1
					+ 1.0/8.0*c1*c2 + 1.0/24.0*c1*c3 + 5.0/48.0*c2*c2 + 1.0/24.0*c2*c3
					+ 1.0/48.0*c3*c3 - 35.0/352.0*c1*c1*c1 - 45.0/352.0*c1*c1*c2
					- 15.0/352.0*c1*c1*c3 - 45.0/352.0*c1*c2*c2 - 9.0/176.0*c1*c2*c3
					- 9.0/352.0*c1*c3*c3 - 35.0/352.0*c2*c2*c2 - 15.0/352.0*c2*c2*c3
					- 9.0/352.0*c2*c3*c3 - 5.0/352.0*c3*c3*c3;
			}
			std::swap(c1, c3);
			std::swap(c2, c3);
		}
		else
		{
			// two eigenvalues are close together

			scalar* ci = &c1122;

			if (b3_b2 <= eps) {
Info << "two eigenvalues are close together case 1" << endl;
				c1122 = 0.5*(a1 - a2)/b2_b1;
				c1133 = 0.5*(a1 - a3)/b3_b1;
				ci = &c2233;
				std::swap(b2, b3);
				std::swap(b1, b3);
			}
			else {
Info << "two eigenvalues are close together case 2" << endl;
				ci = &c1122;
				c1133 = 0.5*(a1 - a3)/b3_b1;
				c2233 = 0.5*(a2 - a3)/b3_b2;
			}

			scalar b0 = 0.5*(b1 + b2);
			scalar epsilon = 0.5*(b1 - b2);
			scalar I[6];

			if (b0 > b3) {
				I[1] = 2.0/std::sqrt(b0 - b3)*std::acos(std::sqrt(b3/b0));
			}
			else { // b0 < b3
				I[1] = 2.0/std::sqrt(b3 - b0)*::acosh(std::sqrt(b3/b0));
			}

			for (int n = 1; n < 5; n++) {
				I[n+1] = ((2.0*n - 1.0)/(2.0*n*(b0 - b3)))*I[n] - std::sqrt(b3)/(n*pow(b0,n)*(b0 - b3));
			}

			*ci = 1.0/4.0*I[3] + 3.0/8.0*I[5]*epsilon*epsilon;

			if (b3_b2 <= eps) {
				std::swap(b1, b3);
				std::swap(b2, b3);
			}
		}

		// Eqn. 16
		c1111 = 0.5/b1 - c1122 - c1133;
		c2222 = 0.5/b2 - c1122 - c2233;
		c3333 = 0.5/b3 - c1133 - c2233;

		// tensor with C components
		symmTensor C(c1111, c1122, c1133, c2222, c2233, c3333);

Info << "C: " << C << endl;
	
		// MB tilde (Eqn. 21) for C
		symmTensor CT(0.0, c1122, c1133, 0.0, c2233, 0.0);

Info << "CT: " << CT << endl;
	
		// Eqn. 20
		symmTensor D(inv(C));

Info << "D: " << D << endl;

		// MB tilde (Eqn. 21) for D
		symmTensor DT(0.0, 1.0/(4.0*c1122), 1.0/(4.0*c1133), 0.0, 1.0/(4.0*c2233), 0.0);
	
		// -F(B) und G(A) Eqn. 11/12
		// note: Omega and Gamma are do have an factor of 1/2 included
		symmTensor FB(-2.0*symm(B2i & (Omega[i] + lambda*Gamma[i])));
		symmTensor GA(2.0*Ci*magD[i]*(Foam::I - 3.0*A2i));

Info << "GA: " << GA << endl;
Info << "FB: " << FB << endl;

		// C : FB
		symmTensor CFB;
		fec_dyad(C, CT, symm(V.T() & FB & V), CFB);

Info << "CFB: " << CFB << endl;
	
		// D : GA
		symmTensor DGA;
		fec_dyad(D, DT, symm(V.T() & GA & V), DGA);

Info << "DGA: " << DGA << endl;
	
		// compute source terms
		sourceA2[i] = GA - symm(V & CFB & V.T());
		sourceB2[i] = FB - symm(V & DGA & V.T());

Info << "sourceA2: " << sourceA2[i] << endl;
Info << "sourceB2: " << sourceB2[i] << endl;

	}
}



void A_from_B(const vector& b, vector& a)
{
	a[0] = (1.0/3.0)*boost::math::ellint_rj(b[0], b[1], b[2], b[0]);
	a[1] = (1.0/3.0)*boost::math::ellint_rj(b[0], b[1], b[2], b[1]);
	a[2] = 1.0-a[0]-a[1];
}

scalar A4_from_AB_analytical(const vector& a, const vector& b, std::size_t i, std::size_t j)
{
	std::size_t k = 0;
	if (i == k || j == k) k++;
	if (i == k || j == k) k++;

	if (i != j) {
		return 0.5*(a[j] + b[i]*(a[i]-a[j])/((b[i]-b[j])));
	}

	std::size_t l = 0;
	if (i == l || j == l || k == l) l++;
	if (i == l || j == l || k == l) l++;

	return 3.0/2.0*a[j] - 0.5*(1.0 + b[i]*(a[k]-a[i])/(b[k]-b[i]) + b[i]*(a[l]-a[i])/(b[l]-b[i]));
}


void stable_A(vector& a)
{
	scalar eps = 1e-8;

	if (std::abs(a[0]-a[1]) < eps) {
		if (a[2] > 0.5) {
			a[0] += eps;
			a[2] -= eps;
		}
		else {
			a[0] -= eps;
			a[2] += eps;
		}
	}

	if (std::abs(a[0]-a[2]) < eps) {
		if (a[1] > 0.5) {
			a[0] += eps;
			a[1] -= eps;
		}
		else {
			a[0] -= eps;
			a[1] += eps;
		}
	}

	if (std::abs(a[1]-a[2]) < eps) {
		if (a[0] > 0.5) {
			a[1] += eps;
			a[0] -= eps;
		}
		else {
			a[1] -= eps;
			a[0] += eps;
		}
	}
}




void B_from_A_fixedpoint(const vector& ac, vector& b, scalar tol = 1e-10)
{
	scalar eps = 1e-32;

	vector a_log_a, e, ab, da;
	vector a = ac;

	for (std::size_t i = 0; i < 3; i++) {
		a[i] += eps;
		a_log_a[i] = a[i]*std::log(a[i]);
		e[i] = 1.65;
	}

	for (;;)
	{
		for (std::size_t i = 0; i < 3; i++) {
			b[i] = 1.0/(std::pow(a[i], e[i]) + eps);
		}

		scalar b_scale = std::pow(b[0]*b[1]*b[2], -1.0/3.0);
		for (std::size_t i = 0; i < 3; i++) {
			b[i] *= b_scale;
		}

		A_from_B(b, ab);
		da = a - ab;
		scalar err = 0;
		for (std::size_t i = 0; i < 3; i++) {
			err = std::max(std::abs(da[i]), err);
		}

		if (err <= tol) {
			break;
		}

		for (std::size_t i = 0; i < 3; i++) {
			e[i] += std::min(std::max(e[i]*da[i]/a_log_a[i], -1.0), 1.0);
			e[i]  = std::min(std::max(e[i], -64.0), 64.0);
		}
	};
}


// compute A4_iijj, k = {0,1,2} \ {i, j}
// approximation error is about 0.67%
scalar A4_from_A_approx_1(const vector& a, std::size_t i, std::size_t j, std::size_t k)
{
	static const scalar p[4] = {0.50838709, -0.32444078, 0.12406763, 1.02006009};

	return (a[i]*a[j]/(p[3] - a[k]))*(p[2]*a[k]*a[k] + p[1]*a[k] + p[0]);
}

// approximation error is about 0.4%
scalar A4_from_A_approx_2(const vector& a, std::size_t i, std::size_t j, std::size_t k)
{
	return (0.49855692 - 0.88951736*a[2] + 0.50273877*a[2]*a[2] - 0.11278585*a[2]*a[2]*a[2])*std::pow((a[0]*a[1])/(1e-8 + (a[0]+a[1])*(a[0]+a[1])), 1.0009197 - 0.09127384*a[2] - 0.05057483*a[2]*a[2]);
}

// approximation error is about 0.4%
scalar A4_from_A_approx_3(const vector& a, std::size_t i, std::size_t j, std::size_t k)
{
	return 0.5*a[i]*(a[i]+1.0) - 3.0/5.0*a[i]*a[j]*a[k];
}

// approximation error is about 0.4%
scalar A4_from_A_approx(const vector& a, std::size_t i, std::size_t j, std::size_t k)
{
	const scalar p = 0.5;
	return 0.5*a[i]*(a[i]+1.0); // - a[i]*a[j]*a[k]*(p*(1.0/3.0-a[i]) + 3.0/5.0);
}

// approximation error is about 0.4%
scalar A4_from_A_approx_new(const vector& a, std::size_t i, std::size_t j, std::size_t k)
{
	const scalar p = 0.5;
	return 0.5*a[i]*a[j]*(1 + 3*a[k])/(1 + 3*(a[i] + a[j])*a[k]);
}



scalar ORW3_A1111(const vector& a)
{
	scalar a1 = a[0];
	scalar a2 = std::max(a[1], a[2]);
	return -0.1480648093 + 0.8084618453*a1 + 0.3722003446*a1*a1 + 0.7765597096*a2 - 1.3431772379*a2*a2 - 1.7366749542*a1*a2 + 0.889594639*a1*a1*a2 + 1.7367571741*a1*a2*a2 - 0.0324756095*a1*a1*a1 + 0.6631716575*a2*a2*a2;
}

scalar ORW3_A2222(const vector& a)
{
	scalar a1 = std::max(a[0], a[2]);
	scalar a2 = a[1];
	return -0.2106349673 + 0.9092350296*a1 + -1.2840654776*a1*a1 + 1.1104441966*a2 + 0.1260059291*a2*a2 + -2.5375632310*a1*a2 + 1.9988098293*a1*a1*a2 + 1.4863151577*a1*a2*a2 + 0.5856304774*a1*a1*a1 + -0.0756740034*a2*a2*a2;
}

scalar ORW3_A3333(const vector& a)
{
	scalar a1 = a[0];
	scalar a2 = a[1];
	return 0.4868019601 + 0.5776328438*a1 + -2.2462007509*a1*a1 + 0.4605743789*a2 + -1.9088154281*a2*a2 + -4.8900459209*a1*a2 + 4.0544348937*a1*a1*a2 + 3.8542602127*a1*a2*a2 + 1.1817992322*a1*a1*a1 + 0.9512305286*a2*a2*a2;
}


/*
Cmk m=1          m=2           m=3
K=1 -0.1480648093 -0.2106349673 0.4868019601
K=2 0.8084618453 0.9092350296 0.5776328438
K=3 0.3722003446 -1.2840654776 -2.2462007509
K=4 0.7765597096 1.1104441966 0.4605743789
K=5 -1.3431772379 0.1260059291 -1.9088154281
K=6 -1.7366749542 -2.5375632310 -4.8900459209
K=7 0.8895946393 1.9988098293 4.0544348937
K=8 1.7367571741 1.4863151577 3.8542602127
K=9 -0.0324756095 0.5856304774 1.1817992322
K=10 0.6631716575 -0.0756740034 0.9512305286
*/


void ORW3_closure(Foam::fvMesh& mesh, const volSymmTensorField& A2, const volTensorField eigenVectorsA2, const volVectorField& eigenValuesA2, const volSymmTensorField& D, volSymmTensorField& ClosureA4)
{
	forAll(mesh.C(), k)
	{
		// FO tensors
		// const symmTensor& A2i = A2[i];

		// eigenvalues / vectors of A2i
		const vector& a = eigenValuesA2[k];
		const tensor& V = eigenVectorsA2[k];

		// compute six coefficients for A4
		scalar A1111 = ORW3_A1111(a);
		scalar A2222 = ORW3_A2222(a);
		scalar A3333 = ORW3_A3333(a);
		scalar A1122 = 0.5*(a[0] + a[1] - a[2] - A1111 - A2222 + A3333);
		scalar A2233 = 0.5*(a[1] + a[2] - a[0] - A2222 - A3333 + A1111);
		scalar A1133 = 0.5*(a[0] + a[2] - a[1] - A1111 - A3333 + A2222);

		tensor DD(
			A1111/3.0, A1122, A1133,
			A1122, A2222/3.0, A2233,
			A1133, A2233, A3333/3.0);

		tensor X(0, 0, 0, 0, 0, 0, 0, 0, 0);

		for(label i = 0; i < 3; i++) {
			vector vi(V[i+0], V[i+3], V[i+6]);
			for(label j = 0; j < 3; j++) {
				vector vj(V[j+0], V[j+3], V[j+6]);
				X += (DD[3*i + j]*vi)*(vi*(vj & D[k] & vj) + vj*((vi & D[k] & vj) + (vj & D[k] & vi)));
			}
		}

		ClosureA4[k] = symm(X);
	}
}


void approx_exact_new_closure(Foam::fvMesh& mesh, const volSymmTensorField& A2, const volTensorField eigenVectorsA2, const volVectorField& eigenValuesA2, const volSymmTensorField& D, volSymmTensorField& ClosureA4)
{
	forAll(mesh.C(), k)
	{
		// FO tensors
		// const symmTensor& A2i = A2[i];

		// eigenvalues / vectors of A2i
		const vector& a = eigenValuesA2[k];
		const tensor& V = eigenVectorsA2[k];

		// compute six coefficients for A4
#if 0
		scalar A1122 = A4_from_A_approx(a, 0, 1, 2); //0.5*a[0]*a[1];
		scalar A2233 = A4_from_A_approx(a, 1, 2, 0); //0.5*a[1]*a[2];
		scalar A1133 = A4_from_A_approx(a, 0, 2, 1); //0.5*a[0]*a[2];
		scalar A1111 = a[0] - A1122 - A1133;
		scalar A2222 = a[1] - A1122 - A2233;
		scalar A3333 = a[2] - A2233 - A1133;
#else
		scalar A1122 = A4_from_A_approx_new(a, 0, 1, 2);
		scalar A2233 = A4_from_A_approx_new(a, 1, 2, 1);
		scalar A1133 = A4_from_A_approx_new(a, 2, 0, 1);
		scalar A1111 = a[0] - A1122 - A1133;
		scalar A2222 = a[1] - A1122 - A2233;
		scalar A3333 = a[2] - A2233 - A1133;

#endif
		tensor DD(
			A1111/3.0, A1122, A1133,
			A1122, A2222/3.0, A2233,
			A1133, A2233, A3333/3.0);

		tensor X(0, 0, 0, 0, 0, 0, 0, 0, 0);

		for(label i = 0; i < 3; i++) {
			vector vi(V[i+0], V[i+3], V[i+6]);
			for(label j = 0; j < 3; j++) {
				vector vj(V[j+0], V[j+3], V[j+6]);
				X += (DD[3*i + j]*vi)*(vi*(vj & D[k] & vj) + vj*((vi & D[k] & vj) + (vj & D[k] & vi)));
			}
		}

		ClosureA4[k] = symm(X);
	}
}

void approx_exact_closure(Foam::fvMesh& mesh, const volSymmTensorField& A2, const volTensorField eigenVectorsA2, const volVectorField& eigenValuesA2, const volSymmTensorField& D, volSymmTensorField& ClosureA4)
{
	forAll(mesh.C(), k)
	{
		// FO tensors
		// const symmTensor& A2i = A2[i];

		// eigenvalues / vectors of A2i
		const vector& a = eigenValuesA2[k];
		const tensor& V = eigenVectorsA2[k];

		// compute six coefficients for A4
#if 0
		scalar A1122 = A4_from_A_approx(a, 0, 1, 2); //0.5*a[0]*a[1];
		scalar A2233 = A4_from_A_approx(a, 1, 2, 0); //0.5*a[1]*a[2];
		scalar A1133 = A4_from_A_approx(a, 0, 2, 1); //0.5*a[0]*a[2];
		scalar A1111 = a[0] - A1122 - A1133;
		scalar A2222 = a[1] - A1122 - A2233;
		scalar A3333 = a[2] - A2233 - A1133;
#else
		scalar A1111 = A4_from_A_approx(a, 0, 1, 2);
		scalar A2222 = A4_from_A_approx(a, 1, 0, 2);
		scalar A3333 = A4_from_A_approx(a, 2, 0, 1);
		scalar A1122 = 0.5*(a[0] + a[1] - a[2] - A1111 - A2222 + A3333);
		scalar A2233 = 0.5*(a[1] + a[2] - a[0] - A2222 - A3333 + A1111);
		scalar A1133 = 0.5*(a[0] + a[2] - a[1] - A1111 - A3333 + A2222);

#endif
		tensor DD(
			A1111/3.0, A1122, A1133,
			A1122, A2222/3.0, A2233,
			A1133, A2233, A3333/3.0);

		tensor X(0, 0, 0, 0, 0, 0, 0, 0, 0);

		for(label i = 0; i < 3; i++) {
			vector vi(V[i+0], V[i+3], V[i+6]);
			for(label j = 0; j < 3; j++) {
				vector vj(V[j+0], V[j+3], V[j+6]);
				X += (DD[3*i + j]*vi)*(vi*(vj & D[k] & vj) + vj*((vi & D[k] & vj) + (vj & D[k] & vi)));
			}
		}

		ClosureA4[k] = symm(X);
	}
}


void exact_closure(Foam::fvMesh& mesh, const volSymmTensorField& A2, const volTensorField eigenVectorsA2, const volVectorField& eigenValuesA2, const volSymmTensorField& D, volSymmTensorField& ClosureA4)
{
	forAll(mesh.C(), k)
	{
		// FO tensors
		// const symmTensor& A2i = A2[i];

		// eigenvalues / vectors of A2i
		vector a = eigenValuesA2[k];
		stable_A(a);
		const tensor& V = eigenVectorsA2[k];

		// compute b
		vector b;
		B_from_A_fixedpoint(a, b);

		// compute six coefficients for A4
		scalar A1122 = A4_from_AB_analytical(a, b, 0, 1);
		scalar A2233 = A4_from_AB_analytical(a, b, 1, 2);
		scalar A1133 = A4_from_AB_analytical(a, b, 0, 2);
		scalar A1111 = a[0] - A1122 - A1133;
		scalar A2222 = a[1] - A1122 - A2233;
		scalar A3333 = a[2] - A2233 - A1133;
		
		tensor DD(
			A1111/3.0, A1122, A1133,
			A1122, A2222/3.0, A2233,
			A1133, A2233, A3333/3.0);

		tensor X(0, 0, 0, 0, 0, 0, 0, 0, 0);

		for(label i = 0; i < 3; i++) {
			vector vi(V[i+0], V[i+3], V[i+6]);
			for(label j = 0; j < 3; j++) {
				vector vj(V[j+0], V[j+3], V[j+6]);
				X += (DD[3*i + j]*vi)*(vi*(vj & D[k] & vj) + vj*((vi & D[k] & vj) + (vj & D[k] & vi)));
			}
		}

		ClosureA4[k] = symm(X);
	}
}


void orthotropic_smooth_closure(Foam::fvMesh& mesh, const volSymmTensorField& A2, const volTensorField eigenVectorsA2, const volVectorField& eigenValuesA2, const volSymmTensorField& D, volSymmTensorField& ClosureA4)
{
	forAll(mesh.C(), k)
	{
		const symmTensor& Di = D[k];

		// eigenvalues / vectors of A2i
		vector a = eigenValuesA2[k];
		const tensor& V = eigenVectorsA2[k];
		
		scalar lambda1 = a[0], lambda2 = a[1], lambda3 = a[2];
		vector p1(V[0+0], V[0+3], V[0+6]);
		vector p2(V[1+0], V[1+3], V[1+6]);
		//vector p3(V[2+0], V[2+3], V[2+6]);
		
		vector Dp1 = Di & p1;
		vector Dp2 = Di & p2;

		ClosureA4[k] = symm(
			(((lambda1-lambda2)*p1) & Dp1)*(p1*p1) + 
			(lambda3/5.0)*((Di && I)*I + 2.0*Di) + 
			(0.25*(lambda2-lambda3))*(
				(3.0*(p1 & Dp1) + (p2 & Dp2))*(p1*p1) +
				(3.0*(p2 & Dp2) + (p1 & Dp1))*(p2*p2) +
				(p1 & Dp2)*(p1*p2 + p2*p1)
			)
		);
	}
}




