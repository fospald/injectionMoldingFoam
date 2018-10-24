
#include "tensor.H"
#include "mathematicalConstants.H"
#include "eig3.H"


namespace Foam
{

void eig3sym(const tensor& t, tensor& e, vector& x)
{
	double A[3][3];
	double E[3][3];
	double X[3];

	A[0][0] = t.xx();
	A[1][1] = t.yy();
	A[2][2] = t.zz();
	A[0][1] = A[1][0] = t.xy();
	A[0][2] = A[2][0] = t.xz();
	A[1][2] = A[2][1] = t.yz();

	eigen_decomposition(A, E, X);

	e = tensor(E[0][0], E[0][1], E[0][2], E[1][0], E[1][1], E[1][2], E[2][0], E[2][1], E[2][2]);
	x = vector(X[0], X[1], X[2]);
}


// Return eigenvalues in ascending order of absolute values
vector eigenValuesSymmTensor(const tensor& t)
{
    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if
    (
        (
            mag(t.xy()) + mag(t.xz()) + mag(t.yx())
          + mag(t.yz()) + mag(t.zx()) + mag(t.zy())
        )
      < SMALL
    )
    {
        // diagonal matrix
        i = t.xx();
        ii = t.yy();
        iii = t.zz();
    }
    else
    {
        scalar a = -t.xx() - t.yy() - t.zz();

        scalar b = t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
            - t.xy()*t.yx() - t.xz()*t.zx() - t.yz()*t.zy();

        scalar c = - t.xx()*t.yy()*t.zz() - t.xy()*t.yz()*t.zx()
            - t.xz()*t.yx()*t.zy() + t.xz()*t.yy()*t.zx()
            + t.xy()*t.yx()*t.zz() + t.xx()*t.yz()*t.zy();

        // If there is a zero root
        if (mag(c) < 1.0e-100)
        {
            scalar disc = sqr(a) - 4*b;

            if (disc >= -SMALL)
            {
                scalar q = -0.5*sqrt(max(0.0, disc));

                i = 0;
                ii = -0.5*a + q;
                iii = -0.5*a - q;
            }
            else
            {
                FatalErrorIn("eigenValues(const tensor&)")
                    << "zero and complex eigenvalues in tensor: " << t
                    << abort(FatalError);
            }
        }
        else
        {
            scalar Q = (a*a - 3*b)/9;
            scalar R = (2*a*a*a - 9*a*b + 27*c)/54;

            scalar R2 = sqr(R);
            scalar Q3 = pow3(Q);

            // Three different real roots
            if (R2 < Q3)
            {
                scalar sqrtQ = sqrt(Q);
                scalar theta = acos(R/(Q*sqrtQ));

                scalar m2SqrtQ = -2*sqrtQ;
                scalar aBy3 = a/3;

                i = m2SqrtQ*cos(theta/3) - aBy3;
                ii =
                    m2SqrtQ
                   *cos((theta + constant::mathematical::twoPi)/3)
                  - aBy3;
                iii =
                    m2SqrtQ
                   *cos((theta - constant::mathematical::twoPi)/3)
                  - aBy3;
            }
            else
            {
 //               scalar A = cbrt(R + sqrt(R2 - Q3));

                // Three equal real roots
//                if (A < SMALL)
                {
                    scalar root = -a/3;
                    return vector(root, root, root);
                }
 /*               else
                {
                    // Complex roots
                    WarningIn("eigenValues(const tensor&)")
                        << "complex eigenvalues detected for tensor: " << t
                        << endl;

                    return vector::zero;
                }
*/
            }
        }
    }


    // Sort the eigenvalues into ascending order
    if (i > ii)
    {
        Swap(i, ii);
    }

    if (ii > iii)
    {
        Swap(ii, iii);
    }

    if (i > ii)
    {
        Swap(i, ii);
    }

    return vector(i, ii, iii);
}



vector eigenVector(const tensor& t, const scalar lambda, int& flag)
{
    // Construct the matrix for the eigenvector problem
    tensor A(t - lambda*I);

    // Calculate the sub-determinants of the 3 components
    scalar sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    scalar sd1 = A.xx()*A.zz() - A.xz()*A.zx();
    scalar sd2 = A.xx()*A.yy() - A.xy()*A.yx();

    scalar magSd0 = mag(sd0);
    scalar magSd1 = mag(sd1);
    scalar magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 > magSd1 && magSd0 > magSd2 && magSd0 > SMALL && (flag & 1) == 0)
    {
        vector ev
        (
            1,
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );
        ev /= mag(ev);
	flag += 1;

        return ev;
    }
    else if (magSd1 > magSd2 && magSd1 > SMALL && (flag & 2) == 0)
    {
        vector ev
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            1,
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );
        ev /= mag(ev);
	flag += 2;

        return ev;
    }
    else if (magSd2 > SMALL && (flag & 4) == 0)
    {
        vector ev
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );
        ev /= mag(ev);
	flag += 4;
        return ev;
    }

    if ((flag & 1) == 0) {
	flag += 1;
	return vector(1,0,0);
    }

    if ((flag & 2) == 0) {
	flag += 2;
        return vector(0,1,0);
    }

    flag += 4;
    return vector(0,0,1); 
}



tensor eigenVectors(const tensor& t, vector& evals)
{
    int flag = 0;

    tensor evs 
    (   
        eigenVector(t, evals.x(), flag),
        eigenVector(t, evals.y(), flag),
        eigenVector(t, evals.z(), flag)
    );  

    return evs;
}

} // End namespace Foam
