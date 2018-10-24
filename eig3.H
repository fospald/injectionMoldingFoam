

//extern "C" {

void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);
void eigen_decomposition2(double A[3][3], double V[3][3], double d[3]);
void eigen_decomposition(double** A, double** V, double* d);

//}

// Return eigenvalues in ascending order of absolute values
vector eigenValues2(const tensor& t)
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

 //           if (disc >= -SMALL)
            {
                scalar q = -0.5*::sqrt(max(0.0, disc));

                i = 0;
                ii = -0.5*a + q;
                iii = -0.5*a - q;
            }
  /*          else
            {
                FatalErrorIn("eigenValues(const tensor&)")
                    << "zero and complex eigenvalues in tensor: " << t
                    << abort(FatalError);
            }
*/
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
                scalar sqrtQ = ::sqrt(Q);
                scalar theta = ::acos(R/(Q*sqrtQ));

                scalar m2SqrtQ = -2*sqrtQ;
                scalar aBy3 = a/3;

                i = m2SqrtQ*::cos(theta/3) - aBy3;
                ii =
                    m2SqrtQ
                   *::cos((theta + constant::mathematical::twoPi)/3)
                  - aBy3;
                iii =
                    m2SqrtQ
                   *::cos((theta - constant::mathematical::twoPi)/3)
                  - aBy3;
            }
            else
            {
//                scalar A = ::cbrt(R + ::sqrt(R2 - Q3));

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


