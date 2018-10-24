
/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
   domain Java Matrix library JAMA. */

#include <math.h>

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

/*
#define n 3

static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  int i,j,k;
  double f,g,h,hh;
  for (j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (i = n-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      f = d[i-1];
      g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      hh = f / (h + h);
      for (j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    h = d[i+1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (j = 0; j <= i; j++) {
        g = 0.0;
        for (k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  int i,j,m,l,k;
  double g,p,r,dl1,h,f,tst1,eps;
  double c,c2,c3,el1,s,s2;

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  f = 0.0;
  tst1 = 0.0;
  eps = pow(2.0,-52.0);
  for (l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        g = d[l];
        p = (d[l+1] - g) / (2.0 * e[l]);
        r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        dl1 = d[l+1];
        h = g - d[l];
        for (i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        c = 1.0;
        c2 = c;
        c3 = c;
        el1 = e[l+1];
        s = 0.0;
        s2 = 0.0;
        for (i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (i = 0; i < n-1; i++) {
    k = i;
    p = d[i];
    for (j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}

void eigen_decomposition2(double A[n][n], double V[n][n], double d[n])
{
  int i,j;
  double e[n];
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);
}

*/




#define EIG3_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau)

#define EIG3_MAX_ROTATIONS 20


// Jacobi iteration for the solution of eigenvectors/eigenvalues of a nxn
// real symmetric matrix. Square nxn matrix a; size of matrix in n;
// output eigenvalues in w; and output eigenvectors in v. Resulting
// eigenvalues/vectors are sorted in decreasing order; eigenvectors are
// normalized.
template<class T>
int JacobiN(T **a, int n, T *w, T **v)
{
  int i, j, k, iq, ip, numPos;
  T tresh, theta, tau, t, sm, s, h, g, c, tmp;
  T bspace[4], zspace[4];
  T *b = bspace;
  T *z = zspace;

  // only allocate memory if the matrix is large
  if (n > 4)
    {
    b = new T[n];
    z = new T[n];
    }

  // initialize
  for (ip=0; ip<n; ip++)
    {
    for (iq=0; iq<n; iq++)
      {
      v[ip][iq] = 0.0;
      }
    v[ip][ip] = 1.0;
    }
  for (ip=0; ip<n; ip++)
    {
    b[ip] = w[ip] = a[ip][ip];
    z[ip] = 0.0;
    }

  // begin rotation sequence
  for (i=0; i<EIG3_MAX_ROTATIONS; i++)
    {
    sm = 0.0;
    for (ip=0; ip<n-1; ip++)
      {
      for (iq=ip+1; iq<n; iq++)
        {
        sm += fabs(a[ip][iq]);
        }
      }
    if (sm == 0.0)
      {
      break;
      }

    if (i < 3)                                // first 3 sweeps
      {
      tresh = 0.2*sm/(n*n);
      }
    else
      {
      tresh = 0.0;
      }

    for (ip=0; ip<n-1; ip++)
      {
      for (iq=ip+1; iq<n; iq++)
        {
        g = 100.0*fabs(a[ip][iq]);

        // after 4 sweeps
        if (i > 3 && (fabs(w[ip])+g) == fabs(w[ip])
        && (fabs(w[iq])+g) == fabs(w[iq]))
          {
          a[ip][iq] = 0.0;
          }
        else if (fabs(a[ip][iq]) > tresh)
          {
          h = w[iq] - w[ip];
          if ( (fabs(h)+g) == fabs(h))
            {
            t = (a[ip][iq]) / h;
            }
          else
            {
            theta = 0.5*h / (a[ip][iq]);
            t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0)
              {
              t = -t;
              }
            }
          c = 1.0 / sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          w[ip] -= h;
          w[iq] += h;
          a[ip][iq]=0.0;

          // ip already shifted left by 1 unit
          for (j = 0;j <= ip-1;j++)
            {
            EIG3_ROTATE(a,j,ip,j,iq);
            }
          // ip and iq already shifted left by 1 unit
          for (j = ip+1;j <= iq-1;j++)
            {
            EIG3_ROTATE(a,ip,j,j,iq);
            }
          // iq already shifted left by 1 unit
          for (j=iq+1; j<n; j++)
            {
            EIG3_ROTATE(a,ip,j,iq,j);
            }
          for (j=0; j<n; j++)
            {
            EIG3_ROTATE(v,j,ip,j,iq);
            }
          }
        }
      }

    for (ip=0; ip<n; ip++)
      {
      b[ip] += z[ip];
      w[ip] = b[ip];
      z[ip] = 0.0;
      }
    }

  //// this is NEVER called
  if ( i >= EIG3_MAX_ROTATIONS )
    {
    // "Error extracting eigenfunctions"
    return 0;
    }

  // sort eigenfunctions                 these changes do not affect accuracy
  for (j=0; j<n-1; j++)                  // boundary incorrect
    {
    k = j;
    tmp = w[k];
    for (i=j+1; i<n; i++)                // boundary incorrect, shifted already
      {
      if (w[i] >= tmp)                   // why exchage if same?
        {
        k = i;
        tmp = w[k];
        }
      }
    if (k != j)
      {
      w[k] = w[j];
      w[j] = tmp;
      for (i=0; i<n; i++)
        {
        tmp = v[i][j];
        v[i][j] = v[i][k];
        v[i][k] = tmp;
        }
      }
    }
  // insure eigenvector consistency (i.e., Jacobi can compute vectors that
  // are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
  // reek havoc in hyperstreamline/other stuff. We will select the most
  // positive eigenvector.
  int ceil_half_n = (n >> 1) + (n & 1);
  for (j=0; j<n; j++)
    {
    for (numPos=0, i=0; i<n; i++)
      {
      if ( v[i][j] >= 0.0 )
        {
        numPos++;
        }
      }
//    if ( numPos < ceil(double(n)/double(2.0)) )
    if ( numPos < ceil_half_n)
      {
      for(i=0; i<n; i++)
        {
        v[i][j] *= -1.0;
        }
      }
    }

  if (n > 4)
    {
    delete [] b;
    delete [] z;
    }
  return 1;
}


//----------------------------------------------------------------------------
template<class T1, class T2>
inline void Transpose3x3(const T1 A[3][3], T2 AT[3][3])
{
  T2 tmp;
  tmp = A[1][0];
  AT[1][0] = A[0][1];
  AT[0][1] = tmp;
  tmp = A[2][0];
  AT[2][0] = A[0][2];
  AT[0][2] = tmp;
  tmp = A[2][1];
  AT[2][1] = A[1][2];
  AT[1][2] = tmp;

  AT[0][0] = A[0][0];
  AT[1][1] = A[1][1];
  AT[2][2] = A[2][2];
}

template<class T>
inline double Determinant3x3(T A[3][3])
{
  return A[0][0] * A[1][1] * A[2][2] + A[1][0] * A[2][1] * A[0][2] +
         A[2][0] * A[0][1] * A[1][2] - A[0][0] * A[2][1] * A[1][2] -
         A[1][0] * A[0][1] * A[2][2] - A[2][0] * A[1][1] * A[0][2];
}

template<class T>
 inline void Identity3x3(T A[3][3])
{
  for (int i = 0; i < 3; i++)
    {
    A[i][0] = A[i][1] = A[i][2] = T(0.0);
    A[i][i] = 1.0;
    }
}

// helper function, swap two 3-vectors
template<class T>
inline void SwapVectors3(T v1[3], T v2[3])
{
  for (int i = 0; i < 3; i++)
    {
    T tmp = v1[i];
    v1[i] = v2[i];
    v2[i] = tmp;
    }
}

// Cross product of two 3-vectors. Result (a x b) is stored in z[3].
inline void Cross(const double x[3], const double y[3], double z[3])
{
  double Zx = x[1] * y[2] - x[2] * y[1];
  double Zy = x[2] * y[0] - x[0] * y[2];
  double Zz = x[0] * y[1] - x[1] * y[0];
  z[0] = Zx; z[1] = Zy; z[2] = Zz;
}

// Compute the norm of 3-vector (double-precision version).
double Norm(const double x[3]) {
    return sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
}

inline double Normalize(double x[3])
{
  double den;
  if ( ( den = Norm( x ) ) != 0.0 )
    {
    for (int i=0; i < 3; i++)
      {
      x[i] /= den;
      }
    }
  return den;
}


// Extract the eigenvalues and eigenvectors from a 3x3 matrix.
// The eigenvectors (the columns of V) will be normalized.
// The eigenvectors are aligned optimally with the x, y, and z
// axes respectively.
template <class T1, class T2>
inline void Diagonalize3x3(const T1 A[3][3], T2 w[3], T2 V[3][3])
{
  int i,j,k,maxI;
  T2 tmp, maxVal;

  // do the matrix[3][3] to **matrix conversion for Jacobi
  T2 C[3][3];
  T2 *ATemp[3],*VTemp[3];
  for (i = 0; i < 3; i++)
    {
    C[i][0] = A[i][0];
    C[i][1] = A[i][1];
    C[i][2] = A[i][2];
    ATemp[i] = C[i];
    VTemp[i] = V[i];
    }

  // diagonalize using Jacobi
  JacobiN(ATemp,3,w,VTemp);

  // if all the eigenvalues are the same, return identity matrix
  if (w[0] == w[1] && w[0] == w[2])
    {
    Identity3x3(V);
    return;
    }

  // transpose temporarily, it makes it easier to sort the eigenvectors
  Transpose3x3(V,V);

  // if two eigenvalues are the same, re-orthogonalize to optimally line
  // up the eigenvectors with the x, y, and z axes
  for (i = 0; i < 3; i++)
    {
    if (w[(i+1)%3] == w[(i+2)%3]) // two eigenvalues are the same
      {
      // find maximum element of the independant eigenvector
      maxVal = fabs(V[i][0]);
      maxI = 0;
      for (j = 1; j < 3; j++)
        {
        if (maxVal < (tmp = fabs(V[i][j])))
          {
          maxVal = tmp;
          maxI = j;
          }
        }
      // swap the eigenvector into its proper position
      if (maxI != i)
        {
        tmp = w[maxI];
        w[maxI] = w[i];
        w[i] = tmp;
        SwapVectors3(V[i],V[maxI]);
        }
      // maximum element of eigenvector should be positive
      if (V[maxI][maxI] < 0)
        {
        V[maxI][0] = -V[maxI][0];
        V[maxI][1] = -V[maxI][1];
        V[maxI][2] = -V[maxI][2];
        }

      // re-orthogonalize the other two eigenvectors
      j = (maxI+1)%3;
      k = (maxI+2)%3;

      V[j][0] = 0.0;
      V[j][1] = 0.0;
      V[j][2] = 0.0;
      V[j][j] = 1.0;
      Cross(V[maxI],V[j],V[k]);
      Normalize(V[k]);
      Cross(V[k],V[maxI],V[j]);

      // transpose vectors back to columns
      Transpose3x3(V,V);
      return;
      }
    }

  // the three eigenvalues are different, just sort the eigenvectors
  // to align them with the x, y, and z axes

  // find the vector with the largest x element, make that vector
  // the first vector
  maxVal = fabs(V[0][0]);
  maxI = 0;
  for (i = 1; i < 3; i++)
    {
    if (maxVal < (tmp = fabs(V[i][0])))
      {
      maxVal = tmp;
      maxI = i;
      }
    }
  // swap eigenvalue and eigenvector
  if (maxI != 0)
    {
    tmp = w[maxI];
    w[maxI] = w[0];
    w[0] = tmp;
    SwapVectors3(V[maxI],V[0]);
    }
  // do the same for the y element
  if (fabs(V[1][1]) < fabs(V[2][1]))
    {
    tmp = w[2];
    w[2] = w[1];
    w[1] = tmp;
    SwapVectors3(V[2],V[1]);
    }

  // ensure that the sign of the eigenvectors is correct
  for (i = 0; i < 2; i++)
    {
    if (V[i][i] < 0)
      {
      V[i][0] = -V[i][0];
      V[i][1] = -V[i][1];
      V[i][2] = -V[i][2];
      }
    }
  // set sign of final eigenvector to ensure that determinant is positive
  if (Determinant3x3(V) < 0)
    {
    V[2][0] = -V[2][0];
    V[2][1] = -V[2][1];
    V[2][2] = -V[2][2];
    }

  // transpose the eigenvectors back again
  Transpose3x3(V,V);
}



void eigen_decomposition(double A[3][3], double V[3][3], double d[3])
{
	Diagonalize3x3(A, d, V);
}


// of symmetric matrix
int eigen_decomposition(double** A, double** V, double* d)
{
	return JacobiN(A, 3, d, V);
}


