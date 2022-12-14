#! /usr/bin/env python3
#
def bernstein_matrix_determinant ( n ):

#*****************************************************************************80
#
## bernstein_matrix_determinant() returns the determinant of the Bernstein matrix.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 March 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#  Output:
#
#    real VALUE, the determinant.
#
  from scipy.special import comb

  value = 1.0
  for i in range ( 0, n ):
    value = value * comb ( n - 1, i )

  return value

def bernstein_matrix_determinant_test ( ):

#*****************************************************************************80
#
## bernstein_matrix_determinant_test() tests bernstein_matrix_determinant().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_matrix_determinant_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_matrix_determinant() computes the determinant of' )
  print ( '  the Bernstein matrix.' )
  print ( '' )
  print ( '    N     ||A||            det(A)      np.linalg.det(A)' )
  print ( '' )

  for n in range ( 5, 16 ):

    a = bernstein_matrix ( n )
    a_norm_frobenius = r8mat_norm_fro ( n, n, a )

    d1 = bernstein_matrix_determinant ( n )
    d2 = np.linalg.det ( a )

    print ( '  %4d  %14g  %14g  %14g' % \
      ( n, a_norm_frobenius, d1, d2 ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_matrix_determinant_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_matrix_inverse ( n ):

#*****************************************************************************80
#
## bernstein_matrix_inverse() returns the inverse of the Bernstein matrix.
#
#  Discussion:
#
#    The inverse Bernstein matrix of order N is an NxN matrix A which can 
#    be used to transform a vector of Bernstein basis coefficients B
#    representing a polynomial P(X) to a corresponding power basis 
#    coefficient vector C:
#
#      C = A * B
#
#    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
#    Bernstein basis vectors as ((1-X)^(N-1), X*(1-X)^(N-2),...,X^(N-1)).
#
#  Example:
#
#    N = 5
#
#   1.0000    1.0000    1.0000    1.0000    1.0000
#        0    0.2500    0.5000    0.7500    1.0000
#        0         0    0.1667    0.5000    1.0000
#        0         0         0    0.2500    1.0000
#        0         0         0         0    1.0000
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#  Output:
#
#    real A(N,N), the inverse Bernstein matrix.
#
  from scipy.special import comb
  import numpy as np

  a = np.zeros ( ( n, n ) )

  for j in range ( 0, n ):
    for i in range ( 0, j + 1 ):
      a[i,j] = comb ( j, i ) / comb ( n - 1, i )

  return a

def bernstein_matrix_inverse_test ( ):

#*****************************************************************************80
#
## bernstein_matrix_inverse_test() tests bernstein_matrix_inverse().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_matrix_inverse_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_matrix returns a matrix A which transforms a' )
  print ( '  polynomial coefficient vector from the power basis to' )
  print ( '  the Bernstein basis.' )
  print ( '  bernstein_matrix_inverse computes the inverse B.' )
  print ( '' )
  print ( '    N     ||A||            ||B||      ||I-A*B||' )
  print ( '' )

  for n in range ( 5, 16 ):

    a = bernstein_matrix ( n )
    a_norm_frobenius = r8mat_norm_fro ( n, n, a )

    b = bernstein_matrix_inverse ( n )
    b_norm_frobenius = r8mat_norm_fro ( n, n, b )

    c = np.dot ( a, b )
    error_norm_frobenius = r8mat_is_identity ( n, c )

    print ( '  %4d  %14g  %14g  %14g' % \
      ( n, a_norm_frobenius, b_norm_frobenius, error_norm_frobenius ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_matrix_inverse_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_matrix ( n ):

#*****************************************************************************80
#
## bernstein_matrix() returns the Bernstein matrix.
#
#  Discussion:
#
#    The Bernstein matrix of order N is an NxN matrix A which can be used to
#    transform a vector of power basis coefficients C representing a polynomial 
#    P(X) to a corresponding Bernstein basis coefficient vector B:
#
#      B = A * C
#
#    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
#    Bernstein basis vectors as ((1-X)^(N-1), X*(1-X)^(N-2),...,X^(N-1)).
#
#  Example:
#
#    N = 5
#
#    1    -4     6    -4     1
#    0     4   -12    12    -4
#    0     0     6   -12     6
#    0     0     0     4    -4
#    0     0     0     0     1
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 March 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#  Output:
#
#    real A(N,N), the Bernstein matrix.
#
  from scipy.special import comb
  import numpy as np

  a = np.zeros ( ( n, n ) )

  for j in range ( 0, n ):
    for i in range ( 0, j + 1 ):
      a[i,j] = r8_mop ( j - i ) * comb ( n - 1 - i, j - i ) * comb ( n - 1, i )

  return a

def bernstein_matrix_test ( ):

#*****************************************************************************80
#
## bernstein_matrix_test() tests bernstein_matrix().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'bernstein_matrix_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_matrix() computes the Bernstein matrix.' )

  m = 5
  n = m

  a = bernstein_matrix ( n )
 
  r8mat_print ( m, n, a, '  Bernstein matrix:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_matrix_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_matrix_test2 ( ):

#*****************************************************************************80
#
## bernstein_matrix_test2() uses bernstein_matrix() to describe Bernstein polynomials.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_matrix_test2():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_matrix() returns a matrix which' )
  print ( '  transforms a polynomial coefficient vector' )
  print ( '  from the the Bernstein basis to the power basis.' )
  print ( '  We can use this to get explicit values of the' )
  print ( '  4-th degree Bernstein polynomial coefficients as' )
  print ( '' )
  print ( '    B(4,K)(X) = C4 * x^4' )
  print ( '              + C3 * x^3' )
  print ( '              + C2 * x^2' )
  print ( '              + C1 * x' )
  print ( '              + C0 * 1' )

  n = 5
  print ( '' )
  print ( '     K             C4             C3             C2            C1             C0' )
  print ( '' )

  a = bernstein_matrix ( n )

  for k in range ( 0, n ):

    x = np.zeros ( n )
    x[k] = 1.0

    ax = np.dot ( a, x )

    print ( '  %4d' % ( k ), end = '' )
    for i in range ( 0, n ):
      print ( '%14g' % ( ax[i] ), end = '' )
    print ( '' )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_matrix_test2():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_01_matrix ( m, n, x ):

#*****************************************************************************80
#
## bernstein_poly_01() evaluates the Bernstein polynomials defined on [0,1].
#
#  Discussion:
#
#    The Bernstein polynomials are assumed to be based on [0,1].
#
#  Formula:
#
#    B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
#
#  First values:
#
#    B(0,0)(X) = 1
#
#    B(1,0)(X) =      1-X
#    B(1,1)(X) =                X
#
#    B(2,0)(X) =     (1-X)^2
#    B(2,1)(X) = 2 * (1-X)   * X
#    B(2,2)(X) =                X^2
#
#    B(3,0)(X) =     (1-X)^3
#    B(3,1)(X) = 3 * (1-X)^2 * X
#    B(3,2)(X) = 3 * (1-X)   * X^2
#    B(3,3)(X) =               X^3
#
#    B(4,0)(X) =     (1-X)^4
#    B(4,1)(X) = 4 * (1-X)^3 * X
#    B(4,2)(X) = 6 * (1-X)^2 * X^2
#    B(4,3)(X) = 4 * (1-X)   * X^3
#    B(4,4)(X) =               X^4
#
#  Special values:
#
#    B(N,I)(X) has a unique maximum value at X = I/N.
#
#    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
#
#    B(N,I)(1/2) = C(N,K) / 2^N
#
#    For a fixed X and N, the polynomials add up to 1:
#
#      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 January 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the number of evaluation points.
#
#    integer N, the degree of the Bernstein polynomials to be
#    used.  For any N, there is a set of N+1 Bernstein polynomials,
#    each of degree N, which form a basis for polynomials on [0,1].
#
#    real X[M], the evaluation points.
#
#  Output:
#
#    real B[M,N+1], the values of the N+1 Bernstein polynomials
#    at the evaluation points.
#
  import numpy as np

  b = np.zeros ( [ m, n + 1 ] )

  for i in range ( 0, m ):

    if ( n == 0 ):
 
      b[i,0] = 1.0
 
    elif ( 0 < n ):
 
      b[i,0] = 1.0 - x[i]
      b[i,1] = x[i]
 
      for j in range ( 2, n + 1 ):
        b[i,j] = x[i] * b[i,j-1]
        for k in range ( j - 1, 0, -1 ):
          b[i,k] = x[i] * b[i,k-1] + ( 1.0 - x[i] ) * b[i,k]
        b[i,0] = ( 1.0 - x[i] ) * b[i,0]

  return b

def bernstein_poly_01_matrix_test ( ):

#*****************************************************************************80
#
## bernstein_poly_01_matrix_test() tests bernstein_poly_01_matrix().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 January 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_poly_01_matrix_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_poly_01_matrix() is given M data values X,' )
  print ( '  and a degree N, and returns an Mx(N+1) matrix B such that' )
  print ( '  B(i,j) is the j-th Bernstein polynomial evaluated at the' )
  print ( '  i-th data value.' )

  m = 5
  x = np.linspace ( 0.0, 1.0, m )
  n = 1
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(5,1+1):' )

  m = 5
  x = np.linspace ( 0.0, 1.0, m )
  n = 4
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(5,4+1):' )

  m = 10
  x = np.linspace ( 0.0, 1.0, m )
  n = 4
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(10,4+1):' )

  m = 3
  x = np.linspace ( 0.0, 1.0, m )
  n = 5
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(3,5+1):' )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_poly_01_matrix_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_01 ( n, x ):

#*****************************************************************************80
#
## bernstein_poly_01() evaluates the Bernstein polynomials defined on [0,1].
#
#  Discussion:
#
#    The Bernstein polynomials are assumed to be based on [0,1].
#
#  Formula:
#
#    B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
#
#  First values:
#
#    B(0,0)(X) = 1
#
#    B(1,0)(X) =      1-X
#    B(1,1)(X) =                X
#
#    B(2,0)(X) =     (1-X)^2
#    B(2,1)(X) = 2 * (1-X)   * X
#    B(2,2)(X) =                X^2
#
#    B(3,0)(X) =     (1-X)^3
#    B(3,1)(X) = 3 * (1-X)^2 * X
#    B(3,2)(X) = 3 * (1-X)   * X^2
#    B(3,3)(X) =               X^3
#
#    B(4,0)(X) =     (1-X)^4
#    B(4,1)(X) = 4 * (1-X)^3 * X
#    B(4,2)(X) = 6 * (1-X)^2 * X^2
#    B(4,3)(X) = 4 * (1-X)   * X^3
#    B(4,4)(X) =               X^4
#
#  Special values:
#
#    B(N,I)(X) has a unique maximum value at X = I/N.
#
#    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
#
#    B(N,I)(1/2) = C(N,K) / 2^N
#
#    For a fixed X and N, the polynomials add up to 1:
#
#      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the degree of the Bernstein polynomials to be
#    used.  For any N, there is a set of N+1 Bernstein polynomials,
#    each of degree N, which form a basis for polynomials on [0,1].
#
#    real X, the evaluation point.
#
#  Output:
#
#    real B(1:N+1), the values of the N+1 Bernstein polynomials at X.
#
  import numpy as np

  b = np.zeros ( n + 1 )

  if ( n == 0 ):
 
    b[0] = 1.0
 
  elif ( 0 < n ):
 
    b[0] = 1.0 - x
    b[1] = x
 
    for i in range ( 2, n + 1 ):
      b[i] = x * b[i-1]
      for j in range ( i - 1, 0, -1 ):
        b[j] = x * b[j-1] + ( 1.0 - x ) * b[j]
      b[0] = ( 1.0 - x ) * b[0]

  return b

def bernstein_poly_01_test ( ):

#*****************************************************************************80
#
## bernstein_poly_01_test() tests bernstein_poly_01().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'bernstein_poly_01_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_poly_01() evaluates Bernstein polynomials.' )
  print ( '' )
  print ( '       N       K             X                 F                          F' )
  print ( '                                               tabulated                  computed' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, n, k, x, f1 = bernstein_poly_01_values ( n_data )

    if ( n_data == 0 ):
      break

    f = bernstein_poly_01 ( n, x )
    f2 = f[k]

    print ( '  %6d  %6d  %12f  %24.16g  %24.16g' % ( n, k, x, f1, f2 ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_poly_01_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_01_test2 ( ):

#*****************************************************************************80
#
## bernstein_poly_01_test2() tests the Partition-of-Unity property.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_poly_01_test2():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_poly_01() evaluates the Bernstein polynomials' )
  print ( '  based on the interval [0,1].' )
  print ( '' )
  print ( '  Here we test the partition of unity property.' )
  print ( '' )
  print ( '     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)' )
  print ( '' )

  for n in range ( 0, 11 ):

    x = np.random.rand ( )

    bvec = bernstein_poly_01 ( n, x )

    print ( '  %4d  %7.4f  %14.6g' % ( n, x, np.sum ( bvec ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_poly_01_test2():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_01_values ( n_data ):

#*****************************************************************************80
#
## bernstein_poly_01_values() returns some values of the Bernstein polynomials.
#
#  Discussion:
#
#    The Bernstein polynomials are assumed to be based on [0,1].
#
#    The formula for the Bernstein polynomials is
#
#      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
#
#    In Mathematica, the function can be evaluated by:
#
#      Binomial[n,i] * (1-x)^(n-i) * x^i
#
#  First values:
#
#    B(0,0)(X) = 1
#
#    B(1,0)(X) =      1-X
#    B(1,1)(X) =               X
#
#    B(2,0)(X) =     (1-X)^2
#    B(2,1)(X) = 2 * (1-X)   * X
#    B(2,2)(X) =               X^2
#
#    B(3,0)(X) =     (1-X)^3
#    B(3,1)(X) = 3 * (1-X)^2 * X
#    B(3,2)(X) = 3 * (1-X)   * X^2
#    B(3,3)(X) =               X^3
#
#    B(4,0)(X) =     (1-X)^4
#    B(4,1)(X) = 4 * (1-X)^3 * X
#    B(4,2)(X) = 6 * (1-X)^2 * X^2
#    B(4,3)(X) = 4 * (1-X)   * X^3
#    B(4,4)(X) =               X^4
#
#  Special values:
#
#    B(N,I)(X) has a unique maximum value at X = I/N.
#
#    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
#
#    B(N,I)(1/2) = C(N,K) / 2^N
#
#    For a fixed X and N, the polynomials add up to 1:
#
#      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    26 December 2014
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Stephen Wolfram,
#    The Mathematica Book,
#    Fourth Edition,
#    Wolfram Media / Cambridge University Press, 1999.
#
#  Input:
#
#    integer N_DATA.  The user sets N_DATA to 0 before the first call.
#
#  Output:
#
#    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
#    returns the corresponding data; when there is no more data, the
#    output value of N_DATA will be 0 again.
#
#    integer N, the degree of the polynomial.
#
#    integer K, the index of the polynomial.
#
#    real X, the argument of the polynomial.
#
#    real F, the value of the polynomial B(N,K)(X).
#
  import numpy as np

  n_max = 15

  f_vec = np.array ( ( \
     0.1000000000000000E+01, \
     0.7500000000000000E+00, \
     0.2500000000000000E+00, \
     0.5625000000000000E+00, \
     0.3750000000000000E+00, \
     0.6250000000000000E-01, \
     0.4218750000000000E+00, \
     0.4218750000000000E+00, \
     0.1406250000000000E+00, \
     0.1562500000000000E-01, \
     0.3164062500000000E+00, \
     0.4218750000000000E+00, \
     0.2109375000000000E+00, \
     0.4687500000000000E-01, \
     0.3906250000000000E-02 ) )

  k_vec = np.array ( ( \
    0, \
    0, 1, \
    0, 1, 2, \
    0, 1, 2, 3, \
    0, 1, 2, 3, 4 ))

  n_vec = np.array ( ( \
    0, \
    1, 1, \
    2, 2, 2, \
    3, 3, 3, 3, \
    4, 4, 4, 4, 4 ))

  x_vec = np.array ( ( \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00, \
     0.25E+00  ) )

  if ( n_data < 0 ):
    n_data = 0

  if ( n_max <= n_data ):
    n_data = 0
    n = 0
    k = 0
    x = 0.0
    f = 0.0
  else:
    n = n_vec[n_data]
    k = k_vec[n_data]
    x = x_vec[n_data]
    f = f_vec[n_data]
    n_data = n_data + 1

  return n_data, n, k, x, f

def bernstein_poly_01_values_test ( ):

#*****************************************************************************80
#
## bernstein_poly_01_values_test() tests bernstein_poly_01_values().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    26 December 2014
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'bernstein_poly_01_values_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_poly_01_values() stores values of Bernstein polynomials.' )
  print ( '' )
  print ( '      N       K            X            F' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, n, k, x, f = bernstein_poly_01_values ( n_data )

    if ( n_data == 0 ):
      break

    print ( '  %6d  %6d  %12f  %24.16g' % ( n, k, x, f ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_poly_01_values_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_ab_approx ( n, a, b, ydata, nval, xval ):

#*****************************************************************************80
#
## bernstein_poly_ab_approx(): Bernstein polynomial approximant to F(X) on [A,B].
#
#  Formula:
#
#    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
#
#    where
#
#      X(I) = ( ( N - I ) * A + I * B ) / N
#      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
#
#  Discussion:
#
#    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
#    not an interpolant; in other words, its value is not guaranteed to equal
#    that of F at any particular point.  However, for a fixed interval
#    [A,B], if we let N increase, the Bernstein polynomial converges
#    uniformly to F everywhere in [A,B], provided only that F is continuous.
#    Even if F is not continuous, but is bounded, the polynomial converges
#    pointwise to F(X) at all points of continuity.  On the other hand,
#    the convergence is quite slow compared to other interpolation
#    and approximation schemes.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    David Kahaner, Cleve Moler, Steven Nash,
#    Numerical Methods and Software,
#    Prentice Hall, 1989,
#    ISBN: 0-13-627258-4,
#    LC: TA345.K34.
#
#  Input:
#
#    integer N, the degree of the Bernstein polynomial
#    to be used.  N must be at least 0.
#
#    real A, B, the endpoints of the interval on which the
#    approximant is based.  A and B should not be equal.
#
#    real YDATA(N+1), the data values at N+1 equally
#    spaced points in [A,B].  If N = 0, then the evaluation point should
#    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
#    ( (N-I)*A + I*B ) / N ).
#
#    integer NVAL, the number of points at which the
#    approximant is to be evaluated.
#
#    real XVAL(NVAL), the point at which the Bernstein 
#    polynomial approximant is to be evaluated.  The entries of XVAL do not 
#    have to lie in the interval [A,B].
#
#  Output:
#
#    real YVAL(NVAL), the values of the Bernstein 
#    polynomial approximant for F, based in [A,B], evaluated at XVAL.
#
  import numpy as np

  yval = np.zeros ( nval )

  for i in range ( 0, nval ):
#
#  Evaluate the Bernstein basis polynomials at XVAL.
#
    bvec = bernstein_poly_ab ( n, a, b, xval[i] )
#
#  Now compute the sum of YDATA(I) * BVEC(I).
#
    yval[i] = np.dot ( ydata, bvec )

  return yval

def bernstein_poly_ab_approx_test ( ):

#*****************************************************************************80
#
## bernstein_poly_ab_approx_test() tests bernstein_poly_ab_approx().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_poly_ab_approx_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_poly_ab_approx() evaluates the Bernstein polynomial' )
  print ( '  approximant to a function F(X) defined over [A,B].' )

  a = 1.0
  b = 3.0

  print ( '' )
  print ( '     N      Max Error' )
  print ( '' )

  for degree in range ( 0, 21 ):
#
#  Generate data values.
#
    xdata = np.zeros ( degree + 1 )
    ydata = np.zeros ( degree + 1 )

    for i in range ( 0, degree + 1 ):

      if ( degree == 0 ):
        xdata[i] = 0.5 * ( a + b );
      else:
        xdata[i] = ( float ( degree - i ) * a   \
                   + float (          i ) * b ) \
                   / float ( degree     )

      ydata[i] = np.sin ( xdata[i] )
#
#  Compare the true function and the approximant.
#
    nval = 501

    xval = np.linspace ( a, b, nval )

    yval = bernstein_poly_ab_approx ( degree, a, b, ydata, nval, xval )

    error_max = max ( abs ( yval - np.sin ( xval ) ) )

    print ( '  %4d  %14.6g' % ( degree, error_max ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_poly_ab_approx_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_ab ( n, a, b, x ):

#*****************************************************************************80
#
## bernstein_poly_ab() evaluates at X the Bernstein polynomials based in [A,B].
#
#  Formula:
#
#    BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
#
#  First values:
#
#    B(0,0)(X) =   1
#
#    B(1,0)(X) = (      B-X                ) / (B-A)
#    B(1,1)(X) = (                 X-A     ) / (B-A)
#
#    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
#    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
#    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
#
#    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
#    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
#    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
#    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
#
#    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
#    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
#    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4 
#    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4 
#    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4 
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the degree of the Bernstein polynomials to be used.
#    For any N, there is a set of N+1 Bernstein polynomials, each of
#    degree N, which form a basis for polynomials on [A,B].
#
#    real A, B, the endpoints of the interval on which the
#    polynomials are to be based.  A and B should not be equal.
#
#    real X, the point at which the polynomials are to be evaluated.
#
#  Output:
#
#    real P(N+1), the values of the N+1 Bernstein polynomials at X.
#
  import numpy as np

  if ( b == a ):
    print ( '' )
    print ( 'bernstein_poly_ab - Fatal error!' )
    print ( '  A = B = %g' % ( a ) )
    raise Exception ( 'bernstein_poly_ab - Fatal error!' )

  p = np.zeros ( n + 1 )

  if ( n == 0 ):
 
    p[0] = 1.0
 
  elif ( 0 < n ):
 
    p[0] = ( b - x ) / ( b - a )
    p[1] = ( x - a ) / ( b - a )
 
    for i in range ( 2, n + 1 ):
      p[i] = ( x - a ) * p[i-1] / ( b - a )
      for j in range ( i - 1, 0, -1 ):
        p[j] = ( ( b - x ) * p[j] + ( x - a ) * p[j-1] ) / ( b - a )
      p[0] = ( b - x ) * p[0] / ( b - a )
 
  return p

def bernstein_poly_ab_test ( ):

#*****************************************************************************80
#
## bernstein_poly_ab_test() tests bernstein_poly_ab().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  n = 10

  print ( '' )
  print ( 'bernstein_poly_ab_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_poly_ab() evaluates Bernstein polynomials over an' )
  print ( '  arbitrary interval [A,B].' )
  print ( '' )
  print ( '  Here, we demonstrate that ' )
  print ( '    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)' )
  print ( '  provided only that' )
  print ( '    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).' )

  x = 0.3
  a = 0.0
  b = 1.0
  p = bernstein_poly_ab ( n, a, b, x )
 
  print ( '' )
  print ( '     N     K     A        B        X       BPAB(N,K,A,B)(X)' )
  print ( '' )
  for k in range ( 0, n + 1 ):
    print ( '  %4d  %4d  %7.4f  %7.4f  %7.4f  %14.6g' % ( n, k, a, b, x, p[k] ) )
 
  x = 1.3
  a = 1.0
  b = 2.0
  p = bernstein_poly_ab ( n, a, b, x )
 
  print ( '' )
  print ( '     N     K     A        B        X       BPAB(N,K,A,B)(X)' )
  print ( '' )
  for k in range ( 0, n + 1 ):
    print ( '  %4d  %4d  %7.4f  %7.4f  %7.4f  %14.6g' % ( n, k, a, b, x, p[k] ) )

  x = 2.6
  a = 2.0
  b = 4.0
  p = bernstein_poly_ab ( n, a, b, x )
 
  print ( '' )
  print ( '     N     K     A        B        X       BPAB(N,K,A,B)(X)' )
  print ( '' )
  for k in range ( 0, n + 1 ):
    print ( '  %4d  %4d  %7.4f  %7.4f  %7.4f  %14.6g' % ( n, k, a, b, x, p[k] ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_poly_ab_test():' )
  print ( '  Normal end of execution.' )
  return

def bernstein_polynomial_test ( ):

#*****************************************************************************80
#
## bernstein_polynomial_test() tests bernstein_polynomial().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    16 March 2016
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'bernstein_polynomial_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test bernstein_polynomial().' )

  bernstein_matrix_test ( )
  bernstein_matrix_test2 ( )
  bernstein_matrix_determinant_test ( )
  bernstein_matrix_inverse_test ( )

  bernstein_poly_01_test ( )
  bernstein_poly_01_test2 ( )
  bernstein_poly_01_matrix_test ( )
  bernstein_poly_01_values_test ( )

  bernstein_poly_ab_test ( )
  bernstein_poly_ab_approx_test ( )

  bernstein_to_legendre_test ( )
  bernstein_to_power_test ( )

  bernstein_vandermonde_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_polynomial_test():' )
  print ( '  Normal end of execution.' )
  print ( '' )
  return

def bernstein_to_legendre ( n ):

#*****************************************************************************80
#
## bernstein_to_legendre() returns the Bernstein-to-Legendre matrix.
#
#  Discussion:
#
#    The Legendre polynomials are often defined on [-1,+1], while the
#    Bernstein polynomials are defined on [0,1].  For this function,
#    the Legendre polynomials have been shifted to share the [0,1]
#    interval of definition.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the maximum degree of the polynomials.
#
#  Output:
#
#    real A(N+1,N+1), the Bernstein-to-Legendre matrix.
#
  from scipy.special import comb
  import numpy as np

  a = np.zeros ( [ n + 1, n + 1 ] )

  for i in range ( 0, n + 1 ):
    for j in range ( 0, n + 1 ):
      for k in range ( 0, i + 1 ):
        a[i,j] = a[i,j] \
          + r8_mop ( i + k ) * ( comb ( i, k ) ) ** 2 / comb ( n + i, j + k )
      a[i,j] = a[i,j] * comb ( n, j ) \
        * ( 2 * i + 1 ) / ( n + i + 1 )

  return a

def bernstein_to_legendre_test ( ):

#*****************************************************************************80
#
## bernstein_to_legendre_test() tests bernstein_to_legendre().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_to_legendre_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_to_legendre() returns the matrix A which maps' )
  print ( '  polynomial coefficients from Bernstein to Legendre form.' )

  n = 5
  a = bernstein_to_legendre ( n )
  r8mat_print ( n + 1, n + 1, a, '  A = bernstein_to_legendre(5):' )

  b = legendre_to_bernstein ( n )
  r8mat_print ( n + 1, n + 1, b, '  B = legendre_to_bernstein(5):' )

  c = np.dot ( a, b )
  e = r8mat_is_identity ( n + 1, c )
  print ( '' )
  print ( '  ||A*B-I|| = %g' % ( e ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_to_legendre_test():' )
  print ( '  Normal end of execution.' )
  return

def legendre_to_bernstein ( n ):

#*****************************************************************************80
#
## legendre_to_bernstein() returns the Legendre-to-Bernstein matrix.
#
#  Discussion:
#
#    The Legendre polynomials are often defined on [-1,+1], while the
#    Bernstein polynomials are defined on [0,1].  For this function,
#    the Legendre polynomials have been shifted to share the [0,1]
#    interval of definition.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the maximum degree of the polynomials.
#
#  Output:
#
#    real A(N+1,N+1), the Legendre-to-Bernstein matrix.
#
  from scipy.special import comb
  import numpy as np

  a = np.zeros ( [ n + 1, n + 1 ] )

  for i in range ( 0, n + 1 ):
    for j in range ( 0, n + 1 ):
      for k in range ( max ( 0, i + j - n ), min ( i, j ) + 1 ):
        a[i,j] = a[i,j] \
          + r8_mop ( j + k ) * ( comb ( j, k ) ) ** 2 * comb ( n - j, i - k )
      a[i,j] = a[i,j] / comb ( n, i )

  return a

def bernstein_to_power ( n ):

#*****************************************************************************80
#
## bernstein_to_power() returns the Bernstein-to-Power matrix.
#
#  Discussion:
#
#    The Bernstein-to-Power matrix of degree N is an N+1xN+1 matrix A which can 
#    be used to transform the N+1 coefficients of a polynomial of degree N
#    from a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
#    to a vector P of coefficients of the power basis (1,x,x^2,...,x^n).
#
#    If we are using N=4-th degree polynomials, the matrix has the form:
#
#      1   0   0   0  0
#     -4   4   0   0  0
#      6 -12   6   0  0
#     -4  12 -12   4  0
#      1  -4   6  -4  1
#
#   and a polynomial with the Bernstein basis representation
#     p(x) = 3/4 * b(4,1) + 1/2 b(4,2)
#   whose Bernstein coefficient vector is
#     B = ( 0, 3/4, 1/2, 0, 0 )
#   will have the Bernstein basis coefficients 
#     P = A * B = ( 0, 3, -6, 3, 0 ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    16 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the degree of the polynomials.
#
#  Output:
#
#    real A(N+1,N+1), the Bernstein-to-Power matrix.
#
  from scipy.special import comb
  import numpy as np

  a = np.zeros ( [ n + 1, n + 1 ] )

  for j in range ( 0, n + 1 ):
    for i in range ( 0, j + 1 ):
      a[n-i,n-j] = r8_mop ( j - i ) * comb ( n - i, j - i ) * comb ( n, i )

  return a

def bernstein_to_power_test ( ):

#*****************************************************************************80
#
## bernstein_to_power_test() tests bernstein_to_power().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    16 March 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'bernstein_to_power_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_to_power() returns the matrix A which maps' )
  print ( '  polynomial coefficients from Bernstein to Power form.' )

  n = 5
  a = bernstein_to_power ( n )
  r8mat_print ( n + 1, n + 1, a, '  A = bernstein_to_power(5):' )

  b = power_to_bernstein ( n )
  r8mat_print ( n + 1, n + 1, b, '  B = power_to_bernstein(5):' )

  c = np.dot ( a, b )
  e = r8mat_is_identity ( n + 1, c )
  print ( '' )
  print ( '  ||A*B-I|| = %g' % ( e ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_to_power_test():' )
  print ( '  Normal end of execution.' )
  return

def power_to_bernstein ( n ):

#*****************************************************************************80
#
## power_to_bernstein() returns the Power-to-Bernstein matrix.
#
#  Discussion:
#
#    The Power-to-Bernstein matrix of degree N is an N+1xN+1 matrix A which can 
#    be used to transform the N+1 coefficients of a polynomial of degree N
#    from a vector P of coefficients of the power basis (1,x,x^2,...,x^n)
#    to a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
#
#    If we are using N=4-th degree polynomials, the matrix has the form:
#
#          1   0    0    0   0
#          1  1/4   0    0   0
#      A = 1  1/2  1/6   0   0
#          1  3/4  1/2  1/4  1
#          1   1    1    1   1
#
#   and a polynomial 
#     p(x) = 3x - 6x^2 + 3x^3
#   whose power coefficient vector is
#     P = ( 0, 3, -6, 3, 0 )
#   will have the Bernstein basis coefficients 
#     B = A * P = ( 0, 3/4, 1/2, 0, 0 ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    16 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the degree of the polynomials.
#
#  Output:
#
#    real A[0:N,0:N], the Power-to-Bernstein matrix.
#
  from scipy.special import comb
  import numpy as np

  a = np.zeros ( [ n + 1, n + 1 ] )

  for j in range ( 0, n + 1 ):
    for i in range ( 0, j + 1 ):
      a[n-i,n-j] = comb ( j, i ) / comb ( n, i )

  return a

def bernstein_vandermonde ( n ):

#*****************************************************************************80
#
## bernstein_vandermonde() returns the Bernstein Vandermonde matrix.
#
#  Discussion:
#
#    The Bernstein Vandermonde matrix of order N is constructed by
#    evaluating the N Bernstein polynomials of degree N-1 at N equally
#    spaced points between 0 and 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#  Output:
#
#    real A(N,N), the Bernstein Vandermonde matrix.
#
  import numpy as np

  v = np.zeros ( [ n, n ] )

  if ( n == 1 ):
    v[0,0] = 1.0
    return v

  for i in range ( 0, n ):
    x = float ( i ) / float ( n - 1 )
    b = bernstein_poly_01 ( n - 1, x )
    for j in range ( 0, n ):
      v[i,j] = b[j]

  return v

def bernstein_vandermonde_test ( ):

#*****************************************************************************80
#
## bernstein_vandermonde_test() tests bernstein_vandermonde().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'bernstein_vandermonde_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bernstein_vandermonde() returns an NxN matrix whose (I,J) entry' )
  print ( '  is the value of the J-th Bernstein polynomial of degree N-1' )
  print ( '  evaluated at the I-th equally spaced point in [0,1].' )

  n = 8
  a = bernstein_vandermonde ( n )
  r8mat_print ( n, n, a, '  Bernstein vandermonde ( 8 ):' )
#
#  Terminate.
#
  print ( '' )
  print ( 'bernstein_vandermonde_test():' )
  print ( '  Normal end of execution.' )
  return

def r8mat_is_identity ( n, a ):

#*****************************************************************************80
#
## r8mat_is_identity() determines if a matrix is the identity.
#
#  Discussion:
#
#    The routine returns the Frobenius norm of A - I.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A(N,N), the matrix.
#
#  Output:
#
#    real ERROR_FROBENIUS, the Frobenius norm
#    of the difference matrix A - I, which would be exactly zero
#    if A were the identity matrix.
#
  import numpy as np

  error_frobenius = 0.0

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      if ( i == j ):
        error_frobenius = error_frobenius + ( a[i,j] - 1.0 ) ** 2
      else:
        error_frobenius = error_frobenius + a[i,j] ** 2

  error_frobenius = np.sqrt ( error_frobenius );

  return error_frobenius

def r8mat_is_identity_test ( ):

#*****************************************************************************80
#
## r8mat_is_identity_test() tests r8mat_is_identity().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'r8mat_is_identity_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8mat_is_identity() reports the Frobenius norm difference' )
  print ( '  between a given matrix A and the identity matrix.' )

  n = 4
  a = np.zeros ( [ n, n ] )
  r8mat_print ( n, n, a, '  Zero matrix:' )
  e = r8mat_is_identity ( n, a )
  print ( '' )
  print ( '  Difference is %g' % ( e ) )

  for i in range ( 0, n ):
    a[i,i] = 1.0
  r8mat_print ( n, n, a, '  Identity matrix:' )
  e = r8mat_is_identity ( n, a )
  print ( '' )
  print ( '  Difference is %g' % ( e ) )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      a[i,j] = a[i,j] + float ( i * j ) / 1000
  r8mat_print ( n, n, a, '  Almost identity matrix:' )
  e = r8mat_is_identity ( n, a )
  print ( '' )
  print ( '  Difference is %g' % ( e ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8mat_is_identity_test():' )
  print ( '  Normal end of execution.' )
  return

def r8mat_norm_fro ( m, n, a ):

#*****************************************************************************80
#
## r8mat_norm_fro() returns the Frobenius norm of an R8MAT.
#
#  Discussion:
#
#    The Frobenius norm is defined as
#
#      value = sqrt ( sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
#
#    The matrix Frobenius norm is not derived from a vector norm, but
#    is compatible with the vector L2 norm, so that:
#
#      vec_norm_l2 ( A * x ) <= mat_norm_fro ( A ) * vec_norm_l2 ( x ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 December 2014
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the number of rows in A.
#
#    integer N, the number of columns in A.
#
#    real A(M,N), the matrix whose Frobenius
#    norm is desired.
#
#  Output:
#
#    real VALUE, the Frobenius norm of A.
#
  import numpy as np
 
  value = np.sqrt ( sum ( sum ( a ** 2 ) ) )

  return value

def r8mat_norm_fro_test ( ):

#*****************************************************************************80
#
## r8mat_norm_fro_test() tests r8mat_norm_fro().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 December 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  m = 5
  n = 4

  a = np.zeros ( ( m, n ) )

  k = 0
  t1 = 0.0

  for i in range ( 0, m ):
    for j in range ( 0, n ):
      k = k + 1
      a[i,j] = k
      t1 = t1 + k * k

  t1 = np.sqrt ( t1 )

  print ( '' )
  print ( 'r8mat_norm_fro_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8mat_norm_fro() computes the Frobenius norm of an R8MAT;' )

  t2 = r8mat_norm_fro ( m, n, a )

  r8mat_print ( m, n, a, '  A:' )
  print ( '' )
  print ( '  Expected Frobenius norm = %g' % ( t1 ) )
  print ( '  Computed Frobenius norm = %g' % ( t2 ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8mat_norm_fro_test():' )
  print ( '  Normal end of execution.' )
  return

def r8mat_print ( m, n, a, title ):

#*****************************************************************************80
#
## r8mat_print() prints an R8MAT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 August 2014
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the number of rows in A.
#
#    integer N, the number of columns in A.
#
#    real A(M,N), the matrix.
#
#    string TITLE, a title.
#
  r8mat_print_some ( m, n, a, 0, 0, m - 1, n - 1, title )

  return

def r8mat_print_test ( ):

#*****************************************************************************80
#
## r8mat_print_test() tests r8mat_print().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    10 February 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'r8mat_print_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8mat_print() prints an R8MAT.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ], 
    [ 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 ], 
    [ 31.0, 32.0, 33.0, 34.0, 35.0, 36.0 ], 
    [ 41.0, 42.0, 43.0, 44.0, 45.0, 46.0 ] ], dtype = np.float64 )
  r8mat_print ( m, n, v, '  Here is an R8MAT:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8mat_print_test:' )
  print ( '  Normal end of execution.' )
  return

def r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title ):

#*****************************************************************************80
#
## r8mat_print_some() prints out a portion of an R8MAT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    10 February 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, N, the number of rows and columns of the matrix.
#
#    real A(M,N), an M by N matrix to be printed.
#
#    integer ILO, JLO, the first row and column to print.
#
#    integer IHI, JHI, the last row and column to print.
#
#    string TITLE, a title.
#
  incx = 5

  print ( '' )
  print ( title )

  if ( m <= 0 or n <= 0 ):
    print ( '' )
    print ( '  (None)' )
    return

  for j2lo in range ( max ( jlo, 0 ), min ( jhi + 1, n ), incx ):

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )
    
    print ( '' )
    print ( '  Col: ', end = '' )

    for j in range ( j2lo, j2hi + 1 ):
      print ( '%7d       ' % ( j ), end = '' )

    print ( '' )
    print ( '  Row' )

    i2lo = max ( ilo, 0 )
    i2hi = min ( ihi, m )

    for i in range ( i2lo, i2hi + 1 ):

      print ( '%7d :' % ( i ), end = '' )
      
      for j in range ( j2lo, j2hi + 1 ):
        print ( '%12g  ' % ( a[i,j] ), end = '' )

      print ( '' )

  return

def r8mat_print_some_test ( ):

#*****************************************************************************80
#
## r8mat_print_some_test() tests r8mat_print_some().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 October 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'r8mat_print_some_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8mat_print_some() prints some of an R8MAT.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ], 
    [ 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 ], 
    [ 31.0, 32.0, 33.0, 34.0, 35.0, 36.0 ], 
    [ 41.0, 42.0, 43.0, 44.0, 45.0, 46.0 ] ], dtype = np.float64 )
  r8mat_print_some ( m, n, v, 0, 3, 2, 5, '  Here is an R8MAT:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8mat_print_some_test():' )
  print ( '  Normal end of execution.' )
  return

def r8_mop ( i ):

#*****************************************************************************80
#
## r8_mop() returns the I-th power of -1 as an R8 value.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    01 June 2013
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer I, the power of -1.
#
#  Output:
#
#    real r8_mop, the I-th power of -1.
#
  if ( ( i % 2 ) == 0 ):
    value = + 1.0
  else:
    value = - 1.0

  return value

def timestamp ( ):

#*****************************************************************************80
#
## timestamp() prints the date as a timestamp.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    06 April 2013
#
#  Author:
#
#    John Burkardt
#
  import time

  t = time.time ( )
  print ( time.ctime ( t ) )

  return None

if ( __name__ == '__main__' ):
  timestamp ( )
  bernstein_polynomial_test ( )
  timestamp ( )
