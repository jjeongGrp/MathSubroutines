#! /usr/bin/env python3
#
def ball01_monomial_integral ( e ):

#*****************************************************************************80
#
## ball01_monomial_integral() returns monomial integrals in the unit ball.
#
#  Discussion:
#
#    The integration region is 
#
#      X^2 + Y^2 + Z^2 <= 1.
#
#    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 September 2016
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Gerald Folland,
#    How to Integrate a Polynomial Over a Sphere,
#    American Mathematical Monthly,
#    Volume 108, Number 5, May 2001, pages 446-448.
#
#  Input:
#
#    integer E(3), the exponents of X, Y and Z in the 
#    monomial.  Each exponent must be nonnegative.
#
#  Output:
#
#    real INTEGRAL, the integral.
#
  import numpy as np
  import scipy.special

  if ( any ( e < 0 ) ):
    print ( '' )
    print ( 'ball01_monomial_integral - Fatal error!' )
    print ( '  All exponents must be nonnegative.' )
    raise Exception ( 'ball01_monomial_integral - Fatal error!' )
#
#  Integrate over the surface.
#
  if ( all ( e == 0 ) ):

    integral = 2.0 * np.sqrt ( np.pi ** 3 ) / scipy.special.gamma ( 1.5 )

  elif ( any ( ( e % 2 ) == 1 ) ):

    integral = 0.0

  else:

    integral = 2.0

    for i in range ( 0, 3 ):
      integral = integral * scipy.special.gamma ( 0.5 * ( e[i] + 1 ) )

    integral = integral / scipy.special.gamma ( 0.5 * ( np.sum ( e + 1 ) ) )
#
#  The surface integral is now adjusted to give the volume integral.
#
  r = 1.0
  s = np.sum ( e ) + 3

  integral = integral * r ** s / float ( s )

  return integral

def ball01_sample ( n ):

#*****************************************************************************80
#
## ball01_sample() uniformly samples the unit ball.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 September 2016
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Russell Cheng,
#    Random Variate Generation,
#    in Handbook of Simulation,
#    edited by Jerry Banks,
#    Wiley, 1998, pages 168.
#
#    Reuven Rubinstein,
#    Monte Carlo Optimization, Simulation, and Sensitivity 
#    of Queueing Networks,
#    Krieger, 1992,
#    ISBN: 0894647644,
#    LC: QA298.R79.
#
#  Input:
#
#    integer N, the number of points.
#
#  Output:
#
#    real X(3,N), the points.
#
  import numpy as np

  x = np.random.normal ( 0.0, 1.0, [ 3, n ] )
  for j in range ( 0, n ):
    norm = np.sqrt ( x[0,j] ** 2 + x[1,j] ** 2 + x[2,j] ** 2 )
    for i in range ( 0, 3 ):
      x[i,j] = x[i,j] / norm

  for j in range ( 0, n ):
    r = np.random.random ( 1 )
    x[:,j] = x[:,j] * r ** ( 1.0 / 3.0 )

  return x

def ball01_volume ( ):

#*****************************************************************************80
#
## ball01_volume() returns the volume of the unit ball.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 September 2016
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VOLUME, the volume of the unit ball.
#
  import numpy as np

  r = 1.0
  volume = 4.0 * np.pi * r ** 3 / 3.0

  return volume

def ball_monte_carlo_test ( ):

#*****************************************************************************80
#
## ball_monte_carlo_test() uses ball01_sample() with an increasing number of points.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 September 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  e_test = np.array ( [ \
    [ 0, 0, 0 ], \
    [ 2, 0, 0 ], \
    [ 0, 2, 0 ], \
    [ 0, 0, 2 ], \
    [ 4, 0, 0 ], \
    [ 2, 2, 0 ], \
    [ 0, 0, 4 ] ] )

  print ( '' )
  print ( 'ball_monte_carlo_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Estimate integrals over the interior of the unit ball' )
  print ( '  using the Monte Carlo method.' )

  print ( '' )
  print ( '         N        1              X^2             Y^2             Z^2             X^4           X^2Y^2           Z^4' )
  print ( '' )

  n = 1

  while ( n <= 65536 ):

    x = ball01_sample ( n )

    print ( '  %8d' % ( n ), end = '' )
    for j in range ( 0, 7 ):
      e = e_test[j,:]
      value = monomial_value ( 3, n, e, x )
      result = ball01_volume ( ) * np.sum ( value ) / float ( n )
      print ( '  %14.6g' % ( result ), end = '' )
    print ( '' )

    n = 2 * n

  print ( '' )
  print ( '     Exact', end = '' )
  for j in range ( 0, 7 ):
    e = e_test[j,:]
    result = ball01_monomial_integral ( e )
    print ( '  %14.6g' % ( result ), end = '' )
  print ( '' )

  return

def monomial_value ( m, n, e, x ):

#*****************************************************************************80
#
## monomial_value() evaluates a monomial.
#
#  Discussion:
#
#    This routine evaluates a monomial of the form
#
#      product ( 1 <= i <= m ) x(i)^e(i)
#
#    The combination 0.0^0, if encountered, is treated as 1.0.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    07 April 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the spatial dimension.
#
#    integer N, the number of evaluation points.
#
#    integer E(M), the exponents.
#
#    real X(M,N), the point coordinates.
#
#  Output:
#
#    real V(N), the monomial values.
#
  import numpy as np

  v = np.ones ( n )

  for i in range ( 0, m ):
    if ( 0 != e[i] ):
      for j in range ( 0, n ):
        v[j] = v[j] * x[i,j] ** e[i]

  return v

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

  return

if ( __name__ == '__main__' ):
  timestamp ( )
  ball_monte_carlo_test ( )
  timestamp ( )
