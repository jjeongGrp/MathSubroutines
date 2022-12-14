#! /usr/bin/env python3
#
def newton_test ( ):

#*****************************************************************************80
#
## newton_test tests newton().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 March 2021
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'newton_test:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test newton()' )

  newton1_test ( )
  newton2_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'newton_test:' )
  print ( '  Normal end of execution.' )
  return
def newton1 ( f, x, dfdx, epsilon = 1.0E-7, n_max = 100 ):

#*****************************************************************************80
#
## newton1() implements Newton's method.
#
#  Discussion:
#
#    A solution of the nonlinear equation F(X)=0 is desired.
#
#    The user supplies functions to evaluate F and F', an initial guess for X,
#    an error tolerance, and a maximum number of steps.
#
#    The reference describes this as a "first draft" of an implementation of
#    Newton's method, and suggests some improvements.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 May 2015
#
#  Author:
#
#    Original Python version by Hans Petter Langtangen,
#    Modifications by John Burkardt.
#
#  Reference:
#
#    Hans Petter Langtangen,
#    A primer on scientific programming with Python,
#    Springer, 2012,
#    ISBN13: 978-3-642-30293-0.
#
#  Input:
#
#    function F(X), returns the value of the function at X.
#
#    real X, the initial estimate for the solution.
#
#    function dfdx(X), returns the value of the derivative at X.
#
#    real EPSILON, an absolute error tolerance for F(X).
#    If omitted, the default value is 1.0E-7.
#
#    integer N_MAX, the maximum number of iterations.
#    If omitted, the default value is 100.
#
#  Output:
#
#    real X, the estimate for the solution.
#
#    integer N, the number of iterations taken.
#
#    real FX, the value of the function at X.
#
  n = 0
  while epsilon < abs ( f ( x ) ) and n <= n_max:
    x = x - f ( x ) / dfdx ( x )
    n = n + 1

  return x, n, f ( x )

def newton1_test ( ):

#*****************************************************************************80
#
## newton1_test() tests newton1().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 May 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'NEWTON1_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  NEWTON1 is a simple implementation of Newton\'s method' )
  print ( '  for estimating a solution of F(X)=0.' )

  x0 = 1.5
  epsilon = 1.0E-7
  n_max = 100

  x, n, fx = newton1 ( f1, x0, dfdx1, epsilon, n_max )

  print ( '' )
  print ( '  F1 = sin ( x ) - 0.5 * x' )
  print ( '  X0 = %g, F(X0) = %g' % ( x0, f1 ( x0 ) ) )
  print ( '  X* = %g, F(X*) = %g' % ( x,  f1 ( x ) ) )
  print ( '  Number of steps was %d' % ( n ) )

  x0 = 0.1
  epsilon = 1.0E-7
  n_max = 100

  x, n, fx = newton1 ( f2, x0, dfdx2, epsilon, n_max )

  print ( '' )
  print ( '  F2 = 2.0 * X - EXP ( - X )' )
  print ( '  X0 = %g, F(X0) = %g' % ( x0, f2 ( x0 ) ) )
  print ( '  X* = %g, F(X*) = %g' % ( x,  f2 ( x ) ) )
  print ( '  Number of steps was %d' % ( n ) )

  x0 = 0.5
  epsilon = 1.0E-7
  n_max = 100

  x, n, fx = newton1 ( f3, x0, dfdx3, epsilon, n_max )

  print ( '' )
  print ( '  F3 = X * EXP ( - X )' )
  print ( '  X0 = %g, F(X0) = %g' % ( x0, f3 ( x0 ) ) )
  print ( '  X* = %g, F(X*) = %g' % ( x,  f3 ( x ) ) )
  print ( '  Number of steps was %d' % ( n ) )

  x0 = 0.03
  epsilon = 1.0E-7
  n_max = 100

  x, n, fx = newton1 ( f4, x0, dfdx4, epsilon, n_max )

  print ( '' )
  print ( '  F4 = EXP ( X ) - 1 / ( 100 * X^2 )' )
  print ( '  X0 = %g, F(X0) = %g' % ( x0, f4 ( x0 ) ) )
  print ( '  X* = %g, F(X*) = %g' % ( x,  f4 ( x ) ) )
  print ( '  Number of steps was %d' % ( n ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'NEWTON1_TEST:' )
  print ( '  Normal end of execution.' )
  return

def newton2 ( f, x, dfdx, epsilon = 1.0E-7, n_max = 100, store = False ):

#*****************************************************************************80
#
## newton2() implements Newton's method.
#
#  Discussion:
#
#    A solution of the nonlinear equation F(X)=0 is desired.
#
#    The user supplies functions to evaluate F and F', an initial guess for X,
#    an error tolerance, and a maximum number of steps.
#
#    This version includes some improvements and additions relative to
#    the simpler version in newton1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 May 2015
#
#  Author:
#
#    Original Python version by Hans Petter Langtangen,
#    Modifications by John Burkardt.
#
#  Reference:
#
#    Hans Petter Langtangen,
#    A primer on scientific programming with Python,
#    Springer, 2012,
#    ISBN13: 978-3-642-30293-0.
#
#  Input:
#
#    function F(X), returns the value of the function at X.
#
#    real X, the initial estimate for the solution.
#
#    function dfdx(X), returns the value of the derivative at X.
#
#    real EPSILON, an absolute error tolerance for F(X).
#    If omitted, the default value is 1.0E-7.
#
#    integer N_MAX, the maximum number of iterations.
#    If omitted, the default value is 100.
#
#    bool STORE, indicates whether the user this procedure to store
#    the successive iterates ( x, f(x) ) in an array to be returned upon
#    completion.
#
#  Output:
#
#    real X, the estimate for the solution.
#
#    integer N, if STORE is FALSE, the number of iterations taken.
#
#    real FX, if STORE is FALSE, the value of the function at X.
#
#    array INFO, if STORE is TRUE, contains the pairs of values ( x, f(x) ) 
#    of the iterates, including the starting point.
#
  n = 0
  f_value = f ( x ) 

  if store:
    info = [ ( x, f_value ) ]

  while epsilon < abs ( f_value ) and n <= n_max:

    dfdx_value = float ( dfdx ( x ) )

    if abs ( dfdx_value ) < 1.0E-14:
      raise ValueError ( "Newton2: f'(%g)=%g" % ( x, dfdx_value ) )

    x = x - f_value / dfdx_value
    n = n + 1

    f_value = f ( x )
    if store:
      info.append ( ( x, f_value ) )

  if store:
    return x, info
  else:
    return x, n, f_value

def newton2_test ( ):

#*****************************************************************************80
#
## newton2_test() tests newton2().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 May 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'NEWTON2_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  NEWTON2 is a more robust implementation of Newton\'s method' )
  print ( '  for estimating a solution of F(X)=0.' )
#
#  Example 1, don't store the data.
#
  x0 = 1.5
  epsilon = 1.0E-7
  n_max = 100
  store = False

  x, n, fx = newton2 ( f1, x0, dfdx1, epsilon, n_max, store )

  print ( '' )
  print ( '  F1 = sin ( x ) - 0.5 * x' )
  print ( '  X0 = %g, F(X0) = %g' % ( x0, f1 ( x0 ) ) )
  print ( '  X* = %g, F(X*) = %g' % ( x,  f1 ( x ) ) )
  print ( '  Number of steps was %d' % ( n ) )
#
#  Example 2, store the data, and print it on return.
#
  x0 = 0.1
  epsilon = 1.0E-7
  n_max = 100
  store = True

  x, info = newton2 ( f2, x0, dfdx2, epsilon, n_max, store )

  print ( '' )
  print ( '  F2 = 2.0 * X - EXP ( - X )' )
  print ( '  Root X = %g' % ( x ) )
  for i in range ( len ( info ) ):
    print ( '  %3d: f(%g) = %g' % ( i, info[i][0], info[i][1] ) )
#
#  Example 3, store the data, and print it on return.
#
  x0 = 0.5
  epsilon = 1.0E-7
  n_max = 100
  store = True

  x, info = newton2 ( f3, x0, dfdx3, epsilon, n_max, store )

  print ( '' )
  print ( '  F3 = X * EXP ( - X )' )
  print ( '  Root X = %g' % ( x ) )
  for i in range ( len ( info ) ):
    print ( '  %3d: f(%g) = %g' % ( i, info[i][0], info[i][1] ) )
#
#  Example 4, store the data, and print it on return.
#
  x0 = 0.03
  epsilon = 1.0E-7
  n_max = 100
  store = True

  x, info = newton2 ( f4, x0, dfdx4, epsilon, n_max, store )

  print ( '' )
  print ( '  F4 = EXP ( X ) - 1 / ( 100 * X^2 )' )
  print ( '  Root X = %g' % ( x ) )
  for i in range ( len ( info ) ):
    print ( '  %3d: f(%g) = %g' % ( i, info[i][0], info[i][1] ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'NEWTON2_TEST:' )
  print ( '  Normal end of execution.' )
  return

def f1 ( x ):

#*****************************************************************************80
#
## f1() defines the function for example 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = np.sin ( x ) - 0.5 * x

  return value

def dfdx1 ( x ):

#*****************************************************************************80
#
## dfdx1() defines the derivative for example 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = np.cos ( x ) - 0.5

  return value

def f2 ( x ):

#*****************************************************************************80
#
## f2() defines the function for example 2.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = 2.0 * x - np.exp ( - x )

  return value

def dfdx2 ( x ):

#*****************************************************************************80
#
## dfdx2() defines the derivative for example 2.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = 2.0 + np.exp ( - x )

  return value

def f3 ( x ):

#*****************************************************************************80
#
## f3() defines the function for example 3.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = x * np.exp ( - x )

  return value

def dfdx3 ( x ):

#*****************************************************************************80
#
## dfdx3() defines the derivative for example 3.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = ( 1.0 - x ) * np.exp ( - x )

  return value

def f4 ( x ):

#*****************************************************************************80
#
## f4() defines the function for example 4.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = np.exp ( x ) - 1.0 / ( 100.0 * x ** 2 )

  return value

def dfdx4 ( x ):

#*****************************************************************************80
#
## dfdx4() defines the derivative for example 4.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    05 October 2021
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the evaluation point.
#
#  Output:
#
#    real VALUE, the function value.
#
  import numpy as np

  value = np.exp ( x ) + 2.0 / ( 100.0 * x ** 3 )

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
  newton_test ( )
  timestamp ( )
