#! /usr/bin/env python3
#
def besselj_test ( ):

#*****************************************************************************80
#
## besselj_test() tests the besselj().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    13 January 2016
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'besselj_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test besselj().' )

  bessel_jx_values_test ( )
  scipy_jv_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'besselj_test():' )
  print ( '  Normal end of execution.' )
  return

def bessel_jx_values ( n_data ):

#*****************************************************************************80
#
## bessel_jx_values() returns some values of the Jx Bessel function.
#
#  Discussion:
#
#    This set of data considers the less common case in which the
#    index of the Bessel function Jn is actually not an integer.
#    We may suggest this case by occasionally replacing the symbol
#    "Jn" by "Jx".
#
#    In Mathematica, the function can be evaluated by:
#
#      besselj[n,x]
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    11 January 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Milton Abramowitz, Irene Stegun,
#    Handbook of Mathematical Functions,
#    National Bureau of Standards, 1964,
#    ISBN: 0-486-61272-4,
#    LC: QA47.A34.
#
#    Stephen Wolfram,
#    The Mathematica Book,
#    Fourth Edition,
#    Cambridge University Press, 1999,
#    ISBN: 0-521-64314-7,
#    LC: QA76.95.W65.
#
#  Input:
#
#    integer N_DATA, the user sets N_DATA to 0 before the first call.
#
#  Output:
#
#    integer N_DATA.  On each call, the routine increments N_DATA by 1, and
#    returns the corresponding data; when there is no more data, the
#    output value of N_DATA will be 0 again.
#
#    real NU, the order of the function.
#
#    real X, the argument of the function.
#
#    real FX, the value of the function.
#
  import numpy as np

  n_max = 28

  fx_vec = np.array ( ( \
        0.3544507442114011E+00, \
        0.6713967071418031E+00, \
        0.5130161365618278E+00, \
        0.3020049060623657E+00, \
        0.06500818287737578E+00, \
       -0.3421679847981618E+00, \
       -0.1372637357550505E+00, \
        0.1628807638550299E+00, \
        0.2402978391234270E+00, \
        0.4912937786871623E+00, \
       -0.1696513061447408E+00, \
        0.1979824927558931E+00, \
       -0.1094768729883180E+00, \
        0.04949681022847794E+00, \
        0.2239245314689158E+00, \
        0.2403772011113174E+00, \
        0.1966584835818184E+00, \
        0.02303721950962553E+00, \
        0.3314145508558904E+00, \
        0.5461734240402840E+00, \
       -0.2616584152094124E+00, \
        0.1296035513791289E+00, \
       -0.1117432171933552E+00, \
        0.03142623570527935E+00, \
        0.1717922192746527E+00, \
        0.3126634069544786E+00, \
        0.1340289119304364E+00, \
        0.06235967135106445E+00 ) )

  nu_vec = np.array ( ( \
    0.50E+00, \
    0.50E+00, \
    0.50E+00, \
    0.50E+00, \
    0.50E+00, \
    0.50E+00, \
    0.50E+00, \
    0.50E+00, \
    1.50E+00, \
    1.50E+00, \
    1.50E+00, \
    1.50E+00, \
    1.50E+00, \
    2.50E+00, \
    2.50E+00, \
    2.50E+00, \
    2.50E+00, \
    2.50E+00, \
    1.25E+00, \
    1.25E+00, \
    1.25E+00, \
    1.25E+00, \
    1.25E+00, \
    2.75E+00, \
    2.75E+00, \
    2.75E+00, \
    2.75E+00, \
    2.75E+00 ))

  x_vec = np.array ( ( \
      0.2E+00, \
      1.0E+00, \
      2.0E+00, \
      2.5E+00, \
      3.0E+00, \
      5.0E+00, \
     10.0E+00, \
     20.0E+00, \
      1.0E+00, \
      2.0E+00, \
      5.0E+00, \
     10.0E+00, \
     50.0E+00, \
      1.0E+00, \
      2.0E+00, \
      5.0E+00, \
     10.0E+00, \
     50.0E+00, \
      1.0E+00, \
      2.0E+00, \
      5.0E+00, \
     10.0E+00, \
     50.0E+00, \
      1.0E+00, \
      2.0E+00, \
      5.0E+00, \
     10.0E+00, \
     50.0E+00  ) )

  if ( n_data < 0 ):
    n_data = 0

  if ( n_max <= n_data ):
    n_data = 0
    nu = 0
    x = 0.0
    fx = 0.0
  else:
    nu = nu_vec[n_data]
    x = x_vec[n_data]
    fx = fx_vec[n_data]
    n_data = n_data + 1

  return n_data, nu, x, fx

def bessel_jx_values_test ( ):

#*****************************************************************************80
#
## bessel_jx_values_test() tests bessel_jx_values().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    11 January 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'bessel_jx_values_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  bessel_jx_values() stores values of the Bessel J function. of real order NU.' )
  print ( '' )
  print ( '      NU          X           J(NU,X)' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, nu, x, fx = bessel_jx_values ( n_data )

    if ( n_data == 0 ):
      break

    print ( '  %12f  %12f  %24.16g' % ( nu, x, fx ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'bessel_jx_values_test():' )
  print ( '  Normal end of execution.' )
  return

def scipy_jv_test ( ):

#*****************************************************************************80
#
## scipy_jv_test() checks some Bessel J functions computed by scipy jv()
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 January 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  from scipy.special import jv

  print ( '' )
  print ( 'scipy_jv_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  scipy.special.jv(nu,x) evaluates the Bessel J function. of real order NU.' )
  print ( '' )
  print ( '      NU          X           jv(nu,x)        Tabulated' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, nu, x, fx2 = bessel_jx_values ( n_data )

    if ( n_data == 0 ):
      break

    fx1 = jv ( nu, x )

    print ( '  %12f  %12f  %24.16g  %24.16g' % ( nu, x, fx1, fx2 ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'scipy_jv_test():' )
  print ( '  Normal end of execution.' )
  return

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
  besselj_test ( )
  timestamp ( )
