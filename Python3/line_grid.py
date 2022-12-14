#! /usr/bin/env python3
#
def line_grid ( n, a, b, c ):

#*****************************************************************************80
#
## line_grid(): grid points over the interior of a line segment in 1D.
#
#  Discussion:
#
#    In 1D, a grid is created using N points.
#
#    Over the interval [A,B], we have 5 choices for grid centering:
#      1: 0,   1/3, 2/3, 1
#      2: 1/5, 2/5, 3/5, 4/5
#      3: 0,   1/4, 2/4, 3/4
#      4: 1/4, 2/4, 3/4, 1
#      5: 1/8, 3/8, 5/8, 7/8
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    24 April 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of points.
#
#    real A, B, the endpoints for each dimension.
#
#    integer C, the grid centering for each dimension.
#    1 <= C <= 5.
#
#  Output:
#
#    real X(N), the points.
#
  import numpy as np
#
#  Create the 1D grids in each dimension.
#
  x = np.zeros ( n )

  for j in range ( 1, n + 1 ):
    jm1 = j - 1

    if ( c == 1 ):

      if ( n == 1 ):
        x[jm1] = 0.5 * ( a + b )
      else:
        x[jm1] = (   float ( n - j     ) * a   \
                   + float (     j - 1 ) * b ) \
                   / float ( n     - 1 )
    elif ( c == 2 ):
      x[jm1] = (   float ( n - j + 1 ) * a   \
                 + float (     j     ) * b ) \
                 / float ( n     + 1 )
    elif ( c == 3 ):
      x[jm1] = (   float ( n - j + 1 ) * a   \
                 + float (     j - 1 ) * b ) \
                 / float ( n         )
    elif ( c == 4 ):
      x[jm1] = (   float ( n - j ) * a   \
                 + float (     j ) * b ) \
                 / float ( n     )
    elif ( c == 5 ):
      x[jm1] = (   float ( 2 * n - 2 * j + 1 ) * a   \
                 + float (         2 * j - 1 ) * b ) \
                 / float ( 2 * n             )

  return x

def line_grid_test01 ( ):

#*****************************************************************************80
#
## line_grid_test01() uses simple parameters.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    24 April 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  n = 11
  a = -1.0
  b = +1.0
  c = 1

  print ( '' )
  print ( 'line_grid_test01' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Create a grid using line_grid.' )
  print ( '  Use simple parameters.' )
  print ( '  Number of grid points N = %d' % ( n ) )
  print ( '' )
  print ( '     N     C      A         B' )
  print ( '' )
  print ( '  %4d  %4d  %8.4f  %8.4f' % ( n, c, a, b ) )

  x = line_grid ( n, a, b, c )

  r8vec_print ( n, x, '  Grid points:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'line_grid_test01:' )
  print ( '  Normal end of execution.' )
  return

def line_grid_test02 ( ):

#*****************************************************************************80
#
## line_grid_test02() tries an increasing number of points.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    24 April 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  a = 0.0
  b = 1.0
  c = 2

  print ( '' )
  print ( 'line_grid_test02' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Create a grid using line_grid.' )
  print ( '  Try an increasing number of points.' )

  n = 4

  for test in range ( 0, 3 ):

    n = 2 * n + 1
    print ( '' )
    print ( '     N     C      A         B' )
    print ( '' )
    print ( '  %4d  %4d  %8.4f  %8.4f' % ( n, c, a, b ) )

    x = line_grid ( n, a, b, c )

    r8vec_print ( n, x, '  Grid points:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'line_grid_test02:' )
  print ( '  Normal end of execution.' )
  return

def line_grid_test03 ( ):

#*****************************************************************************80
#
## line_grid_test03() tries all the centering options.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    24 April 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  n = 11
  a =    0.0
  b = +100.0

  print ( '' )
  print ( 'line_grid_test03' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Try the different centering options.' )
  print ( '  Number of grid points N = %d' % ( n ) )

  for c in range ( 1, 6 ):

    print ( '' )
    print ( '     N     C      A         B' )
    print ( '' )
    print ( '  %4d  %4d  %8.4f  %8.4f' % ( n, c, a, b ) )

    x = line_grid ( n, a, b, c )

    r8vec_print ( n, x, '  Grid points:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'line_grid_test03:' )
  print ( '  Normal end of execution.' )
  return

def line_grid_test ( ):

#*****************************************************************************80
#
## line_grid_test() tests line_grid().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    24 April 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'line_grid_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test line_grid()' )
#
#  Library.
#
  line_grid_test01 ( )
  line_grid_test02 ( )
  line_grid_test03 ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'line_grid_test:' )
  print ( '  Normal end of execution.' )
  return

def r8vec_print ( n, a, title ):

#*****************************************************************************80
#
## r8vec_print() prints an R8VEC.
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
#    integer N, the dimension of the vector.
#
#    real A(N), the vector to be printed.
#
#    string TITLE, a title.
#
  print ( '' )
  print ( title )
  print ( '' )
  for i in range ( 0, n ):
    print ( '%6d:  %12g' % ( i, a[i] ) )

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
  line_grid_test ( )
  timestamp ( )