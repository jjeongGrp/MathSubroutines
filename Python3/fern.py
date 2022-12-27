#! /usr/bin/env python3
#
def fern ( n = 10000 ):

#*****************************************************************************80
#
## fern() displays the Barnsley fractal fern.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 June 2022
#
#  Author:
#
#    This program was modified by John Burkardt from an original
#    program by Cleve Moler.
#
#  Reference:
#
#    Michael Barnsley,
#    Fractals Everywhere,
#    Academic Press, 1988,
#    ISBN: 0120790696,
#    LC: QA614.86.B37.
#
#    Cleve Moler,
#    Experiments with MATLAB,
#    ebook: https://www.mathworks.com/moler/exm/index.html
#
#  Input:
#
#    integer N, the number of points to display.
#    A value of 10,000 or so is enough to see the fern.
#    For values of 500 or less, larger dots are displayed to suggest
#    how the plot is drawn.
#
  import matplotlib.pyplot as plt
  import numpy as np

  prob = np.array ( [ 0.85, 0.92, 0.99, 1.00 ] )
#
#  Compute the points.
#
  p = np.zeros ( [ n, 2 ] )

  p[0,0:2] = np.random.rand ( 2 )

  for i in range ( 1, n ):

    r = np.random.rand ( 1 )

    if ( r < prob[0] ):
      p[i,0] =   0.85 * p[i-1,0] + 0.04 * p[i-1,1] + 0.0
      p[i,1] = - 0.04 * p[i-1,0] + 0.85 * p[i-1,1] + 1.6
    elif ( r < prob[1] ):
      p[i,0] =   0.20 * p[i-1,0] - 0.26 * p[i-1,1] + 0.0
      p[i,1] =   0.23 * p[i-1,0] + 0.22 * p[i-1,1] + 1.6
    elif ( r < prob[2] ):
      p[i,0] = - 0.15 * p[i-1,0] + 0.28 * p[i-1,1] + 0.0
      p[i,1] =   0.26 * p[i-1,0] + 0.24 * p[i-1,1] + 0.44
    else:
      p[i,0] =   0.00 * p[i-1,0] + 0.00 * p[i-1,1] + 0.0
      p[i,1] =   0.00 * p[i-1,0] + 0.16 * p[i-1,1] + 0.0
#
#  Plot the points as small green dots.
#
  if ( n <= 500 ):
    dot_size = 5
  else:
    dot_size = 1

  plt.clf ( )
  plt.plot ( p[:,0], p[:,1], 'g.', markersize = dot_size )
  plt.axis ( 'equal' )
  plt.title ( 'Fractal Fern, N = ' + str ( n ) )
  filename = 'fern_' + str ( n ) + '.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "' + filename + '"' )
  plt.show ( block = False )
  plt.close ( )

  return

def fern_test ( ):

#*****************************************************************************80
#
## fern_test() tests fern().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 June 2022
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'fern_test():' )
  print ( '  Python version: ' + platform.python_version ( ) )
  print ( '  fern() makes an n-point image of the Barnsley fractal fern.' )

  n = 10000
  fern ( n )
#
#  Terminate.
#
  print ( '' )
  print ( 'fern_test():' )
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
#    21 August 2019
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
  fern_test ( )
  timestamp ( )
