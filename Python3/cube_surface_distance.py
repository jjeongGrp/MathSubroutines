#! /usr/bin/env python3
#
def cube_surface_distance_histogram ( n ):

#*****************************************************************************80
#
## cube_surface_distance_histogram() histograms cube surface distance statistics.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 March 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of samples to use.
#
  import matplotlib.pyplot as plt
  import numpy as np

  p1 = cube_surface_sample ( n )
  p2 = cube_surface_sample ( n )
  t = np.zeros ( n )
  for i in range ( 0, n ):
    t[i] = np.linalg.norm ( p1[i] - p2[i] )

  bins = 40
  plt.hist ( t, bins = bins, rwidth = 0.95, \
    range = np.array ( [ 0.0, np.sqrt ( 3.0 ) ] ), density = True )
  plt.grid ( True )
  plt.xlabel ( '<-- Distance -->' )
  plt.ylabel ( '<-- Frequency -->' )
  plt.title ( 'Distance between random points on unit cube surface' )

  return

def cube_surface_distance_stats ( n ):

#*****************************************************************************80
#
## cube_surface_distance_stats() estimates cube surface distance statistics.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 March 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of sample points to use.
#
#  Output:
#
#    real DMU, DVAR, the estimated mean and variance of the
#    distance between two random points on the surface of the unit cube.
#
  import numpy as np

  p1 = cube_surface_sample ( n )
  p2 = cube_surface_sample ( n )
  t = np.zeros ( n )
  for i in range ( 0, n ):
    t[i] = np.linalg.norm ( p1[i] - p2[i] )

  dmu = np.mean ( t )
  dvar = np.var ( t )

  return dmu, dvar

def cube_surface_distance_test ( ):

#*****************************************************************************80
#
## cube_surface_distance_test() tests cube_surface_distance().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 March 2022
#
#  Author:
#
#    John Burkardt
#
  import matplotlib.pyplot as plt
  import numpy as np
  import platform

  print ( '' )
  print ( 'cube_surface_distance_test():' )
  print ( '  Python version: ' + platform.python_version ( ) )
  print ( '  Test cube_surface_distance()' )
  print ( '' )
  print ( '  The mean distance between two random points on the surface' )
  print ( '  of a cube can be computed as:' )
  print ( '    mu = ( 5 mu_diff + mu_same ) / 6' )
  print ( '  where:' )
  print ( '    mu_diff = mean distance, points are on different faces' )
  print ( '    mu_same = mean distance, points are on same face.' )

  n = 10000
  dmu, dvar = cube_surface_distance_stats ( n )
  dmu_diff_exact = 0.9263900551740467
  dmu_same_exact = ( 2.0 + np.sqrt ( 2.0 ) \
    + 5.0 * np.log ( 1.0 + np.sqrt ( 2.0 ) ) ) / 15.0
  dmu_exact = ( 5.0 * dmu_diff_exact + dmu_same_exact ) / 6.0
  print ( '' )
  print ( '  Using N =', n, 'sample points,' )
  print ( '  Exact mean distance (diff) =', dmu_diff_exact )
  print ( '  Exact mean distance (same) =', dmu_same_exact )
  print ( '  Exact mean distance        =', dmu_exact )
  print ( '  Estimated mean distance    =', dmu )
  print ( '  Estimated variance         =', dvar )

  n = 10000
  cube_surface_distance_histogram ( n )
  filename = 'cube_surface_distance_histogram.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "' + filename + '"' )
  plt.show ( block = False )
  plt.close ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'cube_surface_distance_test():' )
  print ( '  Normal end of execution.' )

  return

def cube_surface_sample ( n ):

#*****************************************************************************80
#
## cube_surface_sample() randomly selects points on the surface of a cube.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 March 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of points to sample.
#
#  Output:
#
#    real P(N,3), N points selected uniformly at random from
#    the surface of the unit cube.
#
  import numpy as np

  p = np.random.rand ( n, 3 )

  for i in range ( 0, n ):
    j = np.random.random_integers ( 0, 2 )
    s = np.random.random_integers ( 0, 1 )
    p[i,j] = float ( s )

  return p

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
  cube_surface_distance_test ( )
  timestamp ( )
