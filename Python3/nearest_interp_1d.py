#! /usr/bin/env python3
#
def nearest_interp_1d ( nd, xd, yd, ni, xi ):

#*****************************************************************************80
#
## nearest_interp_1d() evaluates the nearest neighbor interpolant.
#
#  Discussion:
#
#    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
#    constant function which interpolates the data (XD(I),YD(I)) for I = 1
#    to ND.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer ND, the number of data points.
#    ND must be at least 1.
#
#    real XD(ND), the data points.
#
#    real YD(ND), the data values.
#
#    integer NI, the number of interpolation points.
#
#    real XI(NI), the interpolation points.
#
#  Output:
#
#    real YI(NI), the interpolated values.
#
  import numpy as np

  yi = np.zeros ( ni )

  for i in range ( 0, ni ):

    d = float ( 'Inf' )
    k = -1

    for k2 in range ( 0, nd ):

      d2 = abs ( xi[i] - xd[k2] )

      if ( d2 < d ):
        k = k2
        d = d2

    yi[i] = yd[k]

  return yi

def nearest_interp_1d_test01 ( prob, ni ):

#*****************************************************************************80
#
## nearest_interp_1d_test01() tests nearest_interp_1d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer PROB, the index of the problem.
#
#    integer NI, the number of interpolation points.
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'nearest_interp_1d_test01' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Sample the nearest neighbor interpolant for problem #%d' % ( prob ) )

  nd = p00_data_num ( prob )
  d = p00_data ( prob, 2, nd )
  xd = d[0,:]
  yd = d[1,:]
  xd_min = min ( xd )
  xd_max = max ( xd )

  xi = np.linspace ( xd_min, xd_max, ni )
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi )

  title = 'X, Y for problem ' + str ( prob )
  r8vec2_print ( ni, xi, yi, title )
#
#  Terminate.
#
  print ( '' )
  print ( 'nearest_interp_1d_test01:' )
  print ( '  Normal end of execution.' )
  return

def nearest_interp_1d_test02 ( prob ):

#*****************************************************************************80
#
## nearest_interp_1d_test02() tests nearest_interp_1d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer PROB, the index of the problem.
#
  import matplotlib.pyplot as plt
  import numpy as np
  import platform

  print ( '' )
  print ( 'nearest_interp_1d_test02:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Interpolate data from TEST_INTERP problem #%d.' % ( prob ) )

  nd = p00_data_num ( prob )
  print ( '  Number of data points = %d' % ( nd ) )

  xy = p00_data ( prob, 2, nd )
  
  r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

  xd = xy[0,:]
  yd = xy[1,:]
#
#  #1:  Does interpolant match function at interpolation points?
#
  ni = nd
  xi = xd
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi )

  int_error = np.linalg.norm ( yi - yd ) / float ( ni )

  print ( '' )
  print ( '  L2 interpolation error averaged per interpolant node = %g' % ( int_error ) )
#
#  #3: Plot the data.
#
  plt.plot ( xd, yd, 'b-', linewidth = 3.0 )
  plt.plot ( xd, yd, 'ro', markersize = 10 )
  t = 'p0' + str ( prob ) + ' Piecewise Linear Interpolant'
  plt.title ( t )
  plt.grid ( True )
  plt.xlabel ( '<---X--->' )
  plt.ylabel ( '<---Y--->' )
  filename = 'p0' + str ( prob ) + '_pwl.png'
  plt.savefig ( filename )
  plt.show ( block = False )
  plt.close ( )

  print ( '' )
  print ( '  Created plot file "%s"' % ( filename ) )
#
#  #4: Plot the piecewise linear and nearest neighbor interpolants.
#
  ni = 501
  xmin = min ( xd )
  xmax = max ( xd )
  xi = np.linspace ( xmin, xmax, ni )
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi )

  plt.plot ( xi, yi, 'b-', linewidth = 3.0 )
  plt.plot ( xd, yd, 'ro', markersize = 10 )
  t = 'p0' + str ( prob ) + ' Nearest Neighbor Interpolant'
  plt.title ( t )
  plt.grid ( True )
  plt.xlabel ( '<---X--->' )
  plt.ylabel ( '<---Y--->' )
  filename = 'p0' + str ( prob ) + '_near.png'
  plt.savefig ( filename )
  plt.show ( block = False )
  plt.close ( )

  print ( '  Created plot file "%s".' % ( filename ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'nearest_interp_1d_test02:' )
  print ( '  Normal end of execution.' )
  return

def nearest_interp_1d_test ( ):

#*****************************************************************************80
#
## nearest_interp_1d_test() tests nearest_interp_1d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'nearest_interp_1d_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test the nearest_interp_1d library.' )

  prob_num = p00_prob_num ( )

  for prob in range ( 1, prob_num + 1 ):
    nearest_interp_1d_test01 ( prob, 11 )

  for prob in range ( 1, prob_num + 1 ):
    nearest_interp_1d_test02 ( prob )
#
#  Terminate.
#
  print ( '' )
  print ( 'nearest_interp_1d_test:' )
  print ( '  Normal end of execution.' )
  return

def p00_data_num ( prob ):

#*****************************************************************************80
#
## p00_data_num() returns the number of data points for any problem.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer PROB, the problem index.
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  if ( prob == 1 ):
    data_num = p01_data_num ( )
  elif ( prob == 2 ):
    data_num = p02_data_num ( )
  elif ( prob == 3 ):
    data_num = p03_data_num ( )
  elif ( prob == 4 ):
    data_num = p04_data_num ( )
  elif ( prob == 5 ):
    data_num = p05_data_num ( )
  elif ( prob == 6 ):
    data_num = p06_data_num ( )
  elif ( prob == 7 ):
    data_num = p07_data_num ( )
  elif ( prob == 8 ):
    data_num = p08_data_num ( )
  else:
    print ( '' )
    print ( 'P00_data_num - Fatal error!' )
    print ( '  Unexpected input value of PROB.' )
    raise Exception ( 'P00_data_num - Fatal error!' )

  return data_num

def p01_data_num ( ):

#*****************************************************************************80
#
## p01_data_num() returns the number of data points for problem p01.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 18

  return data_num

def p02_data_num ( ):

#*****************************************************************************80
#
## p02_data_num() returns the number of data points for problem p02.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 18

  return data_num

def p03_data_num ( ):

#*****************************************************************************80
#
## p03_data_num() returns the number of data points for problem p03.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 11

  return data_num

def p04_data_num ( ):

#*****************************************************************************80
#
## p04_data_num() returns the number of data points for problem p04.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 8

  return data_num

def p05_data_num ( ):

#*****************************************************************************80
#
## p05_data_num() returns the number of data points for problem p05.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 9

  return data_num

def p06_data_num ( ):

#*****************************************************************************80
#
## p06_data_num() returns the number of data points for problem p06.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 49

  return data_num

def p07_data_num ( ):

#*****************************************************************************80
#
## p07_data_num() returns the number of data points for problem p07.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 4

  return data_num

def p08_data_num ( ):

#*****************************************************************************80
#
## p08_data_num() returns the number of data points for problem p08.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer DATA_NUM, the number of data points.
#
  data_num = 12

  return data_num

def p00_data_num_test ( ):

#*****************************************************************************80
#
## p00_data_num_test() tests p00_data_num.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'p00_data_num_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  p00_data_num returns the number of data points for any problem.' )
  print ( '' )
  print ( '  Problem   Data Num' )
  print ( '' )

  prob_num = p00_prob_num ( )

  for prob in range ( 1, prob_num + 1 ):

    data_num = p00_data_num ( prob )

    print ( '  %7d  %9d' % ( prob, data_num ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'p00_data_num_test:' )
  print ( '  Normal end of execution.' )
  return

def p00_data ( prob, dim_num, data_num ):

#*****************************************************************************80
#
## p00_data() returns the data for any problem.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer PROB, the problem index.
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  if ( prob == 1 ):
    p_data = p01_data ( dim_num, data_num )
  elif ( prob == 2 ):
    p_data = p02_data ( dim_num, data_num )
  elif ( prob == 3 ):
    p_data = p03_data ( dim_num, data_num )
  elif ( prob == 4 ):
    p_data = p04_data ( dim_num, data_num )
  elif ( prob == 5 ):
    p_data = p05_data ( dim_num, data_num )
  elif ( prob == 6 ):
    p_data = p06_data ( dim_num, data_num )
  elif ( prob == 7 ):
    p_data = p07_data ( dim_num, data_num )
  elif ( prob == 8 ):
    p_data = p08_data ( dim_num, data_num )
  else:
    print ( '' )
    print ( 'P00_data - Fatal error!' )
    print ( '  Unexpected input value of PROB.' )
    raise Exception ( 'P00_data - Fatal error!' )

  return p_data

def p01_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p01_data() returns the data for problem p01.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
    [  0.0, 4.0 ], \
    [  1.0, 5.0 ], \
    [  2.0, 6.0 ], \
    [  4.0, 6.0 ], \
    [  5.0, 5.0 ], \
    [  6.0, 3.0 ], \
    [  7.0, 1.0 ], \
    [  8.0, 1.0 ], \
    [  9.0, 1.0 ], \
    [ 10.0, 3.0 ], \
    [ 11.0, 4.0 ], \
    [ 12.0, 4.0 ], \
    [ 13.0, 3.0 ], \
    [ 14.0, 3.0 ], \
    [ 15.0, 4.0 ], \
    [ 16.0, 4.0 ], \
    [ 17.0, 3.0 ], \
    [ 18.0, 0.0 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p02_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p02_data() returns the data for problem p02.
#
#  Discussion:
#
#    Two pairs of identical X values have now been slightly separated.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
   [  0.00,   0.00 ], \
   [  1.34,   5.00 ], \
   [  5.00,   8.66 ], \
   [ 10.00,  10.00 ], \
   [ 10.60,  10.40 ], \
   [ 10.70,  12.00 ], \
   [ 10.705, 28.60 ], \
   [ 10.80,  30.20 ], \
   [ 11.40,  30.60 ], \
   [ 19.60,  30.60 ], \
   [ 20.20,  30.20 ], \
   [ 20.295, 28.60 ], \
   [ 20.30,  12.00 ], \
   [ 20.40,  10.40 ], \
   [ 21.00,  10.00 ], \
   [ 26.00,   8.66 ], \
   [ 29.66,   5.00 ], \
   [ 31.00,   0.00 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p03_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p03_data() returns the data for problem p03.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
    [  0.0,  0.0 ], \
    [  2.0, 10.0 ], \
    [  3.0, 10.0 ], \
    [  5.0, 10.0 ], \
    [  6.0, 10.0 ], \
    [  8.0, 10.0 ], \
    [  9.0, 10.5 ], \
    [ 11.0, 15.0 ], \
    [ 12.0, 50.0 ], \
    [ 14.0, 60.0 ], \
    [ 15.0, 85.0 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p04_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p04_data() returns the data for problem p04.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
   [ 0.00,  0.00 ], \
   [ 0.05,  0.70 ], \
   [ 0.10,  1.00 ], \
   [ 0.20,  1.00 ], \
   [ 0.80,  0.30 ], \
   [ 0.85,  0.05 ], \
   [ 0.90,  0.10 ], \
   [ 1.00,  1.00 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p05_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p05_data() returns the data for problem p05.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real p_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
    [ 0.00, 0.00 ], \
    [ 0.10, 0.90 ], \
    [ 0.20, 0.95 ], \
    [ 0.30, 0.90 ], \
    [ 0.40, 0.10 ], \
    [ 0.50, 0.05 ], \
    [ 0.60, 0.05 ], \
    [ 0.80, 0.20 ], \
    [ 1.00, 1.00 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p06_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p06_data() returns the data for problem p06.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
   [ 595.0, 0.644 ], \
   [ 605.0, 0.622 ], \
   [ 615.0, 0.638 ], \
   [ 625.0, 0.649 ], \
   [ 635.0, 0.652 ], \
   [ 645.0, 0.639 ], \
   [ 655.0, 0.646 ], \
   [ 665.0, 0.657 ], \
   [ 675.0, 0.652 ], \
   [ 685.0, 0.655 ], \
   [ 695.0, 0.644 ], \
   [ 705.0, 0.663 ], \
   [ 715.0, 0.663 ], \
   [ 725.0, 0.668 ], \
   [ 735.0, 0.676 ], \
   [ 745.0, 0.676 ], \
   [ 755.0, 0.686 ], \
   [ 765.0, 0.679 ], \
   [ 775.0, 0.678 ], \
   [ 785.0, 0.683 ], \
   [ 795.0, 0.694 ], \
   [ 805.0, 0.699 ], \
   [ 815.0, 0.710 ], \
   [ 825.0, 0.730 ], \
   [ 835.0, 0.763 ], \
   [ 845.0, 0.812 ], \
   [ 855.0, 0.907 ], \
   [ 865.0, 1.044 ], \
   [ 875.0, 1.336 ], \
   [ 885.0, 1.881 ], \
   [ 895.0, 2.169 ], \
   [ 905.0, 2.075 ], \
   [ 915.0, 1.598 ], \
   [ 925.0, 1.211 ], \
   [ 935.0, 0.916 ], \
   [ 945.0, 0.746 ], \
   [ 955.0, 0.672 ], \
   [ 965.0, 0.627 ], \
   [ 975.0, 0.615 ], \
   [ 985.0, 0.607 ], \
   [ 995.0, 0.606 ], \
   [1005.0, 0.609 ], \
   [1015.0, 0.603 ], \
   [1025.0, 0.601 ], \
   [1035.0, 0.603 ], \
   [1045.0, 0.601 ], \
   [1055.0, 0.611 ], \
   [1065.0, 0.601 ], \
   [1075.0, 0.608 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p07_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p07_data() returns the data for problem p07.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent
#    variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
   [ 0.0, 1.0 ], \
   [ 1.0, 2.0 ], \
   [ 4.0, 2.0 ], \
   [ 5.0, 1.0 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p08_data ( dim_num, data_num ):

#*****************************************************************************80
#
## p08_data() returns the data for problem p08.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer DIM_NUM, the spatial dimension of the dependent variables.
#
#    integer DATA_NUM, the number of data points.
#
#  Output:
#
#    real P_data(DIM_NUM,DATA_NUM), the data.
#
  import numpy as np

  p_data = np.array ( [ \
   [ -1.0, 1.00 ], \
   [ -0.8, 0.64 ], \
   [ -0.6, 0.36 ], \
   [ -0.4, 0.16 ], \
   [ -0.2, 0.04 ], \
   [  0.0, 0.00 ], \
   [  0.2, 0.04 ], \
   [  0.20001, 0.05 ], \
   [  0.4, 0.16 ], \
   [  0.6, 0.36 ], \
   [  0.8, 0.64 ], \
   [  1.0, 1.00 ] ] )

  p_data = np.transpose ( p_data )

  return p_data

def p00_data_test ( ):

#*****************************************************************************80
#
## p00_data_test tests p00_data.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 June 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'p00_data_test tests p00_data' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  p00_data returns the actual (MxN) data for any problem.' )

  prob_num = p00_prob_num ( )

  for prob in range ( 1, prob_num + 1 ):

    print ( '' )
    print ( '  Problem %d' % ( prob ) )

    data_num = p00_data_num ( prob )
    print ( '  DATA_NUM = %d' % ( data_num ) )

    dim_num = p00_dim_num ( prob )
    print ( '  DIM_NUM  = %d' % ( dim_num ) )

    p = p00_data ( prob, dim_num, data_num )

    r8mat_transpose_print ( dim_num, data_num, p, '  Data array:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'p00_data_test:' )
  print ( '  Normal end of execution.' )
  return

def p00_prob_num ( ):

#*****************************************************************************80
#
## p00_prob_num() returns the number of problems.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    integer VALUE, the number of problems.
#
  value = 8

  return value

def p00_prob_num_test ( ):

#*****************************************************************************80
#
## p00_prob_num_test() tests p00_prob_num().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 June 2015
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'p00_prob_num_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  p00_prob_num returns the number of test problems.' )

  num = p00_prob_num ( )
  print ( '' )
  print ( '  TEST_INTERP includes %d test problems.' % ( num ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'p00_prob_num_test:' )
  print ( '  Normal end of execution.' )
  return

def r8mat_transpose_print ( m, n, a, title ):

#*****************************************************************************80
#
## r8mat_transpose_print() prints an R8MAT, transposed.
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
  r8mat_transpose_print_some ( m, n, a, 0, 0, m - 1, n - 1, title )

  return

def r8mat_transpose_print_test ( ):

#*****************************************************************************80
#
## r8mat_transpose_print_test() tests r8mat_transpose_print().
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
  print ( 'r8mat_transpose_print_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8mat_transpose_print prints an R8MAT.' )

  m = 4
  n = 3
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0 ], 
    [ 21.0, 22.0, 23.0 ], 
    [ 31.0, 32.0, 33.0 ], 
    [ 41.0, 42.0, 43.0 ] ], dtype = np.float64 )
  r8mat_transpose_print ( m, n, v, '  Here is an R8MAT, transposed:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8mat_transpose_print_test:' )
  print ( '  Normal end of execution.' )
  return

def r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title ):

#*****************************************************************************80
#
## r8mat_transpose_print_some() prints a portion of an R8MAT, transposed.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 November 2014
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

  for i2lo in range ( max ( ilo, 0 ), min ( ihi, m - 1 ), incx ):

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m - 1 )
    i2hi = min ( i2hi, ihi )
    
    print ( '' )
    print ( '  Row: ', end = '' )

    for i in range ( i2lo, i2hi + 1 ):
      print ( '%7d       ' % ( i ), end = '' )

    print ( '' )
    print ( '  Col' )

    j2lo = max ( jlo, 0 )
    j2hi = min ( jhi, n - 1 )

    for j in range ( j2lo, j2hi + 1 ):

      print ( '%7d :' % ( j ), end = '' )
      
      for i in range ( i2lo, i2hi + 1 ):
        print ( '%12g  ' % ( a[i,j] ), end = '' )

      print ( '' )

  return

def r8mat_transpose_print_some_test ( ):

#*****************************************************************************80
#
## r8mat_transpose_print_some_test() tests r8mat_transpose_print_some().
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
  print ( 'r8mat_transpose_print_some_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8mat_transpose_print_some prints some of an R8MAT, transposed.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ], 
    [ 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 ], 
    [ 31.0, 32.0, 33.0, 34.0, 35.0, 36.0 ], 
    [ 41.0, 42.0, 43.0, 44.0, 45.0, 46.0 ] ], dtype = np.float64 )
  r8mat_transpose_print_some ( m, n, v, 0, 3, 2, 5, '  R8MAT, rows 0:2, cols 3:5:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8mat_transpose_print_some_test:' )
  print ( '  Normal end of execution.' )
  return

def r8vec2_print ( n, a1, a2, title ):

#*****************************************************************************80
#
## r8vec2_print() prints an R8VEC2.
#
#  Discussion:
#
#    An R8VEC2 is a dataset consisting of N pairs of real values, stored
#    as two separate vectors A1 and A2.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of components of the vector.
#
#    real A1(N), A2(N), the vectors to be printed.
#
#    string TITLE, a title.
#
  print ( '' )
  print ( title )
  print ( '' )
  for i in range ( 0, n ):
    print ( '  %6d:   %12g  %12g' % ( i, a1[i], a2[i] ) )

  return

def r8vec2_print_test ( ):

#*****************************************************************************80
#
## r8vec2_print_test() tests r8vec2_print().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 June 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  print ( '' )
  print ( 'r8vec2_print_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  r8vec2_print prints a pair of R8VEC\'s.' )

  n = 6
  v = np.array ( [ 0.0, 0.20, 0.40, 0.60, 0.80, 1.0 ], dtype = np.float64 )
  w = np.array ( [ 0.0, 0.04, 0.16, 0.36, 0.64, 1.0 ], dtype = np.float64 )
  r8vec2_print ( n, v, w, '  Print a pair of R8VEC\'s:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'r8vec2_print_test:' )
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

  return None

if ( __name__ == '__main__' ):
  timestamp ( )
  nearest_interp_1d_test ( )
  timestamp ( )
