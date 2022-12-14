#! /usr/bin/env python3
#
def angle_degree ( x1, y1, x2, y2, x3, y3 ):

#*****************************************************************************80
#
## angle_degree() returns the degree angle defined by three points.
#
#  Discussion:
#
#        P1
#        /
#       /
#      /
#     /
#    P2--------->P3
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X1, Y1, X2, Y2, X3, Y3, the coordinates of the points
#    P1, P2, P3.
#
#  Output:
#
#    real VALUE, the angle swept out by the rays, measured
#    in degrees.  0 <= VALUE < 360.  If either ray has zero length,
#    then VALUE is set to 0.
#
  import numpy as np

  x = ( x3 - x2 ) * ( x1 - x2 ) + ( y3 - y2 ) * ( y1 - y2 )

  y = ( x3 - x2 ) * ( y1 - y2 ) - ( y3 - y2 ) * ( x1 - x2 )

  if ( x == 0.0 and y == 0.0 ):
    value = 0.0
    return value

  value = np.arctan2 ( y, x )

  if ( value < 0.0 ):
    value = value + 2.0 * np.pi

  value = 180.0 * value / np.pi

  return value

def angle_degree_test ( ):

#*****************************************************************************80
#
## angle_degree_test() tests angle_degree().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_angle = 12

  print ( '' )
  print ( 'angle_degree_test():' )
  print ( '  angle_degree() computes an angle.' )
  print ( '' )
  print ( '           X           Y       Theta  atan2(y,x)  angle_degree' )
  print ( '' )

  x2 = 0.0
  y2 = 0.0
  x3 = 1.0
  y3 = 0.0

  for i in range ( 0, n_angle + 1 ):

    thetad = float ( i ) * 360.0       / float ( n_angle )
    thetar = float ( i ) * 2.0 * np.pi / float ( n_angle )

    x1 = np.cos ( thetar )
    y1 = np.sin ( thetar )

    t1 = np.arctan2 ( y1, x1 ) * 180.0 / np.pi

    t2 = angle_degree ( x1, y1, x2, y2, x3, y3 )

    print ( '  %10f  %10f  %10f  %10f  %10f' \
      % ( x1, y1, thetad, t1, t2 ) )

  return

def angle_half ( x1, y1, x2, y2, x3, y3 ):

#*****************************************************************************80
#
## angle_half() finds half an angle.
#
#  Discussion:
#
#    The original angle is defined by the sequence of points P1, P2 and P3.
#
#    The point P4 is calculated so that:
#
#      (P1,P2,P4) = (P1,P2,P3) / 2
#
#        P1
#        /
#       /   P4
#      /  .
#     / .
#    P2--------->P3
#
#    Thanks to Cesar Fraga Bobis for pointing out a typographical error in
#    a previous version of this routine.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X1, Y1, X2, Y2, X2, Y3, points defining the angle.
#
#  Output:
#
#    real P4(2), a point defining the half angle.
#    The vector P4 - P2 will have unit norm.
#
  import numpy as np

  p4 = np.zeros ( 2 )

  norm12 = np.sqrt ( ( x1 - x2 ) ** 2 + ( y1 - y2 ) ** 2 )
  norm32 = np.sqrt ( ( x3 - x2 ) ** 2 + ( y3 - y2 ) ** 2 ) 

  x4 = 0.5 * ( ( x1 - x2 ) / norm12 + ( x3 - x2 ) / norm32 )

  y4 = 0.5 * ( ( y1 - y2 ) / norm12 + ( y3 - y2 ) / norm32 )

  norm = np.sqrt ( x4 ** 2 + y4 ** 2 )

  p4[0] = x2 + x4 / norm
  p4[1] = y2 + y4 / norm

  return p4

def angle_half_test ( ):

#*****************************************************************************80
#
## angle_half_test() tests angle_half().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_angle = 12

  print ( '' )
  print ( 'angle_half_test():' )
  print ( '  angle_half() is given P1, P2, P3, forming an angle.' )
  print ( '  It finds P4 so P4, P2, P3 is half the angle.' )
  print ( '' )
  print ( '  Original Angle    Half Angle' )
  print ( '' )

  x2 = 0.0
  y2 = 0.0
  x3 = 1.0
  y3 = 0.0

  for i in range ( 0, n_angle + 1 ):

    thetad = float ( i ) * 360.0       / float ( n_angle )
    thetar = float ( i ) * 2.0 * np.pi / float ( n_angle )

    x1 = np.cos ( thetar )
    y1 = np.sin ( thetar )

    t1 = np.arctan2 ( y1, x1 )

    t2 = angle_radian ( x1, y1, x2, y2, x3, y3 )

    p4 = angle_half ( x1, y1, x2, y2, x3, y3 )
    x4 = p4[0]
    y4 = p4[1]

    t3 = angle_radian ( x4, y4, x2, y2, x3, y3 )

    print ( '  %10f  %10f' % ( t2, t3 ) )

  return

def angle_radian ( x1, y1, x2, y2, x3, y3 ):

#*****************************************************************************80
#
## angle_radian() returns the angle swept out between two rays.
#
#  Discussion:
#
#    Except for the zero angle case, it should be true that
#
#      angle_radian(P1,P2,P3) + angle_radian(P3,P2,P1) = 2 * PI
#
#        P1
#        /
#       /
#      /
#     /
#    P2--------->P3
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X1, Y1, X2, Y2, X3, Y3, define the rays
#    P1 - P2 and P3 - P2 which in turn define the angle.
#
#  Output:
#
#    real VALUE, the angle swept out by the rays, measured
#    in radians.  0 <= VALUE < 2*PI.  If either ray has zero length,
#    then VALUE is set to 0.
#
  import numpy as np

  p = np.zeros ( 2 )

  p[0] = ( x3 - x2 ) * ( x1 - x2 ) \
       + ( y3 - y2 ) * ( y1 - y2 )

  p[1] = ( x3 - x2 ) * ( y1 - y2 ) \
       - ( y3 - y2 ) * ( x1 - x2 )

  if ( p[0] == 0.0 and p[1] == 0.0 ):
    value = 0.0
    return value

  value = np.arctan2 ( p[1], p[0] )

  if ( value < 0.0 ):
    value = value + 2.0 * np.pi

  return value

def angle_radian_test ( ):

#*****************************************************************************80
#
## angle_radian_test() tests angle_radian().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_angle = 12

  print ( '' )
  print ( 'angle_radian_test():' )
  print ( '  angle_radian() computes an angle in radians.' )
  print ( '' )
  print ( '           X           Y       Theta  atan2(y,x)  angle_radian' )
  print ( '' )

  x2 = 0.0
  y2 = 0.0
  x3 = 1.0
  y3 = 0.0

  for i in range ( 0, n_angle + 1 ):

    thetad = float ( i ) * 360.0       / float ( n_angle )
    thetar = float ( i ) * 2.0 * np.pi / float ( n_angle )

    x1 = np.cos ( thetar )
    y1 = np.sin ( thetar )

    t1 = np.arctan2 ( y1, x1 )

    t2 = angle_radian ( x1, y1, x2, y2, x3, y3 )

    print ( '  %10f  %10f  %10f  %10f  %10f' \
      % ( x1, y1, thetad, t1, t2 ) )

  return

def annulus_area ( r1, r2 ):

#*****************************************************************************80
#
## annulus_area() computes the area of a circular annulus in 2D.
#
#  Discussion:
#
#    A circular annulus with center (XC,YC), inner radius R1 and
#    outer radius R2, is the set of points (X,Y) so that
#
#      R1^2 <= (X-XC)^2 + (Y-YC)^2 <= R2^2
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R1, R2, the inner and outer radii of the annulus.
#
#  Output:
#
#    real AREA, the area of the annulus.
#
  import numpy as np

  area = np.pi * ( r2 + r1 ) * ( r2 - r1 )

  return area

def annulus_area_test ( ):

#*****************************************************************************80
#
## annulus_area_test() tests annulus_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  pc = np.array ( [ 5.0, 3.0 ] )
  r1 = 2.0
  r2 = 3.0

  print ( '' )
  print ( 'annulus_area_test():' )
  print ( '  annulus_area() computes the centroid of a circular annulus.' )
  print ( '' )
  print ( '  The circle has center        %f  %f' % ( pc[0], pc[1] ) )
  print ( '  The inner radius is R1 =     %f' % ( r1 ) )
  print ( '  The outer radius is R2 =     %f' % ( r2 ) )

  area = annulus_area ( r1, r2 )

  print ( '' )
  print ( '  Area: %f' % ( area ) )

  return

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
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VALUE, the volume of the unit ball.
#
  import numpy as np

  r = 1.0
  value = 4.0 * np.pi * r ** 3 / 3.0

  return value

def ball01_volume_test ( ) :

#*****************************************************************************80
#
## ball01_volume_test() tests ball01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'ball01_volume_test():' )
  print ( '  ball01_volume() returns the volume of the unit ball.' )

  value = ball01_volume ( )

  print ( '' )
  print ( '  ball01_volume() = %g' % ( value ) )

  return

def between ( xa, ya, xb, yb, xc, yc ):

#*****************************************************************************80
#
## between() is TRUE if vertex C is between vertices A and B.
#
#  Discussion:
#
#    The points must be (numerically) collinear.
#
#    Given that condition, we take the greater of XA - XB and YA - YB
#    as a "scale" and check where C's value lies.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    real XA, YA, XB, YB, XC, YC, the coordinates of 
#    the vertices.
#
#  Output:
#
#    bool VALUE, is TRUE if C is between A and B.
#
  if ( not collinear ( xa, ya, xb, yb, xc, yc ) ):
    value = False
  elif ( abs ( ya - yb ) < abs ( xa - xb ) ):
    xmax = max ( xa, xb )
    xmin = min ( xa, xb )
    value = ( xmin <= xc and xc <= xmax )
  else:
    ymax = max ( ya, yb )
    ymin = min ( ya, yb )
    value = ( ymin <= yc and yc <= ymax )

  return value
def circle01_length ( ):

#*****************************************************************************80
#
## circle01_length(): length of the circumference of the unit circle in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VALUE, the length.
#
  import numpy as np

  r = 1.0
  value = 2.0 * np.pi * r

  return value

def circle01_length_test ( ) :

#*****************************************************************************80
#
## circle01_length_test() tests circle01_length().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'circle01_length_test():' )
  print ( '  circle01_length() returns the length of the unit circle.' )

  value = circle01_length ( )

  print ( '' )
  print ( '  circle01_length() = %g' % ( value ) )

  return

def circle_area ( r ):

#*****************************************************************************80
#
## circle_area() returns the area of a circle.
#
#  Integration region:
#
#    Points (X,Y) such that
#
#      ( X - XC )^2 + ( Y - YC )^2 <= R^2
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#  Output:
#
#    real VALUE, the area of the circle.
#
  import numpy as np

  value = np.pi * r * r

  return value

def circle_area_test ( ):

#*****************************************************************************80
#
## circle_area_test() tests circle_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'circle_area_test():' )
  print ( '  circle_area() computes the area of a circle of radius R.' )
  print ( '' )
  print ( '      R            Area' )
  print ( '' )

  r = 1.0
  for i in range ( 0, 4 ):
    area = circle_area ( r )
    print ( '  %10f  %10f' % ( r, area ) )
    r = r * 2

  return

def circle_dia2imp_2d ( p1, p2 ):

#*****************************************************************************80
#
## circle_dia2imp_2d() converts a diameter to an implicit circle in 2D.
#
#  Discussion:
#
#    The diameter form of a circle is:
#
#      P1 and P2 are the endpoints of a diameter.
#
#    The implicit form of a circle in 2D is:
#
#      ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 = R^2
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 December 2020
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(2), P2(2), two points that are the
#    endpoints of a diameter of the circle.
#
#  Output:
#
#    real R, the radius of the circle.
#
#    real CENTER(2), the center of the circle.
#
  import numpy as np

  r = 0.5 * np.sqrt ( np.sum ( ( p2 - p1 )**2 ) )

  center = 0.5 * ( p1 + p2 )

  return r, center

def circle_dia2imp_2d_test ( ):

#*****************************************************************************80
#
## circle_dia2imp_2d_test() tests circle_dia2imp_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 September 2020
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'circle_dia2imp_2d_test():' )
  print ( '  circle_dia2imp_2d() converts a diameter to an' )
  print ( '  implicit circle in 2D.' )

  theta = 2.0

  p1 = np.zeros ( 2 )
  p1[0] = 2.0 + 5.0 * np.cos ( theta )
  p1[1] = 3.0 + 5.0 * np.sin ( theta )

  p2 = np.zeros ( 2 )
  p2[0] = 2.0 - 5.0 * np.cos ( theta )
  p2[1] = 3.0 - 5.0 * np.sin ( theta )

  print ( '' )
  print ( '  point P1:' )
  print ( p1 )
  print ( '' )
  print ( '  point P2:' )
  print ( p2 )

  r, center = circle_dia2imp_2d ( p1, p2 )

  circle_imp_print_2d ( r, center, '  The implicit circle:' )

  return

def circle_imp_point_dist_2d ( r, center, p ):

#*****************************************************************************80
#
## circle_imp_point_dist_2d(): distance ( implicit circle, point ) in 2D.
#
#  Discussion:
#
#    The distance is zero if the point is on the circle.
#
#    An implicit circle in 2D satisfies the equation:
#
#      ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 = R^2
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real CENTER(2,1), the center of the circle.
#
#    real P(2,1), the point to be checked.
#
#  Output:
#
#    real DIST, the distance of the point to the circle.
#
  import numpy as np

  r2 = np.sqrt ( ( p[0] - center[0] ) ** 2 + ( p[1] - center[1] ) ** 2 )

  dist = abs ( r2 - r )

  return dist

def circle_imp_point_dist_2d_test ( ):

#*****************************************************************************80
#
## circle_imp_point_dist_2d_test() tests circle_imp_point_dist_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 5
  center = np.array ( [ 0.0, 0.0 ] )
  r = 5.0

  print ( '' )
  print ( 'circle_imp_point_dist_2d_test():' )
  print ( '  circle_imp_point_dist_2d() checks, by finding the' )
  print ( '  distance D from a point (X,Y) to a circle.' )

  print ( '' )
  print ( '  Circle has center (%f,%f) and radius %f' % ( center[0], center[1], r ) )

  print ( '' )
  print ( '       X       Y       D' )
  print ( '' )

  for i in range ( 0, 10 ):

    p = -10.0 + 20.0 * np.random.rand ( 2 )
    d = circle_imp_point_dist_2d ( r, center, p )
    print ( '  %8.4f  %8.4f  %8.4f' % ( p[0], p[1], d ) )

  return

def circle_imp_print_2d ( r, center, title ):

#*****************************************************************************80
#
## circle_imp_print_2d() prints an implicit circle in 2D.
#
#  Discussion:
#
#    An implicit circle in 2D satisfies:
#
#      ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 = R^2
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real CENTER(2), the center of the circle.
#
#    string TITLE, a title.
#
  print ( '' )
  print ( '%s' % ( title ) )
  print ( '' )
  print ( '  Radius = %f' % ( r ) )
  print ( '  Center = ( %f,  %f )' % ( center[0], center[1] ) )

  return

def circle_imp_print_2d_test ( ):

#*****************************************************************************80
#
## circle_imp_print_2d_test() tests circle_imp_print_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  center = np.array ( [ 5.0, -2.0 ] )
  r = 2.0
 
  print ( '' )
  print ( 'circle_imp_print_2d_test():' )
  print ( '  circle_imp_print_2d() prints a circle definition.' )

  circle_imp_print_2d ( r, center, '  An example circle:' )

  return

def circle_lune_angle_by_height_2d ( r, h ):

#*****************************************************************************80
#
## circle_lune_angle_by_height_2d() computes the angle of a circular lune.
#
#  Discussion:
#
#    Draw the chord connecting two points on the circumference of a circle.
#    The region between the chord and the circumference is a "lune".
#    We wish to know the angle subtended by the lune.
#
#    The distance from the center of the circle to the midpoint of the chord
#    is the "height" H of the lune.  It is natural to expect 0 <= H <= R.
#    However, if we allow -R <= H < 0 as well, this allows us to include
#    lunes which involve more than half the circle's area.
#
#    If H < -R or R < H, then no lune is formed, and we return a zero angle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real H, the height of the lune.
#
#  Output:
#
#    real ANGLE, the angle of the lune.
#
  import numpy as np

  if ( -r <= h and h <= r ):
    angle = 2.0 * np.arccos ( h / r )
  else:
    angle = 0.0

  return angle

def circle_lune_angle_by_height_2d_test ( ):

#*****************************************************************************80
#
## circle_lune_angle_by_height_2d_test() tests circle_lune_angle_by_height_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  n_test = 6

  r = 2.0

  print ( '' )
  print ( 'circle_lune_angle_by_height_2d_test():' )
  print ( '  circle_lune_angle_by_height_2d() computes the angle of a' )
  print ( '  circular lune based on the "height" of the circular triangle.' )
  print ( '' )
  print ( '      R            H        Angle' )
  print ( '' )

  for i in range ( - n_test, n_test + 1 ):

    h = i * r / n_test

    angle = circle_lune_angle_by_height_2d ( r, h )

    print ( '  %10f  %10f  %10f' % ( r, h, angle ) )

  return

def circle_lune_area_by_angle_2d ( r, center, theta1, theta2 ):

#*****************************************************************************80
#
## circle_lune_area_by_angle_2d() returns the area of a circular lune in 2D.
#
#  Discussion:
#
#    A lune is formed by drawing a circular arc, and joining its endpoints.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real CENTER(2,1), the center of the circle.
#
#    real THETA1, THETA2, the angles of the rays
#    that begin and end the arc.
#
#  Output:
#
#    real AREA, the area of the lune.
#
  sector = circle_sector_area_2d ( r, center, theta1, theta2 )
  triangle = circle_triangle_area_2d ( r, center, theta1, theta2 )
  area = sector - triangle

  return area

def circle_lune_area_by_angle_2d_test ( ):

#*****************************************************************************80
#
## circle_lune_area_by_angle_2d_test() tests circle_lune_area_by_angle_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_test = 12

  center = np.array ( [ 0.0, 0.0 ] )
  r = 2.0

  print ( '' )
  print ( 'circle_lune_area_by_angle_2d_test():' )
  print ( '  circle_lune_area_by_angle_2d() computes the area of a' )
  print ( '  circular lune, defined by joining the endpoints' )
  print ( '  of a circular arc.' )
  print ( '' )
  print ( '      R            Theta1      Theta2        Area' )
  print ( '' )

  for i in range ( 0, n_test + 1 ): 

    theta1 = 0.0
    theta2 = i * 2.0 * np.pi / n_test

    area = circle_lune_area_by_angle_2d ( r, center, theta1, theta2 )

    print ( '  %10f  %10f  %10f  %10f' % ( r, theta1, theta2, area ) )

  return

def circle_lune_area_by_height_2d ( r, h ):

#*****************************************************************************80
#
## circle_lune_area_by_height_2d() computes the area of a circular lune.
#
#  Discussion:
#
#    Draw the chord connecting two points on the circumference of a circle.
#    The region between the chord and the circumference is a "lune".
#    We wish to know the area of this region.
#
#    The distance from the center of the circle to the midpoint of the chord
#    is the "height" H of the lune.  It is natural to expect 0 <= H <= R.
#    However, if we allow -R <= H < 0 as well, this allows us to include
#    lunes which involve more than half the circle's area.
#
#    If H < -R or R < H, then no lune is formed and we have zero area.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real H, the height of the lune.
#
#  Output:
#
#    real AREA, the area of the lune.
#
  import numpy as np

  if ( -r <= h and h <= r ):
    area = r ** 2 * np.arccos ( h / r ) - h * np.sqrt ( r ** 2 - h ** 2 )
  else:
    area = 0.0

  return area

def circle_lune_area_by_height_2d_test ( ):

#*****************************************************************************80
#
## circle_lune_area_by_height_2d_test() tests circle_lune_area_by_height_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  n_test = 6

  r = 2.0

  print ( '' )
  print ( 'circle_lune_area_by_height_2d_test():' )
  print ( '  circle_lune_area_by_height_2d() computes the area of a' )
  print ( '  circular lune based on the "height" of the circular triangle.' )
  print ( '' )
  print ( '      R            H        Area' )
  print ( '' )

  for i in range ( - n_test, n_test + 1 ):

    h = i * r / n_test

    area = circle_lune_area_by_height_2d ( r, h )

    print ( '  %10f  %10f  %10f' % ( r, h, area ) )

  return

def circle_lune_height_by_angle_2d ( r, angle ):

#*****************************************************************************80
#
## circle_lune_height_by_angle_2d() computes the height of a circular lune.
#
#  Discussion:
#
#    Draw the chord connecting two points on the circumference of a circle.
#    The region between the chord and the circumference is a "lune".
#    The lune subtends a given angle between 0 and 2 pi.
#
#    The distance from the center of the circle to the midpoint of the chord
#    is the "height" H of the lune and we wish to determine this value.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real ANGLE, the angle subtended by the lune.
#
#  Output:
#
#    real HEIGHT, the height of the lune
#
  import numpy as np

  height = r * np.cos ( angle / 2.0 )

  return height

def circle_lune_height_by_angle_2d_test ( ):

#*****************************************************************************80
#
## circle_lune_height_by_angle_2d_test() tests circle_lune_height_by_angle_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_test = 12

  r = 2.0

  print (  '' )
  print ( 'circle_lune_height_by_angle_2d_test():' )
  print ( '  circle_lune_height_by_angle_2d() computes the height of' )
  print ( '  the triangle of a circular lune, given the subtended angle.' )
  print ( '' )
  print ( '      R            Angle        Height' )
  print ( '' )

  for i in range ( 0, n_test + 1 ):

    angle = i * 2.0 * np.pi / n_test

    height = circle_lune_height_by_angle_2d ( r, angle )

    print ( '  %10f  %10f  %10f' % ( r, angle, height ) )

  return

def circle_sector_area_2d ( r, center, theta1, theta2 ):

#*****************************************************************************80
#
## circle_sector_area_2d() returns the area of a circular sector in 2D.
#
#  Discussion:
#
#    A sector is contained within a circular arc and the lines joining each
#    endpoint of the arc to the center of the circle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real CENTER(2,1), the center of the circle.
#
#    real THETA1, THETA2, the angles of the rays
#    that delimit the sector.
#
#  Output:
#
#    real AREA, the area of the sector.
#
  area = 0.5 * r * r * ( theta2 - theta1 )

  return area

def circle_sector_area_2d_test ( ):

#*****************************************************************************80
#
## circle_sector_area_2d_test() tests circle_sector_area_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_test = 12

  center = np.array ( [ 0.0, 0.0 ] )
  r = 2.0

  print ( '' )
  print ( 'circle_sector_area_2d_test():' )
  print ( '  circle_sector_area_2d() computes the area of a' )
  print ( '  circular sector, defined by joining the endpoints' )
  print ( '  of a circular arc.' )
  print ( '' )
  print ( '      R            Theta1      Theta2        Area' )
  print ( '' )

  for i in range ( 0, n_test + 1 ):

    theta1 = 0.0
    theta2 = i * 2.0 * np.pi / n_test

    area = circle_sector_area_2d ( r, center, theta1, theta2 )

    print ( '  %10f  %10f  %10f  %10f' % ( r, theta1, theta2, area ) )

  return

def circles_intersect_area_2d ( r1, r2, d ):

#*****************************************************************************80
#
## circles_intersect_area_2d(): area of the intersection of two circles.
#
#  Discussion:
#
#    Circles of radius R1 and R2 are D units apart.  What is the area of
#    intersection?
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R1, R2, the radiuses of the circles.
#    R1 and R2 should be positive.
#
#    real D, the distance between the circular centers.
#    D must be positive, and should be no greater than R1 + R2.
#
#  Output:
#
#    real AREA, the area of the intersection.
#
  import numpy as np

  if ( r1 + r2 < d ):
    area = 0.0
  elif ( d == 0.0 ):
    area = np.pi * ( min ( r1, r2 ) ) ** 2
  else:
    h1 = 0.5 * ( d ** 2 + r1 ** 2 - r2 ** 2 ) / d
    area1 = circle_lune_area_by_height_2d ( r1, h1 )
    h2 = 0.5 * ( d ** 2 - r1 ** 2 + r2 ** 2 ) / d
    area2 = circle_lune_area_by_height_2d ( r2, h2 )
    area = area1 + area2

  return area

def circles_intersect_area_2d_test ( ):

#*****************************************************************************80
#
## circles_intersect_area_2d_test() tests circles_intersect_area_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 6
  r1_test = np. array ( [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ] )
  r2_test = np. array ( [ 0.5, 0.5, 0.5, 1.0, 1.0, 1.0 ] )
  d_test =  np. array ( [ 1.5, 1.0, 0.5, 1.5, 1.0, 0.0 ] )

  print ( '' )
  print ( 'circles_intersect_area_2d_test():' )
  print ( '  circles_intersect_area_2d() determines the area of the' )
  print ( '  intersection of two circes of radius R1 and R2,' )
  print ( '  with a distance D between the centers.' )
  print ( '' )
  print ( '      R1      R2       D    Area' )
  print ( '' )

  for i in range ( 0, ntest ):

    r1 = r1_test[i]
    r2 = r2_test[i]
    d = d_test[i]
    area = circles_intersect_area_2d ( r1, r2, d )

    print ( '  %6f  %6f  %6f  %6f' % ( r1, r2, d, area ) )

  return

def circles_intersect_points_2d ( r1, center1, r2, center2 ):

#*****************************************************************************80
#
## circles_intersect_points_2d(): intersection points of two circles in 2D.
#
#  Discussion:
#
#    Two circles can intersect in 0, 1, 2 or infinitely many points.
#
#    The 0 and 2 intersection cases are numerically robust the 1 and
#    infinite intersection cases are numerically fragile.  The routine
#    uses a tolerance to try to detect the 1 and infinite cases.
#
#    An implicit circle in 2D satisfies the equation:
#
#      ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 = R^2
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R1, the radius of the first circle.
#
#    real CENTER1(2), the center of the first circle.
#
#    real R2, the radius of the second circle.
#
#    real CENTER2(2), the center of the second circle.
#
#  Output:
#
#    integer NUM_int, the number of intersecting points found.
#    NUM_int will be 0, 1, 2 or 3.  3 indicates that there are an infinite
#    number of intersection points.
#
#    real P(2,2), if NUM_int is 1 or 2,
#    the coordinates of the intersecting points.
#
  import numpy as np

  tol = np.finfo ( float ).eps

  p = np.zeros ( [ 2, 2 ] )
#
#  Take care of the case in which the circles have the same center.
#
  t1 = ( abs ( center1[0] - center2[0] ) \
       + abs ( center1[1] - center2[1] ) ) / 2.0

  t2 = ( abs ( center1[0] ) + abs ( center2[0] ) \
       + abs ( center1[1] ) + abs ( center2[1] ) + 1.0 ) / 5.0

  if ( t1 <= tol * t2 ):

    t1 = abs ( r1 - r2 )
    t2 = ( abs ( r1 ) + abs ( r2 ) + 1.0 ) / 3.0

    if ( t1 <= tol * t2 ):
      num_int = 3
    else:
      num_int = 0

    return num_int, p

  distsq = ( center1[0] - center2[0] ) ** 2 + ( center1[1] - center2[1] ) ** 2

  root = 2.0 * ( r1 * r1 + r2 * r2 ) * distsq - distsq * distsq \
    - ( r1 - r2 ) * ( r1 - r2 ) * ( r1 + r2 ) * ( r1 + r2 )

  if ( root < -tol ):
    num_int = 0
    return num_int, p

  sc1 = ( distsq - ( r2 * r2 - r1 * r1 ) ) / distsq

  if ( root < tol ):
    num_int = 1
    for i in range ( 0, 2 ):
      p[i,0] = center1[i] + 0.5 * sc1 * ( center2[i] - center1[i] )
    return num_int, p

  sc2 = np.sqrt ( root ) / distsq

  num_int = 2

  p[0,0] = center1[0] + 0.5 * sc1 * ( center2[0] - center1[0] ) \
                      - 0.5 * sc2 * ( center2[1] - center1[1] )
  p[1,0] = center1[1] + 0.5 * sc1 * ( center2[1] - center1[1] ) \
                      + 0.5 * sc2 * ( center2[0] - center1[0] )

  p[0,1] = center1[0] + 0.5 * sc1 * ( center2[0] - center1[0] ) \
                      + 0.5 * sc2 * ( center2[1] - center1[1] )
  p[1,1] = center1[1] + 0.5 * sc1 * ( center2[1] - center1[1] ) \
                      - 0.5 * sc2 * ( center2[0] - center1[0] )

  return num_int, p

def circles_intersect_points_2d_test ( ):

#*****************************************************************************80
#
## circles_intersect_points_2d_test() tests circles_intersect_points_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 5
  center1 = np.array ( [ 0.0, 0.0 ] )
  r1 = 5.0
  r2_test = np.array ( [ 0.5, 5.0, 3.0, 3.0, 5.0 ] )
  xc2_test = np.array ( [ 5.0, 7.0710678, 4.0, 6.0, 0.0 ] )
  yc2_test = np.array ( [ 5.0, 7.0710678, 0.0, 0.0, 0.0 ] )

  print ( '' )
  print ( 'circles_intersect_points_2d_test():' )
  print ( '  circles_intersect_points_2d() determines the intersections of' )
  print ( '  two circles in 2D.' )

  circle_imp_print_2d ( r1, center1, '  The first circle:' )

  for i in range ( 0, ntest ):

    r2 = r2_test[i]
    center2 = np.array ( [ xc2_test[i], yc2_test[i] ] )

    circle_imp_print_2d ( r2, center2, '  The second circle:' )

    num_int, x = circles_intersect_points_2d ( r1, center1, r2, center2 )

    if ( num_int == 0 ):

      print ( '' )
      print ( '  The circles do not intersect.' )

    elif ( num_int == 1 ):

      print ( '' )
      print ( '  The circles intersect at one point:' )
      print ( '' )
      print ( '    X       Y' )
      print ( '' )
      print ( '  %6f  %6f' % ( x[0,0], x[1,0] ) )

    elif ( num_int == 2 ):

      print ( '' )
      print ( '  The circles intersect at two points:' )
      print ( '' )
      print ( '    X       Y' )
      print ( '' )
      print ( '  %6f  %6f' % ( x[0,0], x[1,0] ) )
      print ( '  %6f  %6f' % ( x[0,1], x[1,1] ) )

    elif ( num_int == 3 ):

      print ( '' )
      print ( '  The circles coincide (infinite intersection).' )

  return

def circle_triangle_area_2d ( r, center, theta1, theta2 ):

#*****************************************************************************80
#
## circle_triangle_area_2d() returns the area of a circle triangle in 2D.
#
#  Discussion:
#
#    A circle triangle is formed by drawing a circular arc, and considering
#    the triangle formed by the endpoints of the arc plus the center of
#    the circle.
#
#    Note that for angles greater than PI, the triangle will actually
#    have NEGATIVE area.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    05 May 2005
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the circle.
#
#    real CENTER(2,1), the center of the circle.
#
#    real THETA1, THETA2, the angles of the rays that
#    delimit the arc.
#
#  Output:
#
#    real AREA, the (signed) area of the triangle.
#
  import numpy as np

  area = 0.5 * r * r * np.sin ( theta2 - theta1 )

  return area

def circle_triangle_area_2d_test ( ):

#*****************************************************************************80
#
## circle_triangle_area_2d_test() tests circle_triangle_area_2d().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n_test = 12

  center = np.array ( [ 0.0, 0.0 ] )
  r = 2.0

  print ( '' )
  print ( 'circle_triangle_area_2d_test():' )
  print ( '  circle_triangle_area_2d() computes the area of a' )
  print ( '  circular triangle.' )
  print ( '' )
  print ( '      R            Theta1      Theta2        Area' )
  print ( '' )

  for i in range ( 0, n_test + 1 ):

    theta1 = 0.0
    theta2 = i * 2.0 * np.pi / n_test

    area = circle_triangle_area_2d ( r, center, theta1, theta2 )

    print ( '  %10f  %10f  %10f  %10f' % ( r, theta1, theta2, area ) )

  return

def collinear ( xa, ya, xb, yb, xc, yc ):

#*****************************************************************************80
#
## collinear() returns a measure of collinearity for three points.
#
#  Discussion:
#
#    In order to deal with collinear points whose coordinates are not
#    numerically exact, we compare the area of the largest square
#    that can be created by the line segment between two of the points
#    to (twice) the area of the triangle formed by the points.
#
#    If the points are collinear, their triangle has zero area.
#    If the points are close to collinear, then the area of this triangle
#    will be small relative to the square of the longest segment.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    10 September 2016
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    real XA, YA, XB, YB, XC, YC, the coordinates of 
#    the vertices.
#
#  Output:
#
#    bool VALULE, is TRUE if the points are judged 
#    to be collinear.
#
  r8_eps = 2.220446049250313E-016

  area = triangle_area_2 ( xa, ya, xb, yb, xc, yc )

  side_ab_sq = ( xa - xb ) ** 2 + ( ya - yb ) ** 2
  side_bc_sq = ( xb - xc ) ** 2 + ( yb - yc ) ** 2
  side_ca_sq = ( xc - xa ) ** 2 + ( yc - ya ) ** 2

  side_max_sq = max ( side_ab_sq, max ( side_bc_sq, side_ca_sq ) )

  if ( side_max_sq <= r8_eps ):
    value = True
  elif ( 2.0 * abs ( area ) <= r8_eps * side_max_sq ):
    value = True
  else:
    value = False

  return value

def cone_volume ( r, h ):

#*****************************************************************************80
#
## cone_volume() returns the volume of a cone in 3D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 May 2004
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the base of the cone.
#
#    real H, the height of the cone.
#
#  Output:
#
#    real VALUE, the volume of the cone.
#
  import numpy as np

  value = ( np.pi / 3.0 ) * h * r * r

  return value

def cone_volume_test ( ):

#*****************************************************************************80
#
## cone_volume_test() tests cone_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 August 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'cone_volume_test():' )
  print ( '  cone_volume() computes the volume of a cone.' )
  print ( '' )
  print ( '        R        H        ConeVolume' )
  print ( '' )

  r = 1.0
  h = 1.0
  for i in range ( 0, 5 ):
    print ( '  %14.8f  %14.8f  %14.8f' % ( r, h, cone_volume ( r, h ) ) )
    h = h * 2.0

  print ( '' )

  r = 1.0
  h = 1.0
  for i in range ( 0, 5 ):
    print ( '  %14.8f  %14.8f  %14.8f' % ( r, h, cone_volume ( r, h ) ) )
    r = r * 2.0

  return

def cube01_volume ( ):

#*****************************************************************************80
#
## cube01_volume() returns the volume of the unit cube in 3D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VALUE, the volume.
#
  value = 1.0

  return value

def cube01_volume_test ( ) :

#*****************************************************************************80
#
## cube01_volume() tests cube01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'cube01_volume_test():' )
  print ( '  cube01_volume() returns the volume of the unit cube.' )

  value = cube01_volume ( )

  print ( '' )
  print ( '  cube01_volume() = %g' % ( value ) )

  return

def cylinder_volume ( p1, p2, r ):

#*****************************************************************************80
#
## cylinder_volume() returns the volume of a cylinder.
#
#  Discussion:
#
#    A (right) (finite) cylinder in 3D is the set of points
#    contained on or inside a circle of radius R, whose center
#    lies along the line segment from point P1 to P2, and whose
#    plane is perpendicular to that line segment.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    08 March 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(3), P2(3), the first and last points
#    on the axis line of the cylinder.
#
#    real R, the radius of the cylinder.
#
#  Output:
#
#    real VOLUME, the volume of the cylinder.
#
  import numpy as np

  h = np.linalg.norm ( p1 - p2 )

  volume = np.pi * r * r * h

  return volume

def cylinder_volume_test ( ) :

#*****************************************************************************80
#
## cylinder_volume() tests cylinder_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 March 2022
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'cylinder_volume_test():' )
  print ( '  cylinder_volume() returns the volume of a cylinder.' )

  p1 = np.array ( [ 1.0, 2.0, 3.0 ] )
  p2 = np.array ( [ 5.0, 6.0, 5.0 ] )
  r = 5.0

  value = cylinder_volume ( p1, p2, r )

  print ( '' )
  print ( '  cylinder_volume() = ', value )
  print ( '  Exact volume      = ', np.pi * 150.0 )

  return

def degrees_to_radians ( degrees ):

#*****************************************************************************80
#
## degrees_to_radians() converts an angle from degrees to radians.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real DEGREES, an angle in degrees.
#
#  Output:
#
#    real VALUE, the equivalent angle in radians.
#
  import numpy as np

  value = ( degrees / 180.0 ) * np.pi

  return value

def degrees_to_radians_test ( ):

#*****************************************************************************80
#
## degrees_to_radians_test() tests degrees_to_radians().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'degrees_to_radians_test():' )
  print ( '  degrees_to_radians() converts an angle from degrees' )
  print ( '  to radians.' )
  print ( '' )
  print ( '  Degrees     Radians     Degrees' )
  print ( '' )

  for i in range ( -2, 15 ):
    angle_deg = float ( 30 * i )
    angle_rad = degrees_to_radians ( angle_deg )
    angle_deg2 = radians_to_degrees ( angle_rad )
    print ( '  %10f  %10f  %10f' % ( angle_deg, angle_rad, angle_deg2 ) )

  return

def diagonal ( im1, ip1, n, prev_node, next_node, x, y ):

#*****************************************************************************80
#
## diagonal(): VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    integer IM1, IP1, the indices of two vertices.
#
#    integer N, the number of vertices.
#
#    integer PREV_NODE(N), the previous neighbor of each vertex.
#
#    integer NEXT_NODE(N), the next neighbor of each vertex.
#
#    real X(N), Y(N), the coordinates of each vertex.
#
#  Output:
#
#    bool VALUE, the value of the test.
#  
  value1 = in_cone ( im1, ip1, n, prev_node, next_node, x, y )
  value2 = in_cone ( ip1, im1, n, prev_node, next_node, x, y )
  value3 = diagonalie ( im1, ip1, n, next_node, x, y )

  value = ( value1 and value2 and value3 )

  return value

def diagonalie ( im1, ip1, n, next_node, x, y ):

#*****************************************************************************80
#
## diagonalie() is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    integer IM1, IP1, the indices of two vertices.
#
#    integer N, the number of vertices.
#
#    integer NEXT_NODE(N), the next neighbor of each vertex.
#
#    real X(N), Y(N), the coordinates of each vertex.
#
#  Output:
#
#    bool VALUE, the value of the test.
#
  first = im1
  j = first
  jp1 = next_node[first]

  value = True
#
#  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
#
  while ( True ):
#
#  Skip any edge that includes vertex IM1 or IP1.
#
    if ( j == im1 or j == ip1 or jp1 == im1 or jp1 == ip1 ):
      pass
    else:

      value2 = intersect ( \
        x[im1], y[im1], x[ip1], y[ip1], x[j], y[j], x[jp1], y[jp1] )

      if ( value2 ):
        value = False
        break

    j = jp1
    jp1 = next_node[j]

    if ( j == first ):
      break

  return value

def disk01_area ( ):

#*****************************************************************************80
#
## disk01_area() returns the area of the unit disk.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real AREA, the area of the unit disk.
#
  import numpy as np

  r = 1.0
  value = np.pi * r * r

  return value

def disk01_area_test ( ) :

#*****************************************************************************80
#
## disk01_area() tests disk01_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'disk01_area_test():' )
  print ( '  disk01_area() returns the area of the unit disk.' )

  value = disk01_area ( )

  print ( '' )
  print ( '  disk01_area() = %g' % ( value ) )

  return

def disk01_quarter_area ( ):

#*****************************************************************************80
#
## disk01_quarter_area() returns the area of the unit quarter_disk.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    06 May 2016
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real AREA, the area.
#
  import numpy as np

  value = np.pi / 4.0

  return value

def disk01_quarter_area_test ( ) :

#*****************************************************************************80
#
## disk01_quarter_area() tests disk01_quarter_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    06 May 2016
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'disk01_quarter_area_test():' )
  print ( '  disk01_quarter_area() returns the area of the unit quarter disk.' )

  value = disk01_quarter_area ( )

  print ( '' )
  print ( '  disk01_quarter_area() = %g' % ( value ) )

  return

def geometry_test ( ):

#*****************************************************************************80
#
## geometry_test() tests geometry().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    21 July 2021
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'geometry_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test geometry().' )

  angle_degree_test ( )
  angle_half_test ( )
  angle_radian_test ( )

  annulus_area_test ( )

  ball01_volume_test ( )

  circle_area_test ( )
  circle_dia2imp_2d_test ( )
  circle_imp_point_dist_2d_test ( )
  circle_imp_print_2d_test ( )
  circle_lune_angle_by_height_2d_test ( )
  circle_lune_area_by_angle_2d_test ( )
  circle_lune_area_by_height_2d_test ( )
  circle_lune_height_by_angle_2d_test ( )
  circle_sector_area_2d_test ( )
  circle_triangle_area_2d_test ( )

  circle01_length_test ( )

  circles_intersect_area_2d_test ( )
  circles_intersect_points_2d_test ( )

  cone_volume_test ( )

  cube01_volume_test ( )

  cylinder_volume_test ( )

  degrees_to_radians_test ( )

  disk01_area_test ( )

  disk01_quarter_area_test ( )

  hexagon01_area_test ( )

  hyperball01_volume_test ( )

  hypercube01_volume_test ( )

  hypersphere01_area_test ( )

  i4_ceiling_test ( )
  i4_log_10_test ( )
  i4_modp_test ( )
  i4_wrap_test ( )

  line_exp2imp_test ( )
  line_exp_perp_test ( )

  lines_exp_int_test ( )
  lines_imp_int_test ( )

  polygon_area_test ( )

  pyramid_volume_test ( )
  pyramid01_volume_test ( )

  r8_acos_test ( )
  r8_sign_test ( )

  r8mat_det_4d_test ( )
  r8mat_solve_test ( )

  r8vec_indicator1_test ( )

  radians_to_degrees_test ( )

  segment_point_dist_test ( )
  segment_point_near_test ( )

  simplex01_volume_test ( )

  sphere_triangle_sides_to_angles_test ( )

  sphere01_area_test ( )
  sphere01_area_values_test ( )
  sphere01_volume_values_test ( )

  tetrahedron_barycentric_test ( )
  tetrahedron_centroid_test ( )
  tetrahedron_sample_test ( )
  tetrahedron_volume_test ( )

  triangle_angles_test ( )
  triangle_area_test ( )
  triangle_barycentric_test ( )
  triangle_centroid_test ( )
  triangle_circumcircle_test ( )
  triangle_contains_point_test ( )
  triangle_contains_point_1_test ( )
  triangle_diameter_test ( )
  triangle_edge_length_test ( )
  triangle_incircle_test ( )
  triangle_orientation_test ( )
  triangle_orthocenter_test ( )
  triangle_point_dist_test ( )
  triangle_point_near_test ( )
  triangle_quality_test ( )
  triangle_reference_sample_test ( )
  triangle_sample_test ( )
  triangle_xsi_to_xy_test ( )
  triangle_xy_to_xsi_test ( )

  triangle01_area_test ( )
  triangle01_sample_test ( )

  triangle_3d_area_cross_product_test ( )
  triangle_3d_area_heron_test ( )

  wedge01_volume_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'geometry_test():' )
  print ( '  Normal end of execution.' )
  return

def hexagon01_area ( ):

#*****************************************************************************80
#
## hexagon01_area() returns the area of a unit regular hexagon in 2D.
#
#  Discussion:
#
#    The definition is given in terms of THETA, the angle in degrees of the
#    vector (X,Y).  The following six conditions apply, respectively,
#    between the bracketing values of THETA of 0, 60, 120, 180, 240,
#    300, and 360.
#
#                              0 <= Y <= - SQRT(3) * X + SQRT(3)
#                              0 <= Y <=                 SQRT(3)/2
#                              0 <= Y <=   SQRT(3) * X + SQRT(3)
#      - SQRT(3) * X - SQRT(3)   <= Y <= 0
#                    - SQRT(3)/2 <= Y <= 0
#        SQRT(3) * X - SQRT(3)   <= Y <= 0
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2022
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VALUE, the area of the hexagon.
#
  import numpy as np

  value = 3.0 * np.sqrt ( 3.0 ) / 2.0

  return value

def hexagon01_area_test ( ) :

#*****************************************************************************80
#
## hexagon01_area_test() tests hexagon01_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2022
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'hexagon01_area_test():' )
  print ( '  hexagon01_area() returns the area of the unit hexagon.' )
  print ( '' )

  value = hexagon01_area ( )

  print ( '  hexagon01_area() = ', value )

  return

def hyperball01_volume ( m ):

#*****************************************************************************80
#
## hyperball01_volume() returns the volume of the unit hyperball in M dimensions.
#
#  Discussion:
#
#     M  Volume
#
#     1    2
#     2    1        * PI
#     3  ( 4 /   3) * PI
#     4  ( 1 /   2) * PI^2
#     5  ( 8 /  15) * PI^2
#     6  ( 1 /   6) * PI^3
#     7  (16 / 105) * PI^3
#     8  ( 1 /  24) * PI^4
#     9  (32 / 945) * PI^4
#    10  ( 1 / 120) * PI^5
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the spatial dimension.
#
#  Output:
#
#    real VOLUME, the volume of the unit ball.
#
  import numpy as np

  if ( ( m % 2 ) == 0 ):
    m_half = ( m // 2 )
    volume = np.pi ** m_half
    for i in range ( 1, m_half + 1 ):
      volume = volume / float ( i )
  else:
    m_half = ( ( m - 1 ) // 2 )
    volume = np.pi ** m_half * 2.0 ** m
    for i in range ( m_half + 1, 2 * m_half + 2 ):
      volume = volume / float ( i )

  return volume

def hyperball01_volume_test ( ) :

#*****************************************************************************80
#
## hyperball01_volume() tests hyperball01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'hyperball01_volume_test():' )
  print ( '  hyperball01_volume() returns the volume of the unit hyperball' )
  print ( '  in M dimensions.' )
  print ( '' )
  print ( '   M  Volume' )
  print ( '' )

  for m in range ( 1, 11 ):
    value = hyperball01_volume ( m )
    print ( '  %2d  %g' % ( m, value ) )

  return

def hypercube01_volume ( m ):

#*****************************************************************************80
#
## hypercube01_volume() returns the volume of the unit hypercube in M dimensions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the spatial dimension.
#
#  Output:
#
#    real VALUE, the volume.
#
  value = 1.0

  return value

def hypercube01_volume_test ( ) :

#*****************************************************************************80
#
## hypercube01_volume() tests hypercube01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'hypercube01_volume_test():' )
  print ( '  hypercube01_volume() returns the volume of the unit hypercube' )
  print ( '  in M dimensions.' )

  m = 3

  value = hypercube01_volume ( m )

  print ( '' )
  print ( '  hypercube01_volume(%d) = %g' % ( m, value ) )

  return

def hypersphere01_area ( m ):

#*****************************************************************************80
#
## hypersphere01_area() returns the surface area of the unit hypersphere.
#
#  Discussion:
#
#     M   Area
#
#     2    2        * PI
#     3    4        * PI
#     4  ( 2 /   1) * PI^2
#     5  ( 8 /   3) * PI^2
#     6  ( 1 /   1) * PI^3
#     7  (16 /  15) * PI^3
#     8  ( 1 /   3) * PI^4
#     9  (32 / 105) * PI^4
#    10  ( 1 /  12) * PI^5
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    04 January 2014
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the spatial dimension.
#
#  Output:
#
#    real VALUE, the area of the unit hypersphere.
#
  import numpy as np

  if ( ( m % 2 ) == 0 ):
    m_half = ( m // 2 )
    value = 2.0 * np.pi ** m_half
    for i in range (  1, m_half ):
      value = value / float ( i )
  else:
    m_half = ( ( m - 1 ) // 2 )
    value = np.pi ** m_half * 2.0 ** m
    for i in range ( m_half + 1, 2 * m_half + 1 ):
      value = value / float ( i )

  return value

def hypersphere01_area_test ( ) :

#*****************************************************************************80
#
## hypersphere01_area_test() tests hypersphere01_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'hypersphere01_area_test():' )
  print ( '  hypersphere01_area() returns the volume of the unit hypersphere.' )
  print ( '' )
  print ( '   M  Area' )
  print ( '' )

  for m in range ( 1, 11 ):
    value = hypersphere01_area ( m )
    print ( '  %2d  %g' % ( m, value ) )

  return

def i4_ceiling ( x ) :

#*****************************************************************************80
#
## i4_ceiling() rounds an R8 up to the next I4.
#
#  Example:
#
#    X         Value
#
#   -1.1      -1
#   -1.0      -1
#   -0.9       0
#   -0.1       0
#    0.0       0
#    0.1       1
#    0.9       1
#    1.0       1
#    1.1       2
#    2.9       3
#    3.0       3
#    3.14159   4
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    09 May 2013
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the number to be rounded up.
#
#    integer VALUE, the rounded value of X.
#
  import numpy as np

  value = int ( np.ceil ( x ) )

  return value

def i4_ceiling_test ( ):

#*****************************************************************************80
#
## i4_ceiling_test() tests i4_ceiling().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 September 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  r8_lo = -100.0
  r8_hi =  100.0
 
  print ( '' )
  print ( 'i4_ceiling_test():' )
  print ( '  i4_ceiling() evaluates the "ceiling" of a real number.' )
  print ( ' ' )
  print ( '      R8    i4_ceiling(R8)' )
  print ( ' ' )

  for i in range ( 0, 10 ):
    r8 = r8_lo + ( r8_hi - r8_lo ) * np.random.rand ( )
    i4 = i4_ceiling ( r8 )
    print ( '  %8.4f            %4d' % ( r8, i4 ) )

  return

def i4_log_10 ( i ):

#*****************************************************************************80
#
## i4_log_10() returns the integer part of the logarithm base 10 of ABS(X).
#
#  Example:
#
#        I  VALUE
#    -----  --------
#        0    0
#        1    0
#        2    0
#        9    0
#       10    1
#       11    1
#       99    1
#      100    2
#      101    2
#      999    2
#     1000    3
#     1001    3
#     9999    3
#    10000    4
#
#  Discussion:
#
#    i4_log_10 ( I ) + 1 is the number of decimal digits in I.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 May 2013
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer I, the number whose logarithm base 10 is desired.
#
#  Output:
#
#    integer VALUE, the integer part of the logarithm base 10 of
#    the absolute value of X.
#
  import numpy as np

  i = np.floor ( i )

  if ( i == 0 ):

    value = 0

  else:

    value = 0
    ten_pow = 10

    i_abs = np.abs ( i )

    while ( ten_pow <= i_abs ):
      value = value + 1
      ten_pow = ten_pow * 10

  return value

def i4_log_10_test ( ) :

#*****************************************************************************80
#
## i4_log_10_test() tests i4_log_10().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 May 2013
#
#  Author:
#
#    John Burkardt
#
  n = 13

  x = [ 0, 1, 2, 3, 9, 10, 11, 99, 101, -1, -2, -3, -9 ]

  print ( '' )
  print ( 'i4_log_10_test():' )
  print ( '  i4_log_10(): whole part of log base 10,' )
  print ( '' )
  print ( '  X, i4_log_10' )
  print ( '' )

  for i in range ( 0, n ):
    j = i4_log_10 ( x[i] )
    print ( '%6d  %12d' % ( x[i], j ) )

  return

def i4_modp ( i, j ):

#*****************************************************************************80
#
## i4_modp() returns the nonnegative remainder of I4 division.
#
#  Discussion:
#
#    If
#      NREM = i4_modp ( I, J )
#      NMULT = ( I - NREM ) / J
#    then
#      I = J * NMULT + NREM
#    where NREM is always nonnegative.
#
#    The MOD function computes a result with the same sign as the
#    quantity being divided.  Thus, suppose you had an angle A,
#    and you wanted to ensure that it was between 0 and 360.
#    Then mod(A,360) would do, if A was positive, but if A
#    was negative, your result would be between -360 and 0.
#
#    On the other hand, i4_modp(A,360) is between 0 and 360, always.
#
#  Example:
#
#        I     J     MOD  i4_modp    Factorization
#
#      107    50       7       7    107 =  2 *  50 + 7
#      107   -50       7       7    107 = -2 * -50 + 7
#     -107    50      -7      43   -107 = -3 *  50 + 43
#     -107   -50      -7      43   -107 =  3 * -50 + 43
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 May 2013
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer I, the number to be divided.
#
#    integer J, the number that divides I.
#
#  Output:
#
#    integer VALUE, the nonnegative remainder when I is
#    divided by J.
#
  import numpy as np

  if ( j == 0 ):
    print ( '' )
    print ( 'i4_modp - Fatal error!' )
    print ( '  Illegal divisor J = %d' % ( j ) )
    raise Exception ( 'i4_modp - Fatal error!' )

  value = i % j

  if ( value < 0 ):
    value = value + np.abs ( j )

  return value

def i4_modp_test ( ):

#*****************************************************************************80
#
## i4_modp_test() tests i4_modp().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 September 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  test_num = 4

  n_vec = np.array ( ( 107, 107, -107, -107 ) )
  d_vec = np.array ( ( 50, -50, 50, -50 ) )

  print ( '' )
  print ( 'i4_modp_test():' )
  print ( '  i4_modp() factors a number' )
  print ( '  into a multiple M and a positive remainder R.' )
  print ( '' )
  print ( '    Number   Divisor  Multiple Remainder' )
  print ( '' )

  for test in range ( 0, test_num ):
    n = n_vec[test]
    d = d_vec[test]
    r = i4_modp ( n, d )
    m = ( n - r ) // d
    print ( '  %8d  %8d  %8d  %8d' % ( n, d, m, r ) )

  print ( '' )
  print ( '  Repeat using Python % Operator:' )
  print ( '' )

  for test in range ( 0, test_num ):
    n = n_vec[test]
    d = d_vec[test]
    m = n // d
    r = n % d
    print ( '  %8d  %8d  %8d  %8d' % ( n, d, m, r ) )

  return

def i4_wrap ( ival, ilo, ihi ):

#*****************************************************************************80
#
## i4_wrap() forces an integer to lie between given limits by wrapping.
#
#  Example:
#
#    ILO = 4, IHI = 8
#
#    I   Value
#
#    -2     8
#    -1     4
#     0     5
#     1     6
#     2     7
#     3     8
#     4     4
#     5     5
#     6     6
#     7     7
#     8     8
#     9     4
#    10     5
#    11     6
#    12     7
#    13     8
#    14     4
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 May 2013
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer IVAL, an integer value.
#
#    integer ILO, IHI, the desired bounds for the integer value.
#
#  Output:
#
#    integer VALUE, a "wrapped" version of IVAL.
#
  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ):
    value = jlo
  else:
    value = jlo + i4_modp ( ival - jlo, wide )

  return value

def i4_wrap_test ( ):

#*****************************************************************************80
#
## i4_wrap_test() tests i4_wrap().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 May 2013
#
#  Author:
#
#    John Burkardt
#
  ilo = 4
  ihi = 8

  print ( '' )
  print ( 'i4_wrap_test():' )
  print ( '  i4_wrap() forces an integer to lie within given limits.' )
  print ( '' )
  print ( '  ILO = %d' % ( ilo ) )
  print ( '  IHI = %d' % ( ihi ) )
  print ( '' )
  print ( '     I  i4_wrap(I)' )
  print ( '' )

  for i in range ( -10, 21 ):
    j = i4_wrap ( i, ilo, ihi )
    print ( '  %6d  %6d' % ( i, j ) )

  return

def i4mat_print ( m, n, a, title ):

#*****************************************************************************80
#
## i4mat_print() prints an I4MAT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 October 2014
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
#    integer A(M,N), the matrix.
#
#    string TITLE, a title.
#
  i4mat_print_some ( m, n, a, 0, 0, m - 1, n - 1, title )

def i4mat_print_test ( ):

#*****************************************************************************80
#
## i4mat_print_test() tests i4mat_print().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 May 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'i4mat_print_test():' )
  print ( '  i4mat_print() prints an I4MAT.' )

  m = 5
  n = 6
  a = np.array ( ( \
    [ 11, 12, 13, 14, 15, 16 ], \
    [ 21, 22, 23, 24, 25, 26 ], \
    [ 31, 32, 33, 34, 35, 36 ], \
    [ 41, 42, 43, 44, 45, 46 ], \
    [ 51, 52, 53, 54, 55, 56 ] ) )
  title = '  A 5 x 6 integer matrix:'
  i4mat_print ( m, n, a, title )

  return

def i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title ):

#*****************************************************************************80
#
## i4mat_print_some() prints out a portion of an I4MAT.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 May 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, N, the number of rows and columns of the matrix.
#
#    integer A(M,N), an M by N matrix to be printed.
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
      print ( '%7d  ' % ( j ), end = '' )

    print ( '' )
    print ( '  Row' )

    i2lo = max ( ilo, 0 )
    i2hi = min ( ihi, m )

    for i in range ( i2lo, i2hi + 1 ):

      print ( ' %4d: ' % ( i ), end = '' )
      
      for j in range ( j2lo, j2hi + 1 ):
        print ( '%7d  ' % ( a[i,j] ), end = '' )

      print ( '' )

  return

def i4mat_print_some_test ( ):

#*****************************************************************************80
#
## i4mat_print_some_test() tests i4mat_print_some().
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

  print ( '' )
  print ( 'i4mat_print_some_test():' )
  print ( '  i4mat_print_some() prints some of an I4MAT.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11, 12, 13, 14, 15, 16 ], 
    [ 21, 22, 23, 24, 25, 26 ], 
    [ 31, 32, 33, 34, 35, 36 ], 
    [ 41, 42, 43, 44, 45, 46 ] ], dtype = np.int32 )
  i4mat_print_some ( m, n, v, 0, 3, 2, 5, '  Here is I4MAT, rows 0:2, cols 3:5:' )

  return

def i4vec_print ( n, a, title ):

#*****************************************************************************80
#
## i4vec_print() prints an I4VEC.
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
#    integer A(N), the vector to be printed.
#
#    string TITLE, a title.
#
  print ( '' )
  print ( title )
  print ( '' )
  for i in range ( 0, n ):
    print ( '%6d  %6d' % ( i, a[i] ) )

  return

def i4vec_print_test ( ):

#*****************************************************************************80
#
## i4vec_print_test() tests i4vec_print().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 October 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'i4vec_print_test():' )
  print ( '  i4vec_print() prints an I4VEC.' )

  n = 4
  v = np.array ( [ 91, 92, 93, 94 ], dtype = np.int32 )
  i4vec_print ( n, v, '  Here is an I4VEC:' )

  return

def in_cone ( im1, ip1, n, prev_node, next_node, x, y ):

#*****************************************************************************80
#
## in_cone() is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    integer IM1, IP1, the indices of two vertices.
#
#    integer N, the number of vertices.
#
#    integer PREV_NODE(N), the previous neighbor of each vertex.
#
#    integer NEXT_NODE(N), the next neighbor of each vertex.
#
#    real X(N), Y(N), the coordinates of each vertex.
#
#  Output:
#
#    bool VALUE, the value of the test.
#
  im2 = prev_node[im1]
  i = next_node[im1]

  t1 = triangle_area_2 ( x[im1], y[im1], x[i], y[i], x[im2], y[im2] )

  if ( 0.0 <= t1 ):

    t2 = triangle_area_2 ( x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] )
    t3 = triangle_area_2 ( x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] )
    value = ( ( 0.0 < t2 ) and ( 0.0 < t3 ) )

  else:

    t4 = triangle_area_2 ( x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] )
    t5 = triangle_area_2 ( x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] )
    value = not ( ( 0.0 <= t4 ) and ( 0.0 <= t5 ) )

  return value

def intersect ( xa, ya, xb, yb, xc, yc, xd, yd ):

#*****************************************************************************80
#
## intersect() is true if lines VA:VB and VC:VD intersect.
#
#  Discussion:
#
#    Thanks to Gene Dial for correcting the call to intersect_prop(),
#    08 September 2016.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 September 2016
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    real XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
#    coordinates of the four vertices.
#
#  Output:
#
#    bool VALUE, the value of the test.
#
  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ) ):
    value = True
  elif ( between ( xa, ya, xb, yb, xc, yc ) ):
    value = True
  elif ( between ( xa, ya, xb, yb, xd, yd ) ):
    value = True
  elif ( between ( xc, yc, xd, yd, xa, ya ) ):
    value = True
  elif ( between ( xc, yc, xd, yd, xb, yb ) ):
    value = True
  else:
    value = False

  return value

def intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ):

#*****************************************************************************80
#
## intersect_prop() is TRUE if lines VA:VB and VC:VD have a proper intersection.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#    Original C version by Joseph ORourke.
#    FORTRAN90 version by John Burkardt.
#
#  Reference:
#
#    Joseph ORourke,
#    Computational Geometry in C,
#    Cambridge, 1998,
#    ISBN: 0521649765,
#    LC: QA448.D38.
#
#  Input:
#
#    real XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
#    coordinates of the four vertices.
#
#  Output:
#
#    bool VALUE, the result of the test.
#
  if ( collinear ( xa, ya, xb, yb, xc, yc ) ):
    value = False
  elif ( collinear ( xa, ya, xb, yb, xd, yd ) ):
    value = False
  elif ( collinear ( xc, yc, xd, yd, xa, ya ) ):
    value = False
  elif ( collinear ( xc, yc, xd, yd, xb, yb ) ):
    value = False
  else:
    t1 = triangle_area_2 ( xa, ya, xb, yb, xc, yc )
    t2 = triangle_area_2 ( xa, ya, xb, yb, xd, yd )
    t3 = triangle_area_2 ( xc, yc, xd, yd, xa, ya )
    t4 = triangle_area_2 ( xc, yc, xd, yd, xb, yb )

    value1 = ( 0.0 < t1 )
    value2 = ( 0.0 < t2 )
    value3 = ( 0.0 < t3 )
    value4 = ( 0.0 < t4 )

    value = ( l4_xor ( value1, value2 ) ) and ( l4_xor ( value3, value4 ) )
 
  return value

def l4_xor ( l1, l2 ):

#*****************************************************************************80
#
## l4_xor() returns the exclusive OR of two L4's.
#
#  Discussion:
#
#    An L4 is a bool value.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#   John Burkardt
#
#  Input:
#
#    bool L1, L2, two values whose exclusive OR 
#    is needed.
#
#  Output:
#
#    bool VALUE, the exclusive OR of L1 and L2.
#
  value1 = (       l1   and ( not l2 ) )
  value2 = ( ( not l1 ) and       l2   )

  value = ( value1 or value2 )

  return value

def line_exp2imp ( p1, p2 ):

#*****************************************************************************80
#
## line_exp2imp() converts an explicit line to implicit form in 2D.
#
#  Discussion:
#
#    The explicit form of a line in 2D is:
#
#      ( P1, P2 ) = ( (X1,Y1), (X2,Y2) ).
#
#    The implicit form of a line in 2D is:
#
#      A * X + B * Y + C = 0
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(2), P2(2), two points on the line.
#
#    Output, real A, B, C, the implicit form of the line.
#

#
#  Take care of degenerate cases.
#
  if ( p1[0] == p2[0] and p1[1] == p2[1] ):
    print ( '' )
    print ( 'line_exp2imp - Fatal error!' )
    print ( '  P1 = P2' )
    print ( '  P1 = %g  %g' % ( p1[0], p1[1] ) )
    print ( '  P2 = %g  %g' % ( p2[0], p2[1] ) )
    raise Exception ( 'line_exp2imp - Fatal error!' )

  a = p2[1] - p1[1]
  b = p1[0] - p2[0]
  c = p2[0] * p1[1] - p1[0] * p2[1]

  norm = a * a + b * b + c * c

  if ( 0.0 < norm ):
    a = a / norm
    b = b / norm
    c = c / norm

  if ( a < 0.0 ):
    a = -a
    b = -b
    c = -c

  return a, b, c

def line_exp2imp_test ( ):

#*****************************************************************************80
#
## line_exp2imp_test() tests line_exp2imp().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'line_exp2imp_test():' )
  print ( '  line_exp2imp() converts explicit to implicit lines.' )

  a = 1.0
  b = 2.0
  c = 3.0

  print ( '' )
  print ( '  Implicit line A, B, C = %f  %f  %f' % ( a, b, c ) )

  p1, p2 = line_imp2exp ( a, b, c )

  print ( '' )
  print ( '  point P1:' )
  print ( p1 )
  print ( '' )
  print ( '  point P2:' )
  print ( p2 )

  a, b, c = line_exp2imp ( p1, p2 )

  print ( '  Recovered A, B, C =  %f  %f  %f' % ( a, b, c ) )

  return

def line_exp_perp ( p1, p2, p3 ):

#*****************************************************************************80
#
## line_exp_perp() computes a line perpendicular to a line and through a point.
#
#  Discussion:
#
#    The explicit form of a line in 2D is:
#
#      ( P1, P2 ) = ( (X1,Y1), (X2,Y2) ).
#
#    The input point P3 should NOT lie on the line (P1,P2).  If it
#    does, then the output value P4 will equal P3.
#
#    P1-----P4-----------P2
#            |
#            |
#           P3
#
#    P4 is also the nearest point on the line (P1,P2) to the point P3.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(2), P2(2), two points on the line.
#
#    real P3(2), a point (presumably not on the
#    line (P1,P2)), through which the perpendicular must pass.
#
#  Output:
#
#    real P4(2), a point on the line (P1,P2),
#    such that the line (P3,P4) is perpendicular to the line (P1,P2).
#
#    bool FLAG, is TRUE if the value could not be computed.
#
  import numpy as np

  bot = ( p2[0] - p1[0] ) ** 2 + ( p2[1] - p1[1] ) **2

  p4 = np.zeros ( 2 )

  if ( bot == 0.0 ):
    p4[0] = float ( 'inf' )
    p4[1] = float ( 'inf' )
    flag = 1
#
#  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
#
#  (P3-P1) dot (P2-P1) / Norm(P3-P1)**2 = normalized coordinate T
#  of the projection of (P3-P1) onto (P2-P1).
#
  t = ( ( p1[0] - p3[0] ) * ( p1[0] - p2[0] ) \
      + ( p1[1] - p3[1] ) * ( p1[1] - p2[1] ) ) / bot

  p4[0] = p1[0] + t * ( p2[0] - p1[0] )
  p4[1] = p1[1] + t * ( p2[1] - p1[1] )

  flag = 0

  return p4, flag

def line_exp_perp_test ( ):

#*****************************************************************************80
#
## line_exp_perp_test() tests line_exp_perp().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 3

  print ( '' )
  print ( 'line_exp_perp_test():' )
  print ( '  line_exp_perp() is given an explicit line (P1,P2),' )
  print ( '  and another point P3.  It then finds a point' )
  print ( '  P4 on (P1,P2) so that (P1,P2) is perpendicular' )
  print ( '  to (P3,P4).' )

  p1 = np.array ( [ 1.0, 3.0 ] )
  p2 = np.array ( [ 4.0, 0.0 ] )

  p3test = np.array ( [ \
    [ 0.0,  5.0, 5.0 ], \
    [ 0.0, -1.0, 3.0 ] ] )

  print ( '' )
  print ( '  point P1:' )
  print ( p1 )
  print ( '' )
  print ( '  point P2:' )
  print ( p2 )

  p3 = np.zeros ( 2 )

  for j in range ( 0, ntest ):

    p3[0] = p3test[0,j]
    p3[1] = p3test[1,j]
    print ( '' )
    print ( '  point P3:' )
    print ( p3 )

    p4, flag = line_exp_perp ( p1, p2, p3 )
    print ( '' )
    print ( '  point P4:' )
    print ( p4 )

  return

def line_imp2exp ( a, b, c ):

#*****************************************************************************80
#
## line_imp2exp() converts an implicit line to explicit form in 2D.
#
#  Discussion:
#
#    The implicit form of line in 2D is:
#
#      A * X + B * Y + C = 0
#
#    The explicit form of a line in 2D is:
#
#      ( P1, P2 ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer and John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real A, B, C, the implicit line parameters.
#
#  Output:
#
#    real P1(2), P2(2), two points on the line.
#
  import numpy as np

  normsq = a * a + b * b

  if ( normsq == 0.0 ):
    print ( '' )
    print ( 'line_imp2exp - Fatal error!' )
    print ( '  A * A + B * B = 0.' )
    raise Exception ( 'line_imp2exp - Fatal error!' )

  p1 = np.zeros ( 2 )
  p2 = np.zeros ( 2 )

  p1[0] = - a * c / normsq
  p1[1] = - b * c / normsq

  if ( abs ( b ) < abs ( a ) ):
    p2[0] = - ( a - b / a ) * c / normsq
    p2[1] = - ( b + 1.0 ) * c / normsq
  else:
    p2[0] = - ( a + 1.0 ) * c / normsq
    p2[1] = - ( b - a / b ) * c / normsq

  return p1, p2

def line_imp2exp_test ( ):

#*****************************************************************************80
#
## line_imp2exp_test() tests line_imp2exp().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'line_imp2exp_test():' )
  print ( '  line_imp2exp() converts implicit to explicit lines.' )

  a = 1.0
  b = 2.0
  c = 3.0

  print ( '' )
  print ( '  Implicit line A, B, C = %f  %f  %f' % ( a, b, c ) )

  p1, p2 = line_imp2exp ( a, b, c )

  print ( '' )
  print ( '  point P1:' )
  print ( p1 )
  print ( '' )
  print ( '  point P2:' )
  print ( p2 )

  a, b, c = line_exp2imp ( p1, p2 )

  print ( '  Recovered A, B, C =  %f  %f  %f' % ( a, b, c ) )

  return

def lines_exp_int ( p1, p2, q1, q2 ):

#*****************************************************************************80
#
## lines_exp_int() determines where two explicit lines intersect in 2D.
#
#  Discussion:
#
#    The explicit form of a line in 2D is:
#
#      the line through the points P1, P2.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(2,1), P2(2,1), two points on the first line.
#
#    real Q1(2,1), Q2(2,1), two points on the second line.
#
#  Output:
#
#    integer IVAL, reports on the intersection:
#    0, no intersection, the lines may be parallel or degenerate.
#    1, one intersection point, returned in P.
#    2, infinitely many intersections, the lines are identical.
#
#    real P(2,1), if IVAl = 1, P is
#    the intersection point.  Otherwise, P = 0.
#
  import numpy as np

  ival = False
  p = np.zeros ( 2 )
#
#  Check whether either line is a point.
#
  if ( p1[0] == p2[0] and p1[1] == p2[1] ):
    point_1 = True
  else:
    point_1 = False

  if ( q1[0] == q2[0] and q1[1] == q2[1] ):
    point_2 = True
  else:
    point_2 = False
#
#  Convert the lines to ABC format.
#
  if ( not point_1 ):
    a1, b1, c1 = line_exp2imp ( p1, p2 )

  if ( not point_2 ):
    a2, b2, c2 = line_exp2imp ( q1, q2 )
#
#  Search for intersection of the lines.
#
  if ( point_1 and point_2 ):
    if ( p1[0] == q1[0] and p1[1] == q1[1] ):
      ival = True
      p[0] = p1[0]
      p[1] = p1[1]
  elif ( point_1 ):
    if ( a2 * p1[0] + b2 * p1[1] == c2 ):
      ival = True
      p[0] = p1[0]
      p[1] = p1[1]
  elif ( point_2 ):
    if ( a1 * q1[0] + b1 * q1[1] == c1 ):
      ival = True
      p[0] = q1[0]
      p[1] = q1[1]
  else:
    ival, p = lines_imp_int ( a1, b1, c1, a2, b2, c2 )

  return ival, p

def lines_exp_int_test ( ):

#*****************************************************************************80
#
## lines_exp_int_test() tests lines_exp_int().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'lines_exp_int_test():' )
  print ( '  lines_exp_int() finds intersections of' )
  print ( '  two explicit lines in 2D.' )
  print ( '' )

  for test in range ( 0, 3 ):
#
#  x + 2y - 4 = 0
#  x - y - 1 = 0
#
    if ( test == 0 ):

      p1 = np.array ( [ 0.0,  2.0 ] )
      p2 = np.array ( [ 4.0,  0.0 ] )
      q1 = np.array ( [ 0.0, -1.0 ] )
      q2 = np.array ( [ 1.0,  0.0 ] )
#
#  x + 2y - 4 = 0
#  2x + 4y - 1 = 0
#
    elif ( test == 1 ):

      p1 = np.array ( [ 0.00, 2.00 ] )
      p2 = np.array ( [ 4.00, 0.00 ] )
      q1 = np.array ( [ 0.00, 0.25 ] )
      q2 = np.array ( [ 0.50, 0.00 ] )
#
#  x + 2y - 4 = 0
#  -3x - 6y +12 = 0
#
    elif ( test == 2 ):

      p1 = np.array ( [ 0.0, 2.0 ] )
      p2 = np.array ( [ 4.0, 0.0 ] )
      q1 = np.array ( [ 0.0, 2.0 ] )
      q2 = np.array ( [ 4.0, 0.0 ] )

    print ( '' )
    print ( '  P1  %8f  %8f' % ( p1[0], p1[1] ) )
    print ( '  P2  %8f  %8f' % ( p2[0], p2[1] ) )
    print ( '' )
    print ( '  Q1  %8f  %8f' % ( q1[0], q1[1] ) )
    print ( '  Q2  %8f  %8f' % ( q2[0], q2[1] ) )

    ival, p = lines_exp_int ( p1, p2, q1, q2 )

    if ( ival == 1 ):
      print ( '  Intersection at %8f  %8f' % ( p[0], p[1] ) )
    elif ( ival == 0 ):
      print ( '  Lines are parallel, no intersection.' )
    elif ( ival == 2 ):
      print ( '  Lines are coincident.' )
    else:
      print ( '  Unknown return value of IVAL = %d' % ( ival ) )

  return

def lines_imp_int ( a1, b1, c1, a2, b2, c2 ):

#*****************************************************************************80
#
## lines_imp_int() determines where two implicit lines intersect in 2D.
#
#  Discussion:
#
#    The implicit form of a line in 2D is:
#
#      A * X + B * Y + C = 0
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 December 2010
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real A1, B1, C1, define the first line.
#    At least one of A1 and B1 must be nonzero.
#
#    real A2, B2, C2, define the second line.
#    At least one of A2 and B2 must be nonzero.
#
#  Output:
#
#    integer IVAL, reports on the intersection.
#    -1, both A1 and B1 were zero.
#    -2, both A2 and B2 were zero.
#     0, no intersection, the lines are parallel.
#     1, one intersection point, returned in X, Y.
#     2, infinitely many intersections, the lines are identical.
#
#    real P(2,1), if IVAL = 1, then P is
#    the intersection point.  if IVAL = 2, then P is one of the
#    points of intersection.  Otherwise, P = [].
#
  import numpy as np

  p = np.zeros ( 2 )
#
#  Refuse to handle degenerate lines.
#
  if ( a1 == 0.0 and b1 == 0.0 ):
    ival = -1
    return ival, p
  elif ( a2 == 0.0 and b2 == 0.0 ):
    ival = -2
    return ival, p
#
#  Set up and solve a linear system.
#
  a = np.zeros ( [ 2, 3 ] )

  a[0,0] = a1
  a[0,1] = b1
  a[0,2] = - c1

  a[1,0] = a2
  a[1,1] = b2
  a[1,2] = - c2

  a, info = r8mat_solve ( 2, 1, a )
#
#  If the inverse exists, then the lines intersect at the solution point.
#
  if ( info == 0 ):

    ival = 1
    p[0] = a[0,2]
    p[1] = a[1,2]
#
#  If the inverse does not exist, then the lines are parallel
#  or coincident.  Check for parallelism by seeing if the
#  C entries are in the same ratio as the A or B entries.
#
  else:

    ival = 0

    if ( a1 == 0.0 ):
      if ( b2 * c1 == c2 * b1 ):
        ival = 2
        p[0] = 0.0
        p[1] = - c1 / b1
    else:
      if ( a2 * c1 == c2 * a1 ):
        ival = 2
        if ( abs ( a1 ) < abs ( b1 ) ):
          p[0] = 0.0
          p[1] = - c1 / b1
        else:
          p[0] = - c1 / a1
          p[1] = 0.0

  return ival, p

def lines_imp_int_test ( ):

#*****************************************************************************80
#
## lines_imp_int_test() tests lines_imp_int().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'lines_imp_int_test():' )
  print ( '  lines_imp_int() finds the intersection of' )
  print ( '  two lines written in implicit form.' )
#
#  x + 2y - 4 = 0
#
  a1 =  1.0
  b1 =  2.0
  c1 = -4.0
  print ( '' )
  print ( '  Line 1 coefficients:  %8f  %8f  %8f' % ( a1, b1, c1 ) )
#
#  x - y - 1 = 0
#
  a2 =  1.0
  b2 = -1.0
  c2 = -1.0
  print ( '  Line 2 coefficients:  %8f  %8f  %8f' % (  a2, b2, c2 ) )

  ival, p = lines_imp_int ( a1, b1, c1, a2, b2, c2 )

  if ( ival == 1 ):
    print ( '  Intersection at %8f  %8f' % ( p[0], p[1] ) )
  elif ( ival == 0 ):
    print ( '  Lines are parallel, no intersection.' )
  elif ( ival == 2 ):
    print ( '  Lines are coincident.' )
  else:
    print ( '  Unknown return value of ival = %d' % ( ival ) )
#
#  2x + 4y - 1 = 0
#
  print ( '' )
  print ( '  Line 1 coefficients:  %8f  %8f  %8f' % ( a1, b1, c1 ) )

  a2 =  2.0
  b2 = +4.0
  c2 = -1.0
  print ( '  Line 2 coefficients:  %8f  %8f  %8f' % ( a2, b2, c2 ) )

  ival, p = lines_imp_int (a1, b1, c1, a2, b2, c2 )

  if ( ival == 1 ):
    print ( '  Intersection at %8f  %8f' % ( p[0], p[1] ) )
  elif ( ival == 0 ):
    print ( '  Lines are parallel, no intersection.' )
  elif ( ival == 2 ):
    print ( '  Lines are coincident.' )
  else:
    print ( '  Unknown return value of ival = %d' % ( ival ) )
#
#  -3x - 6y +12 = 0
#
  print ( '' )
  print ( '  Line 1 coefficients:  %8f  %8f  %8f' % ( a1, b1, c1 ) )

  a2 =  -3.0
  b2 =  -6.0
  c2 = +12.0
  print ( '  Line 2 coefficients:  %8f  %8f  %8f' % ( a2, b2, c2 ) )
 
  ival, p = lines_imp_int  ( a1, b1, c1, a2, b2, c2 )

  if ( ival == 1 ):
    print ( '  Intersection at %8f  %8f' % ( p[0], p[1] ) )
  elif ( ival == 0 ):
    print ( '  Lines are parallel, no intersection.' )
  elif ( ival == 2 ):
    print ( '  Lines are coincident.' )
  else:
    print ( '  Unknown return value of ival = %d' % ( ival ) )

  return

def polygon_area ( n, v ):

#*****************************************************************************80
#
## polygon_area() computes the area of a polygon.
#
#  Discussion:
#
#    AREA = 1/2 * abs ( sum ( 1 <= I <= N ) X(I) * ( Y(I+1) - Y(I-1) ) )
#    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
#
#    If the vertices are given in counterclockwise order, the area
#    will be positive.  If the vertices are given in clockwise order,
#    the area will be negative.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of vertices of the polygon.
#
#    real V(2,N), the vertices.
#
#  Output:
#
#    real AREA, the area of the polygon.
#
  area = 0.0

  for i in range ( 0, n ):

    im1 = i4_wrap ( i - 1, 0, n - 1 )
    ip1 = i4_wrap ( i + 1, 0, n - 1 )

    area = area + v[0,i] * ( v[1,ip1] - v[1,im1] )

  area = 0.5 * area

  return area

def polygon_area_test ( ):

#*****************************************************************************80
#
## polygon_area_test() tests polygon_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  test_num = 2
  area_exact_test = np.array ( [ 2.0, 6.0 ] )
  n_test = np.array ( [ 4, 8 ] )

  print ( '' )
  print ( 'polygon_area_test():' )
  print ( '  polygon_area() computes the area of a polygon.' )

  for test in range ( 0, test_num ):

    n = n_test[test]
    area_exact = area_exact_test[test]

    if ( test == 0 ):

      v = np.array ( [ \
        [ 1.0, 2.0, 1.0, 0.0 ], \
        [ 0.0, 1.0, 2.0, 1.0 ] ] )

    elif ( test == 1 ):

      v = np.array ( [ \
        [ 0.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0, 0.0 ], \
        [ 0.0, 0.0, 3.0, 3.0, 1.0, 1.0, 2.0, 2.0 ] ] )

    print ( '' )
    print ( '  Number of polygonal vertices = %d' % ( n ) )

    r8mat_transpose_print ( 2, n, v, '  The polygon vertices:' )

    area = polygon_area ( n, v )

    print ( '' )
    print ( '  Exact area is        %g' % ( area_exact ) )
    print ( '  The computed area is %g' % ( area ) )

  return

def polygon_area_3d ( n, v ):

#*****************************************************************************80
#
## polygon_area_3d() computes the area of a polygon in 3D.
#
#  Discussion:
#
#    The computation is not valid unless the vertices of the polygon
#    lie in a plane, so that the polygon that is defined is "flat".
#
#    The polygon does not have to be "regular", that is, neither its
#    sides nor its angles need to be equal.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 May 2022
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Allen Van Gelder,
#    Efficient Computation of Polygon Area and Polyhedron Volume,
#    Graphics Gems V,
#    edited by Alan Paeth,
#    AP Professional, 1995.
#
#  Input:
#
#    integer N, the number of vertices.
#
#    real V(3,N), the coordinates of the vertices.
#    The vertices should be listed in neighboring order.
#
#  Output:
#
#    real AREA, the area of the polygon.
#
#    real NORMAL[3], the unit normal vector to the polygon.
#
  import numpy as np

  normal = np.zeros ( 3 )

  for i in range ( 0, n ):
#
#  When i = n - 1, ip1 = 0
#
    ip1 = ( i + 1 % n )
#
#  Compute the cross product vector.
#
    cross[0] = v[i,1] * v[ip1,2] - v[i,2] * v[ip1,1]
    cross[1] = v[i,2] * v[ip1,0] - v[i,0] * v[ip1,2]
    cross[2] = v[i,0] * v[ip1,1] - v[i,1] * v[ip1,0]

    normal = normal + cross

  area = np.linalg.norm ( normal )

  if ( area != 0.0 ):
    normal = normal / area
  else:
    normal = 1.0 / np.sqrt ( 3.0 )

  area = 0.5 * area

  return area, normal

def triangle_area_2 ( xa, ya, xb, yb, xc, yc ):

#*****************************************************************************80
#
## triangle_area_2() computes the signed area of a triangle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real XA, YA, XB, YB, XC, YC, the coordinates of
#    the vertices of the triangle, given in counterclockwise order.
#
#  Output:
#
#    real VALUE, the signed area of the triangle.
#
  value = 0.5 * ( ( xb - xa ) * ( yc - ya ) - ( xc - xa ) * ( yb - ya ) )

  return value

def triangle_3d_area_cross_product ( v ):

#*****************************************************************************80
#
## triangle_3d_area_cross_product() computes the area of a triangle in 3D.
#
#  Discussion:
#
#    This routine uses the fact that the norm of the cross product vector
#    is the area of the parallelogram they form.  The triangle they
#    form has half that area.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    30 January 2005
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer and John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real V(3,3), the triangle vertices.  
#    The first row contains the coordinates of the first vertex.
#
#  Output:
#
#    real AREA, the area of the triangle.
#
  import numpy as np

  cross = np.zeros ( 3 )

  cross[0] = ( v[1,1] - v[0,1] ) * ( v[2,2] - v[0,2] ) \
           - ( v[1,2] - v[0,2] ) * ( v[2,1] - v[0,1] )

  cross[1] = ( v[1,2] - v[0,2] ) * ( v[2,0] - v[0,0] ) \
           - ( v[1,0] - v[0,0] ) * ( v[2,2] - v[0,2] )

  cross[2] = ( v[1,0] - v[0,0] ) * ( v[2,1] - v[0,1] ) \
           - ( v[1,1] - v[0,1] ) * ( v[2,0] - v[0,0] )

  area = 0.5 * np.linalg.norm ( cross )

  return area

def triangle_3d_area_cross_product_test ( ):

#*****************************************************************************80
#
## triangle_3d_area_cross_product_test() tests triangle_3d_area_cross_product();
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2022
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 1.0,       2.0,       3.0 ], \
    [ 2.4142137, 3.4142137, 3.0 ], \
    [ 1.7071068, 2.7071068, 4.0 ] ] )

  print ( '' )
  print ( 'triangle_3d_area_cross_product_test():' )
  print ( '  triangle_3d_area_cross_product() computes the area' )
  print ( '  of a triangle in 3D.' )
  print ( '' )
  print ( '  Triangle vertices:' )
  print ( t )

  area = triangle_3d_area_cross_product ( t )

  print ( '' )
  print ( '  Computed area is ', area )

  return

def triangle_3d_area_heron ( v ):

#*****************************************************************************80
#
## triangle_3d_area_heron() computes the area of a triangle in 3D.
#
#  Discussion:
#
#    This routine uses Heron's formula.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 February 2005
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer and John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real V(3,3), the triangle vertices.  
#    The first row contains the coordinates of the first vertex.
#
#  Output:
#
#    real AREA, the area of the triangle.
#
  import numpy as np

  s = np.zeros ( 3 )

  for j in range ( 0, 3 ):
    s[j] = np.linalg.norm ( v[j-1] - v[j] )

  area = (   s[0] + s[1] + s[2] ) \
       * ( - s[0] + s[1] + s[2] ) \
       * (   s[0] - s[1] + s[2] ) \
       * (   s[0] + s[1] - s[2] )

  area = 0.25 * np.sqrt ( area )

  return area

def triangle_3d_area_heron_test ( ):

#*****************************************************************************80
#
## triangle_3d_area_heron_test() tests triangle_3d_area_heron().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2022
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 1.0,       2.0,       3.0 ], \
    [ 2.4142137, 3.4142137, 3.0 ], \
    [ 1.7071068, 2.7071068, 4.0 ] ] )

  print ( '' )
  print ( 'triangle_3d_area_heron_test():' )
  print ( '  triangle_3d_area_heron() computes the area' )
  print ( '  of a triangle in 3D.' )
  print ( '' )
  print ( '  Triangle vertices:' )
  print ( t )

  area = triangle_3d_area_heron ( t )
  print ( '' )
  print ( '  Computed area is ', area )

  return

def pyramid01_volume ( ):

#*****************************************************************************80
#
## pyramid01_volume() returns the volume of a unit pyramid.
#
#  Discussion:
#
#    A pyramid with square base can be regarded as the upper half of a
#    3D octahedron.
#
#    The integration region:
#
#      - ( 1 - Z ) <= X <= 1 - Z
#      - ( 1 - Z ) <= Y <= 1 - Z
#                0 <= Z <= 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    Output, real VALUE, the volume of the pyramid.
#
  value = 4.0 / 3.0

  return value

def pyramid01_volume_test ( ) :

#*****************************************************************************80
#
## pyramid01_volume_test() tests pyramid01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'pyramid01_volume_test():' )
  print ( '  pyramid01_volume() returns the volume of the unit pyramid.' )

  value = pyramid01_volume ( )

  print ( '' )
  print ( '  pyramid01_volume() = %g' % ( value ) )

  return

def pyramid_volume ( r, h ):

#*****************************************************************************80
#
## pyramid_volume() returns the volume of a pyramid.
#
#  Discussion:
#
#    A pyramid with square base can be regarded as the upper half of a
#    3D octahedron.
#
#    The integration region:
#
#      - ( R - Z ) <= X <= R - Z
#      - ( R - Z ) <= Y <= R - Z
#                0 <= Z <= H.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the "radius" of the pyramid, that is, half the
#    length of one of the sides of the square base.
#
#    real H, the height of the pyramid.
#
#  Output:
#
#    real VALUE, the volume of the pyramid.
#
  value = ( 4.0 / 3.0 ) * h * r * r

  return value

def pyramid_volume_test ( ) :

#*****************************************************************************80
#
## pyramid_volume_test() tests pyramid_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 June 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'pyramid_volume_test():' )
  print ( '  pyramid_volume() returns the volume of a pyramid.' )

  print ( '' )
  print ( '     Radius     Height     Volume' )
  print ( '' )

  for i in range ( 0, 5 ):
    r = 1.0 + 9.0 * np.random.rand ( )
    h = 1.0 + 9.0 * np.random.rand ( )
    volume = pyramid_volume ( r, h )
    print ( '  %8.4f  %8.4f  %8.4f' % ( r, h, volume ) )

  return

def r8_acos ( c ):

#*****************************************************************************80
#
## r8_acos() computes the arc cosine function, with argument truncation.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 July 2014
#
#  Author:
#
#    Original FORTRAN77 version by Wayne Fullerton.
#    Python version by John Burkardt.
#
#  Input:
#
#    real C, the argument.
#
#  Output:
#
#    real VALUE, the arc-cosine of C.
#
  import numpy as np

  c2 = max ( c,  - 1.0 )
  c2 = min ( c2, +1.0 )
  
  value = np.arccos ( c2 )

  return value

def r8_acos_test ( ):

#*****************************************************************************80
#
## r8_acos_test() tests r8_acos().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    08 November 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
 
  print ( '' )
  print ( 'r8_acos_test():' )
  print ( '  r8_acos() computes the arc-cosine of an angle.' )
  print ( '' )
  print ( '       C            r8_acos(C)        ACOS(C)' )
  print ( '' )

  for test in range ( -1, 14 ):

    c = float ( test - 6 ) / 6.0

    if ( -1.0 <= c and c <= 1.0 ):
      print ( '  %14.6g  %14.6g  %14.6g' % ( c, r8_acos ( c ), np.arccos ( c ) ) )
    else:
      print ( '  %14.6g  %14.6g' % ( c, r8_acos ( c ) ) )

  return

def r8mat_det_4d ( a ):

#*****************************************************************************80
#
## r8mat_det_4d() computes the determinant of a 4 by 4 matrix.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real A(4,4), the matrix whose determinant is desired.
#
#  Output:
#
#    real VALUE, the determinant of the matrix.
#
  value = \
      a[0,0] * ( \
        a[1,1] * ( a[2,2] * a[3,3] - a[2,3] * a[3,2] ) \
      - a[1,2] * ( a[2,1] * a[3,3] - a[2,3] * a[3,1] ) \
      + a[1,3] * ( a[2,1] * a[3,2] - a[2,2] * a[3,1] ) ) \
    - a[0,1] * ( \
        a[1,0] * ( a[2,2] * a[3,3] - a[2,3] * a[3,2] ) \
      - a[1,2] * ( a[2,0] * a[3,3] - a[2,3] * a[3,0] ) \
      + a[1,3] * ( a[2,0] * a[3,2] - a[2,2] * a[3,0] ) ) \
    + a[0,2] * ( \
        a[1,0] * ( a[2,1] * a[3,3] - a[2,3] * a[3,1] ) \
      - a[1,1] * ( a[2,0] * a[3,3] - a[2,3] * a[3,0] ) \
      + a[1,3] * ( a[2,0] * a[3,1] - a[2,1] * a[3,0] ) ) \
    - a[0,3] * ( \
        a[1,0] * ( a[2,1] * a[3,2] - a[2,2] * a[3,1] ) \
      - a[1,1] * ( a[2,0] * a[3,2] - a[2,2] * a[3,0] ) \
      + a[1,2] * ( a[2,0] * a[3,1] - a[2,1] * a[3,0] ) )

  return value

def r8mat_det_4d_test ( ):

#*****************************************************************************80
#
## r8mat_det_4d_test() tests r8mat_det_4d()
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 April 2009
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 4

  x = np.array ( [ 1.0, 10.0, 4.0, 2.0 ] )

  print ( '' )
  print ( 'r8mat_det_4d_test():' )
  print ( '  r8mat_det_4d() computes the determinant of a 4 by 4 matrix' )

  a = r8mat_vand2 ( n, n, x )
  det = r8mat_det_4d ( a )

  r8mat_print ( n, n, a, '  Matrix:' )

  print ( '' )
  print ( '  r8mat_det_4d computes determinant: %g' % ( det ) )
#
#  Special formula for the determinant of a Vandermonde matrix:
#
  det = 1.0
  for i in range ( 0, n ):
    for j in range ( 0, i ):
      det = det * ( x[i] - x[j] )

  print ( '  Exact determinant is %g' % ( det ) )

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
#    08 May 2020
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
#    08 May 2020
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

def r8mat_solve ( n, nrhs, a ):

#*****************************************************************************80
#
## r8mat_solve() uses Gauss-Jordan elimination to solve an N by N linear system.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    integer NRHS, the number of right hand sides.  NRHS
#    must be at least 0.
#
#    real A(N,N+NRHS), contains in rows and
#    columns 1 to N the coefficient matrix, and in columns N+1 through
#    N+NRHS, the right hand sides.
#
#  Output:
#
#    real A(N,N+NRHS), the coefficient matrix
#    area has been destroyed, while the right hand sides have
#    been overwritten with the corresponding solutions.
#
#    integer INFO, singularity flag.
#    0, the matrix was not singular, the solutions were computed.
#    J, factorization failed on step J, and the solutions could not
#    be computed.
#
  info = 0

  for j in range ( 0, n ):
#
#  Choose a pivot row IPIVOT.
#
    ipivot = j
    apivot = a[j,j]

    for i in range ( j + 1, n ):
      if ( abs ( apivot ) < abs ( a[i,j] ) ):
        apivot = a[i,j]
        ipivot = i

    if ( apivot == 0.0 ):
      info = j
      return a, info
#
#  Interchange.
#
    for k in range ( 0, n + nrhs ):
      temp        = a[ipivot,k]
      a[ipivot,k] = a[j,k]
      a[j,k]      = temp
#
#  A(J,J) becomes 1.
#
    a[j,j] = 1.0
    for k in range ( j + 1, n + nrhs ):
      a[j,k] = a[j,k] / apivot
#
#  A(I,J) becomes 0.
#
    for i in range ( 0, n ):

      if ( i != j ):

        factor = a[i,j]
        a[i,j] = 0.0
        for k in range ( j + 1, n + nrhs ):
          a[i,k] = a[i,k] - factor * a[j,k]

  return a, info

def r8mat_solve_test ( ):

#*****************************************************************************80
#
## r8mat_solve_test() tests r8mat_solve().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 February 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 3
  rhs_num = 2
#
#  Each row of this definition is a COLUMN of the matrix.
#
  a = np.array ( [ \
    [ 1.0, 2.0, 3.0, 14.0,  7.0 ], \
    [ 4.0, 5.0, 6.0, 32.0, 16.0 ], \
    [ 7.0, 8.0, 0.0, 23.0,  7.0 ] ] )

  print ( '' )
  print ( 'r8mat_solve_test():' )
  print ( '  r8mat_solve() solves linear systems.' )
#
#  Print out the matrix to be inverted.
#
  r8mat_print ( n, n + rhs_num, a, '  The linear system:' )
#
#  Solve the systems.
#
  a, info = r8mat_solve ( n, rhs_num, a )
 
  if ( info != 0 ):
    print ( '' )
    print ( '  The input matrix was singular.' )
    print ( '  The solutions could not be computed.' )
    return

  r8mat_print ( n, n + rhs_num, a, '  Factored matrix and solutions:' )

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

  print ( '' )
  print ( 'r8mat_transpose_print_test():' )
  print ( '  r8mat_transpose_print() prints an R8MAT.' )

  m = 4
  n = 3
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0 ], 
    [ 21.0, 22.0, 23.0 ], 
    [ 31.0, 32.0, 33.0 ], 
    [ 41.0, 42.0, 43.0 ] ], dtype = np.float64 )
  r8mat_transpose_print ( m, n, v, '  Here is an R8MAT, transposed:' )

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

  print ( '' )
  print ( 'r8mat_transpose_print_some_test():' )
  print ( '  r8mat_transpose_print_some() prints some of an R8MAT, transposed.' )

  m = 4
  n = 6
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ], 
    [ 21.0, 22.0, 23.0, 24.0, 25.0, 26.0 ], 
    [ 31.0, 32.0, 33.0, 34.0, 35.0, 36.0 ], 
    [ 41.0, 42.0, 43.0, 44.0, 45.0, 46.0 ] ], dtype = np.float64 )
  r8mat_transpose_print_some ( m, n, v, 0, 3, 2, 5, '  R8MAT, rows 0:2, cols 3:5:' )

  return

def r8mat_vand2 ( m, n, x ):

#*****************************************************************************80
#
## r8mat_vand2() returns the N by N row Vandermonde matrix A.
#
#  Discussion:
#
#    The row Vandermonde matrix returned by this routine reads "across"
#    rather than down.  In particular, each row begins with a 1, followed by
#    some value X, followed by successive powers of X.
#
#  Formula:
#
#    A(I,J) = X(I)^(J-1)
#
#  Properties:
#
#    A is nonsingular if, and only if, the X values are distinct.
#
#    The determinant of A is
#
#      det(A) = product ( 2 <= I <= N ) (
#        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
#
#    The matrix A is generally ill-conditioned.
#
#  Example:
#
#    N = 5, X = (2, 3, 4, 5, 6)
#
#    1 2  4   8   16
#    1 3  9  27   81
#    1 4 16  64  256
#    1 5 25 125  625
#    1 6 36 216 1296
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 February 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, N, the order of the matrix desired.
#
#    real X(M), the values that define A.
#
#  Output:
#
#    real A(M,N), the M by N row Vandermonde matrix.
#
  import numpy as np

  a = np.zeros ( [ m, n ] )
  
  for i in range ( 0, m ):
    a[i,0] = 1.0
    for j in range ( 1, n ):
      a[i,j] = a[i,j-1] * x[i]

  return a

def r8mat_vand2_test ( ):

#*****************************************************************************80
#
## r8mat_vand2_test() tests r8mat_vand2().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 February 2015
#
#  Author:
#
#    John Burkardt
#
  m = 5
  n = 4

  print ( '' )
  print ( 'r8mat_vand2_test():' )
  print ( '  r8mat_vand2() returns a row Vandermonde matrix.' )

  x = r8vec_indicator1 ( m )
  print ( '' )
  print ( '  factor vector x:' )
  print ( x )
 
  a = r8mat_vand2 ( m, n, x )
  r8mat_print ( m, n, a, '  The row Vandermonde matrix:' )

  return

def r8po_det ( n, a_lu ):

#*****************************************************************************80
#
## r8po_det() computes the determinant of a matrix factored by r8po_fa.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A_LU(N,N), the factor from r8po_fa.
#
#    real VALUE, the determinant of A.
#
  value = 1.0

  for i in range ( 0, n ):
    value = value * a_lu[i,i] ** 2

  return value

def r8po_det_test ( ):

#*****************************************************************************80
#
## r8po_det_test() tests r8po_det().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 5

  print ( '' )
  print ( 'r8po_det_test():' )
  print ( '  r8po_det() finds the determinant of a positive definite symmetric' )
  print ( '  matrix after it has been factored.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )

  a = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      a[i,j] = min ( i, j ) + 1

  r8po_print ( n, a, '  The matrix A:' )
#
#  Get R, the Cholesky factor of A.
#
  r = r8po_fa ( n, a )
#
#  Get the determinant of A.
#
  value = r8po_det ( n, r )
  print ( '' )
  print ( '  Determinant of A = %g' % ( value ) )

  return

def r8po_dif2 ( n ):

#*****************************************************************************80
#
## r8po_dif2() returns the second difference matrix in R8PO format.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of rows and columns of the matrix.
#    N must be positive.
#
#    real A(N,N), the matrix.
#
  a = r8po_zeros ( n )

  for i in range ( 0, n ):
    a[i,i] = 2.0
  
  for i in range ( 0, n - 1 ):
    a[i,i+1] = -1.0

  return a

def r8po_dif2_test ( ):

#*****************************************************************************80
#
## r8po_dif2_test() tests r8po_dif2().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2015
#
#  Author:
#
#    John Burkardt
#
  n = 5

  print ( '' )
  print ( 'r8po_dif2_test():' )
  print ( '  r8po_dif2() returns the second difference matrix in R8PO format.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )
#
#  Set the matrix.
#
  a = r8po_dif2 ( n )

  r8po_print ( n, a, '  The matrix:' )

  return

def r8po_fa ( n, a ):

#*****************************************************************************80
#
## r8po_fa() factors a R8PO matrix.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#    The positive definite symmetric matrix A has a Cholesky factorization
#    of the form:
#
#      A = R' * R
#
#    where R is an upper triangular matrix with positive elements on
#    its diagonal.  This routine overwrites the matrix A with its
#    factor R.
#
#    This function failed miserably when I wrote "r = a", because of a
#    disastrously misconceived feature of Python, which does not copy
#    one matrix to another, but makes them both point to the same object.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2015
#
#  Author:
#
#    John Burkardt.
#
#  Reference:
#
#    Dongarra, Bunch, Moler, Stewart,
#    LINPACK User's Guide,
#    SIAM, 1979.
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A(N,N), the matrix in R8PO storage.
#
#    real R(N,N), the Cholesky factor R in R8GE storage.
#
#    integer INFO, error flag.
#    0, normal return.
#    K, error condition.  The principal minor of order K is not
#    positive definite, and the factorization was not completed.
#
  import numpy as np

  r = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( i, n ):
      r[i,j] = a[i,j]

  for j in range ( 0, n ):

    for k in range ( 0, j ):
      t = 0.0
      for i in range ( 0, k ):
        t = t + r[i,k] * r[i,j]
      r[k,j] = ( r[k,j] - t ) / r[k,k]

    t = 0.0
    for i in range ( 0, j ):
      t = t + r[i,j] ** 2

    s = r[j,j] - t

    if ( s <= 0.0 ):
      print ( '' )
      print ( 'r8po_fa - Fatal error!' )
      print ( '  Factorization fails on column %d.' % ( j ) )
      raise Exception ( 'r8po_fa - Fatal error!' )

    r[j,j] = np.sqrt ( s )
#
#  Since the Cholesky factor is stored in R8GE format, be sure to
#  zero out the lower triangle.
#
  for i in range ( 0, n ):
    for j in range ( 0, i ):
      r[i,j] = 0.0

  return r

def r8po_fa_test ( ):

#*****************************************************************************80
#
## r8po_fa_test() tests r8po_fa().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 5

  print ( '' )
  print ( 'r8po_fa_test():' )
  print ( '  r8po_fa() factors a positive definite symmetric' )
  print ( '  linear system,' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )

  a = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      a[i,j] = min ( i, j ) + 1

  r8po_print ( n, a, '  The matrix A:' )
#
#  Factor the matrix.
#
  r = r8po_fa ( n, a )

  r8ut_print ( n, n, r, '  The factor R (a R8UT matrix):' )
#
#  Compute the product R' * R.
#
  rtr = r8ut_mtm ( n, r, r )

  r8ge_print ( n, n, rtr, '  The product R\' * R:' )

  return

def r8po_indicator ( n ):

#*****************************************************************************80
#
## r8po_indicator() sets up a R8PO indicator matrix.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of rows and columns of the matrix.
#    N must be positive.
#
#    real A(N,N), the R8PO matrix.
#
  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  a = r8po_zeros ( n )

  for i in range ( 0, n ):
    for j in range ( i, n ):
      a[i,j] = float ( fac * ( i + 1 ) + ( j + 1 ) )

  return a

def r8po_indicator_test ( ):

#*****************************************************************************80
#
## r8po_indicator_test() tests r8po_indicator().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
  n = 5

  print ( '' )
  print ( 'r8po_indicator_test():' )
  print ( '  r8po_indicator() sets up an R8PO indicator matrix.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )
#
#  Set the matrix.
#
  a = r8po_indicator ( n )

  r8po_print ( n, a, '  The R8PO indicator matrix:' )

  return

def r8po_inverse ( n, r ):

#*****************************************************************************80
#
## r8po_inverse() computes the inverse of a matrix factored by r8po_fa.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real R(N,N), the Cholesky factor, in R8GE storage, returned by r8po_fa.
#
#    real B(N,N), the inverse matrix, in R8PO storage.
#
  import numpy as np

  b = r8po_zeros ( n )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      b[i,j] = r[i,j]
#
#  Compute Inverse ( R ).
#
  for k in range ( 0, n ):

    b[k,k] = 1.0 / b[k,k]
    for i in range ( 0, k ):
      b[i,k] = - b[i,k] * b[k,k]

    for j in range ( k + 1, n ):
      t = b[k,j]
      b[k,j] = 0.0
      for i in range ( 0, k + 1 ):
        b[i,j] = b[i,j] + t * b[i,k]
#
#  Compute Inverse ( R ) * ( Inverse ( R ) )'.
#
  for j in range ( 0, n ):

    for k in range ( 0, j ):
      t = b[k,j]
      for i in range ( 0, k + 1 ):
        b[i,k] = b[i,k] + t * b[i,j]

    t = b[j,j]
    for k in range ( 0, j + 1 ):
      b[k,j] = b[k,j] * t

  return b

def r8po_inverse_test ( ):

#*****************************************************************************80
#
## r8po_inverse_test() tests r8po_inverse().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 July 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 4

  print ( '' )
  print ( 'r8po_inverse_test():' )
  print ( '  r8po_inverse() computes the inverse of' )
  print ( '  a symmetric positive definite matrix' )
  print ( '  factored by r8po_fa.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )
#
#  Set the matrix.
#
  a = r8po_zeros ( n )

  for i in range ( 0, n ):
    for j in range ( i, n ):
      a[i,j] = min ( i, j ) + 1

  r8po_print ( n, a, '  Matrix A:' )
#
#  Factor the matrix.
#
  r = r8po_fa ( n, a )
#
#  Compute the inverse.
#
  b = r8po_inverse ( n, r )

  r8po_print ( n, b, '  Inverse matrix B:' )
#
#  Check.
#
  c = r8po_mm ( n, a, b )

  r8po_print ( n, c, '  Product A * B:' )

  return

def r8po_ml ( n, r, x ):

#*****************************************************************************80
#
## r8po_ml() computes A * x = b after A has been factored by r8po_fa.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real R(N,N), the Cholesky factor, in R8GE format, returned by r8po_fa.
#
#    real X(N), the vector to be multiplied by A.
#
#    real B(N), the product A * x.
#
  import numpy as np

  b = np.zeros ( n )
#
#  Compute R * x = y.
#
  for i in range ( 0, n ):
    b[i] = r[i,i] * x[i]
    for j in range ( i + 1, n ):
      b[i] = b[i] + r[i,j] * x[j]
#
#  Compute R' * y = b.
#
  for i in range ( n - 1, -1, -1 ):
    b[i] = r[i,i] * b[i]
    for j in range ( 0, i ):
      b[i] = b[i] + b[j] * r[j,i]

  return b

def r8po_ml_test ( ):

#*****************************************************************************80
#
## r8po_ml_test() tests r8po_ml().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 August 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 10

  print ( '' )
  print ( 'r8po_ml_test():' )
  print ( '  r8po_ml() computes A*x for an R8PO matrix A' )
  print ( '  even after it has been factored by r8po_fa.' )

  a = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      a[i,j] = min ( i, j ) + 1
#
#  Set the desired solution.
#
  x = r8vec_indicator1 ( n )
#
#  Compute the corresponding right hand side.
#
  b = r8po_mv ( n, a, x )
#
#  Factor the matrix.
#
  a_lu = r8po_fa ( n, a )
#
#  Solve the linear system.
#
  x = r8po_sl ( n, a_lu, b )
  print ( '' )
  print ( '  solution x:' )
  print ( x )
#
#  Set the desired solution.
#
  for i in range ( 0, n ):
    x[i] = 1.0
#
#  Compute the corresponding right hand side, using the factored matrix.
#
  b = r8po_ml ( n, a, x )
#
#  Solve the linear system.
#
  x = r8po_sl ( n, a, b )
 
  print ( '' )
  print ( '  solution x:' )
  print ( x )

  return

def r8po_mm ( n, a, b ):

#*****************************************************************************80
#
## r8po_mm() multiplies two R8PO matrices.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrices.
#    N must be positive.
#
#    real A(N,N), B(N,N), the factors.
#
#    real C(N,N), the product.
#
  import numpy as np

  c = r8po_zeros ( n )
  
  for i in range ( 0, n ):

    for j in range ( i, n ):
      for k in range ( 0, n ):

        if ( i <= k ):
          aik = a[i,k]
        else:
          aik = a[k,i]

        if ( k <= j ):
          bkj = b[k,j]
        else:
          bkj = b[j,k]

        c[i,j] = c[i,j] + aik * bkj

  return c

def r8po_mm_test ( ):

#*****************************************************************************80
#
## r8po_mm_test() tests r8po_mm().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 July 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 5

  print ( '' )
  print ( 'r8po_mm_test():' )
  print ( '  r8po_mm() computes the product of two' )
  print ( '  symmetric positive definite matrices.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )
#
#  Set (the upper half of) matrix A.
#
  a = r8po_zeros ( n )

  for i in range ( 0, n ):
    a[i,i] = 2.0
  for i in range ( 0, n - 1 ):
    a[i,i+1] = -1.0

  r8po_print ( n, a, '  Matrix A:' )
#
#  Set (the upper half of) matrix B.
#
  b = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    b[i,i] = float ( i + i + 1 )
  for i in range ( 0, n - 1 ):
    b[i,i+1] = float ( i + i + 1 + 1 )

  r8po_print ( n, b, '  Matrix B:' )
#
#  Compute the product.
#
  c = r8po_mm ( n, a, b )

  r8po_print ( n, c, '  Product matrix C = A * B:' )

  return

def r8po_mv ( n, a, x ):

#*****************************************************************************80
#
## r8po_mv() multiplies a R8PO matrix times a vector.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A(N,N), the R8PO matrix.
#
#    real X(N), the vector to be multiplied by A.
#
#    real B(N), the product A * x.
#
  import numpy as np

  b = np.zeros ( n )

  for i in range ( 0, n ):
    for j in range ( 0, i ):
      b[i] = b[i] + a[j,i] * x[j]
    for j in range ( i, n ):
      b[i] = b[i] + a[i,j] * x[j]

  return b

def r8po_mv_test ( ):

#*****************************************************************************80
#
## r8po_mv_test() tests r8po_mv().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
  n = 5

  print ( '' )
  print ( 'r8po_mv_test():' )
  print ( '  r8po_mv() computes the product of an R8PO matrix and a vector.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )
#
#  Set (the upper half of) matrix A.
#
  a = r8po_zeros ( n )

  for i in range ( 0, n ):
    a[i,i] = 2.0
  for i in range ( 0, n - 1 ):
    a[i,i+1] = -1.0

  r8po_print ( n, a, '  Matrix A:' )
#
#  Set the vector V.
#
  v = r8vec_indicator1 ( n )

  print ( '' )
  print ( '  vector v:' )
  print ( v )
#
#  Compute the product.
#
  w = r8po_mv ( n, a, v )

  print ( '' )
  print ( '  product w = A*v:' )
  print ( w )

  return

def r8po_print ( n, a, title ):

#*****************************************************************************80
#
## r8po_print() prints a R8PO matrix.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of an SPO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A(N,N), the R8PO matrix.
#
#    string TITLE, a title to be printed.
#
  r8po_print_some ( n, a, 0, 0, n - 1, n - 1, title )

  return

def r8po_print_test ( ):

#*****************************************************************************80
#
## r8po_print_test() tests r8po_print().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'r8po_print_test():' )
  print ( '  r8po_print() prints an R8PO matrix.' )

  n = 5
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0 ], 
    [ 12.0, 22.0, 23.0, 24.0, 25.0, ], 
    [ 13.0, 23.0, 33.0, 34.0, 35.0 ], 
    [ 14.0, 24.0, 34.0, 44.0, 45.0 ],
    [ 14.0, 25.0, 35.0, 45.0, 55.0 ] ], dtype = np.float64 )

  r8po_print ( n, v, '  Here is an R8PO matrix:' )

  return

def r8po_print_some ( n, a, ilo, jlo, ihi, jhi, title ):

#*****************************************************************************80
#
## r8po_print_some() prints some of a R8PO matrix.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A(N,N), the R8PO matrix.
#
#    integer ILO, JLO, IHI, JHI, the first row and
#    column, and the last row and column to be printed.
#
#    string TITLE, a title.
#
  incx = 5

  print ( '' )
  print ( title )

  if ( n <= 0 ):
    print ( '' )
    print ( '  (None)' )
    return
#
#  Print the columns of the matrix, in strips of 5.
#
  for j2lo in range ( max ( jlo, 0 ), min ( jhi + 1, n ), incx ):

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n - 1 )
    j2hi = min ( j2hi, jhi )
    
    print ( '' )
    print ( '  Col: ', end = '' )

    for j in range ( j2lo, j2hi + 1 ):
      print ( '%7d       ' % ( j ), end = '' )

    print ( '' )
    print ( '  Row' )

    i2lo = max ( ilo, 0 )
    i2hi = min ( ihi, n - 1 )

    for i in range ( i2lo, i2hi + 1 ):

      print ( '%7d :' % ( i ), end = '' )
      
      for j in range ( j2lo, j2hi + 1 ):
        if ( i <= j ):
          print ( '%12g  ' % ( a[i,j] ), end = '' )
        else:
          print ( '%12g  ' % ( a[j,i] ), end = '' )

      print ( '' )

  return

def r8po_print_some_test ( ):

#*****************************************************************************80
#
## r8po_print_some_test() tests r8po_print_some().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 July 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'r8po_print_some_test():' )
  print ( '  r8po_print_some() prints some of an R8PO matrix.' )

  n = 5
  v = np.array ( [ \
    [ 11.0, 12.0, 13.0, 14.0, 15.0 ], 
    [ 12.0, 22.0, 23.0, 24.0, 25.0, ], 
    [ 13.0, 23.0, 33.0, 34.0, 35.0 ], 
    [ 14.0, 24.0, 34.0, 44.0, 45.0 ],
    [ 14.0, 25.0, 35.0, 45.0, 55.0 ] ], dtype = np.float64 )

  r8po_print_some ( n, v, 0, 3, 3, 4, '  Here is an R8PO matrix:' )

  return

def r8po_random ( n ):

#*****************************************************************************80
#
## r8po_random() randomizes a R8PO matrix.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#    The matrix computed here is not simply a set of random numbers in
#    the nonzero slots of the R8PO array.  It is also a positive definite
#    matrix.  It is computed by setting a "random" upper triangular
#    Cholesky factor R, and then computing A = R'*R.
#    The randomness is limited by the fact that all the entries of
#    R will be between 0 and 1.  A truly random R is only required
#    to have positive entries on the diagonal.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#    N must be positive.
#
#  Output:
#
#    real A(N,N), the R8PO matrix.
#
  import numpy as np

  a = r8po_zeros ( n )

  for i in range ( n - 1, -1, -1 ):
#
#  Set row I of R.
#
    for j in range ( i, n ):
      r = np.random.rand ( )
      a[i,j] = r
#
#  Consider element J of row I, last to first.
#
    for j in range ( n - 1, i - 1, -1 ):
#
#  Add multiples of row I to lower elements of column J.
#
      for i2 in range ( i + 1, j + 1 ):
        a[i2,j] = a[i2,j] + a[i,i2] * a[i,j]
#
#  Reset element J.
#
      a[i,j] = a[i,i] * a[i,j]

  return a

def r8po_random_test ( ):

#*****************************************************************************80
#
## r8po_random_test() tests r8po_random().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
  n = 5

  print ( '' )
  print ( 'r8po_random_test():' )
  print ( '  r8po_random() computes a random positive definite' )
  print ( '  symmetric matrix.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )

  a = r8po_random ( n )

  r8po_print ( n, a, '  The random R8PO matrix:' )

  return

def r8po_sl ( n, r, b ):

#*****************************************************************************80
#
## r8po_sl() solves a R8PO system factored by r8po_fa.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2015
#
#  Author:
#
#    John Burkardt.
#
#  Reference:
#
#    Dongarra, Bunch, Moler, Stewart,
#    LINPACK User's Guide,
#    SIAM, 1979.
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real R(N,N), the Cholesky factor, in R8GE storage, 
#    returned by r8po_fa.
#
#    real B(N), the right hand side.
#
#    real X(N), the solution vector.
#
  import numpy as np

  x = np.zeros ( n )

  for i in range ( 0, n ):
    x[i] = b[i]
#
#  Solve R' * y = b.
#
  for k in range ( 0, n ):
    t = 0.0
    for i in range ( 0, k ):
      t = t + x[i] * r[i,k]
    x[k] = ( x[k] - t ) / r[k,k]
#
#  Solve R * x = y.
#
  for k in range ( n - 1, -1, -1 ):
    x[k] = x[k] / r[k,k]
    for i in range ( 0, k ):
      x[i] = x[i] - r[i,k] * x[k]

  return x

def r8po_sl_test ( ):

#*****************************************************************************80
#
## r8po_sl_test() tests r8po_sl().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  n = 5

  print ( '' )
  print ( 'r8po_sl_test():' )
  print ( '  r8po_sl() solves a linear system with an R8PO matrix' )
  print ( '  after it has been factored by r8po_fa.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )
#
#  Set (the upper half of) matrix A.
#
  a = r8po_zeros ( n )

  for i in range ( 0, n ):
    a[i,i] = 2.0
  for i in range ( 0, n - 1 ):
    a[i,i+1] = -1.0

  r8po_print ( n, a, '  Matrix A:' )
#
#  Factor the matrix.
#
  r = r8po_fa ( n, a )
#
#  Set the right hand side.
#
  b = np.zeros ( n )
  b[n-1] = float ( n + 1 )
  print ( '' )
  print ( '  right hand side b:' )
  print ( b )
#
#  Solve the linear system.
#
  x = r8po_sl ( n, r, b )
  print ( '' )
  print ( '  solution x:' )
  print ( x )

  return

def r8po_to_r8ge ( n, a ):

#*****************************************************************************80
#
## r8po_to_r8ge() copies a R8PO matrix to a R8GE matrix.
#
#  Discussion:
#
#    The R8PO storage format is appropriate for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#
#    real A(N,N), the R8PO matrix.
#
#  Output:
#
#    real B(N,N), the R8GE matrix.
#
  import numpy as np

  b = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      if ( i <= j ):
        b[i,j] = a[i,j]
      else:
        b[i,j] = a[j,i]

  return b

def r8po_to_r8ge_test ( ):

#*****************************************************************************80
#
## r8po_to_r8ge_test() tests r8po_to_r8ge().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
  n = 5

  print ( '' )
  print ( 'r8po_to_r8ge_test():' )
  print ( '  r8po_to_r8ge() converts a R8PO matrix to R8GE format.' )
  print ( '' )
  print ( '  Matrix order N =', n )

  a = r8po_random ( n )

  r8po_print ( n, a, '  The random R8PO matrix:' )
 
  r8mat_print ( n, n, a, '  The random R8PO matrix (printed by R8GE_print):' )

  b = r8po_to_r8ge ( n, a )

  r8mat_print ( n, n, b, '  The random R8GE matrix (printed by R8GE_print):' )

  return

def r8po_zeros ( n ):

#*****************************************************************************80
#
## r8po_zeros() zeroes an R8PO matrix.
#
#  Discussion:
#
#    The R8PO storage format is used for a symmetric positive definite 
#    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
#    upper triangular matrix, so it will be in R8GE storage format.)
#
#    Only the diagonal and upper triangle of the square array are used.
#    This same storage scheme is used when the matrix is factored by
#    r8po_fa, or inverted by r8po_inverse.  For clarity, the lower triangle
#    is set to zero.
#
#    R8PO storage is used by LINPACK and LAPACK.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the order of the matrix.
#    N must be positive.
#
#    real A(N,N), the R8PO matrix.
#
  import numpy as np

  a = np.zeros ( [ n, n ] )

  return a

def r8po_zeros_test ( ):

#*****************************************************************************80
#
## r8po_zeros_test() tests r8po_zeros().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2015
#
#  Author:
#
#    John Burkardt
#
  n = 5

  print ( '' )
  print ( 'r8po_zeros_test():' )
  print ( '  r8po_zeros() zeros out space for a' )
  print ( '  symmetric positive definite matrix.' )
  print ( '' )
  print ( '  Matrix order N = %d' % ( n ) )

  a = r8po_zeros ( n )

  r8po_print ( n, a, '  Matrix A:' )

  return

def r8_sign ( x ):

#*****************************************************************************80
#
## r8_sign() returns the sign of an R8.
#
#  Discussion:
#
#    The value is +1 if the number is positive or zero, and it is -1 otherwise.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 June 2013
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real X, the number whose sign is desired.
#
#  Output:
#
#    real VALUE, the sign of X.
#
  if ( x < 0.0 ):
    value = -1.0
  else:
    value = +1.0
 
  return value

def r8_sign_test ( ):

#*****************************************************************************80
#
## r8_sign_test() tests r8_sign().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 September 2014
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  test_num = 5

  r8_test = np.array ( [ -1.25, -0.25, 0.0, +0.5, +9.0 ] )

  print ( '' )
  print ( 'r8_sign_test():' )
  print ( '  r8_sign() returns the sign of an R8.' )
  print ( '' )
  print ( '     R8     r8_sign(R8)' )
  print ( '' )

  for test in range ( 0, test_num ):
    r8 = r8_test[test]
    s = r8_sign ( r8 )
    print ( '  %8.4f  %8.0f' % ( r8, s ) )

  return

def r8vec_angle_3d ( u, v ):

#*****************************************************************************80
#
## r8vec_angle_3d() computes the angle between two vectors in 3D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    07 May 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real U(1:3), V(1:3), the vectors.
#
#  Output:
#
#    real ANGLE, the angle between the two vectors.
#
  import numpy as np

  uv_dot = np.dot ( u, v )
  u_norm = np.linalg.norm ( u )
  v_norm = np.linalg.norm ( v )

  angle_cos = uv_dot / u_norm / v_norm

  angle = np.arccos ( angle_cos )

  return angle

def r8vec_cross_product_3d ( v1, v2 ):

#*****************************************************************************80
#
## r8vec_cross_product_3d() computes the cross product of two R8VEC's in 3D.
#
#  Discussion:
#
#    The cross product in 3D can be regarded as the determinant of the
#    symbolic matrix:
#
#          [  i  j  k ]
#      det [ x1 y1 z1 ]
#          [ x2 y2 z2 ]
#
#      = ( y1 * z2 - z1 * y2 ) * i
#      + ( z1 * x2 - x1 * z2 ) * j
#      + ( x1 * y2 - y1 * x2 ) * k
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 May 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real V1[3], V2[3], the two vectors.
#
#  Output:
#
#    real V3[3], the cross product vector.
#
  import numpy as np

  v3 = np.zeros ( 3 )

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1]
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2]
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0]

  return v3

def r8vec2_print ( a1, a2, title ):

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
#    20 June 2020
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
  n = len ( a1 )

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
#    20 June 2020
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'r8vec2_print_test():' )
  print ( '  r8vec2_print() prints a pair of R8VEC\'s.' )

  n = 6
  v = np.array ( [ 0.0, 0.20, 0.40, 0.60, 0.80, 1.0 ], dtype = np.float64 )
  w = np.array ( [ 0.0, 0.04, 0.16, 0.36, 0.64, 1.0 ], dtype = np.float64 )
  r8vec2_print ( v, w, '  Print a pair of R8VEC\'s:' )

  return

def r8vec3_print ( n, a1, a2, a3, title ):

#*****************************************************************************80
#
## r8vec3_print() prints an R8VEC3.
#
#  Discussion:
#
#    An R8VEC3 is a dataset consisting of 3 vectors of N real values.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 September 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of components of the vector.
#
#    real A1(N), A2(N), A3(N), the vectors to be printed.
#
#    string TITLE, a title.
#
  print ( '' )
  print ( title )
  print ( '' )
  for i in range ( 0, n ):
    print ( '  %6d:   %12g  %12g  %12g' % ( i, a1[i], a2[i], a3[i] ) )

  return

def r8vec3_print_test ( ):

#*****************************************************************************80
#
## r8vec3_print_test() tests r8vec3_print().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 September 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'r8vec3_print_test():' )
  print ( '  r8vec3_print() prints three R8VECs.' )

  n = 6

  t = np.array ( [ 0.0, 0.20, 0.40, 0.60, 0.80, 1.0 ], dtype = np.float64 )
  u = np.array ( [ 0.0, 0.04, 0.16, 0.36, 0.64, 1.0 ], dtype = np.float64 )
  v = np.array ( [ 0.0, 0.24, 0.56, 0.96, 1.44, 2.0 ], dtype = np.float64 )

  r8vec3_print ( n, t, u, v, '  X, X^2, X+X^2\'s:' )

  return

def r8vec_indicator1 ( n ):

#*****************************************************************************80
#
## r8vec_indicator1() sets an R8VEC to the indicator vector (1,2,3,...).
#
#  Discussion:
#
#    An R8VEC is a vector of R8's.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 September 2014
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of elements of the vector.
#
#  Output:
#
#    real A(N), the indicator array.
#
  import numpy as np

  a = np.zeros ( n )

  for i in range ( 0, n ):
    a[i] = i + 1

  return a

def r8vec_indicator1_test ( ):

#*****************************************************************************80
#
## r8vec_indicator1_test() tests r8vec_indicator1().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 February 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'r8vec_indicator1_test():' )
  print ( '  r8vec_indicator1() returns the 1-based indicator matrix.' )

  n = 10
  a = r8vec_indicator1 ( n )
  print ( '' )
  print ( '  the 1-based indicator vector:' )
  print ( a )

  return

def radians_to_degrees ( radians ):

#*****************************************************************************80
#
## radians_to_degrees() converts an angle from radians to degrees.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real RADIANS, an angle in radians.
#
#  Output:
#
#    real VALUE, the equivalent angle in degrees.
#
  import numpy as np

  value = ( radians / np.pi ) * 180.0

  return value

def radians_to_degrees_test ( ):

#*****************************************************************************80
#
## radians_to_degrees_test() tests radians_to_degrees().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'radians_to_degrees_test():' )
  print ( '  radians_to_degrees() converts an angle from radians' )
  print ( '  to degrees.' )
  print ( '' )
  print ( '  Degrees     Radians     Degrees' )
  print ( '' )

  for i in range ( -2, 15 ):
    angle_deg = float ( 30 * i )
    angle_rad = degrees_to_radians ( angle_deg )
    angle_deg2 = radians_to_degrees ( angle_rad )
    print ( '  %10f  %10f  %10f' % ( angle_deg, angle_rad, angle_deg2 ) )

  return

def segment_point_dist ( p1, p2, p ):

#*****************************************************************************80
#
## segment_point_dist(): distance ( line segment, point ) in 2D.
#
#  Discussion:
#
#    A line segment is the finite portion of a line that lies between
#    two points.
#
#    The nearest point will satisfy the condition
#
#      PN = (1-T) * P1 + T * P2.
#
#    T will always be between 0 and 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(2), P2(2), the endpoints of the line segment.
#
#    real P(2), the point whose nearest neighbor on the line
#    segment is to be determined.
#
#  Output:
#
#    real DIST, the distance from the point to the line segment.
#
  import numpy as np
#
#  If the line segment is actually a point, then the answer is easy.
#
  if ( p1[0] == p2[0] and p1[1] == p2[1] ):

    t = 0.0

  else:

    bot = ( p2[0] - p1[0] ) ** 2 + ( p2[1] - p1[1] ) ** 2

    t = ( ( p[0] - p1[0] ) * ( p2[0] - p1[0] ) \
        + ( p[1] - p1[1] ) * ( p2[1] - p1[1] ) ) / bot

    t = max ( t, 0.0 )
    t = min ( t, 1.0 )

  pn = np.zeros ( 2 )

  pn[0] = p1[0] + t * ( p2[0] - p1[0] )
  pn[1] = p1[1] + t * ( p2[1] - p1[1] )

  dist = np.sqrt ( ( pn[0] - p[0] ) ** 2 + ( pn[1] - p[1] ) ** 2 )

  return dist

def segment_point_dist_test ( ):

#*****************************************************************************80
#
## segment_point_dist_test() tests segment_point_dist().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'segment_point_dist_test():' )
  print ( '  segment_point_dist() computes the distance from a point to a line segment.' )

  p1 = np.array ( [ 1.0, 2.0 ] )
  p2 = np.array ( [ 3.0, 4.0 ] )

  print ( '' )
  print ( '  segment endpoint p1:' )
  print ( p1 )
  print ( '' )
  print ( '  segment endpoint p2:' )
  print ( p2 )

  p = np.array ( [ 2.0, 3.0 ] )
  dist = segment_point_dist ( p1, p2, p )
  print ( '' )
  print ( '  test point p:' )
  print ( p )

  print ( '' )
  print ( '  Distance to segment = %g' % ( dist ) )

  p = np.array ( [ 4.0, 5.0 ] )
  dist = segment_point_dist ( p1, p2, p )
  print ( '' )
  print ( '  test point p:' )
  print ( p )

  print ( '' )
  print ( '  Distance to segment = %g' % ( dist ) )

  p = np.array ( [ 1.0, 4.0 ] )
  dist = segment_point_dist ( p1, p2, p )
  print ( '' )
  print ( '  test point p:' )
  print ( p )

  print ( '' )
  print ( '  Distance to segment = %g' % ( dist ) )

  p = np.array ( [ 0.0, 0.0 ] )
  dist = segment_point_dist ( p1, p2, p )
  print ( '' )
  print ( '  test point p:' )
  print ( p )

  print ( '' )
  print ( '  Distance to segment = %g' % ( dist ) )

  return

def segment_point_near ( p1, p2, p ):

#*****************************************************************************80
#
## segment_point_near() finds the line segment point nearest a point in 2D.
#
#  Discussion:
#
#    A line segment is the finite portion of a line that lies between
#    two points.
#
#    The nearest point will satisfy the condition
#
#      PN = (1-T) * P1 + T * P2.
#
#    T will always be between 0 and 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real P1(2,1), P2(2,1), the endpoints of the line segment.
#
#    real P(2,1), the point whose nearest neighbor
#    on the line segment is to be determined.
#
#  Output:
#
#    real PN(2,1), the point on the line segment which is
#    nearest the point (X,Y).
#
#    real DIST, the distance from the point to the
#    nearest point on the line segment.
#
#    real T, the relative position of the point (XN,YN)
#    to the points (X1,Y1) and (X2,Y2).
#
  import numpy as np
#
#  If the line segment is actually a point, then the answer is easy.
#
  if ( p1[0] == p2[0] and p1[1] == p2[1] ):

    t = 0.0

  else:

    bot = ( p2[0] - p1[0] ) ** 2 + ( p2[1] - p1[1] ) ** 2

    t = ( ( p[0] - p1[0] ) * ( p2[0] - p1[0] ) \
        + ( p[1] - p1[1] ) * ( p2[1] - p1[1] ) ) / bot

    t = max ( t, 0.0 )
    t = min ( t, 1.0 )

  pn = np.zeros ( 2 )

  pn[0] = p1[0] + t * ( p2[0] - p1[0] )
  pn[1] = p1[1] + t * ( p2[1] - p1[1] )

  dist = np.sqrt ( ( pn[0] - p[0] ) ** 2 + ( pn[1] - p[1] ) ** 2 )

  return pn, dist, t

def segment_point_near_test ( ):

#*****************************************************************************80
#
## segment_point_near_test() tests segment_point_near().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  test_num = 3

  print ( '' )
  print ( 'segment_point_near_test():' )
  print ( '  segment_point_near() computes the nearest point' )
  print ( '  from a line segment to a point in 2D.' )

  for test in range ( 0, test_num ):

    p1 = np.random.rand ( 2 )
    p2 = np.random.rand ( 2 )
    p = np.random.rand ( 2 )

    pn, dist, t = segment_point_near ( p1, p2, p )

    print ( '' )
    print ( '  TEST = %d' % ( test ) )
    print ( '  P1 =   %12f  %12f' % ( p1[0], p1[1] ) )
    print ( '  P2 =   %12f  %12f' % ( p2[0], p2[1] ) )
    print ( '  P =    %12f  %12f' % ( p[0], p[1] ) )
    print ( '  PN =   %12f  %12f' % ( pn[0], pn[1] ) )
    print ( '  DIST = %12f' % ( dist ) )
    print ( '  T =    %12f' % ( t ) )

  return

def simplex01_volume ( m ):

#*****************************************************************************80
#
## simplex01_volume() returns the volume of the unit simplex in M dimensions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 January 2014
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer M, the spatial dimension.
#
#  Output:
#
#    real VALUE, the volume.
#
  value = 1.0
  for i in range ( 1, m + 1 ):
    value = value / float ( i )

  return value

def simplex01_volume_test ( ) :

#*****************************************************************************80
#
## simplex01_volume_test() tests simplex01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'simplex01_volume_test():' )
  print ( '  simplex01_volume() returns the volume of the unit simplex' )
  print ( '  in M dimensions.' )
  print ( '' )
  print ( '   M   Volume' )
  print ( '' )

  for m in range ( 1, 10 ):
    value = simplex01_volume ( m )
    print ( '  %2d  %g' % ( m, value ) )

  return

def sphere01_area ( ):

#*****************************************************************************80
#
## sphere01_area() returns the area of the surface of the unit sphere in 3D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    24 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VALUE, the area.
#
  import numpy as np

  r = 1.0
  value = 4.0 * np.pi * r * r

  return value

def sphere01_area_test ( ) :

#*****************************************************************************80
#
## sphere01_area_test() tests sphere01_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'sphere01_area_test():' )
  print ( '  sphere01_area() returns the area of the unit sphere.' )
  print ( '' )

  value = sphere01_area ( )

  print ( '  sphere01_area() =  %g' % ( value ) )

  return

def sphere01_area_values ( n_data ):

#*****************************************************************************80
#
## sphere01_area_values() returns some areas of the unit sphere in ND.
#
#  Discussion:
#
#    The formula for the surface area of the unit sphere in N dimensions is:
#
#      sphere_Unit_area ( N ) = 2 * PI^(N/2) / Gamma ( N / 2 )
#
#    Some values of the function include:
#
#       N   Area
#
#       2    2        * PI
#       3  ( 4 /    ) * PI
#       4  ( 2 /   1) * PI^2
#       5  ( 8 /   3) * PI^2
#       6  ( 1 /   1) * PI^3
#       7  (16 /  15) * PI^3
#       8  ( 1 /   3) * PI^4
#       9  (32 / 105) * PI^4
#      10  ( 1 /  12) * PI^5
#
#    For the unit sphere, Area(N) = N * Volume(N)
#
#    In Mathematica, the function can be evaluated by:
#
#      2 * Pi^(n/2) / Gamma[n/2]
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 February 2015
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
#    integer N_data.  The user sets N_data to 0 before the first call.
#
#  Output:
#
#    integer N_data. On each subsequent
#    call, N_data is incremented and that test data is returned.  When
#    there is no more test data, N_data is set to 0.
#
#    integer N, the spatial dimension.
#
#    real AREA, the area of the unit sphere 
#    in that dimension.
#
  import numpy as np

  n_max = 20

  area_vec = np.array ( ( \
     0.2000000000000000E+01, \
     0.6283185307179586E+01, \
     0.1256637061435917E+02, \
     0.1973920880217872E+02, \
     0.2631894506957162E+02, \
     0.3100627668029982E+02, \
     0.3307336179231981E+02, \
     0.3246969701133415E+02, \
     0.2968658012464836E+02, \
     0.2550164039877345E+02, \
     0.2072514267328890E+02, \
     0.1602315322625507E+02, \
     0.1183817381218268E+02, \
     0.8389703410491089E+01, \
     0.5721649212349567E+01, \
     0.3765290085742291E+01, \
     0.2396678817591364E+01, \
     0.1478625959000308E+01, \
     0.8858104195716824E+00, \
     0.5161378278002812E+00 ))

  n_vec = np.array ( ( \
     1, \
     2, \
     3, \
     4, \
     5, \
     6, \
     7, \
     8, \
     9, \
    10, \
    11, \
    12, \
    13, \
    14, \
    15, \
    16, \
    17, \
    18, \
    19, \
    20 ))

  if ( n_data < 0 ):
    n_data = 0

  if ( n_max <= n_data ):
    n_data = 0
    n = 0.0
    area = 0.0
  else:
    n = n_vec[n_data]
    area = area_vec[n_data]
    n_data = n_data + 1

  return n_data, n, area

def sphere01_area_values_test ( ):

#*****************************************************************************80
#
## sphere01_area_values_test() tests sphere01_area_values().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 February 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'sphere01_area_values_test():' )
  print ( '  sphere01_area_values() stores areas of the unit sphere in N dimensions.' )
  print ( '' )
  print ( '      N         sphere01_area(N)' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, n, area = sphere01_area_values ( n_data )

    if ( n_data == 0 ):
      break

    print ( '  %12d  %24.16f' % ( n, area ) )

  return

def sphere01_volume_values ( n_data ):

#*****************************************************************************80
#
## sphere01_volume_values() returns some volumes of the unit sphere in ND.
#
#  Discussion:
#
#    The formula for the volume of the unit sphere in N dimensions is
#
#      Volume(N) = 2 * PI^(N/2) / ( N * Gamma ( N / 2 ) )
#
#    This function satisfies the relationships:
#
#      Volume(N) = 2 * PI * Volume(N-2) / N
#      Volume(N) = Area(N) / N
#
#    Some values of the function include:
#
#       N  Volume
#
#       1    1
#       2    1        * PI
#       3  ( 4 /   3) * PI
#       4  ( 1 /   2) * PI^2
#       5  ( 8 /  15) * PI^2
#       6  ( 1 /   6) * PI^3
#       7  (16 / 105) * PI^3
#       8  ( 1 /  24) * PI^4
#       9  (32 / 945) * PI^4
#      10  ( 1 / 120) * PI^5
#
#    In Mathematica, the function can be evaluated by:
#
#      2 * Pi^(n/2) / ( n * Gamma[n/2] )
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 February 2015
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
#    integer N_data.  The user sets N_data to 0 before the first call.
#
#  Output:
#
#    integer N_data.  On each subsequent
#    call, N_data is incremented and that test data is returned.  When
#    there is no more test data, N_data is set to 0.
#
#    integer N, the spatial dimension.
#
#    real VOLUME, the volume of the unit 
#    sphere in that dimension.
#
  import numpy as np

  n_max = 20

  n_vec = np.array ( ( \
     1, \
     2, \
     3, \
     4, \
     5, \
     6, \
     7, \
     8, \
     9, \
    10, \
    11, \
    12, \
    13, \
    14, \
    15, \
    16, \
    17, \
    18, \
    19, \
    20 ))

  volume_vec = np.array ( ( \
     0.2000000000000000E+01, \
     0.3141592653589793E+01, \
     0.4188790204786391E+01, \
     0.4934802200544679E+01, \
     0.5263789013914325E+01, \
     0.5167712780049970E+01, \
     0.4724765970331401E+01, \
     0.4058712126416768E+01, \
     0.3298508902738707E+01, \
     0.2550164039877345E+01, \
     0.1884103879389900E+01, \
     0.1335262768854589E+01, \
     0.9106287547832831E+00, \
     0.5992645293207921E+00, \
     0.3814432808233045E+00, \
     0.2353306303588932E+00, \
     0.1409811069171390E+00, \
     0.8214588661112823E-01, \
     0.4662160103008855E-01, \
     0.2580689139001406E-01  ))

  n_vec = np.array ( ( \
     1, \
     2, \
     3, \
     4, \
     5, \
     6, \
     7, \
     8, \
     9, \
    10, \
    11, \
    12, \
    13, \
    14, \
    15, \
    16, \
    17, \
    18, \
    19, \
    20 ))

  if ( n_data < 0 ):
    n_data = 0

  if ( n_max <= n_data ):
    n_data = 0
    n = 0.0
    volume = 0.0
  else:
    n = n_vec[n_data]
    volume = volume_vec[n_data]
    n_data = n_data + 1

  return n_data, n, volume

def sphere01_volume_values_test ( ):

#*****************************************************************************80
#
## sphere01_volume_values_test() tests sphere01_volume_values().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 February 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'sphere01_volume_values_test():' )
  print ( '  sphere01_volume_values() stores values of the sphere01_volume function.' )
  print ( '' )
  print ( '      N         sphere01_volume(X)' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, n, volume = sphere01_volume_values ( n_data )

    if ( n_data == 0 ):
      break

    print ( '  %12d  %24.16f' % ( n, volume ) )

  return

def sphere_triangle_angles_to_area ( r, a, b, c ):

#*****************************************************************************80
#
## sphere_triangle_angles_to_area() computes the area of a spherical triangle.
#
#  Discussion:
#
#    A sphere centered at 0 in 3D satisfies the equation:
#
#      X*X + Y*Y + Z*Z = R*R
#
#    A spherical triangle is specified by three points on the surface
#    of the sphere.
#
#    The area formula is known as Girard's formula.
#
#    The area of a spherical triangle is:
#
#      AREA = ( A + B + C - PI ) * R*R
#
#    where A, B and C are the (surface) angles of the triangle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 February 2005
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the sphere.
#
#    real A, B, C, the angles of the triangle.
#
#  Output:
#
#    real AREA, the area of the spherical triangle.
#
  import numpy as np
#
#  Apply Girard's formula.
#
  area = r * r * ( a + b + c - np.pi )

  return area

def sphere_triangle_sides_to_angles ( r, aside, bside, cside ):

#*****************************************************************************80
#
## sphere_triangle_sides_to_angles() computes spherical triangle angles in 3D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the sphere.
#
#    real ASIDE, BSIDE, CSIDE, the (geodesic) length of the 
#    sides of the triangle.
#
#  Output:
#
#    real A, B, C, the spherical angles of the triangle.
#    Angle A is opposite the side of length ASIDE, and so on.
#
  import numpy as np

  asu = aside / r
  bsu = bside / r
  csu = cside / r
  ssu = ( asu + bsu + csu ) / 2.0

  tan_a2 = np.sqrt ( ( np.sin ( ssu - bsu ) * np.sin ( ssu - csu ) ) / \
                     ( np.sin ( ssu )       * np.sin ( ssu - asu ) ) )

  a = 2.0 * np.arctan ( tan_a2 )

  tan_b2 = np.sqrt ( ( np.sin ( ssu - asu ) * np.sin ( ssu - csu ) ) / \
                     ( np.sin ( ssu )       * np.sin ( ssu - bsu ) ) )

  b = 2.0 * np.arctan ( tan_b2 )

  tan_c2 = np.sqrt ( ( np.sin ( ssu - asu ) * np.sin ( ssu - bsu ) ) / \
                     ( np.sin ( ssu )       * np.sin ( ssu - csu ) ) )

  c = 2.0 * np.arctan ( tan_c2 )

  return a, b, c

def sphere_triangle_sides_to_angles_test ( ):

#*****************************************************************************80
#
## sphere_triangle_sides_to_angles_test() tests sphere_triangle_sides_to_angles().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    19 July 2018
#
#  Author:
#
#    John Burkardt
#
  r = 10.0

  print ( '' )
  print ( 'sphere_triangle_sides_to_angles_test():' )
  print ( '  sphere_triangle_sides_to_angles() takes the sides of a' )
  print ( '  spherical triangle and determines the angles.' )

  aside = 121.0 + ( 15.4 / 60.0 )
  bside = 104.0 + ( 54.7 / 60.0 )
  cside =  65.0 + ( 42.5 / 60.0 )

  aside = degrees_to_radians ( aside )
  bside = degrees_to_radians ( bside )
  cside = degrees_to_radians ( cside )

  aside = r * aside
  bside = r * bside
  cside = r * cside
#
#  Get the spherical angles.
#
  a, b, c = sphere_triangle_sides_to_angles ( r, aside, bside, cside )

  print ( '' )
  print ( '  A       = %8f (radians)' % ( a ) )
  a = radians_to_degrees ( a )
  print ( '  A       = %8f (degrees)' % ( a ) )
  a = 117.0 + ( 58.0 / 60.0 )
  print ( '  Correct = %8f (radians)' % ( a ) )

  print ( '' )
  print ( '  B       = %8f (radians)' % ( b ) )
  b = radians_to_degrees ( b )
  print ( '  B       = %8f (degrees)' % ( b ) )
  b = 93.0 + ( 13.8 / 60.0 )
  print ( '  Correct = %8f (radians)' % ( b ) )

  print ( '' )
  print ( '  C       = %8f (radians)' % ( c ) )
  c = radians_to_degrees ( c )
  print ( '  C       = %8f (degrees)' % ( c ) )
  c = 70.0 + ( 20.6 / 60.0 )
  print ( '  Correct = %8f (radians)' % ( c ) )

  return

def sphere_triangle_vertices_to_sides ( r, v1, v2, v3 ):

#*****************************************************************************80
#
## sphere_triangle_vertices_to_sides() computes spherical triangle sides.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 April 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real R, the radius of the sphere.
#
#    real V1(3), V2(3), V3(3), the vertices of the spherical triangle.
#
#  Output:
#
#    real ASIDE, BSIDE, CSIDE, the (geodesic) length of the sides
#    of the triangle.
#
  import numpy as np

  aside = r * np.arccos ( np.dot ( v2, v3 ) / r ** 2 )
  bside = r * np.arccos ( np.dot ( v3, v1 ) / r ** 2 )
  cside = r * np.arccos ( np.dot ( v1, v2 ) / r ** 2 )

  return aside, bside, cside

def tetrahedron_barycentric ( tetra, p ):

#*****************************************************************************80
#
## tetrahedron_barycentric() returns barycentric coordinates of a point in 3D.
#
#  Discussion:
#
#    The barycentric coordinates of a point (X,Y,Z) with respect to
#    a tetrahedron are a set of four values C(1:4), each associated
#    with a vertex of the tetrahedron.  The values must sum to 1.
#    If all the values are between 0 and 1, the point is contained
#    within the tetrahedron.
#
#    The barycentric coordinate of point X related to vertex A can be
#    interpreted as the ratio of the volume of the tetrahedron with 
#    vertex A replaced by vertex X to the volume of the original 
#    tetrahedron.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real TETRA(3,4) the tetrahedron vertices.
#
#    real P(3), the point to be checked.
#
#  Output:
#
#    real C(4), the barycentric coordinates of (X,Y,Z) with
#    respect to the tetrahedron.
#
  import numpy as np
#
#  Set up the linear system
#
#    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
#    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
#    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
#
#  which is satisfied by the barycentric coordinates of (X,Y,Z).
#
  a = np.zeros ( [ 3, 4 ] )

  a[0:3,0:3] = tetra[0:3,1:4]
  for i in range ( 0, 3 ):
    a[i,3] = p[i]

  for i in range ( 0, 3 ):
    a[i,0:4] = a[i,0:4] - tetra[i,0]
#
#  Solve the linear system.
#
  nrhs = 1
  a, info = r8mat_solve ( 3, nrhs, a )

  if ( info != 0 ):
    print ( '' )
    print ( 'tetrahedron_barycentric - Fatal error!' )
    print ( '  The linear system is singular.' )
    print ( '  The input data does not form a proper tetrahedron.' )
    raise Exception ( 'tetrahedron_barycentric - Fatal error!' )

  c = np.zeros ( 4 )

  c[1:4] = a[0:3,3]

  c[0] = 1.0 - np.sum ( c[1:4] )

  return c

def tetrahedron_barycentric_test ( ):

#*****************************************************************************80
#
## tetrahedron_barycentric_test() tests tetrahedron_barycentric().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 1.0, 4.0, 3.0 ], \
    [ 2.0, 4.0, 3.0 ], \
    [ 1.0, 6.0, 3.0 ], \
    [ 1.0, 4.0, 4.0 ] ] )

  t = np.transpose ( t )

  print ( '' )
  print ( 'tetrahedron_barycentric_test():' )
  print ( '  tetrahedron_barycentric() converts XYZ to XSI.' )
  print ( '  We are computing the XSI coordinates just to verify' )
  print ( '  that the points are inside the tetrahedron.' )

  r8mat_transpose_print ( 3, 4, t, '  Tetrahedron vertices' )

  print ( '' )
  print ( '  (X,Y,Z)   (XSI1,XSI2,XSI3,XSI4):' )
  print ( '' )

  for i in range ( 0, 10 ):
    p = tetrahedron_sample ( t, 1 )
    xsi = tetrahedron_barycentric ( t, p )
    print ( '  %8f  %8f  %8f    %8f  %8f  %8f  %8f' \
      % ( p[0], p[1], p[2], xsi[0], xsi[1], xsi[2], xsi[3] ) )

  return

def tetrahedron_centroid ( tetra ):

#*****************************************************************************80
#
## tetrahedron_centroid() computes the centroid of a tetrahedron.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real TETRA(3,4) the tetrahedron vertices.
#
#  Output:
#
#    real CENTROID(3), the coordinates of the centroid.
#
  import numpy as np

  centroid = np.sum ( tetra, 1 ) / 4.0

  return centroid

def tetrahedron_centroid_test ( ):

#*****************************************************************************80
#
## tetrahedron_centroid_test() tests tetrahedron_centroid().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  tetra = np.array ( [ \
    [  0.000000,  0.942809, -0.333333 ], \
    [ -0.816496, -0.816496, -0.333333 ], \
    [  0.816496, -0.816496, -0.333333 ], \
    [  0.000000,  0.000000,  1.000000 ] ] )

  tetra = np.transpose ( tetra )

  print ( '' )
  print ( 'tetrahedron_centroid_test():' )
  print ( '  tetrahedron_centroid() computes the centroid of a tetrahedron' )

  r8mat_transpose_print ( 3, 4, tetra, '  Tetrahedron vertices:' )

  centroid = tetrahedron_centroid ( tetra )

  print ( '' )
  print ( '  tetrahedron centroid:' )
  print ( centroid )

  return

def tetrahedron_sample ( t, n ):

#*****************************************************************************80
#
## tetrahedron_sample() returns random points in a tetrahedron.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(3,4), the tetrahedron vertices.
#
#    integer N, the number of points to generate.
#
#  Output:
#
#    real P(3,N), random points in the tetrahedron.
#
  import numpy as np

  p = np.zeros ( [ 3, n ] )

  for j in range ( 0, n ):

    r = np.random.rand ( )
#
#  Interpret R as a percentage of the tetrahedron's volume.
#
#  Imagine a plane, parallel to face 1, so that the volume between
#  vertex 1 and the plane is R percent of the full tetrahedron volume.
#
#  The plane will intersect sides 12, 13, and 14 at a fraction
#  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
#
    alpha = r ** ( 1.0 / 3.0 )
#
#  Determine the coordinates of the points on sides 12, 13 and 14 intersected
#  by the plane, which form a triangle TR.
#
    tr = np.zeros ( [ 3, 3 ] )

    tr[:,0] = alpha * t[:,0] + ( 1.0 - alpha ) * t[:,1]
    tr[:,1] = alpha * t[:,0] + ( 1.0 - alpha ) * t[:,2]
    tr[:,2] = alpha * t[:,0] + ( 1.0 - alpha ) * t[:,3]
#
#  Now choose, uniformly at random, a point in this triangle.
#
    r = np.random.rand ( )
#
#  Interpret R as a percentage of the triangle's area.
#
#  Imagine a line L, parallel to side 1, so that the area between
#  vertex 1 and line L is R percent of the full triangle's area.
#
#  The line L will intersect sides 2 and 3 at a fraction
#  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
#
    alpha = np.sqrt ( r )
#
#  Determine the coordinates of the points on sides 2 and 3 intersected
#  by line L.
#
    p12 = alpha * tr[:,0] + ( 1.0 - alpha ) * tr[:,1]
    p13 = alpha * tr[:,0] + ( 1.0 - alpha ) * tr[:,2]
#
#  Now choose, uniformly at random, a point on the line L.
#
    beta = np.random.rand ( )

    p[:,j] = beta * p12[:] + ( 1.0 - beta ) * p13[:]

  return p

def tetrahedron_sample_test ( ):

#*****************************************************************************80
#
## tetrahedron_sample_test() tests tetrahedron_sample().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
     [ 1.0, 4.0, 3.0 ], \
     [ 2.0, 4.0, 3.0 ], \
     [ 1.0, 6.0, 3.0 ], \
     [ 1.0, 4.0, 4.0 ] ] )

  t = np.transpose ( t )

  print ( '' )
  print ( 'tetrahedron_sample_test():' )
  print ( '  tetrahedron_sample() samples a tetrahedron.' )
  print ( '  We are computing the XSI coordinates just to verify' )
  print ( '  that the points are inside the tetrahedron.' )

  r8mat_transpose_print ( 3, 4, t, '  Tetrahedron vertices' )

  print ( '' )
  print ( '  (X,Y,Z)   (XSI1,XSI2,XSI3,XSI4):' )
  print ( '' )

  for i in range ( 0, 10 ):
    p = tetrahedron_sample ( t, 1 )
    xsi = tetrahedron_barycentric ( t, p )
    print ( '  %8f  %8f  %8f    %8f  %8f  %8f  %8f' \
      % ( p[0], p[1], p[2], xsi[0], xsi[1], xsi[2], xsi[3] ) )

  return

def tetrahedron_volume ( tetra ):

#*****************************************************************************80
#
## tetrahedron_volume() computes the volume of a tetrahedron.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real TETRA(3,4), the vertices of the tetrahedron.
#
#  Output:
#
#    real VOLUME, the volume of the tetrahedron.
#
  import numpy as np

  a = np.zeros ( [ 4, 4 ] )

  a[0:3,0:4] = tetra.copy ( )
  a[3,0:4] = 1.0

  value = abs ( r8mat_det_4d ( a ) ) / 6.0

  return value

def tetrahedron_volume_test ( ):

#*****************************************************************************80
#
## tetrahedron_volume_test() tests tetrahedron_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    20 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  tetra = np.array ( [ \
    [  0.000000,  0.942809, -0.333333 ], \
    [ -0.816496, -0.816496, -0.333333 ], \
    [  0.816496, -0.816496, -0.333333 ], \
    [  0.000000,  0.000000,  1.000000 ] ] )

  tetra = np.transpose ( tetra )

  print ( '' )
  print ( 'tetrahedron_volume_test():' )
  print ( '  tetrahedron_volume() computes the volume of a tetrahedron' )

  r8mat_transpose_print ( 3, 4, tetra, '  Tetrahedron vertices' )

  volume = tetrahedron_volume ( tetra )

  print ( '' )
  print ( '  Volume = %g' % ( volume ) )

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

def timestamp_test ( ):

#*****************************************************************************80
#
## timestamp_test() tests timestamp().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    03 December 2014
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'timestamp_test():' )
  print ( '  timestamp() prints a timestamp of the current date and time.' )
  print ( '' )

  timestamp ( )

  return

def triangle01_area ( ):

#*****************************************************************************80
#
## triangle01_area() computes the area of the unit triangle in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    13 January 2014
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real AREA, the area.
#
  area = 0.5

  return area

def triangle01_area_test ( ):

#*****************************************************************************80
#
## triangle01_area_test() tests triangle01_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 November 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
 
  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle01_area_test():' )
  print ( '  triangle01_area() computes the area of the unit triangle.' )

  r8mat_print ( 2, 3, t, '  Triangle vertices (columns)' )

  area = triangle01_area ( )

  print ( '' )
  print ( '  Triangle area is %g' % ( area ) )

  return

def triangle01_sample ( n ):

#*****************************************************************************80
#
## triangle01_sample() samples the interior of the unit triangle in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 November 2016
#
#  Author:
#
#    John Burkardt
#
#  Reference:
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
#    real XY(2,N), the points.
#
  import numpy as np

  m = 2

  xy = np.zeros ( [ m, n ] )

  for j in range ( 0, n ):

    e = np.random.rand ( m + 1 )

    e = - np.log ( e )

    d = np.sum ( e )

    xy[0:2,j] = e[0:2] / d

  return xy

def triangle01_sample_test ( ):

#*****************************************************************************80
#
## triangle01_sample_test() tests triangle01_sample().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
 
  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle01_sample_test():' )
  print ( '  triangle01_sample() samples the unit triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '  Sample points (X,Y):' )
  print ( '' )

  n = 10
  xy = triangle01_sample ( n )
  r8mat_transpose_print ( 2, n, xy, '  Sample points:' )

  return

def triangle_angles ( t ):

#*****************************************************************************80
#
## triangle_angles() computes the angles of a triangle in 2D.
#
#  Discussion:
#
#    The law of cosines is used:
#
#      C^2 = A^2 + B^2 - 2 * A * B * COS ( GAMMA )
#
#    where GAMMA is the angle opposite side C.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real ANGLE(3), the angles opposite
#    sides P1-P2, P2-P3 and P3-P1, in radians.
#
  import numpy as np
#
#  Compute the length of each side.
#
  a = np.sqrt ( ( t[0,0] - t[0,1] ) ** 2 + ( t[1,0] - t[1,1] ) ** 2 )
  b = np.sqrt ( ( t[0,1] - t[0,2] ) ** 2 + ( t[1,1] - t[1,2] ) ** 2 )
  c = np.sqrt ( ( t[0,2] - t[0,0] ) ** 2 + ( t[1,2] - t[1,0] ) ** 2 )

  angle = np.zeros ( 3 )
#
#  Take care of unlikely special cases.
#
  if ( a == 0.0 and b == 0.0 and c == 0.0 ):
    for i in range ( 0, 3 ):
      angle[i] = 2.0 * np.pi / 3.0
    return angle

  if ( c == 0.0 or a == 0.0 ):
    angle[0] = np.pi
  else:
    angle[0] = np.arccos ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) )

  if ( a == 0.0 or b == 0.0 ):
    angle[1] = np.pi
  else:
    angle[1] = np.arccos ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) )

  if ( b == 0.0 or c == 0.0 ):
    angle[2] = np.pi
  else:
    angle[2] = np.arccos ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) )

  return angle

def triangle_angles_test ( ):

#*****************************************************************************80
#
## triangle_angles_test() tests triangle_angles().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle_angles_test():' )
  print ( '  triangle_angles() computes the angles of a triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  angle = triangle_angles ( t )

  print ( '' )
  print ( '      Radians      Degrees' )
  print ( '' )
  for i in range ( 0, 3 ):
    print ( '  %12g  %12g' % ( angle[i], angle[i] * 180.0 / np.pi ) )

  return

def triangle_area ( t ):

#*****************************************************************************80
#
## triangle_area() computes the area of a triangle in 2D.
#
#  Discussion:
#
#    If the triangle's vertices are given in counterclockwise order,
#    the area will be positive.  If the triangle's vertices are given
#    in clockwise order, the area will be negative!
#
#    An earlier version of this routine always returned the absolute
#    value of the computed area.  I am convinced now that that is
#    a less useful result!  For instance, by returning the signed 
#    area of a triangle, it is possible to easily compute the area 
#    of a nonconvex polygon as the sum of the (possibly negative) 
#    areas of triangles formed by node 1 and successive pairs of vertices.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real AREA, the area of the triangle.
#
  area = 0.5 * ( \
      t[0,0] * ( t[1,1] - t[1,2] ) \
    + t[0,1] * ( t[1,2] - t[1,0] ) \
    + t[0,2] * ( t[1,0] - t[1,1] ) )

  return area

def triangle_area_test ( ):

#*****************************************************************************80
#
## triangle_area_test() tests triangle_area().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle_area_test():' )
  print ( '  triangle_area() computes the area of a triangle.' )

  r8mat_print ( 2, 3, t, '  Triangle vertices (columns)' )

  area = triangle_area ( t )

  print ( '' )
  print ( '  Triangle area is %g' % ( area ) )

  return

def triangle_barycentric ( t, p ):

#*****************************************************************************80
#
## triangle_barycentric() finds the barycentric coordinates of a point.
#
#  Discussion:
#
#    The barycentric coordinate of point X related to vertex A can be
#    interpreted as the ratio of the area of the triangle with
#    vertex A replaced by vertex X to the area of the original
#    triangle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#    The vertices should be given in counter clockwise order.
#
#    real P(2), the point to be checked.
#
#  Output:
#
#    real XSI(3), the barycentric coordinates of (X,Y)
#    with respect to the triangle.
#
  import numpy as np

  nrhs = 1
#
#  Set up the linear system
#
#    ( X2-X1  X3-X1 ) XSI(1)  = X-X1
#    ( Y2-Y1  Y3-Y1 ) XSI(2)    Y-Y1
#
#  which is satisfied by the barycentric coordinates of (X,Y).
#
  a = np.zeros ( [ 2, 3 ] )

  a[0,0] = t[0,1] - t[0,0]
  a[0,1] = t[0,2] - t[0,0]
  a[0,2] = p[0]   - t[0,0]

  a[1,0] = t[1,1] - t[1,0]
  a[1,1] = t[1,2] - t[1,0]
  a[1,2] = p[1]   - t[1,0]
#
#  Solve the linear system.
#
  a, info = r8mat_solve ( 2, nrhs, a )

  if ( info != 0 ):
    print ( '' )
    print ( 'triangle_barycentric - Fatal error!' )
    print ( '  The linear system is singular.' )
    print ( '  The input data does not form a proper triangle.' )
    raise Exception ( 'triangle_barycentric - Fatal error!' )

  xsi = np.zeros ( 3 )

  xsi[0] = a[0,2]
  xsi[1] = a[1,2]
  xsi[2] = 1.0 - xsi[0] - xsi[1]

  return xsi

def triangle_barycentric_test ( ):

#*****************************************************************************80
#
## triangle_barycentric_test() tests triangle_barycentric().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 July 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 7

  p_test = np.array ( [ \
    [ 0.25, 0.75, 1.00, 11.00, 0.00,   0.50, 0.60 ], \
    [ 0.25, 0.25, 1.00,  0.50, 1.00, -10.00, 0.60 ] ] )

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle_barycentric_test():' )
  print ( '  triangle_barycentric() converts XY coordinates' )
  print ( '  to barycentric XSI coordinates in a triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '   X       Y     XSI' )
  print ( '' )

  for j in range ( 0, ntest ):

    p = p_test[:,j]

    xsi = triangle_barycentric ( t, p )

    print ( '  %10f  %10f    %10f  %10f  %10f' \
      % ( p[0], p[1], xsi[0], xsi[1], xsi[2] ) )

  return

def triangle_centroid ( t ):

#*****************************************************************************80
#
## triangle_centroid() computes the centroid of a triangle in 2D.
#
#  Discussion:
#
#    The centroid of a triangle can also be considered the
#    center of gravity, or center of mass, assuming that the triangle
#    is made of a thin uniform sheet of massy material.
#
#    The centroid of a triangle is the intersection of the medians.
#
#    A median of a triangle is a line connecting a vertex to the
#    midpoint of the opposite side.
#
#    In barycentric coordinates, in which the vertices of the triangle
#    have the coordinates (1,0,0), (0,1,0) and (0,0,1), the centroid
#    has coordinates (1/3,1/3,1/3).
#
#    In geometry, the centroid of a triangle is often symbolized by "G".
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer, John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real CENTROID(2), the coordinates of the centroid.
#
  import numpy as np

  centroid = np.zeros ( 2 )

  for i in range ( 0, 2 ):
    for j in range ( 0, 3 ):
      centroid[i] = centroid[i] + t[i,j]
    centroid[i] = centroid[i] / 3.0

  return centroid

def triangle_centroid_test ( ):

#*****************************************************************************80
#
## triangle_centroid_test() tests triangle_centroid().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 4

  print ( '' )
  print ( 'triangle_centroid_test():' )
  print ( '  triangle_centroid() computes the centroid of a triangle' )

  for i in range ( 0, ntest ):

    if ( i == 0 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.0 ], \
        [ 0.0, 0.0, 1.0 ] ] )
    elif ( i == 1 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 0.86602539 ] ] )
    elif ( i == 2 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 10.0 ] ] )
    elif ( i == 3 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 10.0 ], \
        [ 0.0, 0.0, 2.0 ] ] )

    r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    centroid = triangle_centroid ( t )

    print ( '' )
    print ( '  triangle centroid:' )
    print ( centroid )

  return

def triangle_circumcircle ( t ):

#*****************************************************************************80
#
## triangle_circumcircle() computes the circumcircle of a triangle in 2D.
#
#  Discussion:
#
#    The circumcenter of a triangle is the center of the circumcircle, the
#    circle that passes through the three vertices of the triangle.
#
#    The circumcircle contains the triangle, but it is not necessarily the
#    smallest triangle to do so.
#
#    If all angles of the triangle are no greater than 90 degrees, then
#    the center of the circumscribed circle will lie inside the triangle.
#    Otherwise, the center will lie outside the triangle.
#
#    The circumcenter is the intersection of the perpendicular bisectors
#    of the sides of the triangle.
#
#    In geometry, the circumcenter of a triangle is often symbolized by "O".
#
#    Thanks to Chenguang Zhang for pointing out a mistake in the formula
#    for the center, 02 December 2016.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 December 2016
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real R, CENTER(2,1), the circumradius and circumcenter
#    of the triangle.
#
  import numpy as np

  center = np.zeros ( 2 )
#
#  Circumradius.
#
  a = np.sqrt ( ( t[0,0] - t[0,1] ) ** 2 + ( t[1,0] - t[1,1] ) ** 2 )
  b = np.sqrt ( ( t[0,1] - t[0,2] ) ** 2 + ( t[1,1] - t[1,2] ) ** 2 )
  c = np.sqrt ( ( t[0,2] - t[0,0] ) ** 2 + ( t[1,2] - t[1,0] ) ** 2 )

  bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c ) * (   a + b - c )

  if ( bot <= 0.0 ):
    r = - 1.0
    return r, center

  r = a * b * c / np.sqrt ( bot )
#
#  Circumcenter.
#
  f = np.zeros ( 2 )

  f[0] = ( t[0,1] - t[0,0] ) ** 2 + ( t[1,1] - t[1,0] ) ** 2
  f[1] = ( t[0,2] - t[0,0] ) ** 2 + ( t[1,2] - t[1,0] ) ** 2

  top = np.zeros ( 2 )

  top[0] =    ( t[1,2] - t[1,0] ) * f[0] - ( t[1,1] - t[1,0] ) * f[1]
  top[1] =  - ( t[0,2] - t[0,0] ) * f[0] + ( t[0,1] - t[0,0] ) * f[1]

  det  =    ( t[1,2] - t[1,0] ) * ( t[0,1] - t[0,0] ) \
          - ( t[1,1] - t[1,0] ) * ( t[0,2] - t[0,0] ) 

  center[0] = t[0,0] + 0.5 * top[0] / det
  center[1] = t[1,0] + 0.5 * top[1] / det

  return r, center

def triangle_circumcircle_test ( ):

#*****************************************************************************80
#
## triangle_circumcircle_test() tests triangle_circumcircle().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'triangle_circumcircle_test():' )
  print ( '  triangle_circumcircle() computes the circumcenter of a triangle.' )

  for i in range ( 0, 4 ):

    if ( i == 0 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.0 ], \
        [ 0.0, 0.0, 1.0 ] ] )
    elif ( i == 1 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 0.86602539 ] ] )
    elif ( i == 2 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 10.0 ] ] )
    elif ( i == 3 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 10.0 ], \
        [ 0.0, 0.0, 2.0 ] ] )

    r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    r, center = triangle_circumcircle ( t )

    print ( '' )
    print ( '  triangle_circumcenter:' )
    print ( center )

    print ( '' )
    print ( '  Circumradius: %g' % ( r ) )

  return

def triangle_contains_point_1 ( t, p ):

#*****************************************************************************80
#
## triangle_contains_point_1() finds if a point is inside a triangle.
#
#  Discussion:
#
#    It is conventional to list the triangle vertices in counter clockwise
#    order.  However, this routine does not require a particular order
#    for the vertices.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#    real P(2,1), the point to be checked.
#
#  Output:
#
#    bool INSIDE, is TRUE if the point is inside
#    the triangle or on its boundary.
#
  xsi = triangle_barycentric ( t, p )

  if ( xsi[0] < 0.0 or xsi[1] < 0.0 or xsi[2] < 0.0 ):
    inside = False
  else:
    inside = True

  return inside

def triangle_contains_point_1_test ( ):

#*****************************************************************************80
#
## triangle_contains_point_1_test() tests triangle_contains_point_1().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 June 2006
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  test_num = 7

  p_test = np.array ( [ \
    [ 0.25,   0.75,   1.00,  11.00,   0.00,   0.50,  0.60 ], \
    [ 0.25,   0.25,   1.00,   0.50,   1.00, -10.00,  0.60 ] ] )

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  p = np.zeros ( 2 )

  print ( '' )
  print ( 'triangle_contains_point_1_test()' )
  print ( '  triangle_contains_point_1() reports if a point' )
  print ( '  is inside a triangle' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '     X       Y     Inside' )
  print ( '' )

  for j in range ( 0, test_num ):

    p[0] = p_test[0,j]
    p[1] = p_test[1,j]

    inside = triangle_contains_point_1 ( t, p )

    print ( '  %10g  %10g  %s' % ( p[0], p[1], inside ) )
#
#  Make a copy of the triangle with vertices in reverse order.
#
  print ( '' )
  print ( '  Repeat the test, but reverse the triangle vertex ordering.' )
 
  t2 = np.zeros ( [ 2, 3 ] )
  for j in range ( 0, 3 ):
    for i in range ( 0, 2 ):
      t2[i,j] = t[i,2-j]

  r8mat_transpose_print ( 2, 3, t2, '  Triangle vertices (reversed):' )

  print ( '' )
  print ( '     X       Y     Inside' )
  print ( '' )

  for j in range ( 0, test_num ):

    p[0] = p_test[0,j]
    p[1] = p_test[1,j]

    inside = triangle_contains_point_1 ( t2, p )

    print ( '  %10g  %10g  %s' % ( p[0], p[1], inside ) )

  return

def triangle_contains_point ( t, p ):

#*****************************************************************************80
#
## triangle_contains_point() finds if a point is inside a triangle in 2D.
#
#  Discussion:
#
#    The routine assumes that the vertices are given in counter-clockwise
#    order.  If the triangle vertices are actually given in clockwise
#    order, this routine will behave as though the triangle contains
#    no points whatsoever!
#
#    The routine determines if a point P is "to the right of" each of the lines
#    that bound the triangle.  It does this by computing the cross product
#    of vectors from a vertex to its next vertex, and to P.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#    The vertices should be given in counter clockwise order.
#
#    real P(2,1), the point to be checked.
#
#  Output:
#
#    bool INSIDE, is TRUE if the point is inside
#    the triangle or on its boundary.
#
  inside = True

  for j in range ( 0, 3 ):

    jp1 = i4_wrap ( j + 1, 0, 2 )

    if ( 0.0 < ( p[0] - t[0,j] ) * ( t[1,jp1] - t[1,j] ) \
             - ( p[1] - t[1,j] ) * ( t[0,jp1] - t[0,j] ) ):
      inside = False
      return inside

  return inside

def triangle_contains_point_test ( ):

#*****************************************************************************80
#
## triangle_contains_point_test() tests triangle_contains_point().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 June 2006
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  test_num = 7

  p_test = np.array ( [ \
    [ 0.25,   0.75,   1.00,  11.00,   0.00,   0.50,  0.60 ], \
    [ 0.25,   0.25,   1.00,   0.50,   1.00, -10.00,  0.60 ] ] )

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  p = np.zeros ( 2 )

  print ( '' )
  print ( 'triangle_contains_point_test()' )
  print ( '  triangle_contains_point() reports if a point' )
  print ( '  is inside a triangle' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '     X       Y     Inside' )
  print ( '' )

  for j in range ( 0, test_num ):

    p[0] = p_test[0,j]
    p[1] = p_test[1,j]

    inside = triangle_contains_point ( t, p )

    print ( '  %10g  %10g  %s' % ( p[0], p[1], inside ) )
#
#  Make a copy of the triangle with vertices in reverse order.
#
  print ( '' )
  print ( '  Repeat the test, but reverse the triangle vertex ordering.' )
 
  t2 = np.zeros ( [ 2, 3 ] )
  for j in range ( 0, 3 ):
    for i in range ( 0, 2 ):
      t2[i,j] = t[i,2-j]

  r8mat_transpose_print ( 2, 3, t2, '  Triangle vertices (reversed):' )

  print ( '' )
  print ( '     X       Y     Inside' )
  print ( '' )

  for j in range ( 0, test_num ):

    p[0] = p_test[0,j]
    p[1] = p_test[1,j]

    inside = triangle_contains_point ( t2, p )

    print ( '  %10g  %10g  %s' % ( p[0], p[1], inside ) )

  return

def triangle_diameter ( t ):

#*****************************************************************************80
#
## triangle_diameter() computes the diameter of a triangle in 2D.
#
#  Discussion:
#
#    The diameter of a triangle is the diameter of the smallest circle
#    that can be drawn around the triangle.  At least two of the vertices
#    of the triangle will intersect the circle, but not necessarily
#    all three!
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real DIAMETER, the diameter of the triangle.
#
  import numpy as np
#
#  Compute the squared length of each side.
#
  asq = ( t[0,0] - t[0,1] ) ** 2 + ( t[1,0] - t[1,1] ) ** 2
  bsq = ( t[0,1] - t[0,2] ) ** 2 + ( t[1,1] - t[1,2] ) ** 2
  csq = ( t[0,2] - t[0,0] ) ** 2 + ( t[1,2] - t[1,0] ) ** 2
#
#  Take care of a zero side.
#
  if ( asq == 0.0 ):
    diameter = np.sqrt ( bsq )
    return diameter
  elif ( bsq == 0.0 ):
    diameter = np.sqrt ( csq )
    return diameter
  elif ( csq == 0.0 ):
    diameter = np.sqrt ( asq )
    return diameter
#
#  Make ASQ the largest.
#
  if ( asq < bsq ):
    temp = asq
    asq = bsq
    bsq = temp

  if ( asq < csq ):
    temp = asq
    asq = csq
    csq = temp
#
#  If ASQ is very large...
#
  if ( bsq + csq < asq ):

    diameter = np.sqrt ( asq )

  else:

    a = np.sqrt ( asq )
    b = np.sqrt ( bsq )
    c = np.sqrt ( csq )

    diameter = 2.0 * a * b * c / np.sqrt ( ( a + b + c ) * ( - a + b + c ) \
      * ( a - b + c ) * ( a + b - c ) )

  return diameter

def triangle_diameter_test ( ):

#*****************************************************************************80
#
## triangle_diameter_test() tests triangle_diameter().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
 
  print ( '' )
  print ( 'triangle_diameter_test():' )
  print ( '  triangle_diameter() computes the diameter of' )
  print ( '  the SMALLEST circle around a triangle.' )

  t = np.array ( [ \
    [ 4.0, 1.0, -2.0 ], \
    [ 2.0, 5.0,  2.0 ] ] )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  diameter = triangle_diameter ( t )

  print ( '' )
  print ( '  Diameter = %g' % ( diameter ) )

  t = np.array ( [ \
    [ 4.0, 5.0, 6.0 ], \
    [ 2.0, 4.0, 6.0 ] ] )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  diameter = triangle_diameter ( t )

  print ( '' )
  print ( '  Diameter = %g' % ( diameter ) )

  t = np.array ( [ \
    [ 4.0, 1.0, 4.0 ], \
    [ 2.0, 5.0, 2.0 ] ] )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  diameter = triangle_diameter ( t )

  print ( '' )
  print ( '  Diameter = %g' % ( diameter ) )

  return

def triangle_edge_length ( t ):

#*****************************************************************************80
#
## triangle_edge_length() returns edge lengths of a triangle in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real EDGE_length(3), the length of the edges.
#
  import numpy as np

  edge_length = np.zeros ( 3 )

  for j1 in range ( 0, 3 ):
    j2 = i4_wrap ( j1 + 1, 0, 2 )

    edge_length[j1] = np.sqrt ( ( t[0,j2] - t[0,j1] ) ** 2 \
                              + ( t[1,j2] - t[1,j1] ) ** 2 )

  return edge_length

def triangle_edge_length_test ( ):

#*****************************************************************************80
#
## triangle_edge_length_test() tests triangle_edge_length().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'triangle_edge_length_test():' )
  print ( '  triangle_edge_length() computes the edge lengths' )
  print ( '  of a triangle.' )

  t = np.array ( [ \
    [ 4.0, 1.0, -2.0 ], \
    [ 2.0, 5.0,  2.0 ] ] )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  edge_length = triangle_edge_length ( t )

  print ( '' )
  print ( '  triangle edge lengths:' )
  print ( edge_length )

  t = np.array ( [ \
    [ 4.0, 5.0, 6.0 ], \
    [ 2.0, 4.0, 6.0 ] ] )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  edge_length = triangle_edge_length ( t )

  print ( '' )
  print ( '  triangle edge lengths:' )
  print ( edge_length )

  t = np.array ( [ \
    [ 4.0, 1.0, 4.0 ], \
    [ 2.0, 5.0, 2.0 ] ] )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  edge_length = triangle_edge_length ( t )

  print ( '' )
  print ( '  triangle edge lengths:' )
  print ( edge_length )


  return

def triangle_incircle ( t ):

#*****************************************************************************80
#
## triangle_incircle() computes the inscribed circle of a triangle in 2D.
#
#  Discussion:
#
#    The inscribed circle of a triangle is the largest circle that can
#    be drawn inside the triangle.  It is tangent to all three sides,
#    and the lines from its center to the vertices bisect the angles
#    made by each vertex.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer and John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real R, CENTER(2), the radius and center of the
#    inscribed circle.
#
  import numpy as np

  center = np.zeros ( 2 )
#
#  Compute the length of each side.
#
  a = np.sqrt ( ( t[0,0] - t[0,1] ) ** 2 + ( t[1,0] - t[1,1] ) ** 2 )
  b = np.sqrt ( ( t[0,1] - t[0,2] ) ** 2 + ( t[1,1] - t[1,2] ) ** 2 )
  c = np.sqrt ( ( t[0,2] - t[0,0] ) ** 2 + ( t[1,2] - t[1,0] ) ** 2 )

  perimeter = a + b + c

  if ( perimeter == 0.0 ):
    center[0] = t[0,0]
    center[1] = t[1,0]
    r = 0.0
    return r, center

  center[0] = (  \
      b * t[0,0] \
    + c * t[0,1] \
    + a * t[0,2] ) / perimeter

  center[1] = (  \
      b * t[1,0] \
    + c * t[1,1] \
    + a * t[1,2] ) / perimeter

  r = 0.5 * np.sqrt ( \
      ( - a + b + c ) \
    * ( + a - b + c ) \
    * ( + a + b - c ) / perimeter )

  return r, center

def triangle_incircle_test ( ):

#*****************************************************************************80
#
## triangle_incircle_test() tests triangle_incircle().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle_incircle_test():' )
  print ( '  triangle_incircle() computes the incircle of a triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  r, center = triangle_incircle ( t )

  print ( '' )
  print ( '  triangle incenter:' )
  print ( center )

  print ( '' )
  print ( '  Incircle radius is %g' % ( r ) )

  return

def triangle_orientation ( t ):

#*****************************************************************************80
#
## triangle_orientation() determines the orientation of a triangle in 2D.
#
#  Discussion:
#
#    Three distinct non-colinear points in the plane define a circle.
#    If the points are visited in the order P1, P2, and then
#    P3, this motion defines a clockwise or counterclockwise
#    rotation along the circle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    integer VALUE, reports if the three points lie
#    clockwise on the circle that passes through them.  The possible
#    return values are:
#    0, the points are distinct, noncolinear, and lie counterclockwise
#    on their circle.
#    1, the points are distinct, noncolinear, and lie clockwise
#    on their circle.
#    2, the points are distinct and colinear.
#    3, at least two of the points are identical.
#
  for j in range ( 0, 3 ):
    jp1 = i4_wrap ( j + 1, 0, 2 )
    if ( t[0,j] == t[0,jp1] and t[1,j] == t[1,jp1] ): 
      value = 3
      return value

  det = ( t[0,0] - t[0,2] ) * ( t[1,1] - t[1,2] ) \
      - ( t[0,1] - t[0,2] ) * ( t[1,0] - t[1,2] )

  if ( det == 0.0 ):
    value = 2
  elif ( det < 0.0 ):
    value = 1
  elif ( 0.0 < det ):
    value = 0

  return value

def triangle_orientation_test ( ):

#*****************************************************************************80
#
## triangle_orientation_test() tests triangle_orientation().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
 
  print ( '' )
  print ( 'triangle_orientation_test():' )
  print ( '  triangle_orientation() determines orientation of a triangle.' )

  t = np.array ( [ \
    [ 4.0, 1.0, -2.0 ], \
    [ 2.0, 5.0,  2.0 ] ] )

  i = triangle_orientation ( t )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  if ( i == 0 ):
    print ( '  The points are counterclockwise.' )
  elif ( i == 1 ):
    print ( '  The points are clockwise.' )
  elif ( i == 2 ):
    print ( '  The points are colinear.' )
  elif ( i == 3 ):
    print ( '  The points are not distinct.' )
  else:
    print ( '  The return value makes no sense.' )

  t = np.array ( [ \
    [ 1.0, 4.0,  1.0 ], \
    [ 5.0, 2.0, -1.0 ] ] )

  i = triangle_orientation ( t )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  if ( i == 0 ):
    print ( '  The points are counterclockwise.' )
  elif ( i == 1 ):
    print ( '  The points are clockwise.' )
  elif ( i == 2 ):
    print ( '  The points are colinear.' )
  elif ( i == 3 ):
    print ( '  The points are not distinct.' )
  else:
    print ( '  The return value makes no sense.' )
#
#  Colinear points
#
  t = np.array ( [ \
    [ 1.0, 2.0, 3.0 ], \
    [ 5.0, 7.0, 9.0 ] ] )

  i = triangle_orientation ( t )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  if ( i == 0 ):
    print ( '  The points are counterclockwise.' )
  elif ( i == 1 ):
    print ( '  The points are clockwise.' )
  elif ( i == 2 ):
    print ( '  The points are colinear.' )
  elif ( i == 3 ):
    print ( '  The points are not distinct.' )
  else:
    print ( '  The return value makes no sense.' )
#
#  Nondistinct points.
#
  t = np.array ( [ \
    [ 1.0, 4.0, 1.0 ], \
    [ 5.0, 2.0, 5.0 ] ] )

  i = triangle_orientation ( t )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  if ( i == 0 ):
    print ( '  The points are counterclockwise.' )
  elif ( i == 1 ):
    print ( '  The points are clockwise.' )
  elif ( i == 2 ):
    print ( '  The points are colinear.' )
  elif ( i == 3 ):
    print ( '  The points are not distinct.' )
  else:
    print ( '  The return value makes no sense.' )

  return

def triangle_orthocenter ( t ):

#*****************************************************************************80
#
## triangle_orthocenter() computes the orthocenter of a triangle in 2D.
#
#  Discussion:
#
#    The orthocenter is defined as the intersection of the three altitudes
#    of a triangle.
#
#    An altitude of a triangle is the line through a vertex of the triangle
#    and perpendicular to the opposite side.
#
#    In geometry, the orthocenter of a triangle is often symbolized by "H".
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer and John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real CENTER(2,1), the orthocenter of the triangle.
#
#    bool FLAG, is TRUE if the value could not be computed.
#
  import numpy as np

  p1 = np.zeros ( 2 )
  p1[0] = t[0,0]
  p1[1] = t[1,0]
  p2 = np.zeros ( 2 )
  p2[0] = t[0,1]
  p2[1] = t[1,1]
  p3 = np.zeros ( 2 )
  p3[0] = t[0,2]
  p3[1] = t[1,2]
  center = np.zeros ( 2 )
#
#  Determine a point P23 common to the line (P2,P3) and
#  its perpendicular through P1.
#
  p23, flag = line_exp_perp ( p2, p3, p1 )

  if ( flag ):
    center[0] = float ( 'inf' )
    center[1] = float ( 'inf' )
    return center, flag
#
#  Determine a point P31 common to the line (P3,P1) and
#  its perpendicular through P2.
#
  p31, flag = line_exp_perp ( p3, p1, p2 )

  if ( flag ):
    center[0] = float ( 'inf' )
    center[1] = float ( 'inf' )
    return center, flag
#
#  Determine CENTER, the intersection of the lines (P1,P23) and (P2,P31).
#
  ival, center = lines_exp_int ( p1, p23, p2, p31 )

  if ( ival != 1 ):
    flag = 1
    center[0] = float ( 'inf' )
    center[1] = float ( 'inf' )
    return center, flag

  return center, flag

def triangle_orthocenter_test ( ):

#*****************************************************************************80
#
## triangle_orthocenter_test() tests triangle_orthocenter().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  print ( '' )
  print ( 'triangle_orthocenter_test():' )
  print ( '  triangle_orthocenter() computes the orthocenter of a triangle.' )

  for i in range ( 0, 4 ):

    if ( i == 0 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.0 ], \
        [ 0.0, 0.0, 1.0 ] ] )
    elif ( i == 1 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 0.86602539 ] ] )
    elif ( i == 2 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 10.0 ] ] )
    elif ( i == 3 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 10.0 ], \
        [ 0.0, 0.0, 2.0 ] ] )

    r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    center, flag = triangle_orthocenter ( t )

    print ( '' )
    print ( '  triangle orthocenter:' )
    print ( center )

  return

def triangle_point_dist ( t, p ):

#*****************************************************************************80
#
## triangle_point_dist(): distance ( triangle, point ) in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#    real P(2), the point to be checked.
#
#  Output:
#
#    real DIST, the distance from the point to the triangle.
#
  import numpy as np

  p1 = np.zeros ( 2 )
  p2 = np.zeros ( 2 )
#
#  Find the distance to each of the line segments.
#
  dist = float ( 'inf' )

  for j in range ( 0, 3 ):

    jp1 = i4_wrap ( j + 1, 0, 2 )

    p1[0] = t[0,j]
    p1[1] = t[1,j]

    p2[0] = t[0,jp1]
    p2[1] = t[1,jp1]

    dist2 = segment_point_dist ( p1, p2, p )

    if ( dist2 < dist ):
      dist = dist2

  return dist

def triangle_point_dist_test ( ):

#*****************************************************************************80
#
## triangle_point_dist_test() tests triangle_point_dist().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 7

  ptest = np.array ( [ \
    [ 0.25, 0.75, 1.00, 11.00, 0.00,  0.50,  0.60 ], \
    [ 0.25, 0.25, 1.00,  0.50, 1.00, -10.00, 0.60 ] ] )

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  p = np.zeros ( 2 )

  print ( '' )
  print ( 'triangle_point_dist_test():' )
  print ( '  triangle_point_dist() computes the distance' )
  print ( '  between a point and a triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '       P       DIST' )
  print ( '' )

  for j in range ( 0, ntest ):

    p[0] = ptest[0,j]
    p[1] = ptest[1,j]

    dist = triangle_point_dist ( t, p )

    print ( '  %10g  %10g  %10g' % ( p[0], p[1], dist ) )

  return

def triangle_point_near ( t, p ):

#*****************************************************************************80
#
## triangle_point_near() computes the nearest point on a triangle in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#    real P(2), the point whose nearest triangle point
#    is to be determined.
#
#  Output:
#
#    real PN(2), the nearest point to P.
#
#    real DIST, the distance from the point to the triangle.
#
  import numpy as np

  pn = np.zeros ( 2 )
  p1 = np.zeros ( 2 )
  p2 = np.zeros ( 2 )
#
#  Find the distance to each of the line segments that make up the edges
#  of the triangle.
#
  dist = float ( 'inf' )

  for j in range ( 0, 3 ):

    jp1 = i4_wrap ( j + 1, 0, 2 )

    p1[0] = t[0,j]
    p1[1] = t[1,j]


    p2[0] = t[0,jp1]
    p2[1] = t[1,jp1]

    pn2, dist2, tval = segment_point_near ( p1, p2, p )

    if ( dist2 < dist ):
      dist = dist2
      pn[0] = pn2[0]
      pn[1] = pn2[1]

  return pn, dist

def triangle_point_near_test ( ):

#*****************************************************************************80
#
## triangle_point_near_test() tests triangle_point_near().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 December 2010
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 7

  ptest = np.array ( [ \
    [ 0.25, 0.75, 1.00, 11.00, 0.00,   0.50, 0.60 ], \
    [ 0.25, 0.25, 1.00,  0.50, 1.00, -10.00, 0.60 ] ] )

  t = np.array ( [ \
    [ 0.0, 0.0, 1.0 ], \
    [ 1.0, 0.0, 0.0 ] ] )

  print ( '' )
  print ( 'triangle_point_near_test():' )
  print ( '  triangle_point_near() computes the nearest' )
  print ( '  triangle point to a point.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '           P                PN' )
  print ( '' )

  p = np.zeros ( 2 )

  for j in range ( 0, ntest ):

    p[0] = ptest[0,j]
    p[1] = ptest[1,j]

    [ pn, dist ] = triangle_point_near ( t, p )

    print ( '  %10g  %10g    %10g  %10g' % ( p[0], p[1], pn[0], pn[1] ) )

  return

def triangle_quality ( t ):

#*****************************************************************************80
#
## triangle_quality() computes the "quality" of a triangle in 2D.
#
#  Discussion:
#
#    The quality of a triangle is 2 times the ratio of the radius of the inscribed
#    circle divided by that of the circumscribed circle.  An equilateral
#    triangle achieves the maximum possible quality of 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Adrian Bowyer and John Woodwark,
#    A Programmer's Geometry,
#    Butterworths, 1983.
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#  Output:
#
#    real QUALITY, the quality of the triangle.
#
  import numpy as np
#
#  Compute the length of each side.
#
  a = np.sqrt ( ( t[0,0] - t[0,1] ) ** 2 + ( t[1,0] - t[1,1] ) ** 2 )
  b = np.sqrt ( ( t[0,1] - t[0,2] ) ** 2 + ( t[1,1] - t[1,2] ) ** 2 )
  c = np.sqrt ( ( t[0,2] - t[0,0] ) ** 2 + ( t[1,2] - t[1,0] ) ** 2 )

  if ( a * b * c == 0.0 ):
    value = 0.0
  else:
    value = ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) / ( a * b * c )

  return value

def triangle_quality_test ( ):

#*****************************************************************************80
#
## triangle_quality_test() tests triangle_quality().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  ntest = 4

  print ( '' )
  print ( 'triangle_quality_test():' )
  print ( '  triangle_quality() computes the quality of a triangle.' )

  for i in range ( 0, ntest ):
 
    if ( i == 0 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.0 ], \
        [ 0.0, 0.0, 1.0 ] ] )
    elif ( i == 1 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 0.5 ], \
        [ 0.0, 0.0, 0.86602539 ] ] )
    elif ( i == 2 ):
      t = np.array ( [ \
        [ 0.0, 1.0,  0.5 ], \
        [ 0.0, 0.0, 10.0 ] ] )
    elif ( i == 3 ):
      t = np.array ( [ \
        [ 0.0, 1.0, 10.0 ], \
        [ 0.0, 0.0, 2.0 ] ] )

    r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

    quality = triangle_quality ( t )

    print ( '' )
    print ( '  Quality = %g' % ( quality ) )

  return

def triangle_reference_sample ( n ):

#*****************************************************************************80
#
## triangle_reference_sample() returns random points in the reference triangle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 November 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of points to generate.
#
#  Output:
#
#    real P(2,N), random points in the triangle.
#
  import numpy as np

  alpha = np.random.rand ( n )
#
#  Interpret R as a percentage of the triangle's area.
#
#  Imagine a line L, parallel to side 1, so that the area between
#  vertex 1 and line L is R percent of the full triangle's area.
#
#  The line L will intersect sides 2 and 3 at a fraction
#  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
#
  for i in range ( 0, n ):
    alpha[i] = np.sqrt ( alpha[i] )
#
#  Now choose, uniformly at random, a point on the line L.
#
  beta = np.random.rand ( n )

  p = np.zeros ( [ 2, n ] )

  for j in range ( 0, n ):
    p[0,j] = ( 1.0 - beta[j] ) * alpha[j]
    p[1,j] =         beta[j]   * alpha[j]

  return p

def triangle_reference_sample_test ( ):

#*****************************************************************************80
#
## triangle_reference_sample_test() tests triangle_reference_sample().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 November 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 0.0, 1.0,  0.0 ], \
    [ 0.0, 0.0,  1.0 ] ] )

  print ( '' )
  print ( 'triangle_reference_sample_test():' )
  print ( '  triangle_reference_sample() samples the reference triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:' )
  print ( '' )

  for i in range ( 0, 10 ):

    p = triangle_reference_sample ( 1 )

    xsi = triangle_xy_to_xsi ( t, p )

    print ( '  %10g  %10g    %10g  %10g  %10g' % ( p[0], p[1], xsi[0], xsi[1], xsi[2] ) )

  return

def triangle_sample ( t, n ):

#*****************************************************************************80
#
## triangle_sample() returns random points in a triangle.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#    integer N, the number of points to generate.
#
#  Output:
#
#    real P(2,N), random points in the triangle.
#
  import numpy as np

  alpha = np.random.rand ( n )
#
#  Interpret R as a percentage of the triangle's area.
#
#  Imagine a line L, parallel to side 1, so that the area between
#  vertex 1 and line L is R percent of the full triangle's area.
#
#  The line L will intersect sides 2 and 3 at a fraction
#  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
#
  for i in range ( 0, n ):
    alpha[i] = np.sqrt ( alpha[i] )
#
#  Determine the coordinates of the points on sides 2 and 3 intersected
#  by line L.
#
  p12 = np.zeros ( [ 2, n ] )
  p13 = np.zeros ( [ 2, n ] )

  for i in range ( 0, 2 ):
    for j in range ( 0, n ):
      p12[i,j] = ( 1.0 - alpha[j] ) * t[i,0] \
                       + alpha[j]   * t[i,1]

      p13[i,j] = ( 1.0 - alpha[j] ) * t[i,0] \
                       + alpha[j]   * t[i,2]
#
#  Now choose, uniformly at random, a point on the line L.
#
  alpha = np.random.rand ( n )

  p = np.zeros ( [ 2, n ] )

  for i in range ( 0, 2 ):
    for j in range ( 0, n ):
      p[i,j] = ( 1.0 - alpha[j] ) * p12[i,j] \
                     + alpha[j]   * p13[i,j]

  return p

def triangle_sample_test ( ):

#*****************************************************************************80
#
## triangle_sample_test() tests triangle_sample().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 4.0, 1.0, -2.0 ], \
    [ 2.0, 5.0,  2.0 ] ] )

  print ( '' )
  print ( 'triangle_sample_test():' )
  print ( '  triangle_sample() samples a triangle.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:' )
  print ( '' )

  for i in range ( 0, 10 ):

    p = triangle_sample ( t, 1 )

    xsi = triangle_xy_to_xsi ( t, p )

    print ( '  %10g  %10g    %10g  %10g  %10g' % ( p[0], p[1], xsi[0], xsi[1], xsi[2] ) )

  return

def triangle_xsi_to_xy ( t, xsi ):

#*****************************************************************************80
#
## triangle_xsi_to_xy() converts from barycentric to XY coordinates in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#    real XSI(3,1), the barycentric coordinates of a point.
#    XSI(1) + XSI(2) + XSI(3) should equal 1, but this is not checked.
#
#  Output:
#
#    real P(2,1), the XY coordinates of the point.
#
  import numpy as np

  p = np.zeros ( 2 )

  p[0] = t[0,0] * xsi[0] + t[0,1] * xsi[1] + t[0,2] * xsi[2]
  p[1] = t[1,0] * xsi[0] + t[1,1] * xsi[1] + t[1,2] * xsi[2]

  return p

def triangle_xsi_to_xy_test ( ):

#*****************************************************************************80
#
## triangle_xsi_to_xy_test() tests triangle_xsi_to_xy().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 4.0, 1.0, -2.0 ], \
    [ 2.0, 5.0,  2.0 ] ] )

  print ( '' )
  print ( 'triangle_xsi_to_xy_test():' )
  print ( '  triangle_xsi_to_xy() converts XSI to XY coordinates.' )
  print ( '' )
  print ( '  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)' )
  print ( '  works properly.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '  Sample points:' )
  print ( '' )

  p = np.zeros ( 2 )

  for j in range ( 0, 10 ):

    if ( j == 0 ):
      p[0] = ( t[0,0] + t[0,1] + t[0,2] ) / 3.0
      p[1] = ( t[1,0] + t[1,1] + t[1,2] ) / 3.0
    elif ( j == 1 ):
      p[0] = 3.0
      p[1] = 0.0
    else:
      p = triangle_sample ( t, 1 )

    xsi = triangle_xy_to_xsi ( t, p )

    p2 = triangle_xsi_to_xy ( t, xsi )

    print ( '' )
    print ( '  %8g  %8g    %8g  %8g  %8g' % ( p[0], p[1], xsi[0], xsi[1], xsi[2] ) )
    print ( '  %8g  %8g' % ( p2[0], p2[1] ) )

  return

def triangle_xy_to_xsi ( t, p ):

#*****************************************************************************80
#
## triangle_xy_to_xsi() converts from XY to barycentric in 2D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T(2,3), the triangle vertices.
#
#    real P(2), the XY coordinates of a point.
#
#  Output:
#
#    real XSI(3), the barycentric coordinates of the point.
#    XSI1 + XSI2 + XSI3 should equal 1.
#
  import numpy as np

  xsi = np.zeros ( 3 )

  det = ( t[0,0] - t[0,2] ) * ( t[1,1] - t[1,2] ) \
      - ( t[0,1] - t[0,2] ) * ( t[1,0] - t[1,2] )

  xsi[0] = (   ( t[1,1] - t[1,2] ) * ( p[0] - t[0,2] ) \
             - ( t[0,1] - t[0,2] ) * ( p[1] - t[1,2] ) ) / det

  xsi[1] = ( - ( t[1,0] - t[1,2] ) * ( p[0] - t[0,2] ) \
             + ( t[0,0] - t[0,2] ) * ( p[1] - t[1,2] ) ) / det

  xsi[2] = 1.0 - xsi[0] - xsi[1]

  return xsi

def triangle_xy_to_xsi_test ( ):

#*****************************************************************************80
#
## triangle_xy_to_xsi_test() tests triangle_xy_to_xsi().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 October 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np

  t = np.array ( [ \
    [ 4.0, 1.0, -2.0 ], \
    [ 2.0, 5.0,  2.0 ] ] )

  print ( '' )
  print ( 'triangle_xy_to_xsi_test():' )
  print ( '  triangle_xy_to_xsi() converts XY to XSI coordinates.' )
  print ( '  Verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)' )
  print ( '  works properly.' )

  r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  print ( '' )
  print ( '  Sample points:' )
  print ( '' )

  p = np.zeros ( 2 )

  for j in range ( 0, 10 ):

    if ( j == 0 ):
      p[0] = ( t[0,0] + t[0,1] + t[0,2] ) / 3.0
      p[1] = ( t[1,0] + t[1,1] + t[1,2] ) / 3.0
    elif ( j == 1 ):
      p[0] = 3.0
      p[1] = 0.0
    else:
      p = triangle_sample ( t, 1 )

    xsi = triangle_xy_to_xsi ( t, p )

    p2 = triangle_xsi_to_xy ( t, xsi )

    print ( '' )
    print ( '  %8g  %8g    %8g  %8g  %8g' % ( p[0], p[1], xsi[0], xsi[1], xsi[2] ) )
    print ( '  %8g  %8g' % ( p2[0], p2[1] ) )

  return

def wedge01_volume ( ):

#*****************************************************************************80
#
## wedge01_volume() returns the volume of the unit wedge in 3D.
#
#  Discussion:
#
#    The unit wedge is:
#
#      0 <= X
#      0 <= Y
#      X + Y <= 1
#      -1 <= Z <= 1.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    23 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Output:
#
#    real VALUE, the volume of the unit wedge.
#
  value = 1.0

  return value

def wedge01_volume_test ( ) :

#*****************************************************************************80
#
## wedge01_volume_test() tests wedge01_volume().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 June 2015
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'wedge01_volume_test():' )
  print ( '  wedge01_volume() returns the volume of the unit wedge.' )

  value = wedge01_volume ( )

  print ( '' )
  print ( '  wedge01_volume() = %g' % ( value ) )

  return

if ( __name__ == '__main__' ):
  timestamp ( )
  geometry_test ( )
  timestamp ( )
