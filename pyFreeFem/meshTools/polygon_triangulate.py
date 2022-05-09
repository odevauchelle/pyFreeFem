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
#      /7
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
  import platform

  n_angle = 12

  print ( '' )
  print ( 'angle_degree_test:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  angle_degree() computes an angle;' )
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
#    bool VALUE, is TRUE if the points are judged
#    to be collinear.
#
  r8_eps = 2.220446049250313E-016

  area = triangle_area ( xa, ya, xb, yb, xc, yc )

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
    print ( '  Col: ' ),

    for j in range ( j2lo, j2hi + 1 ):
      print ( '%7d  ' % ( j ) ),

    print ( '' )
    print ( '  Row' )

    i2lo = max ( ilo, 0 )
    i2hi = min ( ihi, m )

    for i in range ( i2lo, i2hi + 1 ):

      print ( ' %4d: ' % ( i ) ),

      for j in range ( j2lo, j2hi + 1 ):
        print ( '%7d  ' % ( a[i,j] ) ),

      print ( '' )

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

  t1 = triangle_area ( x[im1], y[im1], x[i], y[i], x[im2], y[im2] )

  if ( 0.0 <= t1 ):

    t2 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] )
    t3 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] )
    value = ( ( 0.0 < t2 ) and ( 0.0 < t3 ) )

  else:

    t4 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] )
    t5 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] )
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
    t1 = triangle_area ( xa, ya, xb, yb, xc, yc )
    t2 = triangle_area ( xa, ya, xb, yb, xd, yd )
    t3 = triangle_area ( xc, yc, xd, yd, xa, ya )
    t4 = triangle_area ( xc, yc, xd, yd, xb, yb )

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
#    An L4 is a logical value.
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

def polygon_area ( n, x, y ):

#*****************************************************************************80
#
## polygon_area() returns the area of a polygon.
#
#  Discussion:
#
#    The vertices should be listed in counter-clockwise order so that
#    the area will be positive.
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
#    John Burkardt.
#
#  Input:
#
#    integer N, the number of vertices.
#
#    real X(N), Y(N), the vertex coordinates.
#
#  Output:
#
#    real AREA, the area of the polygon.
#
  area = 0.0
  im1 = n - 1

  for i in range ( 0, n ):
    area = area + x[im1] * y[i] - x[i] * y[im1]
    im1 = i

  area = 0.5 * area

  return area

def polygon_triangulate ( n, x, y ):

#******************************************************************************/
#
## polygon_triangulate() determines a triangulation of a polygon.
#
#  Discussion:
#
#    There are N-3 triangles in the triangulation.
#
#    For the first N-2 triangles, the first edge listed is always an
#    internal diagonal.
#
#    Thanks to Gene Dial for pointing out a mistake in the area calculation,
#    10 September 2016.
#
#    Thanks to Gene Dial for suggesting that the next() array should be
#    renamed next_node() to avoid the Python keyword next, 22 September 2016.
#
#    Gene Dial requested an angle tolerance of about 1 millionth radian or
#    5.7E-05 degrees, 26 June 2018.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    26 June 2018
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
#    integer N, the number of vertices.
#
#    real X[N], Y[N], the coordinates of each vertex.
#
#  Output:
#
#    integer TRIANGLES[N-2,3], the triangles.
#
  import numpy as np

  angle_tol = 5.7E-05
#
#  We must have at least 3 vertices.
#
  if ( n < 3 ):
    print ( '' )
    print ( 'polygon_triangulate - Fatal error!' )
    print ( '  N < 3' )
    raise Exception ( 'polygon_triangulate - Fatal error!' )
#
#  Consecutive vertices cannot be equal.
#
  node_m1 = n - 1
  for node in range ( 0, n ):
    if ( x[node_m1] == x[node] and y[node_m1] == y[node] ):
      print ( '' )
      print ( 'polygon_triangulate - Fatal error!' )
      print ( '  Two consecutive nodes are identical.' )
      raise Exception ( 'polygon_triangulate - Fatal error!' )

    node_m1 = node
#
#  No node can be the vertex of an angle less than 1 degree
#  in absolute value.
#
  node1 = n - 1

  for node2 in range ( 0, n ):

    node3 = ( ( node2 + 1 ) % n )

    angle = angle_degree ( \
      x[node1], y[node1], \
      x[node2], y[node2], \
      x[node3], y[node3] )

    if ( abs ( angle ) <= angle_tol ):
      print ( '' )
      print ( 'polygon_triangulate - Fatal error!' )
      print ( '  Polygon has an angle %g smaller than %g degrees.' \
        % ( angle, angle_tol ) )
      print ( '  occurring at node %d' % ( node2 ) )
      raise Exception ( 'polygon_triangulate - Fatal error!' )

    node1 = node2
#
#  Area must be positive.
#
  area = polygon_area ( n, x, y )

  if ( area <= 0.0 ):
    print ( '' )
    print ( 'polygon_triangulate - Fatal error!' )
    print ( '  Polygon has zero or negative area.' )
    raise Exception ( 'polygon_triangulate - Fatal error!' )

  triangles = np.zeros ( [ n - 2,  3 ], dtype = np.int32 )
#
#  PREV_NODE and NEXT_NODE point to the previous and next nodes.
#
  prev_node = np.zeros ( n, dtype = np.int32 )
  next_node = np.zeros ( n, dtype = np.int32 )

  i = 0
  prev_node[i] = n - 1
  next_node[i] = i + 1

  for i in range ( 1, n - 1 ):
    prev_node[i] = i - 1
    next_node[i] = i + 1

  i = n - 1
  prev_node[i] = i - 1
  next_node[i] = 0
#
#  EAR indicates whether the node and its immediate neighbors form an ear
#  that can be sliced off immediately.
#
  ear = np.zeros ( n, dtype = np.bool )
  for i in range ( 0, n ):
    ear[i] = diagonal ( prev_node[i], next_node[i], n, prev_node, \
      next_node, x, y )

  triangle_num = 0

  i2 = 0

  while ( triangle_num < n - 3 ):
#
#  If I2 is an ear, gather information necessary to carry out
#  the slicing operation and subsequent "healing".
#
    if ( ear[i2] ):
      i3 = next_node[i2]
      i4 = next_node[i3]
      i1 = prev_node[i2]
      i0 = prev_node[i1]
#
#  Make vertex I2 disappear.
#
      next_node[i1] = i3
      prev_node[i3] = i1
#
#  Update the earity of I1 and I3, because I2 disappeared.
#
      ear[i1] = diagonal ( i0, i3, n, prev_node, next_node, x, y )
      ear[i3] = diagonal ( i1, i4, n, prev_node, next_node, x, y )
#
#  Add the diagonal [I3, I1, I2] to the list.
#
      triangles[triangle_num,0] = i3
      triangles[triangle_num,1] = i1
      triangles[triangle_num,2] = i2
      triangle_num = triangle_num + 1
#
#  Try the next vertex.
#
    i2 = next_node[i2]
#
#  The last triangle is formed from the three remaining vertices.
#
  i3 = next_node[i2]
  i1 = prev_node[i2]

  triangles[triangle_num,0] = i3
  triangles[triangle_num,1] = i1
  triangles[triangle_num,2] = i2
  triangle_num = triangle_num + 1

  return triangles

def polygon_triangulate_comb_test ( ):

#*****************************************************************************80
#
## polygon_triangulate_comb_test() tests the "comb_10" polygon.
#
#  Discussion:
#
#    There are N-3 triangles in the triangulation.
#
#    For the first N-2 triangles, the first edge listed is always an
#    internal diagonal.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    14 October 2015
#
  import numpy as np
  import platform

  n = 10
  x = np.array ( [ 8.0,  7.0,  6.0,  5.0,  4.0,  3.0,  2.0,  1.0,  0.0,  4.0 ] )
  y = np.array ( [ 0.0, 10.0,  0.0, 10.0,  0.0, 10.0,  0.0, 10.0,  0.0, -2.0 ] )

  print ( '' )
  print ( 'polygon_triangulate_comb_test' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Triangulate the comb_10 polygon.' )

  triangles = polygon_triangulate ( n, x, y )

  i4mat_print ( n - 2, 3, triangles, '  Triangles' )
#
#  Terminate.
#
  print ( '' )
  print ( 'polygon_triangulate_comb_test:' )
  print ( '  Normal end of execution.' )
  return

def polygon_triangulate_test ( ):

#*****************************************************************************80
#
## polygon_triangulate_test() tests polygon_triangulate().
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
  print ( 'polygon_triangulate_test():' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test polygon_triangulate()' )

  angle_degree_test ( )
  polygon_triangulate_comb_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'polygon_triangulate_test:' )
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

def triangle_area ( xa, ya, xb, yb, xc, yc ):

#*****************************************************************************80
#
## triangle_area() computes the signed area of a triangle.
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

if ( __name__ == '__main__' ):
  timestamp ( )
  polygon_triangulate_test ( )
  timestamp ( )
