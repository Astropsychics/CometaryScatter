# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

//****************************************************************************80

int i4_div_rounded ( int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    I4_DIV_ROUNDED computes the rounded result of I4 division.
//
//  Discussion:
//
//    This routine computes C = A / B, where A, B and C are integers
//    and C is the closest integer value to the exact real result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the number to be divided,
//    and the divisor.
//
//    Output, int I4_DIV_ROUNDED, the rounded result
//    of the division.
//
{
  int a_abs;
  int b_abs;
  static int i4_huge = 2147483647;
  int value;

  if ( a == 0 && b == 0 )
  {
    value = i4_huge;
  }
  else if ( a == 0 )
  {
    value = 0;
  }
  else if ( b == 0 )
  {
    if ( a < 0 )
    {
      value = - i4_huge;
    }
    else
    {
      value = + i4_huge;
    }
  }
  else
  {
    a_abs = abs ( a );
    b_abs = abs ( b );

    value = a_abs / b_abs;
//
//  Round the value.
//
    if ( ( 2 * value + 1 ) * b_abs < 2 * a_abs )
    {
      value = value + 1;
    }
//
//  Set the sign.
//
    if ( ( a < 0 && 0 < b ) || ( 0 < a && b < 0 ) )
    {
      value = - value;
    }
  }
  return value;
}
