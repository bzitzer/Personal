

#include "VASlalib.h"
#include <iostream>
#include <cstring> 
#include <sys/time.h>

using namespace std;

double slaDsep(double a1, double b1, double a2, double b2)
/*
**  - - - - - - - -
**   s l a D s e p
**  - - - - - - - -
**
**  Angle between two points on a sphere.
**
**  (double precision)
**
**  Given:
**     a1,b1    double    spherical coordinates of one point
**     a2,b2    double    spherical coordinates of the other point
**
**  (The spherical coordinates are [RA,Dec], [Long,Lat] etc, in radians.)
**
**  The result is the angle, in radians, between the two points.  It
**  is always positive.
**
**  Called:  slaDcs2c, slaDsepv
**
**  Last revision:   7 May 2000
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
	double v1[3], v2[3];
	
	/* Convert coordinates from spherical to Cartesian. */
	slaDcs2c(a1, b1, v1);
	slaDcs2c(a2, b2, v2);
	
	/* Angle between the vectors. */
	return slaDsepv(v1, v2);
}

void slaDcs2c(double a, double b, double v[3])
/*
**  - - - - - - - - -
**   s l a D c s 2 c
**  - - - - - - - - -
**
**  Spherical coordinates to direction cosines.
**
**  (double precision)
**
**  Given:
**     a,b       double      spherical coordinates in radians
**                           (RA,Dec), (long,lat) etc
**
**  Returned:
**     v         double[3]   x,y,z unit vector
**
**  The spherical coordinates are longitude (+ve anticlockwise
**  looking from the +ve latitude pole) and latitude.  The
**  Cartesian coordinates are right handed, with the x axis
**  at zero longitude and latitude, and the z axis at the
**  +ve latitude pole.
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
	double cosb;
	
	cosb = cos(b);
	v[0] = cos(a) * cosb;
	v[1] = sin(a) * cosb;
	v[2] = sin(b);
	return;
}

double slaDsepv(double v1[3], double v2[3])
/*
**  - - - - - - - - -
**   s l a D s e p v
**  - - - - - - - - -
**
**  Angle between two vectors.
**
**  (double precision)
**
**  Given:
**     v1      double[3]    first vector
**     v2      double[3]    second vector
**
**  The result is the angle, in radians, between the two vectors.  It
**  is always positive.
**
**  Notes:
**
**  1  There is no requirement for the vectors to be unit length.
**
**  2  If either vector is null, zero is returned.
**
**  3  The simplest formulation would use dot product alone.  However,
**     this would reduce the accuracy for angles near zero and pi.  The
**     algorithm uses both cross product and dot product, which maintains
**     accuracy for all sizes of angle.
**
**  Called:  slaDvxv, slaDvn, slaDvdv
**
**  Last revision:   7 May 2000
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
	double v1xv2[3], wv[3], s, c;
	
	/* Modulus of cross product = sine multiplied by the two moduli. */
	slaDvxv(v1, v2, v1xv2);
	slaDvn(v1xv2, wv, &s);
	
	/* Dot product = cosine multiplied by the two moduli. */
	c = slaDvdv(v1, v2);
	
	/* Angle between the vectors. */
	return s != 0.0 ? atan2(s, c) : 0.0;
}

void slaDvxv(double va[3], double vb[3], double vc[3])
/*
**  - - - - - - - -
**   s l a D v x v
**  - - - - - - - -
**
**  Vector product of two 3-vectors.
**
**  (double precision)
**
**  Given:
**     va      double[3]     first vector
**     vb      double[3]     second vector
**
**  Returned:
**     vc      double[3]     vector result
**
**  Note:  the same vector may be specified more than once.
**
**  Last revision:   6 November 1999
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
	double vw[3];
	int i;
	
	/* Form the vector product va cross vb */
	vw[0] = va[1] * vb[2] - va[2] * vb[1];
	vw[1] = va[2] * vb[0] - va[0] * vb[2];
	vw[2] = va[0] * vb[1] - va[1] * vb[0];
	
	/* Return the result */
	for(i = 0; i < 3; i++)
	{
		vc[i] = vw[i];
	}
	return;
}

void slaDvn(double v[3], double uv[3], double* vm)
/*
**  - - - - - - -
**   s l a D v n
**  - - - - - - -
**
**  Normalizes a 3-vector also giving the modulus.
**
**  (double precision)
**
**  Given:
**     v       double[3]      vector
**
**  Returned:
**     uv      double[3]      unit vector in direction of v
**     *vm     double         modulus of v
**
**  Note:  v and uv may be the same array.
**
**  If the modulus of v is zero, uv is set to zero as well.
**
**  Last revision:   6 December 2001
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
	int i;
	double w1, w2;
	
	/* Modulus */
	w1 = 0.0;
	for(i = 0; i < 3; i++)
	{
		w2 = v[i];
		w1 += w2 * w2;
	}
	w1 = sqrt(w1);
	*vm = w1;
	
	/* Normalize the vector */
	w1 = (w1 > 0.0) ? w1 : 1.0;
	
	for(i = 0; i < 3; i++)
	{
		uv[i] = v[i] / w1;
	}
	return;
}

double slaDvdv(double va[3], double vb[3])
/*
**  - - - - - - - -
**   s l a D v d v
**  - - - - - - - -
**
**  Scalar product of two 3-vectors.
**
**  (double precision)
**
**
**  Given:
**      va      double(3)     first vector
**      vb      double(3)     second vector
**
**
**  The result is the scalar product va.vb (double precision)
**
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
	return va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
}
