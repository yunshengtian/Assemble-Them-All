///////////////////////////////////////////////////////////////////////////////
// From ODE's box.cpp
///////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <math.h>
#include "odeBoxBox.h"

namespace ode {
/* debugging:
 *   IASSERT  is an internal assertion, i.e. a consistency check. if it fails
 *            we want to know where.
 *   UASSERT  is a user assertion, i.e. if it fails a nice error message
 *            should be printed for the user.
 *   AASSERT  is an arguments assertion, i.e. if it fails "bad argument(s)"
 *            is printed.
 *   DEBUGMSG just prints out a message
 */

#  if defined(__STDC__) && __STDC_VERSION__ >= 199901L
#    define __FUNCTION__ __func__
#  endif
#ifndef dNODEBUG
#  ifdef __GNUC__
#    define dIASSERT(a) { if (!(a)) { dDebug (d_ERR_IASSERT, \
			"assertion \"" #a "\" failed in %s() [%s:%u]",__FUNCTION__,__FILE__,__LINE__); } }
#    define dUASSERT(a,msg) { if (!(a)) { dDebug (d_ERR_UASSERT, \
			msg " in %s()", __FUNCTION__); } }
#    define dDEBUGMSG(msg) { dMessage (d_ERR_UASSERT,				\
	msg " in %s() [%s:%u]", __FUNCTION__,__FILE__,__LINE__); }
#  else // not __GNUC__
#    define dIASSERT(a) { if (!(a)) { dDebug (d_ERR_IASSERT, \
			"assertion \"" #a "\" failed in %s:%u",__FILE__,__LINE__); } }
#    define dUASSERT(a,msg) { if (!(a)) { dDebug (d_ERR_UASSERT, \
			msg " (%s:%u)", __FILE__,__LINE__); } }
#    define dDEBUGMSG(msg) { dMessage (d_ERR_UASSERT, \
			msg " (%s:%u)", __FILE__,__LINE__); }
#  endif
#  define dIVERIFY(a) dIASSERT(a)
#else
#  define dIASSERT(a) ((void)0)
#  define dUASSERT(a,msg) ((void)0)
#  define dDEBUGMSG(msg) ((void)0)
#  define dIVERIFY(a) ((void)(a))
#endif

#  ifdef __GNUC__
#    define dICHECK(a) { if (!(a)) { dDebug (d_ERR_IASSERT, \
			"assertion \"" #a "\" failed in %s() [%s:%u]",__FUNCTION__,__FILE__,__LINE__); *(int *)0 = 0; } }
#  else // not __GNUC__
#    define dICHECK(a) { if (!(a)) { dDebug (d_ERR_IASSERT, \
			"assertion \"" #a "\" failed in %s:%u",__FILE__,__LINE__); *(int *)0 = 0; } }
#  endif

// Argument assert is a special case of user assert
#define dAASSERT(a) dUASSERT(a,"Bad argument(s)")


/* Define the dInfinity macro */
#ifdef INFINITY
	#define dInfinity INFINITY
#elif defined(HUGE_VAL)
	#ifdef dSINGLE
		#ifdef HUGE_VALF
			#define dInfinity HUGE_VALF
		#else
			#define dInfinity ((float)HUGE_VAL)
		#endif
	#else
		#define dInfinity HUGE_VAL
	#endif
#else
	#ifdef dSINGLE
		#define dInfinity ((float)(1.0/0.0))
	#else
		#define dInfinity (1.0/0.0)
	#endif
#endif


typedef double dReal;

/* these types are mainly just used in headers */
typedef dReal dVector3[4];
typedef dReal dVector4[4];
typedef dReal dMatrix3[4*3];
typedef dReal dMatrix4[4*4];
typedef dReal dMatrix6[8*6];
typedef dReal dQuaternion[4];

#define REAL(x) (x)
#define dRecip(x) (1.0/(x))
#define dSqrt(x) sqrt(x)
#define dRecipSqrt(x) (1.0/sqrt(x))
#define dSin(x) sin(x)
#define dCos(x) cos(x)
#define dFabs(x) fabs(x)
#define dAtan2(y,x) atan2((y),(x))
#define dFMod(a,b) (fmod((a),(b)))
#define dFloor(x) floor(x)
#define dCeil(x) ceil(x)
#define dCopySign(a,b) (copysign((a),(b)))
#define dNextAfter(x, y) nextafter(x, y)


static __inline dReal _dCalcVectorDot3(const dReal *a, const dReal *b, unsigned step_a, unsigned step_b)
{
	return a[0] * b[0] + a[step_a] * b[step_b] + a[2 * step_a] * b[2 * step_b];
}

static __inline dReal dCalcVectorDot3    (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,1,1); }
static __inline dReal dCalcVectorDot3_13 (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,1,3); }
static __inline dReal dCalcVectorDot3_31 (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,3,1); }
static __inline dReal dCalcVectorDot3_33 (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,3,3); }
static __inline dReal dCalcVectorDot3_14 (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,1,4); }
static __inline dReal dCalcVectorDot3_41 (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,4,1); }
static __inline dReal dCalcVectorDot3_44 (const dReal *a, const dReal *b) { return _dCalcVectorDot3(a,b,4,4); }

static __inline void dMultiplyHelper1_331(dReal *res, const dReal *a, const dReal *b)
{
	dReal res_0, res_1, res_2;
	res_0 = dCalcVectorDot3_41(a, b);
	res_1 = dCalcVectorDot3_41(a + 1, b);
	res_2 = dCalcVectorDot3_41(a + 2, b);
	// Only assign after all the calculations are over to avoid incurring memory aliasing
	res[0] = res_0; res[1] = res_1; res[2] = res_2;
}

static __inline void dMultiply1_331(dReal *res, const dReal *a, const dReal *b)
{
	dMultiplyHelper1_331(res, a, b);
}

static __inline void dMultiplyHelper0_331(dReal *res, const dReal *a, const dReal *b)
{
	dReal res_0, res_1, res_2;
	res_0 = dCalcVectorDot3(a, b);
	res_1 = dCalcVectorDot3(a + 4, b);
	res_2 = dCalcVectorDot3(a + 8, b);
	// Only assign after all the calculations are over to avoid incurring memory aliasing
	res[0] = res_0; res[1] = res_1; res[2] = res_2;
}

static __inline void dMultiply0_331(dReal *res, const dReal *a, const dReal *b)
{
	dMultiplyHelper0_331(res, a, b);
}


#define NUMC_MASK (0xffff)

#define CONTACTS_UNIMPORTANT			0x80000000


typedef struct dContactGeom {
	dVector3 pos;          ///< contact position
	dVector3 normal;       ///< normal vector
	dReal depth;           ///< penetration depth
} dContactGeom;

//****************************************************************************
// box-box collision utility


// given two lines
//    qa = pa + alpha* ua
//    qb = pb + beta * ub
// where pa,pb are two points, ua,ub are two unit length vectors, and alpha,
// beta go from [-inf,inf], return alpha and beta such that qa and qb are
// as close as possible
void dLineClosestApproach (const dVector3 pa, const dVector3 ua,
							 const dVector3 pb, const dVector3 ub,
							 dReal *alpha, dReal *beta)
{
	dVector3 p;
	p[0] = pb[0] - pa[0];
	p[1] = pb[1] - pa[1];
	p[2] = pb[2] - pa[2];
	dReal uaub = dCalcVectorDot3(ua,ub);
	dReal q1 =  dCalcVectorDot3(ua,p);
	dReal q2 = -dCalcVectorDot3(ub,p);
	dReal d = 1-uaub*uaub;
	if (d <= REAL(0.0001)) {
		// @@@ this needs to be made more robust
		*alpha = 0;
		*beta  = 0;
	}
	else {
		d = dRecip(d);
		*alpha = (q1 + uaub*q2)*d;
		*beta  = (uaub*q1 + q2)*d;
	}
}



// find all the intersection points between the 2D rectangle with vertices
// at (+/-h[0],+/-h[1]) and the 2D quadrilateral with vertices (p[0],p[1]),
// (p[2],p[3]),(p[4],p[5]),(p[6],p[7]).
//
// the intersection points are returned as x,y pairs in the 'ret' array.
// the number of intersection points is returned by the function (this will
// be in the range 0 to 8).

static int intersectRectQuad (dReal h[2], dReal p[8], dReal ret[16])
{
	// q (and r) contain nq (and nr) coordinate points for the current (and
	// chopped) polygons
	int nq=4,nr;
	dReal buffer[16];
	dReal *q = p;
	dReal *r = ret;
	for (int dir=0; dir <= 1; dir++) {
		// direction notation: xy[0] = x axis, xy[1] = y axis
		for (int sign=-1; sign <= 1; sign += 2) {
			// chop q along the line xy[dir] = sign*h[dir]
			dReal *pq = q;
			dReal *pr = r;
			nr = 0;
			for (int i=nq; i > 0; i--) {
		// go through all points in q and all lines between adjacent points
		if (sign*pq[dir] < h[dir]) {
			// this point is inside the chopping line
			pr[0] = pq[0];
			pr[1] = pq[1];
			pr += 2;
			nr++;
			if (nr & 8) {
				q = r;
				goto done;
			}
		}
		dReal *nextq = (i > 1) ? pq+2 : q;
		if ((sign*pq[dir] < h[dir]) ^ (sign*nextq[dir] < h[dir])) {
			// this line crosses the chopping line
			pr[1-dir] = pq[1-dir] + (nextq[1-dir]-pq[1-dir]) /
				(nextq[dir]-pq[dir]) * (sign*h[dir]-pq[dir]);
			pr[dir] = sign*h[dir];
			pr += 2;
			nr++;
			if (nr & 8) {
				q = r;
				goto done;
			}
		}
		pq += 2;
			}
			q = r;
			r = (q==ret) ? buffer : ret;
			nq = nr;
		}
	}
 done:
	if (q != ret) memcpy (ret,q,nr*2*sizeof(dReal));
	return nr;
}


// given n points in the plane (array p, of size 2*n), generate m points that
// best represent the whole set. the definition of 'best' here is not
// predetermined - the idea is to select points that give good box-box
// collision detection behavior. the chosen point indexes are returned in the
// array iret (of size m). 'i0' is always the first entry in the array.
// n must be in the range [1..8]. m must be in the range [1..n]. i0 must be
// in the range [0..n-1].

void cullPoints (int n, dReal p[], int m, int i0, int iret[])
{
	// compute the centroid of the polygon in cx,cy
	int i,j;
	dReal a,cx,cy,q;
	if (n==1) {
		cx = p[0];
		cy = p[1];
	}
	else if (n==2) {
		cx = REAL(0.5)*(p[0] + p[2]);
		cy = REAL(0.5)*(p[1] + p[3]);
	}
	else {
		a = 0;
		cx = 0;
		cy = 0;
		for (i=0; i<(n-1); i++) {
			q = p[i*2]*p[i*2+3] - p[i*2+2]*p[i*2+1];
			a += q;
			cx += q*(p[i*2]+p[i*2+2]);
			cy += q*(p[i*2+1]+p[i*2+3]);
		}
		q = p[n*2-2]*p[1] - p[0]*p[n*2-1];
		a = dRecip(REAL(3.0)*(a+q));
		cx = a*(cx + q*(p[n*2-2]+p[0]));
		cy = a*(cy + q*(p[n*2-1]+p[1]));
	}

	// compute the angle of each point w.r.t. the centroid
	dReal A[8];
	for (i=0; i<n; i++) A[i] = dAtan2(p[i*2+1]-cy,p[i*2]-cx);

	// search for points that have angles closest to A[i0] + i*(2*pi/m).
	int avail[8];
	for (i=0; i<n; i++) avail[i] = 1;
	avail[i0] = 0;
	iret[0] = i0;
	iret++;
	for (j=1; j<m; j++) {
		a = (dReal)(dReal(j)*(2*M_PI/m) + A[i0]);
		if (a > M_PI) a -= (dReal)(2*M_PI);
		dReal maxdiff=1e9,diff;
#ifndef dNODEBUG
		*iret = i0;			// iret is not allowed to keep this value
#endif
		for (i=0; i<n; i++) {
			if (avail[i]) {
		diff = dFabs (A[i]-a);
		if (diff > M_PI) diff = (dReal) (2*M_PI - diff);
		if (diff < maxdiff) {
			maxdiff = diff;
			*iret = i;
		}
			}
		}
		avail[*iret] = 0;
		iret++;
	}
}


// given two boxes (p1,R1,side1) and (p2,R2,side2), collide them together and
// generate contact points. this returns 0 if there is no contact otherwise
// it returns the number of contacts generated.
// `normal' returns the contact normal.
// `depth' returns the maximum penetration depth along that normal.
// `return_code' returns a number indicating the type of contact that was
// detected:
//        1,2,3 = box 2 intersects with a face of box 1
//        4,5,6 = box 1 intersects with a face of box 2
//        7..15 = edge-edge contact
// `maxc' is the maximum number of contacts allowed to be generated, i.e.
// the size of the `contact' array.
// `contact' and `skip' are the contact array information provided to the
// collision functions. this function only fills in the position and depth
// fields.
int dBoxBox (const dVector3 p1, const dMatrix3 R1,
				 const dVector3 side1, const dVector3 p2,
				 const dMatrix3 R2, const dVector3 side2,
				 dVector3 normal, dReal *depth, int *return_code,
				 int flags, dContactGeom *contact)
{
	const dReal fudge_factor = REAL(1.05);
	dVector3 p,pp,normalC={0,0,0};
	const dReal *normalR = 0;
	dReal A[3],B[3],R11,R12,R13,R21,R22,R23,R31,R32,R33,
		Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33,s,s2,l,expr1_val;
	int i,j,invert_normal,code;

	// get vector from centers of box 1 to box 2, relative to box 1
	p[0] = p2[0] - p1[0];
	p[1] = p2[1] - p1[1];
	p[2] = p2[2] - p1[2];
	dMultiply1_331 (pp,R1,p);		// get pp = p relative to body 1

	// get side lengths / 2
	A[0] = side1[0]*REAL(0.5);
	A[1] = side1[1]*REAL(0.5);
	A[2] = side1[2]*REAL(0.5);
	B[0] = side2[0]*REAL(0.5);
	B[1] = side2[1]*REAL(0.5);
	B[2] = side2[2]*REAL(0.5);

	// Rij is R1'*R2, i.e. the relative rotation between R1 and R2
	R11 = dCalcVectorDot3_44(R1+0,R2+0); R12 = dCalcVectorDot3_44(R1+0,R2+1); R13 = dCalcVectorDot3_44(R1+0,R2+2);
	R21 = dCalcVectorDot3_44(R1+1,R2+0); R22 = dCalcVectorDot3_44(R1+1,R2+1); R23 = dCalcVectorDot3_44(R1+1,R2+2);
	R31 = dCalcVectorDot3_44(R1+2,R2+0); R32 = dCalcVectorDot3_44(R1+2,R2+1); R33 = dCalcVectorDot3_44(R1+2,R2+2);

	Q11 = dFabs(R11); Q12 = dFabs(R12); Q13 = dFabs(R13);
	Q21 = dFabs(R21); Q22 = dFabs(R22); Q23 = dFabs(R23);
	Q31 = dFabs(R31); Q32 = dFabs(R32); Q33 = dFabs(R33);

	// for all 15 possible separating axes:
	//   * see if the axis separates the boxes. if so, return 0.
	//   * find the depth of the penetration along the separating axis (s2)
	//   * if this is the largest depth so far, record it.
	// the normal vector will be set to the separating axis with the smallest
	// depth. note: normalR is set to point to a column of R1 or R2 if that is
	// the smallest depth normal so far. otherwise normalR is 0 and normalC is
	// set to a vector relative to body 1. invert_normal is 1 if the sign of
	// the normal should be flipped.

	do {
#define TST(expr1,expr2,norm,cc) \
		expr1_val = (expr1); /* Avoid duplicate evaluation of expr1 */ \
		s2 = dFabs(expr1_val) - (expr2); \
		if (s2 > 0) return 0; \
		if (s2 > s) { \
			s = s2; \
			normalR = norm; \
			invert_normal = ((expr1_val) < 0); \
			code = (cc); \
			if (flags & CONTACTS_UNIMPORTANT) break; \
		}

		s = -dInfinity;
		invert_normal = 0;
		code = 0;

		// separating axis = u1,u2,u3
		TST (pp[0],(A[0] + B[0]*Q11 + B[1]*Q12 + B[2]*Q13),R1+0,1);
		TST (pp[1],(A[1] + B[0]*Q21 + B[1]*Q22 + B[2]*Q23),R1+1,2);
		TST (pp[2],(A[2] + B[0]*Q31 + B[1]*Q32 + B[2]*Q33),R1+2,3);

		// separating axis = v1,v2,v3
		TST (dCalcVectorDot3_41(R2+0,p),(A[0]*Q11 + A[1]*Q21 + A[2]*Q31 + B[0]),R2+0,4);
		TST (dCalcVectorDot3_41(R2+1,p),(A[0]*Q12 + A[1]*Q22 + A[2]*Q32 + B[1]),R2+1,5);
		TST (dCalcVectorDot3_41(R2+2,p),(A[0]*Q13 + A[1]*Q23 + A[2]*Q33 + B[2]),R2+2,6);

		// note: cross product axes need to be scaled when s is computed.
		// normal (n1,n2,n3) is relative to box 1.
#undef TST
#define TST(expr1,expr2,n1,n2,n3,cc) \
		expr1_val = (expr1); /* Avoid duplicate evaluation of expr1 */ \
		s2 = dFabs(expr1_val) - (expr2); \
		if (s2 > 0) return 0; \
		l = dSqrt ((n1)*(n1) + (n2)*(n2) + (n3)*(n3)); \
		if (l > 0) { \
			s2 /= l; \
			if (s2*fudge_factor > s) { \
				s = s2; \
				normalR = 0; \
				normalC[0] = (n1)/l; normalC[1] = (n2)/l; normalC[2] = (n3)/l; \
				invert_normal = ((expr1_val) < 0); \
				code = (cc); \
				if (flags & CONTACTS_UNIMPORTANT) break; \
			} \
		}

		// We only need to check 3 edges per box
		// since parallel edges are equivalent.

		// separating axis = u1 x (v1,v2,v3)
		TST(pp[2]*R21-pp[1]*R31,(A[1]*Q31+A[2]*Q21+B[1]*Q13+B[2]*Q12),0,-R31,R21,7);
		TST(pp[2]*R22-pp[1]*R32,(A[1]*Q32+A[2]*Q22+B[0]*Q13+B[2]*Q11),0,-R32,R22,8);
		TST(pp[2]*R23-pp[1]*R33,(A[1]*Q33+A[2]*Q23+B[0]*Q12+B[1]*Q11),0,-R33,R23,9);

		// separating axis = u2 x (v1,v2,v3)
		TST(pp[0]*R31-pp[2]*R11,(A[0]*Q31+A[2]*Q11+B[1]*Q23+B[2]*Q22),R31,0,-R11,10);
		TST(pp[0]*R32-pp[2]*R12,(A[0]*Q32+A[2]*Q12+B[0]*Q23+B[2]*Q21),R32,0,-R12,11);
		TST(pp[0]*R33-pp[2]*R13,(A[0]*Q33+A[2]*Q13+B[0]*Q22+B[1]*Q21),R33,0,-R13,12);

		// separating axis = u3 x (v1,v2,v3)
		TST(pp[1]*R11-pp[0]*R21,(A[0]*Q21+A[1]*Q11+B[1]*Q33+B[2]*Q32),-R21,R11,0,13);
		TST(pp[1]*R12-pp[0]*R22,(A[0]*Q22+A[1]*Q12+B[0]*Q33+B[2]*Q31),-R22,R12,0,14);
		TST(pp[1]*R13-pp[0]*R23,(A[0]*Q23+A[1]*Q13+B[0]*Q32+B[1]*Q31),-R23,R13,0,15);
#undef TST
	} while (0);

	if (!code) return 0;

	// if we get to this point, the boxes interpenetrate. compute the normal
	// in global coordinates.
	if (normalR) {
		normal[0] = normalR[0];
		normal[1] = normalR[4];
		normal[2] = normalR[8];
	}
	else {
		dMultiply0_331 (normal,R1,normalC);
	}
	if (invert_normal) {
		normal[0] = -normal[0];
		normal[1] = -normal[1];
		normal[2] = -normal[2];
	}
	*depth = -s;

	// compute contact point(s)

	if (code > 6) {
		// An edge from box 1 touches an edge from box 2.
		// find a point pa on the intersecting edge of box 1
		dVector3 pa;
		dReal sign;
		// Copy p1 into pa
		for (i=0; i<3; i++) pa[i] = p1[i]; // why no memcpy?
		// Get world position of p2 into pa
		for (j=0; j<3; j++) {
			sign = (dCalcVectorDot3_14(normal,R1+j) > 0) ? REAL(1.0) : REAL(-1.0);
			for (i=0; i<3; i++) pa[i] += sign * A[j] * R1[i*4+j];
		}

		// find a point pb on the intersecting edge of box 2
		dVector3 pb;
		// Copy p2 into pb
		for (i=0; i<3; i++) pb[i] = p2[i]; // why no memcpy?
		// Get world position of p2 into pb
		for (j=0; j<3; j++) {
			sign = (dCalcVectorDot3_14(normal,R2+j) > 0) ? REAL(-1.0) : REAL(1.0);
			for (i=0; i<3; i++) pb[i] += sign * B[j] * R2[i*4+j];
		}

		dReal alpha,beta;
		dVector3 ua,ub;
		// Get direction of first edge
		for (i=0; i<3; i++) ua[i] = R1[((code)-7)/3 + i*4];
		// Get direction of second edge
		for (i=0; i<3; i++) ub[i] = R2[((code)-7)%3 + i*4];
		// Get closest points between edges (one at each)
		dLineClosestApproach (pa,ua,pb,ub,&alpha,&beta);
		for (i=0; i<3; i++) pa[i] += ua[i]*alpha;
		for (i=0; i<3; i++) pb[i] += ub[i]*beta;
		// Set the contact point as halfway between the 2 closest points
		for (i=0; i<3; i++) contact[0].pos[i] = REAL(0.5)*(pa[i]+pb[i]);
		contact[0].depth = *depth;
		*return_code = code;
		return 1;
	}

	// okay, we have a face-something intersection (because the separating
	// axis is perpendicular to a face). define face 'a' to be the reference
	// face (i.e. the normal vector is perpendicular to this) and face 'b' to be
	// the incident face (the closest face of the other box).
	// Note: Unmodified parameter values are being used here
	const dReal *Ra,*Rb,*pa,*pb,*Sa,*Sb;
	if (code <= 3) { // One of the faces of box 1 is the reference face
		Ra = R1; // Rotation of 'a'
		Rb = R2; // Rotation of 'b'
		pa = p1; // Center (location) of 'a'
		pb = p2; // Center (location) of 'b'
		Sa = A;  // Side Lenght of 'a'
		Sb = B;  // Side Lenght of 'b'
	}
	else { // One of the faces of box 2 is the reference face
		Ra = R2; // Rotation of 'a'
		Rb = R1; // Rotation of 'b'
		pa = p2; // Center (location) of 'a'
		pb = p1; // Center (location) of 'b'
		Sa = B;  // Side Lenght of 'a'
		Sb = A;  // Side Lenght of 'b'
	}

	// nr = normal vector of reference face dotted with axes of incident box.
	// anr = absolute values of nr.
	/*
		The normal is flipped if necessary so it always points outward from box 'a',
		box 'b' is thus always the incident box
	*/
	dVector3 normal2,nr,anr;
	if (code <= 3) {
		normal2[0] = normal[0];
		normal2[1] = normal[1];
		normal2[2] = normal[2];
	}
	else {
		normal2[0] = -normal[0];
		normal2[1] = -normal[1];
		normal2[2] = -normal[2];
	}
	// Rotate normal2 in incident box opposite direction
	dMultiply1_331 (nr,Rb,normal2);
	anr[0] = dFabs (nr[0]);
	anr[1] = dFabs (nr[1]);
	anr[2] = dFabs (nr[2]);

	// find the largest compontent of anr: this corresponds to the normal
	// for the incident face. the other axis numbers of the incident face
	// are stored in a1,a2.
	int lanr,a1,a2;
	if (anr[1] > anr[0]) {
		if (anr[1] > anr[2]) {
			a1 = 0;
			lanr = 1;
			a2 = 2;
		}
		else {
			a1 = 0;
			a2 = 1;
			lanr = 2;
		}
	}
	else {
		if (anr[0] > anr[2]) {
			lanr = 0;
			a1 = 1;
			a2 = 2;
		}
		else {
			a1 = 0;
			a2 = 1;
			lanr = 2;
		}
	}

	// compute center point of incident face, in reference-face coordinates
	dVector3 center;
	if (nr[lanr] < 0) {
		for (i=0; i<3; i++) center[i] = pb[i] - pa[i] + Sb[lanr] * Rb[i*4+lanr];
	}
	else {
		for (i=0; i<3; i++) center[i] = pb[i] - pa[i] - Sb[lanr] * Rb[i*4+lanr];
	}

	// find the normal and non-normal axis numbers of the reference box
	int codeN,code1,code2;
	if (code <= 3) codeN = code-1; else codeN = code-4;
	if (codeN==0) {
		code1 = 1;
		code2 = 2;
	}
	else if (codeN==1) {
		code1 = 0;
		code2 = 2;
	}
	else {
		code1 = 0;
		code2 = 1;
	}

	// find the four corners of the incident face, in reference-face coordinates
	dReal quad[8];	// 2D coordinate of incident face (x,y pairs)
	dReal c1,c2,m11,m12,m21,m22;
	c1 = dCalcVectorDot3_14 (center,Ra+code1);
	c2 = dCalcVectorDot3_14 (center,Ra+code2);
	// optimize this? - we have already computed this data above, but it is not
	// stored in an easy-to-index format. for now it's quicker just to recompute
	// the four dot products.
	m11 = dCalcVectorDot3_44 (Ra+code1,Rb+a1);
	m12 = dCalcVectorDot3_44 (Ra+code1,Rb+a2);
	m21 = dCalcVectorDot3_44 (Ra+code2,Rb+a1);
	m22 = dCalcVectorDot3_44 (Ra+code2,Rb+a2);
	{
		dReal k1 = m11*Sb[a1];
		dReal k2 = m21*Sb[a1];
		dReal k3 = m12*Sb[a2];
		dReal k4 = m22*Sb[a2];
		quad[0] = c1 - k1 - k3;
		quad[1] = c2 - k2 - k4;
		quad[2] = c1 - k1 + k3;
		quad[3] = c2 - k2 + k4;
		quad[4] = c1 + k1 + k3;
		quad[5] = c2 + k2 + k4;
		quad[6] = c1 + k1 - k3;
		quad[7] = c2 + k2 - k4;
	}

	// find the size of the reference face
	dReal rect[2];
	rect[0] = Sa[code1];
	rect[1] = Sa[code2];

	// intersect the incident and reference faces
	dReal ret[16];
	int n = intersectRectQuad (rect,quad,ret);
	if (n < 1) return 0;		// this should never happen

	// convert the intersection points into reference-face coordinates,
	// and compute the contact position and depth for each point. only keep
	// those points that have a positive (penetrating) depth. delete points in
	// the 'ret' array as necessary so that 'point' and 'ret' correspond.
	dReal point[3*8];		// penetrating contact points
	dReal dep[8];			// depths for those points
	dReal det1 = dRecip(m11*m22 - m12*m21);
	m11 *= det1;
	m12 *= det1;
	m21 *= det1;
	m22 *= det1;
	int cnum = 0;			// number of penetrating contact points found
	for (j=0; j < n; j++) {
		dReal k1 =  m22*(ret[j*2]-c1) - m12*(ret[j*2+1]-c2);
		dReal k2 = -m21*(ret[j*2]-c1) + m11*(ret[j*2+1]-c2);
		for (i=0; i<3; i++) point[cnum*3+i] =
							center[i] + k1*Rb[i*4+a1] + k2*Rb[i*4+a2];
		dep[cnum] = Sa[codeN] - dCalcVectorDot3(normal2,point+cnum*3);
		if (dep[cnum] >= 0) {
			ret[cnum*2] = ret[j*2];
			ret[cnum*2+1] = ret[j*2+1];
			cnum++;
			if ((cnum | CONTACTS_UNIMPORTANT) == (flags & (NUMC_MASK | CONTACTS_UNIMPORTANT))) {
					break;
			}
		}
	}
	if (cnum < 1) {
			return 0;	// this should not happen, yet does at times (demo_plane2d single precision).
	}

	// we can't generate more contacts than we actually have
	int maxc = flags & NUMC_MASK;
	if (maxc > cnum) maxc = cnum;
	if (maxc < 1) maxc = 1;	// Even though max count must not be zero this check is kept for backward compatibility as this is a public function

	if (cnum <= maxc) {
		// we have less contacts than we need, so we use them all
		for (j=0; j < cnum; j++) {
			dContactGeom *con = &contact[j];
			for (i=0; i<3; i++) con->pos[i] = point[j*3+i] + pa[i];
			con->depth = dep[j];
		}
	}
	else {
		// we have more contacts than are wanted, some of them must be culled.
		// find the deepest point, it is always the first contact.
		int i1 = 0;
		dReal maxdepth = dep[0];
		for (i=1; i<cnum; i++) {
			if (dep[i] > maxdepth) {
		maxdepth = dep[i];
		i1 = i;
			}
		}

		int iret[8];
		cullPoints (cnum,ret,maxc,i1,iret);

		for (j=0; j < maxc; j++) {
			dContactGeom *con = &contact[j];
			for (i=0; i<3; i++) con->pos[i] = point[iret[j]*3+i] + pa[i];
			con->depth = dep[iret[j]];
		}
		cnum = maxc;
	}

	*return_code = code;
	return cnum;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace Eigen;

Contacts odeBoxBox(
	const Matrix4d &Ma, const Vector3d &scaleA,
	const Matrix4d &Mb, const Vector3d &scaleB)
{
	// dBoxBox expects vec4 and mat43 (transposed)
	double pa4[4], pb4[4], Ra4[12], Rb4[12];
	Ra4[0] = Ma(0,0); Ra4[1] = Ma(0,1); Ra4[2]  = Ma(0,2);
	Ra4[4] = Ma(1,0); Ra4[5] = Ma(1,1); Ra4[6]  = Ma(1,2);
	Ra4[8] = Ma(2,0); Ra4[9] = Ma(2,1); Ra4[10] = Ma(2,2);
	Ra4[3] = 0.0;     Ra4[7] = 0.0;     Ra4[11] = 0.0;
	pa4[0] = Ma(0,3); pa4[1] = Ma(1,3); pa4[2]  = Ma(2,3);
	pa4[3] = 0.0;
	Rb4[0] = Mb(0,0); Rb4[1] = Mb(0,1); Rb4[2]  = Mb(0,2);
	Rb4[4] = Mb(1,0); Rb4[5] = Mb(1,1); Rb4[6]  = Mb(1,2);
	Rb4[8] = Mb(2,0); Rb4[9] = Mb(2,1); Rb4[10] = Mb(2,2);
	Rb4[3] = 0.0;     Rb4[7] = 0.0;     Rb4[11] = 0.0;
	pb4[0] = Mb(0,3); pb4[1] = Mb(1,3); pb4[2]  = Mb(2,3);
	pb4[3] = 0.0;

	const int nContactGeoms = 8;
	dContactGeom contactGeoms[nContactGeoms];
	memset(contactGeoms, 0, nContactGeoms*sizeof(dContactGeom));
	double sidesA[3], sidesB[3];
	sidesA[0] = scaleA(0);
	sidesA[1] = scaleA(1);
	sidesA[2] = scaleA(2);
	sidesB[0] = scaleB(0);
	sidesB[1] = scaleB(1);
	sidesB[2] = scaleB(2);
	int rc = 0;
	double depthMax = 0;
	double normal[3] = {0};
	int hits = dBoxBox(pa4, Ra4, sidesA, pb4, Rb4, sidesB,
						normal, &depthMax, &rc, nContactGeoms, contactGeoms);

	// Structure to be returned
	Contacts contacts;
	contacts.count = hits;
	contacts.depthMax = depthMax;
	contacts.normal << normal[0], normal[1], normal[2];

	for(int k = 0; k < hits; ++k) {
		// Get contact point and normal in world space
		dContactGeom *contactGeom = &contactGeoms[k];
		contacts.positions[k] << contactGeom->pos[0], contactGeom->pos[1], contactGeom->pos[2];
		contacts.depths[k] = contactGeom->depth;
	}
	return contacts;
}

}