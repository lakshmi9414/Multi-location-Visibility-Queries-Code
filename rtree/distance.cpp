/* this file contains the implementation of functions defining
Tinf and Tmininf */

#include "distance.h"
#include "../func/gendef.h"
#include <memory.h>
#include <math.h>
#include<algorithm>
#include<cfloat>
#include <limits>



using namespace std;
//-----------------------------------------------------------
bool chk_in(float _a, float _b1, float _b2)
{
	if (_b1 <= _a && _a <= _b2)
		return true;
	else 
		return false;
}
//-----------------------------------------------------------
bool chk_intr(float _a1, float _a2, float _b1, float _b2)
{
	if (fmax(_a1, _b1) <= fmin(_a2, _b2))
		return true;
	else 
		return false;
}
//-----------------------------------------------------------
bool chk_mbr_contains(float * _m1, float * _m2, int _dim)
{
	for (int i = 0; i < _dim; i ++)
	{
		float m1xl = _m1[2 * i], m1xr = _m1[2 * i + 1];
		float m2xl = _m2[2 * i], m2xr = _m2[2 * i + 1];
		if (!(m1xl <= m2xl && m2xr <= m1xr))
			return false;
	}
	return true;
}
//-----------------------------------------------------------
bool chk_part_intr(float _a1, float _a2, float _b1, float _b2)
{
	if ((_a1 < _b1 && _b1 <= _a2 && _a2 < _b2) || 
		(_b1 < _a1 && _a1 <= _b2 && _b2 < _a2))
		return true;
	else 
		return false;
}
//-----------------------------------------------------------
void intr_period(float *_m1, float *_m2, float *_v1, float *_v2,
				 int _dim, float &_ts, float &_te)
//this function computes the intersection period in the future
//for two rectangles to intersect
//_m1 and _m2 are mbrs of the current time
//_v1 and _v2 are velocities of mbrs
//the returned period is stored in [_ts, _te]
{
	float *tsi = new float [_dim];
	float *tei = new float [_dim];
	  //tsi and tei is the intersect period along each dimension

	for (int i = 0; i < _dim; i ++)
	{
		float m1xl = _m1[2 * i], m1xr = _m1[2 * i + 1];
		float m2xl = _m2[2 * i], m2xr = _m2[2 * i + 1];
		float m1vl = _v1[2 * i], m1vr = _v1[2 * i + 1];
		float m2vl = _v2[2 * i], m2vr = _v2[2 * i + 1];
		float tlr = -MAXREAL, trl = -MAXREAL;
		  //tlr and trl are the intersection time of left-right
		  //and right-left end points respectively
		if (m2vr != m1vl)
			tlr = (m1xl - m2xr) / (m2vr - m1vl);
		if (m2vl != m1vr)
			trl = (m1xr - m2xl) / (m2vl - m1vr);
		if (tlr >= 0 && trl >= 0)
		{
			tsi[i] = fmin(tlr, trl); tei[i] = fmax(tlr, trl);
		}
		else
			if (fmax(tlr, trl) >= 0)
			{
				if (chk_intr(m1xl, m1xr, m2xl, m2xr))
				{
					tsi[i] = 0; tei[i] = fmax(tlr, trl);
				}
				else
				{
					tsi[i] = fmax(tlr, trl); tei[i] = MAXREAL;
				}
			}
			else
			{
				if (chk_intr(m1xl, m1xr, m2xl, m2xr))
				{
					tsi[i] = 0; tei[i] = MAXREAL;
				}
				else
				{
					tsi[i] = MAXREAL; tei[i] = MAXREAL;
				}
			}
	}

	_ts = -MAXREAL; _te = MAXREAL;
	for (int i = 0; i < _dim; i ++)
	{
		_ts = fmax(_ts, tsi[i]); _te = fmin(_te, tei[i]);
	}
	if (_ts > _te)
	{
		_ts = _te = MAXREAL;
	}

	delete [] tei;
	delete [] tsi;
}
//-----------------------------------------------------------
float NNinf(float *_m1, float *_m2, float *_m3, int _dim)
// _m3 is the query point
// _m1 is the NN and _m2 is the object
{
	float *m1 = (float *) _m1, *m2 = (float *) _m2, *m3 = (float *) _m3;
	float *V1 = m1 + 2 * _dim, *V2 = m2 + 2 * _dim, *V3 = m3 + 2 * _dim;
	//the function handles 2D points 
	float x0 = m1[0], y0 = m1[2], u0 = V1[0], v0 = V1[2];
	float x1 = m2[0], y1 = m2[2], u1 = V2[0], v1 = V2[2];
	float x2 = m3[0], y2 = m3[2], u2 = V3[0], v2 = V3[2];

	if ((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2) >=
		(x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))
		return 0;


	float A = x0 - x2, B = u0 - u2, C = y0 - y2, D = v0 - v2;
	float H = x1 - x2, I = u1 - u2, J = y1 - y2, K = v1 - v2;
	float E = B * B + D * D, F = 2 * A * B + 2 * C * D, G = A * A + C * C;
	float L = I * I + K * K, M = 2 * H * I + 2 * J * K, N = H * H + J * J;
	float t1, t2;
	bool ret = solve_quadratic(E - L, F - M, G - N, t1, t2);
	if (!ret)
		// no solution
	{
		if (E - L > 0)
		{
			return 0;
		}
		else
			return MAXREAL;
	}
	if (E - L > 0)
	{
		if (t1 >= 0)
			return 0;
		else
			return fmax(0.0, t2);
	}
	else
	{
		if (t2 >= 0)
			return fmax(0.0, t1);
		else
			return MAXREAL;
	}
}
//-----------------------------------------------------------
/*****************************************************************
this function returns the minimum influence time for a TPNN query
NOTE: THIS FUNCTION HANDLES ONLY 2D STATIC POINTS 
para:
p1: the current NN found 
p2: the object that is being tested (for which we want to return its
  minimum inflence time
p3: the query line direction (4-size array: x0, y0, u0, v0)
NOTE: m1 and m2 and m3 are all points 
dim: the dimensionality

coded by Qiongmao Shen 08/02/02
*****************************************************************/
float NNinf2(float *p1,float *p2,float *p3)
//p1 p2 p3 are the points(p1 is the NN,p2 is point which are examed,p3 is the partition point
{
   float pLength;
   float numerator;
   float denominator;
   //denominator=2*p2[0]*cos(Xangle)+2*p2[1]*sin(Xangle)-2*p1[0]*cos(Xangle)-2*p1[1]*sin(Xangle);
   denominator=2*p2[0]*p3[2]+2*p2[1]*p3[3]-2*p1[0]*p3[2]-2*p1[1]*p3[3];
   if (denominator==0)
   {
	   return MAXREAL;
   }
   else
   {
	  //numerator=2*p1[0]*start[0]-p1[0]*p1[0]-p1[1]*p1[1]+2*p1[1]*start[1]-2*p2[0]*start[0]+
	  //    p2[0]*p2[0]+p2[1]*p2[1]-2*p2[1]*start[1];
      numerator=2*p1[0]*p3[0]-p1[0]*p1[0]-p1[1]*p1[1]+2*p1[1]*p3[1]-2*p2[0]*p3[0]+
		  p2[0]*p2[0]+p2[1]*p2[1]-2*p2[1]*p3[1];
	  pLength=numerator/denominator;
   };
   return pLength;
};

//-----------------------------------------------------------
float NNmininf(float *_m1, float *_m2, float *_m3, int _dim)
// _m3 is the query point
// _m1 is the NN and _m2 is the MBR
{
	float *m1 = (float *) _m1, *m3 = (float *) _m3;
	float *v1 = m1 + 2 * _dim, *v3 = m3 + 2 * _dim;
	float *m2 = new float [4 * _dim], *v2 = m2 + 2 * _dim;
	memcpy(m2, _m2, 2 * _dim * sizeof(float));
	memcpy(v2, _m2 + 2 * _dim, 2 * _dim * sizeof(float));

	//check if _m3 is in _m2
	bool in = true;
	for (int i = 0; i < _dim; i ++)
	{
		if (!(m2[2 * i] <= m3[2 * i] && m3[2 * i] <= m2[2 * i + 1]))
		{
			in = false; i = _dim;
		}
	}
	if (in) 
		return 0;

	int cal_dim = -1;
	  //this is the dimension you want to calculate distance
	float max = -MAXREAL;
	for (int i = 0; i < _dim; i ++)
	{	
		float this_dist=-MAXREAL;
		if (m3[2 * i] < m2[2 * i])
			this_dist = m2[2 * i] - m3[2 * i];
		else
			this_dist = m3[2 * i] - m2[2 * i + 1];
		if (this_dist > max)
		{
			max = this_dist; cal_dim = i;
		}
	}

	// now construct a new m2 to compute Tmininf
	for (int i = 0; i < _dim; i ++)
	{
		if (i == cal_dim)
		{
			if (m3[2 * i] < m2[2 * i])
			{  m2[2 * i + 1] = m2[2 * i]; v2[2 * i + 1] = v2[2 * i];}
			else
				{  m2[2 * i] = m2[2 * i + 1]; v2[2 * i] = v2[2 * i + 1];}

		}
		else
		{
			m2[2 * i] = m2[2 * i + 1] = m3[2 * i];
			v2[2 * i] = v2[2 * i + 1] = v3[2 * i];
		}
	}

	float ret = NNinf(m1, m2, m3, _dim);

	delete [] m2;
	return ret;
}
//-----------------------------------------------------------
float partintr_time(float *_m1, float *_m2, float *_v1, float *_v2, int _dim, 
					float _ts, float _te)
//this function computes the time (return value) of partial intersection 
//between two rectangles
//_m1 and _m2 are mbrs of the current time
//_v1 and _v2 are velocities of mbrs
//[_ts, _te] is the intersection period
{
	float tpar = MAXREAL;
	float *tpari = new float [_dim];
	  //tpari is the ealiest time extents along the ith dimension partially intersect

	for (int i = 0; i < _dim; i ++)
	{
		float m1xl = _m1[2 * i], m1xr = _m1[2 * i + 1];
		float m2xl = _m2[2 * i], m2xr = _m2[2 * i + 1];
		float m1vl = _v1[2 * i], m1vr = _v1[2 * i + 1];
		float m2vl = _v2[2 * i], m2vr = _v2[2 * i + 1];
		float tll = -MAXREAL, trr = -MAXREAL;
		  //tll and trr are the intersection time of left-left
		  //and right-right end points respectively
		if (m2vl != m1vl)
			tll = (m1xl - m2xl) / (m2vl - m1vl);
		if (m2vr != m1vr)
			trr = (m1xr - m2xr) / (m2vr - m1vr);

		float m1xlts = m1xl + m1vl * _ts, m1xrts = m1xr + m1vr * _ts;
		float m2xlts = m2xl + m2vl * _ts, m2xrts = m2xr + m2vr * _ts;

		if (!chk_in(tll, _ts, _te) && !chk_in(trr, _ts, _te))
		{
			if (chk_part_intr(m1xlts, m1xrts, m2xlts, m2xrts))
				tpari[i] = _ts;
			else
				tpari[i] = MAXREAL;
		}
		else
			if (chk_in(tll, _ts, _te) && chk_in(trr, _ts, _te))
			{
				if (chk_part_intr(m1xlts, m1xrts, m2xlts, m2xrts))
					tpari[i] = _ts;
				else
					tpari[i] = fmin(tll, trr);
			}
			else
			{
				if (chk_part_intr(m1xlts, m1xrts, m2xlts, m2xrts))
					tpari[i] = _ts;
				else
				{
					if (chk_in(tll, _ts, _te))
						tpari[i] = tll;
					else 
						tpari[i] = trr;
				}
			}
	}

	tpar = MAXREAL; 
	for (int i = 0; i < _dim; i ++)
	{
		tpar = fmin(tpar, tpari[i]);
	}

	delete [] tpari;
	return tpar;
}
//-----------------------------------------------------------
bool solve_quadratic(float _a, float _b, float _c, 
					  float &_x1, float &_x2)
{
	if (_a == 0 && _b == 0)
	{
		return false;
	}

	if (_a == 0)
	{
		_x1 = _x2 = -_c / _b;
		return true;
	}

	float delta_sqr = _b * _b - 4 * _a * _c;
	if (delta_sqr < 0)
		return false;
	float delta = sqrt(delta_sqr);
	_x1 = (-_b - delta) / 2 / _a;
	_x2 = (-_b + delta) / 2 / _a;
	return true;
}
//-----------------------------------------------------------
float testdist(float *_m1, float *_m2, float *_m3, int _dim)
//query is always the second place
{
	float *m1, *m2;
	m1 = (float *)_m1; m2 = (float *)_m2;
	float *p = new float [_dim];
	for (int i = 0; i < _dim; i ++)
		p[i] = m2[2 * i];
	float min = MINDIST(p, m1, _dim);
	delete [] p;
	return min;
}
//-----------------------------------------------------------
float testmin(float *_m1, float *_m2, float *_m3, int _dim)
//query is always the second place
{
	float *m1, *m2;
	m1 = (float *)_m1; m2 = (float *)_m2;
	float *p = new float [_dim];
	for (int i = 0; i < _dim; i ++)
		p[i] = m2[2 * i];
	float min = MINDIST(p, m1, _dim);
	delete [] p;
	return min;
}
//-----------------------------------------------------------
float WQinf(float *_m1, float *_m2, float *_m3, int _dim)
{
	float *m1 = (float *) _m1, *m2 = (float *) _m2;
	float *v1 = m1 + 2 * _dim, *v2 = m2 + 2 * _dim;
	float ts, te;
	intr_period(m1, m2, v1, v2, _dim, ts, te);
	if (ts == 0) 
		return te;
	else 
		return ts;
}
//-----------------------------------------------------------
float WQmininf(float *_m1, float *_m2, float *_m3, int _dim)
{
	float *m1 = (float *) _m1, *m2 = (float *) _m2;
	float *v1 = m1 + 2 * _dim, *v2 = m2 + 2 * _dim;
	float ts, te;
	intr_period(m1, m2, v1, v2, _dim, ts, te);
	
	if (!chk_mbr_contains(m2, m1, _dim))
		//query should always be put on _m2
		return ts;
	else
		return partintr_time(m1, m2, v1, v2, _dim, ts, te);
}
