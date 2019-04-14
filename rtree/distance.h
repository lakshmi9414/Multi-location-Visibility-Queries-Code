/* this file contains declarations of functions calculating
Tinf and Tmininf */
#ifndef DISTANCE_H
#define DISTANCE_H

bool chk_in(float _a, float _b1, float _b2);
bool chk_intr(float _a1, float _a2, float _b1, float _b2);
bool chk_mbr_contains(float * _m1, float * _m2, int _dim);
bool chk_part_intr(float _a1, float _a2, float _b1, float _b2);
void intr_period(float *_m1, float *_m2, float *_v1, float *_v2,
				 int _dim, float &_ts, float &_te);
float partintr_time(float *_m1, float *_m2, float *_v1, float *_v2, int _dim);
bool solve_quadratic(float _a, float _b, float _c, float &_x1, float &_x2);
//-----------------------------------------------------------
float NNinf(float *_m1, float *_m2, float *_m3, int _dim);
float NNinf2(float *p1,float *p2,float *p3);
float NNmininf(float *_m1, float *_m2, float *_m3, int _dim);
float testdist(float *_m1, float *_m2, float *_m3, int _dim);
float testmin(float *_m1, float *_m2, float *_m3, int _dim);
float WQinf(float *_m1, float *_m2, float *_m3, int _dim);
float WQmininf(float *_m1, float *_m2, float *_m3, int _dim);
//-----------------------------------------------------------
#endif