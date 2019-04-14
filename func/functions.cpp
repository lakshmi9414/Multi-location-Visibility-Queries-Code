#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gendef.h"
//////////////////////////////////////////////////////////////////////////////
// globals
//////////////////////////////////////////////////////////////////////////////

#ifdef UNIX
void strupr(char *_msg)
{
	int dist = 'A' - 'a';
	char *c_ptr = _msg;

	while (*c_ptr)
	{
		if (*c_ptr >= 'a')
			*c_ptr += dist;
		c_ptr ++;
	}
}
#endif

void error(char *t, bool ex)
{
    fprintf(stderr, t);
    if (ex)
        exit(0);
}

float area(int dimension, float *mbr)
// berechnet Flaeche (Volumen) eines mbr der Dimension dimension
{
    int i;
    float sum;

    sum = 1.0;
    for (i = 0; i < dimension; i++)
	sum *= mbr[2*i+1] - mbr[2*i];

    return sum;
}

float margin(int dimension, float *mbr)
// berechnet Umfang eines mbr der Dimension dimension
{
    float *ml, *mu, *m_last, sum;

    sum = 0.0;
    m_last = mbr + 2*dimension;
    ml = mbr;
    mu = ml + 1;
    while (mu < m_last)
    {
	sum += *mu - *ml;
	ml += 2;
	mu += 2;
    }

    return sum;
}

bool inside(float &p, float &lb, float &ub)
// ist ein Skalar in einem Intervall ?
{
    return (p >= lb && p <= ub);
}

bool inside(float *v, float *mbr, int dimension)
// ist ein Vektor in einer Box ?
{
    int i;

    for (i = 0; i < dimension; i++)
	if (!inside(v[i], mbr[2*i], mbr[2*i+1]))
	    return FALSE;

    return TRUE;
}

// calcutales the overlapping rectangle between r1 and r2
// if rects do not overlap returns null
float* overlapRect(int dimension, float *r1, float *r2)
{
        float *overlap = new float[2*dimension];
        for (int i=0; i<dimension; i++)
        {
            if ((r1[i*2]>r2[i*2+1]) || (r1[i*2+1]<r2[i*2])) // non overlapping
	    {
                delete [] overlap;
		return NULL;
	    }
	    overlap[2*i] = max(r1[i*2], r2[i*2]);
            overlap[2*i+1] = min(r1[i*2+1], r2[i*2+1]);
        }

        return overlap;
}

float overlap(int dimension, float *r1, float *r2)
// calcutales the overlapping area of r1 and r2
// calculate overlap in every dimension and multiplicate the values
{
    float sum;
    float *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;

    sum = 1.0;
    r1pos = r1; r2pos = r2;
    r1last = r1 + 2 * dimension;

    while (r1pos < r1last)
    {
	r1_lb = *(r1pos++);
	r1_ub = *(r1pos++);
	r2_lb = *(r2pos++);
	r2_ub = *(r2pos++);

        // calculate overlap in this dimension

        if (inside(r1_ub, r2_lb, r2_ub))
        // upper bound of r1 is inside r2 
	{
            if (inside(r1_lb, r2_lb, r2_ub))
            // and lower bound of r1 is inside
                sum *= (r1_ub - r1_lb);
            else
                sum *= (r1_ub - r2_lb);
	}
	else
	{
            if (inside(r1_lb, r2_lb, r2_ub))
	    // and lower bound of r1 is inside
		sum *= (r2_ub - r1_lb);
	    else 
	    {
		if (inside(r2_lb, r1_lb, r1_ub) &&
		    inside(r2_ub, r1_lb, r1_ub))
	        // r1 contains r2
		    sum *= (r2_ub - r2_lb);
		else
		// r1 and r2 do not overlap
		    sum = 0.0;
	    }
	}
    }

    return sum;
}

void enlarge(int dimension, float **mbr, float *r1, float *r2)
// enlarge r in a way that it contains s
{
    int i;

    *mbr = new float[2*dimension];
    for (i = 0; i < 2*dimension; i += 2)
    {
	(*mbr)[i]   = min(r1[i],   r2[i]);

	(*mbr)[i+1] = max(r1[i+1], r2[i+1]);
    }
}

bool section(int dimension, float *mbr1, float *mbr2)
{
    int i;

    for (i = 0; i < dimension; i++)
    {
	if (mbr1[2*i] > mbr2[2*i + 1] ||
	    mbr1[2*i + 1] < mbr2[2*i]) 
	    return FALSE;
    }
    return TRUE;
}

bool section_c(int dimension, float *mbr1, float *center, float radius)
{
	float r2;

	r2 = radius * radius;

	if ( (r2 - MINDIST(center,mbr1,dimension)) < 1.0e-8)
		return TRUE;
	else
		return FALSE;
	
}

int sort_lower_mbr(const void *d1, const void *d2)
// Vergleichsfunktion fuer qsort, sortiert nach unterem Wert der mbr bzgl.
// der angegebenen Dimension
{
    SortMbr *s1, *s2;
    float erg;
    int dimension;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;
    erg = s1->mbr[2*dimension] - s2->mbr[2*dimension];
    if (erg < 0.0)
	return -1;
    else if (erg == 0.0)
	return 0;
    else 
	return 1;
}

int sort_upper_mbr(const void *d1, const void *d2)
// Vergleichsfunktion fuer qsort, sortiert nach oberem Wert der mbr bzgl.
// der angegebenen Dimension
{
    SortMbr *s1, *s2;
    float erg;
    int dimension;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;
    erg = s1->mbr[2*dimension+1] - s2->mbr[2*dimension+1];
    if (erg < 0.0)
	return -1;
    else if (erg == 0.0)
	return 0;
    else 
	return 1;
}

int sort_center_mbr(const void *d1, const void *d2)
// Vergleichsfunktion fuer qsort, sortiert nach center of mbr 
{
    SortMbr *s1, *s2;
    int i, dimension;
    float d, e1, e2;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;

    e1 = e2 = 0.0;
    for (i = 0; i < dimension; i++)
    {
        d = ((s1->mbr[2*i] + s1->mbr[2*i+1]) / 2.0) - s1->center[i];
        e1 += d*d;
        d = ((s2->mbr[2*i] + s2->mbr[2*i+1]) / 2.0) - s2->center[i];
        e2 += d*d;
    }

    if (e1 < e2)
	return -1;
    else if (e1 == e2)
	return 0;
    else 
	return 1;
}

int sortmindist(const void *element1, const void *element2)
{
    //
    // Vergleichsfunktion fuer die Sortierung der BranchList bei der
    // NearestNarborQuery (Sort, Branch and Prune)
    //

	//The consequence of this sorting is that active branches are sorted
	//in acsending order of their mindist.

    BranchList *e1,*e2;

    e1 = (BranchList *) element1;
    e2 = (BranchList *) element2;
    
    if (e1->mindist < e2->mindist) 
	return(-1);
    else if (e1->mindist > e2->mindist)
	return(1);
    else
	return(0);

}

int pruneBrunchList(float *nearest_distanz, const void *activebrunchList, 
		    int n)
{

    // Schneidet im Array BranchList alle Eintraege ab, deren Distanz groesser
    // ist als die aktuell naeheste Distanz
    //

    BranchList *bl; 
    
    int i,j,k, aktlast;
    
    bl = (BranchList *) activebrunchList;
    
    // 1. Strategie:
    // 
    // Ist MINDIST(P,M1) > MINMAXDIST(P,M2), so kann der 
    // NearestNeighbor auf keinen Fall mehr in M1 liegen!
    //


    aktlast = n;
    

    for( i = 0; i < aktlast ; i++ )
    {
	if (bl[i].minmaxdist < bl[aktlast-1].mindist)
		//bl[aktlast-1] is the maximum mindist of all the branches 
		//if this condition holds, it implies that there exist some branches whose
		//mindists are greater than minmandist of this branch. in the sequel, all these
		//branches will be removed from further consideration. 
	    for( j = 0; (j < aktlast) ; j++ )
		if ((i!=j) && (bl[j].mindist>bl[i].minmaxdist))
		{
		    aktlast = j;

			//checkpoint
			if (j <= i)
				error("Error in pruneBrunchList\n", true);
		    break;
		}
    }

    

    // 2. Strategie:
    //
    // nearest_distanz > MINMAXDIST(P,M)
    // -> nearest_distanz = MIMMAXDIST(P,M)
    //


    for( i = 0; i < aktlast; i++)
	if (*nearest_distanz > bl[i].minmaxdist)
	    *nearest_distanz = bl[i].minmaxdist;

    
    // 3. Strategie:
    // 
    // nearest_distanz < MINDIST(P,M)
    //
    // in M kann der Nearest-Narbor sicher nicht mehr liegen.
    //

    for( i = 0; (i < aktlast) && (*nearest_distanz >= bl[i].mindist) ; i++);
    
    aktlast = i;

    // printf("n: %d aktlast: %d \n",n,aktlast);   

    return (aktlast);
   
}

float objectDIST(float *p1, float *p2, int dim)
{

    //
    // Berechnet das Quadrat der euklid'schen Metrik.
    // (Der tatsaechliche Abstand ist uninteressant, weil
    // die anderen Metriken (MINDIST und MINMAXDIST fuer
    // die NearestNarborQuery nur relativ nie absolut 
    // gebraucht werden; nur Aussagen "<" oder ">" sind
    // relevant.
    //

    float summe = 0;
    int i;
    
    for( i = 0; i < dim; i++)
	summe += pow(p1[i] - p2[i], 2 );

    //return( sqrt(summe) );
    return(summe);
}

float MINDIST(float *p, float *bounces, int dim)
{
    //
    // Berechne die kuerzeste Entfernung zwischen einem Punkt Point
    // und einem MBR bounces (Lotrecht!)
    //

    float summe = 0.0;
    float r;
    int i;

    for(i = 0; i < dim; i++)
    {
		if (p[i] < bounces[2*i])
			r = bounces[2*i];
		else
		{
			if (p[i] > bounces[2*i+1])
			r = bounces[2*i+1];
			else 
			r = p[i];
		}    

		summe += pow(p[i]-r,2);
    }
    return(summe);
    
}

float MINMAXDIST(float *p, float *bounces)
{

    // Berechne den kleinsten maximalen Abstand von einem Punkt Point
    // zu einem MBR bounces.
    // Wird benutzt zur Abschaetzung von Abstaenden bei NearestNarborQuery.
    // Kann als Supremum fuer die aktuell kuerzeste Distanz: 
    // Alle Punkte mit einem Abstand > MINMAXDIST sind keine Kandidaten mehr
    // fuer den NearestNarbor
    // vgl. Literatur: 
    // Nearest Narbor Query v. Roussopoulos, Kelley und Vincent, 
    // University of Maryland
    
    float summe = 0;
    float minimum = 1.0e20;
    float S = 0;

    float rmk, rMi;
    int k,i;

    for( i = 0; i < DIMENSION; i++) 
    { 
	rMi = (	p[i] >= (bounces[2*i]+bounces[2*i+1])/2 )
	    ? bounces[2*i] : bounces[2*i+1];
	S += pow( p[i] - rMi , 2 );
    }

    for( k = 0; k < DIMENSION; k++)
    {  
	
	rmk = ( p[k] <=  (bounces[2*k]+bounces[2*k+1]) / 2 ) ? 
	    bounces[2*k] : bounces[2*k+1];

	summe = pow( p[k] - rmk , 2 );	
	
	rMi = (	p[k] >= (bounces[2*k]+bounces[2*k+1]) / 2 )
	    ? bounces[2*k] : bounces[2*k+1];

	summe += S - pow( p[k] - rMi , 2 );
	
	minimum = min( minimum,summe);
    }

    return(minimum);
}


//---added for validity region---------------------------------

// decide the distance between p0 and p1
float point_point_dist(float *p0, float *p1)
{
	return sqrt((p1[0]-p0[0])*(p1[0]-p0[0])+(p1[1]-p0[1])*(p1[1]-p0[1]));
}

// compute the distance from the point p to the line l0l1
float point_line_dist(float *l0, float *l1, float *p)
{
	return fabs((l1[1]-l0[1])*(p[0]-l0[0])-(l1[0]-l0[0])*(p[1]-l0[1]))
		/ point_point_dist(l0, l1);
}

// compute which direction does the point p lies in with repect to line l0l1
int clockwise_point(float *l0, float *l1, float *p)
{
	int rtn = 0;

	if (l0[0] == l1[0])
	{
		if (p[0] > l0[0])
			rtn = 1;
		else if (p[0] < l0[0])
			rtn = -1;

		if (l0[1] > l1[1])
			rtn = -rtn;
	}
	else if (l0[1] == l1[1])
	{
		if (p[1] > l0[1])
			rtn = 1;
		else if (p[1] < l0[1])
			rtn = -1;

		if (l0[0] > l1[0])
			rtn = -rtn;
	}
	else
	{
		float val = (l1[1]-l0[1])/(l1[0]-l0[0])*(p[0]-l0[0])+l0[1]-p[1];
		if (val > 0.0)
			rtn = 1;
		else if (val < 0.0)
			rtn = -1;

		if (l0[0] < l1[0])
			rtn = -rtn;
	}

	return rtn;
}

float distinit(float *p, float *e)
{
	if (e[0] == p[0] && e[2] == p[1]) // the same point, define the distance as MAXREAL
		return MAXREAL;
	else 
		return MINDIST(p, e, 2);
}

// distance function for the consecutive influence points
float distcont(float *l0, float *l1, float *e)
{
//test
//if (e[0] - 5016 < 0.1 && e[2] == 4920)
//{
//	int j;
//	j =1;
//}
//if (e[0] == 5242 && e[2] == 5384 && l1[0] == 5177.999512 && l1[1] == 5366)
//{
//	int j;
//	j =1;
//}
	
	if (e[0] == e[1] && e[2] == e[3]) // is it a point?
	{
		float p[2];
		p[0] = e[0];
		p[1] = e[2];

		if (clockwise_point(l0, l1, p) <= 0)
			return MAXREAL;
		float sqrl0e = (l0[0] - p[0]) * (l0[0] - p[0]) + (l0[1] - p[1]) * (l0[1] - p[1]);
		float sqrl1e = (l1[0] - p[0]) * (l1[0] - p[0]) + (l1[1] - p[1]) * (l1[1] - p[1]);
		float sqrl0l1 = (l1[0] - l0[0]) * (l1[0] - l0[0]) + (l1[1] - l0[1]) * (l1[1] - l0[1]);
		return (sqrl0e + sqrl1e - sqrl0l1) / (2 * sqrt(sqrl0e) * sqrt(sqrl1e));
	}
	else // it is a rectangle, approximate it as a circle
	{
		float o[2]; // center of the circle
		o[0] = (e[0] + e[1]) / 2;
		o[1] = (e[2] + e[3]) / 2;
		
		float r = sqrt((o[0] - e[0]) * (o[0] - e[0]) + (o[1] - e[2]) * (o[1] - e[2])); // radius of the circle

		if (point_line_dist(l0, l1, o) >= r && clockwise_point(l0, l1, o) == -1) // the circle is not in the clockwise direction
			return MAXREAL;
		
		if (point_point_dist(l0, o) <= r || point_point_dist(l1, o) <= r) // the line segment intersects the circle
			return -MAXREAL;

		float c[2]; // center of the line segment l0l1
		c[0] = (l0[0] + l1[0]) / 2;
		c[1] = (l0[1] + l1[1]) / 2;

		float a = point_point_dist(c, l0); // half length of the line segment

		float f[2]; // the projection of o on the perpendicular bisector of l0l1

		if (l1[0] - l0[0] == 0)
		{
			f[0] = l0[0];
			f[1] = o[1];
		}
		else if (l1[1] - l0[1] == 0)
		{
			f[0] = o[0];
			f[1] = l1[1];
		}
		else 
		{
			float k = (l1[1] - l0[1]) / (l1[0] - l0[0]);
			float kk = -1 / k;

			f[0] = -(k*o[0]-kk*c[0]-o[1]+c[1])/(-k+kk);
			f[1] = -(-k*kk*o[0]+k*kk*c[0]+kk*o[1]-k*c[1])/(k-kk);
		}

		float t = point_point_dist(f, o);
		float l = point_point_dist(f, c);

		// find the two solutions for radius rr
		float rr0 = (pow(a,2)*r - pow(l,2)*r + pow(r,3) - r*pow(t,2) - 
			l*sqrt(pow(a,4) + 2*pow(a,2)*pow(l,2) + pow(l,4) - 2*pow(a,2)*pow(r,2) - 2*pow(l,2)*pow(r,2) + pow(r,4) - 
				2*pow(a,2)*pow(t,2) + 2*pow(l,2)*pow(t,2) - 2*pow(r,2)*pow(t,2) + pow(t,4)))/(2.*(pow(l,2) - pow(r,2)));
		float rr1 = (pow(a,2)*r - pow(l,2)*r + pow(r,3) - r*pow(t,2) + 
			l*sqrt(pow(a,4) + 2*pow(a,2)*pow(l,2) + pow(l,4) - 2*pow(a,2)*pow(r,2) - 2*pow(l,2)*pow(r,2) + pow(r,4) - 
				2*pow(a,2)*pow(t,2) + 2*pow(l,2)*pow(t,2) - 2*pow(r,2)*pow(t,2) + pow(t,4)))/(2.*(pow(l,2) - pow(r,2)));
		
		float rr;
		
		// choose the larger one as the solution
		if (rr0 > rr1)
			rr = rr0;
		else
			rr = rr1;

		// compute the cosine value
		float sine = a / rr;
		float cose = sqrt(1 - sine*sine);
		
		// if the angle is larger than pi/2, make the cosine value negetive
		if (r + a > point_point_dist(o, c))
			cose = -cose;

		return cose;
	}
}

SECTION section_new(int dimension, float *p, float *q) // inside means q is inside p
{
    bool inside;
    bool overlap;

    overlap = TRUE;
    inside = TRUE;

    for (int i = 0; i < dimension; i++)
    {
		if (q[2 * i] > p[2 * i + 1] ||  q[2 * i + 1] < p[2 * i])
			overlap = FALSE;
		if (q[2 * i] < p[2 * i] ||
			q[2 * i + 1] > p[2 * i + 1])
			inside = FALSE;
    }
    if (inside)
		return INSIDE;
    else if (overlap)
		return OVERLAP;
    else
		return S_NONE;
}

//added by Tanzima
float Dist(Pointlocation p, Point2D q)
{
	return sqrt((p.x-q[0])*(p.x-q[0])+(p.y-q[1])*(p.y-q[1]));
	
}

float Dist1(Pointlocation p, Point2D q)
{
	return ((p.x-q[0])*(p.x-q[0])+(p.y-q[1])*(p.y-q[1]));
	
}
float Dist(Point2D p, Point2D q)
{
	return sqrt((p[0]-q[0])*(p[0]-q[0])+(p[1]-q[1])*(p[1]-q[1]));
	
}


float MINRECTDIST(Rectangle1 p, float *bounces)
{	
	/*
	//partial overlap
	if (((p.x1>=bounces[0] && p.x1<=bounces[1]) || (p.x2>=bounces[0] && p.x2<=bounces[1])) && (p.y2 <= bounces[2]))
		return pow(bounces[2]-p.y2,2);
	if (((p.x1>=bounces[0] && p.x1<=bounces[1]) || (p.x2>=bounces[0] && p.x2<=bounces[1])) && (p.y1 >= bounces[3]))
		return pow(bounces[3]-p.y1,2);
	if (((p.y1>=bounces[2] && p.y1<=bounces[3]) || (p.y2>=bounces[2] && p.y2<=bounces[3])) && (p.x2 <= bounces[0]))
		return pow(bounces[0]-p.x2,2);
	if (((p.y1>=bounces[2] && p.y1<=bounces[3]) || (p.y2>=bounces[2] && p.y2<=bounces[3])) && (p.x1 >= bounces[1]))
		return pow(bounces[1]-p.x1,2);
	*/

	//disjoint
	if(p.x2 < bounces[0] && p.y2 < bounces[2])
		return pow(bounces[0]-p.x2,2)+pow(bounces[2]-p.y2,2);
	if(p.x2 < bounces[0] && p.y1 > bounces[3])
		return pow(bounces[0]-p.x2,2)+pow(bounces[3]-p.y1,2);
	if(p.x1 > bounces[1] && p.y2 < bounces[2])
		return pow(bounces[1]-p.x1,2)+pow(bounces[2]-p.y2,2);
	if(p.x1 > bounces[1] && p.y1 > bounces[3])
		return pow(bounces[1]-p.x1,2)+pow(bounces[3]-p.y1,2);
	

	//overlap
	if (((p.x1>=bounces[0] && p.x1<=bounces[1]) || (p.x2>=bounces[0] && p.x2<=bounces[1]) || (p.x1<=bounces[0] && p.x2>=bounces[1])) && ((p.y1>=bounces[2] && p.y1<=bounces[3]) || (p.y2>=bounces[2] && p.y2<=bounces[3]) || (p.y1<=bounces[2] && p.y2>=bounces[3])))
		return 0.0;

	//one axis overlap
	if (p.y2 <= bounces[2])
		return pow(bounces[2]-p.y2,2);
	if (p.y1 >= bounces[3])
		return pow(bounces[3]-p.y1,2);
	if (p.x2 <= bounces[0])
		return pow(bounces[0]-p.x2,2);
	if (p.x1 >= bounces[1])
		return pow(bounces[1]-p.x1,2);

	return 0.0;    
}


float MAXRECTDIST(Rectangle1 p, float *bounces)
{  
	float d1, d2, d3, d4;

	d1 = pow(bounces[1]-p.x1,2)+pow(bounces[3]-p.y1,2);
	d2 = pow(bounces[1]-p.x1,2)+pow(bounces[2]-p.y2,2);
	d3 = pow(bounces[0]-p.x2,2)+pow(bounces[3]-p.y1,2);
	d4 = pow(bounces[0]-p.x2,2)+pow(bounces[2]-p.y2,2);

	if(d1 >= d2 && d1 >= d3 && d1 >= d4)
		return d1;
	else if (d2 >= d1 && d2 >= d3 && d2 >= d4)
		return d2;
	else if (d3 >= d1 && d3 >= d2 && d3 >= d4)
		return d3;
	else return d4;
}

float MINRECTDIST1(Pointlocation p, Rectangle1 bounces)
{	
	float summe = 0.0;
    float r;
    int i;

	if (p.x < bounces.x1)
		r = bounces.x1;
	else
	{
		if (p.x > bounces.x2)
			r = bounces.x2;
		else 
			r = p.x;
	}    

	summe += pow(p.x-r,2);

	if (p.y < bounces.y1)
		r = bounces.y1;
	else
	{
		if (p.y > bounces.y2)
			r = bounces.y2;
		else 
			r = p.y;
	}    

	summe += pow(p.y-r,2);
    
    return(summe);
}

float MAXRECTDIST1(Pointlocation p, Rectangle1 bounces)
{
	float d1, d2, d3, d4;

	d1 = pow(bounces.x2-p.x,2)+pow(bounces.y2-p.y,2);
	d2 = pow(bounces.x2-p.x,2)+pow(bounces.y1-p.y,2);
	d3 = pow(bounces.x1-p.x,2)+pow(bounces.y2-p.y,2);
	d4 = pow(bounces.x1-p.x,2)+pow(bounces.y1-p.y,2);

	if(d1 >= d2 && d1 >= d3 && d1 >= d4)
		return d1;
	else if (d2 >= d1 && d2 >= d3 && d2 >= d4)
		return d2;
	else if (d3 >= d1 && d3 >= d2 && d3 >= d4)
		return d3;
	else return d4;
}


//added by Eunus


float MAXGROUPDIST(Rectangle1 R[], int g_size, Point2D x){
	int i;

	float maxDist = 0.0;

	for(i=0; i<g_size; i++){
		Point2D q;
		q[0] = R[i].x1;
		q[1] = R[i].y1;
			//float dist = Dist(x, q);
		float dqi = Dist(x, q);

		if(dqi > maxDist){
			maxDist = dqi;
		}
	}
	
	return maxDist;
}

void GEOCENTROID(Rectangle1 R[], int g_size, Point2D x){
	
	//Point2D x;
	
	int i;

	x[0] = 0.0;
	x[1] = 0.0;

	for(i=0; i<g_size; i++){
			//Point2D q;
			x[0] += R[i].x1;
			x[1] += R[i].y1;
	}
	x[0] /= g_size;
	x[1] /= g_size;

	//return x;
}


float fDIST(float aggDist, float curDist, int f){
	
	float totalDist;
	
	if(f == 1){//sum
		totalDist = aggDist + curDist;
	}else if(f==2){//max
		
		totalDist = aggDist > curDist ? aggDist : curDist;
	}else if(f==3){//min
		totalDist = aggDist < curDist ? aggDist : curDist;
	}else{
		totalDist = aggDist;
	}

	return totalDist;
}

bool TERMINATE_LB(Rectangle1 R[], int g_size, int sg_size, Point2D x, Rectangle1 p, float maxQDist, int f, double BestAggDist[], bool status[])
{

	bool cont = false;
	float lb[128];
	int i;


	for(i=0; i<g_size-sg_size+1; i++){
		lb[i] = 0.0;
	}
	
	Pointlocation xp;
	xp.x = x[0];
	xp.y = x[1];

	float r = MINRECTDIST1(xp, p);


//maxQDist == MAXGROUPDIST(Rectangle1 R[], int g_size, Point2D x)
	if(r > maxQDist){
		cont = true;
	}else{
		
		int sq[128]; //sorted query index
		float smq[128];//

		//initialize:
		/*
		for(i=0; i<g_size; i++){
			smq[i] = MAXDOUBLE;
			sq[i] = i; 
		}
		*/

		for(i=0; i<g_size; i++){


			Point2D q;
			q[0] = R[i].x1;
			q[1] = R[i].y1;
			//float dist = Dist(x, q);
			float mqi = r-Dist(x, q);
			
			//initialize one row
			smq[i] = MAXREAL;
			sq[i] = i; 
			//end init

			//insert in a sorted list
			for(int l=0; l<i+1; l++)
			{
				if (mqi < smq[l])
				{
					for(int j=i; j>l; j--)
					{
						//maxdist[j]=maxdist[j-1];	
						smq[j] = smq[j-1];
						sq[j] = sq[j-1];
					}
					//maxdist[l]=edist2;
					smq[l] = mqi;
					sq[l] = i;

					break;
				}
			}
			//end sorted list sqi,smqi

		}//end of creating list based on query distance

		i = 1;
		float aggDist = smq[0];
		while(i < g_size){
			
			aggDist = fDIST(aggDist,smq[i],f);
			if(i >= sg_size && aggDist > BestAggDist[i-sg_size+1]){
				//lb[i-sg_size+1] = aggDist;
				status[i-sg_size+1] = true; // final
			}
		}

		for(i=0; i<g_size-sg_size+1; i++){
			if(status[i] == false){ //if one is false we need to keep searching
				cont = true;
			}
		}

		
	}


	//double BestAggDist[g_size];
	//bool status[g_size];		   
	//Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][], Pointlocation rslt[], int *num_of_data

	return cont;
}

