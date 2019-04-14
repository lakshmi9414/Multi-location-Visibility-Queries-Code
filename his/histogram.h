/* The implementation of the histogram file reader
   by Bobby Last modified: Oct.15, 2002 */

#ifndef __HISTOGRAM
#define __HISTOGRAM

class Histogram
{
public:
	int dim;			// dimensionality
	float *bounces;		// bounces of the buckets
	int *num;			// number of objects falling into the buckets
	int totalnum;		// total number of objects
	int num_bucket;		// number of buckets
	int max_bucket;		// maximum number of buckets
	float univarea;			// area of the universe

	Histogram(int dim = 2, int max_bucket = 10000);		// constructor, dim: dimensionality, max_buckets: maximum number of buckets
	~Histogram();										// destructor
	void read(char* hisfile);							// read histogram from a file
	float wq(float *bounces);							// number of objects returned estimation for window query
	float wq_N(float *bounces);							// total number of objects estimated for window query by expanding the density to the universe
	float wintq_N(float *bounces);						// total number of objects estimated for window query (considering only buckets near the boundary of query window) by expanding the density to the universe
	float nnq_N(float *bounces, float radius = 0.0f);	// total number of objects estimated for nearest neighbor query by expanding the density to the universe
};

#endif
