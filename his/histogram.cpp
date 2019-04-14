/* The implementation of the histogram file reader
   by Bobby Last modified: Oct.15, 2002 */

#include <stdio.h>
#include <stdlib.h>
#include "../func/gendef.h"
#include "histogram.h"

Histogram::Histogram(int dim, int max_bucket)
{
	this->dim = dim;
	bounces = new float[2*dim*max_bucket];
	num = new int[max_bucket];
	totalnum = 0;
	this->max_bucket = max_bucket;
	num_bucket = 0;
	univarea = 0.0f;

	int i;
	for (i = 0; i < max_bucket; i++)
	{
		int j;
		for (j = 0; j < 2*dim; j++)
			bounces[2*dim*i+j] = 0.0f;
		num[i] = 0;
	}
}

Histogram::~Histogram()
{
	if (bounces)  delete []bounces;
	if (num)  delete []num;
}

void Histogram::read(char* hisfile)
{
	// read histogram
	FILE *fp = fopen(hisfile, "r");

	if (!fp)
	{
		printf("Cannot open histogram file.\n");
		return;
	}
	
	num_bucket = 0;

	float *min = new float[dim];
	float *max = new float[dim];
	
	int i;
	for (i = 0; i < dim; i++)
	{
		min[i] = 1e30f;
		max[i] = -1e30f;
	}

	while (!feof(fp))
	{
		int id;
		int read = fscanf(fp, "%d", &id);
		if (read != 1)
			break;

		for (i = 0; i < dim; i++)
			fscanf(fp, "%f %f", &(bounces[2*dim*num_bucket+2*i]), &(bounces[2*dim*num_bucket+2*i+1]));

		fscanf(fp, "%d", &(num[num_bucket]));
		totalnum += num[num_bucket];

		fscanf(fp, "\n");

		// check the min & max coordinates of the buckets to determine the universe size
		for (i = 0; i < dim; i++)
		{
			if (bounces[2*dim*num_bucket+2*i] < min[i])
				min[i] = bounces[2*dim*num_bucket+2*i];
			if (bounces[2*dim*num_bucket+2*i+1] > max[i])
				max[i] = bounces[2*dim*num_bucket+2*i+1];
		}

		num_bucket++;

		if (num_bucket >= max_bucket)
		{
			printf("max_bucket is too small when reading histogram.\n");
			exit(1);
		}
	}

	fclose(fp);

	// compute the area of the universe
	univarea = 1.0f;
	for (i = 0; i < dim; i++)
		univarea *= max[i] - min[i];

	printf("%d buckets read from histogram.\n", num_bucket);
}

float Histogram::wq(float *bounces)
{
	float rslt = 0.0f;
	
	int i;
	for (i = 0; i < num_bucket; i++)
	{
		int section = section_new(dim, bounces, &(this->bounces[2*dim*i]));
		
		if (section == INSIDE)
			rslt += num[i];
		else if (section == OVERLAP)
			rslt += num[i]*overlap(dim, bounces, &(this->bounces[2*dim*i]))/area(dim, &(this->bounces[2*dim*i]));
	}

	return rslt;
}

float Histogram::wq_N(float *bounces)
{
	return wq(bounces) / area(dim, bounces) * univarea;
}

float Histogram::wintq_N(float *bounces)
{
	int num_objs = 0;
	float total_area = 0.0f;
	
	int i;
	int count = 0;
	for (i = 0; i < num_bucket; i++)
	{
		int section = section_new(dim, bounces, &(this->bounces[2*dim*i]));
		
		if (section == OVERLAP)
		{
			num_objs += num[i];
			total_area += area(dim, &(this->bounces[2*dim*i]));
			count++;
		}
	}

	// printf("count: %d  num_objs: %d\n", count, num_objs);
	return num_objs / total_area * univarea;
}

float Histogram::nnq_N(float *bounces, float radius)
{
	return 0.0f;
}
