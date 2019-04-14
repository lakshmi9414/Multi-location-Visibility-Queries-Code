/*rtree_cmd.cpp
  this file implements the class RTreeCmdIntrpr*/

#include <stdio.h>
#include<cstring>
#include <string.h>
#include <stdlib.h>
#include "rtree_cmd.h"
#include "rtree.h"
#include "../blockfile/cache.h"
#include "../blockfile/blk_file.h"
#include "../linlist/linlist.h"
using namespace std;
//----------------------------------------------------------
RTreeCmdIntrpr::RTreeCmdIntrpr() : CmdIntrpr()
{
	version();
	tree = NULL;
	cache = NULL;
	o_blen = 4096;
	o_dim = 2;
	o_csize = 0;
	
}
//----------------------------------------------------------
RTreeCmdIntrpr::~RTreeCmdIntrpr()
{
	if (tree)
	{
		if (cnfrm_cmd("The tree has not be freed.  Free it now"))
			free_tree();
	}

	if (cache)
		delete cache;
}
//----------------------------------------------------------
bool RTreeCmdIntrpr::build_tree(char *_tree_fname, char *_data_fname, int _b_len, int _dim, int _csize)
{
	
	if (tree)
	{
		if (cnfrm_cmd("The tree has not be freed.  Free it now"))
		{
			//printf("Freeing the tree (possibly writing to disk)\n");
			delete tree; tree = NULL;
		}
		else
		{
			//printf("No tree was ructed.\n"); return false;
		}
	}

	if (cache)
		delete cache;

	cache = NULL;
	cache = new Cache(_csize, _b_len);
	
	//printf("Building in progress..\n");
	tree = new RTree(_data_fname, _tree_fname, _b_len, cache, _dim);
	//printf("Tree constructed from %s and residing in memory\n", _data_fname);
	return true;
}
//----------------------------------------------------------
bool RTreeCmdIntrpr::chk_file_exist(char *_fname)
{
	FILE *f = fopen(_fname, "r");
	if (f)
	{
		fclose(f); return true;
	}
	else return false;
}
//----------------------------------------------------------
void RTreeCmdIntrpr::free_tree()
{
	if (tree)
	{
		delete tree; tree = NULL; 
		//printf("Tree freed.\n");
	}

	if (cache)
	{
		delete cache; cache = NULL;
	}
}
//----------------------------------------------------------
int RTreeCmdIntrpr::knn_sngle(float *_mbr, int *_io_count)
{
	if (!cache)
	{ //printf("Error: Cache is not available\n"); 
	  *_io_count = 0; return 0;};

	int pgflt_bfr = cache -> page_faults;

	SortedLinList *res_list = new SortedLinList();
	//tree -> NNQuery(_mbr, res_list);
	int res_num = res_list -> get_num();
	delete res_list;

	int pgflt_aft = cache -> page_faults;
	*_io_count = pgflt_aft - pgflt_bfr;

	return res_num;
}
//----------------------------------------------------------
int RTreeCmdIntrpr::qry_sngle(float *_mbr, int *_io_count)
{
	if (!cache)
	{ //printf("Error: Cache is not available\n"); 
	  *_io_count = 0; return 0;};

	int pgflt_bfr = cache -> page_faults;

	SortedLinList *res_list = new SortedLinList();
	tree -> rangeQuery(_mbr, res_list);
	int res_num = res_list -> get_num();
	delete res_list;

	int pgflt_aft = cache -> page_faults;
	*_io_count = pgflt_aft - pgflt_bfr;

	return res_num;
}
//----------------------------------------------------------
void RTreeCmdIntrpr::run()
{
	bool RUNNING = true;
	while (RUNNING)
	{
		char cmd[100];
		get_cmd(">", cmd);
		//strupr(cmd);
		toupper(*cmd); // changed by me
		
		if (strncmp(cmd, "EXIT", strlen(cmd) + 1) == 0)
		{
			if (cnfrm_cmd("Really quit Tree Runner"))
			{
				printf("Quitting...\n");
				RUNNING = false;
			}
		}
		if (strncmp(cmd, "BUILD", strlen(cmd) + 1) == 0)
		{
			bool again = true;
			char tree_fname[100], data_fname[100];
			int blen;
			FILE *f;

			while (again)
			{
				again = false;
				get_cmd("Tree file name: ", tree_fname);
			    if (chk_file_exist(tree_fname))
				{
					printf("The tree file already exists.\n");
					again = true;
				}
			}

			again = true;
			while (again)
			{
				again = false;
				get_cmd("Data file name: ", data_fname);
				if (!chk_file_exist(data_fname)) 
				{
					printf("The data set file does not exist.\n");
					again = true;
				}
			}

			again = true;
			while (again)
			{
				again = false;

				printf("Dimension(%d): ", o_dim);

				char dim_char[100];
				get_cmd("", dim_char);
				if (strlen(dim_char) > 0)
				{
					int new_dim = atoi(dim_char);
					if (new_dim == 0)
						again = true;
					else
						o_dim = new_dim;
				}
			}

			get_blen();
			get_csize();
	
			if (cnfrm_cmd("Proceed"))
				build_tree(tree_fname, data_fname, o_blen, o_dim, o_csize);
		}

		if (strncmp(cmd, "FREE", strlen(cmd) + 1) == 0)
			free_tree();

		if (strncmp(cmd, "LOAD", strlen(cmd) + 1) == 0)
		{
			bool again = true;
			char tree_fname[100];

			while(again)
			{
				again = false;
				get_cmd("Tree file name: ", tree_fname);
				if (!chk_file_exist(tree_fname))
				{
					printf("The tree does not exist!\n");
					again = true;
				}
			}


			get_blen();
			get_csize();
			
			free_tree();
			cache = new Cache(o_csize, o_blen);

			if (cnfrm_cmd("Proceed"))
			{
				tree = new RTree(tree_fname, cache);
				if (o_csize > 0 && tree -> file -> get_blocklength() != o_blen)
					printf("Warning: the block length of cache does not agree with the tree!\n");
				printf("Tree residing in memory.\n");
			}
		}

		if (strncmp(cmd, "QRYSNGL", strlen(cmd) + 1) == 0)
		{
			if (tree)
			{
				int dim = tree -> dimension;
				float *mbr = new float[2 * dim];

				for (int i = 0; i < 2 * dim; i++)
				{
					char num_char[100];
					printf("Enter mbr[%d]: ", i);
					get_cmd("", num_char);
					mbr[i] = atof(num_char);
				}

				int io_count = 0, res_count = 0;
				res_count = qry_sngle(mbr, &io_count);

				printf("%d records retrieved in %d I/O's\n", res_count, io_count);

				delete [] mbr;
			}
			else
				printf("No tree loaded\n");
		}

		if (strncmp(cmd, "NNSNGL", strlen(cmd) + 1) == 0)
		{
			if (tree)
			{
				int dim = tree -> dimension;
				float *mbr = new float[dim];

				for (int i = 0; i < dim; i++)
				{
					char num_char[100];
					printf("Enter mbr[%d]: ", i);
					get_cmd("", num_char);
					mbr[i] = atof(num_char);
				}

				int io_count = 0, res_count = 0;
				res_count = knn_sngle(mbr, &io_count);

				printf("%d records retrieved in %d I/O's\n", res_count, io_count);

				delete [] mbr;
			}
			else
				printf("No tree loaded\n");
		}

		if (strncmp(cmd, "FLUSH", strlen(cmd) + 1) == 0)
		{
			if (cache)
			{
				cache -> flush();
				printf("Cache flushed\n");
			}
			else
				printf("No cache initialized\n");
		}
	}
}
//----------------------------------------------------------
void RTreeCmdIntrpr::version()
{
	//printf("Tree Runner V1.0 for RTree\n");
	//printf("Designed by TAO Yufei, Dec. 2000\n");
}
//=========================================================
int RTreeCmdIntrpr::get_blen()
{
	bool again = true;
	while (again)
	{
		again = false;
		
		printf("Block length(%d): ",o_blen);

		char bln_char[100];
		get_cmd("", bln_char);

		if (strlen(bln_char) > 0)
		{
			int new_bln = atoi(bln_char);
			if (new_bln == 0)
				again = true;
			else
				o_blen = new_bln;
		}
	}

	return o_blen;
}
//=========================================================
int RTreeCmdIntrpr::get_csize()
{
	bool again = true;
	while (again)
	{
		again = false;

		printf("Cache size(%d): ", o_csize);

		char csize_char[100];
		get_cmd("", csize_char);
		if (strlen(csize_char) > 0)
		{
			int new_csize = atoi(csize_char);
			if (new_csize < 0)
				again = true;
			else
				o_csize = new_csize;
		}
	}

	return o_csize;
}
