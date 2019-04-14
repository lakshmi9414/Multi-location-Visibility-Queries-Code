/*rtree_cmd.h
  this file defines the command interpretor RTreeCmdIntrpr*/
#ifndef _RTREE_CMD_H
#define _RTREE_CMD_H

#include <vector>
#include "../func/gendef.h"
//----------------------------------------------------------
class RTree;
class Cache;

class RTreeCmdIntrpr : public CmdIntrpr
{
public:
	RTree *tree;
	Cache *cache;
	int o_blen;
	int o_dim;
	int o_csize;
	

	RTreeCmdIntrpr();
	~RTreeCmdIntrpr();

	//=======================================
	int get_blen();
	int get_csize();
	//=======================================
	bool build_tree(char *_tree_fname, char *_data_fname, int _b_len, int _dim, int _csize);
	bool chk_file_exist(char *_fname);
	void free_tree();
	int knn_sngle(float *_mbr, int *_io_count);
	int qry_sngle(float *_mbr, int *_io_count);
	void run();
	void version();
};

#endif
