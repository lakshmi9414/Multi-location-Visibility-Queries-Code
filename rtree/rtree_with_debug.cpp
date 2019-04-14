/*rtree.cpp
  this file implements the RTree class*/
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "rtree.h"
#include "entry.h"
#include "rtnode.h"
#include "distance.h"
#include "../blockfile/cache.h"
#include "../blockfile/blk_file.h"
#include "../linlist/linlist.h"


//Added by Tanzima
extern int io_access;
extern int disktime;
extern int updatetime;
extern int counttime;
extern int kmintime;
//......................
//------------------------------------------------------------
RTree::RTree(char *fname, int _b_length, Cache *c, int _dimension)
  //use this constructor to build a new tree
{
    file = new BlockFile(fname, _b_length);
    cache = c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
	  //note that when a tree is constructed, the root is automatically created
	  //though at this time there is no entry in the root yet.
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr -> block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
RTree::RTree(char *fname, Cache *c)
  //use this constructor to restore a tree from a file
{
    file = new BlockFile(fname, 0);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    file -> read_header(header);
    read_header(header);
	delete [] header;

    root_ptr = NULL;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
RTree::RTree(char *inpname, char *fname, int _b_length, Cache *c, int _dimension)
  // construct new R-tree from a specified input textfile with rectangles
{
    Entry *d;
    FILE *fp;
    file = new BlockFile(fname, _b_length);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);
	delete [] header;

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr->block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;

	int record_count = 0;

    if((fp = fopen(inpname,"r")) == NULL)
    {
      delete this;
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
      while (!feof(fp))
      {
		record_count ++;

		d = new Entry(dimension, NULL);

    	fscanf(fp, "%d", &(d -> son));
		//printf("ID=%d ",d->son);
		//for (int i = 0; i < 2 * dimension; i ++)
		//{
			//fscanf(fp, " %f", &(d -> bounces[i]));
		    fscanf(fp, " %f %f %f %f", &(d->bounces[0]),
		 	&(d->bounces[1]), &(d->bounces[2]), &(d->bounces[3]));
		//}
		//fscanf(fp, "\n");

		//if(record_count==20000 || record_count==20001)	printf(": %f %f %f %f\n", d->bounces[0], d->bounces[2], d->bounces[1], d->bounces[3]);
    	insert(d);

		  //d will be deleted in insert()

		if (record_count % 100 == 0)
		{
			for (int i = 0; i < 79; i ++)  //clear a line
				printf("\b");

			printf("inserting object %d", record_count);
		}
      }
    }
	//TANZIMA
	printf("inserting object %d", record_count);
	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
RTree::~RTree()
{
	char *header = new char[file -> get_blocklength()];
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL)
    {
        delete root_ptr;
        root_ptr = NULL;
    }

	if (cache)
      cache -> flush();

    delete file;

    delete re_data_cands;
	delete deletelist;

    //printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
	//   num_of_inodes, num_of_dnodes, num_of_data);
}
//------------------------------------------------------------
void RTree::del_root()
{
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
bool RTree::delete_entry(Entry *d)
{
	load_root();

	R_DELETE del_ret;
	del_ret=root_ptr->delete_entry(d);

	if (del_ret == NOTFOUND) return false;
	if (del_ret == ERASED) 
		error("RTree::delete_entry--The root has been deleted\n",true);
 
	if (root_ptr -> level > 0 && root_ptr -> num_entries == 1)
		//there is only one entry in the root but the root
		//is not leaf.  in this case, the child of the root is exhalted to root
	{
		root = root_ptr -> entries[0].son;
		delete root_ptr;
		root_ptr = NULL;
		load_root();
		num_of_inodes--;
	}

	//Now will reinsert the entries
	while (deletelist -> get_num() > 0)
	{
		Linkable *e;
		e = deletelist -> get_first();
		Entry *new_e = new Entry(dimension, NULL);
		new_e -> set_from_Linkable(e);
		deletelist -> erase();
		insert(new_e);
	}

	delete root_ptr;
	root_ptr = NULL;

	return true;
}
//------------------------------------------------------------
void RTree::insert(Entry* d)
{
    int i, j;
    RTNode *sn;
    RTNode *nroot_ptr;
    int nroot;
    Entry *de;
    R_OVERFLOW split_root;
    Entry *dc;
    float *nmbr;

    // load root into memory
    load_root();

    // no overflow occured until now
    re_level = new bool[root_ptr -> level + 1];
    for (i = 0; i <= root_ptr -> level; i++)
        re_level[i] = FALSE;

    // insert d into re_data_cands as the first entry to insert
    // make a copy of d because it should be erased later
    Linkable *new_link;
	new_link = d -> gen_Linkable();
	re_data_cands -> insert(new_link);

	delete d;  //we follow the convention that the entry will be deleted when insertion finishes

    j = -1;
    while (re_data_cands -> get_num() > 0)
    {
        // first try to insert data, then directory entries
	    Linkable *d_cand;
		d_cand = re_data_cands -> get_first();
        if (d_cand != NULL)
        {
            // since "erase" deletes the data itself from the
            // list, we should make a copy of the data before
            // erasing it
			dc = new Entry(dimension, NULL);
            dc -> set_from_Linkable(d_cand);
            re_data_cands -> erase();

            // start recursive insert with root
			split_root = root_ptr -> insert(dc, &sn);
        }
        else
	        error("RTree::insert: inconsistent list re_data_cands", TRUE);

    	if (split_root == SPLIT)
    	// insert has lead to split --> new root-page with two sons (i.e. root and sn)
    	{
    	    nroot_ptr = new RTNode(this);
    	    nroot_ptr -> level = root_ptr -> level + 1;
    	    num_of_inodes++;
    	    nroot = nroot_ptr -> block;

    	    de = new Entry(dimension, this);
    	    nmbr = root_ptr -> get_mbr();
    	    memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de->son = root_ptr->block;
    	    de->son_ptr = root_ptr;
    	    nroot_ptr -> enter(de);

    	    de = new Entry(dimension, this);
    	    nmbr = sn -> get_mbr();
    	    memcpy(de -> bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de -> son = sn -> block;
    	    de -> son_ptr = sn;
    	    nroot_ptr->enter(de);

    	    root = nroot;
            root_ptr = nroot_ptr;

            root_is_data = FALSE;
        }
        j++;
    }

    num_of_data++;

    delete [] re_level;

	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
void RTree::load_root()
{
    if (root_ptr == NULL)
        root_ptr = new RTNode(this, root);
}
//------------------------------------------------------------
void RTree::rangeQuery(float *mbr, SortedLinList *res)
{
    load_root();
	//Added by Tanzima
	io_access++;
	//..............
    root_ptr -> rangeQuery(mbr,res);

	delete root_ptr;
	root_ptr = NULL;	
}
//------------------------------------------------------------
void RTree::read_header(char *buffer)
{
    int i;

    memcpy(&dimension, buffer, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
int RTree::update_rslt(Entry *_e, float _dist, Entry *_rslt, 
						float *_key, int _k)
{
	for (int i = 0; i < _k; i ++)
	{
		if (_dist < _key[i])
		{
			for (int j = _k - 1; j > i; j --)
			{
				_rslt[j] = _rslt[j - 1];
				_key[j] = _key[j - 1];
			}
			_rslt[i] = *_e;
			_key[i] = _dist;
			return i;
		}
	}
	error("Error in update_rslt\n", true);
	return -1;
}
//------------------------------------------------------------
void RTree::write_header(char *buffer)
{
    int i;

    memcpy(buffer, &dimension, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
//---added for valdity region queries-------------------------
// perform a window query mbr to get the query result into in_objs, put the outer objects retrieved into out_objs_so_far
void RTree::rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far)
{
    load_root();

    root_ptr->rect_win_query(mbr, in_objs, out_objs_so_far);

	delete root_ptr;
	root_ptr = NULL;
}

// perform a window query mbr (excluding the window excr) to get the query result into c_inf_objs
void RTree::rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs)
{
    load_root();

    root_ptr->rect_win_query(mbr, exclmbr, c_inf_objs);

	delete root_ptr;
	root_ptr = NULL;
}

void RTree::BFNN(float *_qpt, int _k, Entry *_rslt)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = (float) MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				if (edist < key[_k - 1])
					update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
			}
			else
			{
				if (edist<key[_k - 1]) 
					//meaning that edist is valid and we insert it to heap
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					heap -> insert(he);
					delete he;
				}
			}
		}
	
		delete rtn;

		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->key>key[_k - 1]) //the algorithm terminates
					son = -1;
				else if (he->level == 0) 
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
					 else
						son=he->son1;
			}
		}
		delete he;
		//--------------------------------------------------------
	}

	delete [] key;
	delete heap;
}



//RTree::BFNN --- Modified by Tanzima for Rectangle

float RTree::KMin(Point2D m, int c, float cl, int k, Pointlocation _rslt[], int *num_of_data, float d_safe_corner,DistfromPOI _cornertoPOI[])
{
	float *kmindist = new float[k];
	float tempdist;
	for(int i=0; i<k; i++)
	{
			kmindist[i]=(float) MAXREAL;
	}	
	
	for(int i=0; i<*num_of_data;i++)
	{
		if(cl* _cornertoPOI[i].d[c] <= d_safe_corner)
		{
			tempdist=Dist(_rslt[i],m);			

			for(int j=0; j<k; j++)
			{				

				if(kmindist[j]>tempdist)
				{
					for(int l=k-1;l>j;l--)
					{
						kmindist[l]= kmindist[l-1];
					}
					kmindist[j]= tempdist;
					break;
				}
			}
		}
		
	}
	return kmindist[k-1];
}


void RTree::UpdateCount(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, float d_safe_corner, int *count, DistfromPOI _cornertoPOI[])
{
	Point2D c[4];
	c[0][0] = R.x1; c[0][1]=R.y1;
	c[1][0] = R.x1; c[1][1]=R.y2;
	c[2][0] = R.x2; c[2][1]=R.y2;
	c[3][0] = R.x2; c[3][1]=R.y1;
	
		
	_cornertoPOI[*num_of_data-1].d[0]= Dist(_rslt[*num_of_data-1],c[0]);
	_cornertoPOI[*num_of_data-1].d[1]= Dist(_rslt[*num_of_data-1],c[1]);
	_cornertoPOI[*num_of_data-1].d[2]= Dist(_rslt[*num_of_data-1],c[2]);
	_cornertoPOI[*num_of_data-1].d[3]= Dist(_rslt[*num_of_data-1],c[3]);
	
	for(int j=0; j<4; j++)
	{
		if(count[j]<k)
		{
			count[j]=0;
			for(int i=0; i<*num_of_data;i++)
			{
				if(cl*_cornertoPOI[i].d[j] <= d_safe_corner)	count[j]++;
				
			}
		}
	}
	/*
	for(int i=0; i<_rslt.size();i++)
	{
		if(i<_cornertoPOI.size())
		{
			if(cl*_cornertoPOI[i].d[0] <= d_safe_corner)	count[0]++;
			if(cl*_cornertoPOI[i].d[1] <= d_safe_corner)	count[1]++;
			if(cl*_cornertoPOI[i].d[2] <= d_safe_corner)	count[2]++;
			if(cl*_cornertoPOI[i].d[3] <= d_safe_corner)	count[3]++;
		
		}
		else
		{
			Pointlocation p = _rslt[i];
			dist_c_P.d[0]= Dist(p,c[0]);
			if(cl*dist_c_P.d[0] <= d_safe_corner)	count[0]++;
			dist_c_P.d[1]= Dist(p,c[1]);
			if(cl*dist_c_P.d[1] <= d_safe_corner)	count[1]++;
			dist_c_P.d[2]= Dist(p,c[2]);
			if(cl*dist_c_P.d[2] <= d_safe_corner)	count[2]++;
			dist_c_P.d[3]= Dist(p,c[3]);
			if(cl*dist_c_P.d[3] <= d_safe_corner)	count[3]++;
			_cornertoPOI.push_back(dist_c_P);

		}
	}
	*/
}

int RTree::UpdateStatus(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, int *count, DistfromPOI _cornertoPOI[])
{
	Point2D o;
	o[0] = (R.x1+R.x2)/2;
	o[1] = (R.y1+R.y2)/2;
	float r = Dist(_rslt[*num_of_data-1],o);
	
	Point2D c[4];
	c[0][0] = R.x1; c[0][1]=R.y1;
	c[1][0] = R.x1; c[1][1]=R.y2;
	c[2][0] = R.x2; c[2][1]=R.y2;
	c[3][0] = R.x2; c[3][1]=R.y1;

	float d1=Dist(c[0],c[1]);
	float d2=Dist(c[1],c[2]);
	float d_safe_corner = r-0.5*sqrt(d1*d1+d2*d2);

	

	
	UpdateCount(R,cl,k,_rslt,num_of_data,d_safe_corner,count,_cornertoPOI);
	

	for(int i=0; i<4; i++)
	{
		if (count[i]<k)	return 0;
	}

	

	float d_i, d_j,d_max_m,d_max = 0;
	int i,j;
	Point2D m;
	for (i=0; i<4; i++)
	{
		j=(i+1)%4;
		
		m[0]= (c[i][0]+c[j][0])/2;
		m[1]= (c[i][1]+c[j][1])/2;
		d_i= KMin(m,i,cl,k,_rslt,num_of_data,d_safe_corner,_cornertoPOI);
		d_j= KMin(m,j,cl,k,_rslt,num_of_data,d_safe_corner,_cornertoPOI);
		if(d_i>d_j) d_max_m = d_i;
		else d_max_m = d_j;
		if(d_max_m > d_max)
			d_max = d_max_m;
		
	}
	float c_dmax=d1;
	if(d2>d1)	c_dmax=d2;
	float d_safe = r- 0.5 * c_dmax;	

	if(cl*d_max>d_safe)
		return ceil(r+cl*(d_max-d_safe));
	else
		return -1;
}
void RTree::Rect_kNNQ(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data)
{
	int io=0;
	//init status
	int count[4];
	for(int i=0; i<4; i++)	count[i]=0;
	int status=0;
	Point2D o;
	o[0] = (R.x1+R.x2)/2;
	o[1] = (R.y1+R.y2)/2;
	DistfromPOI cornertoPOI[10000];	

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	while (son != -1 && status != -1)
	{
		
		
		RTNode *rtn = new RTNode(this, son);
		//Experiment
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			float o1[2];
			o1[0]=(float)o[0];
			o1[1]=(float)o[1];
			float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;

		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{	
				if (he->level == 0) //p is an object 
				{
					if(*num_of_data==10000)
						//printf("\nGreater than 10000\n");
						error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
					
					_rslt[*num_of_data].x=he->x1;
					_rslt[*num_of_data].y=he->y1;
					*num_of_data=*num_of_data+1;
					//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);
						
					
					
					if(status==0)
					{												
						status = UpdateStatus(R,cl,k,_rslt,num_of_data,count,cornertoPOI);
												
						if(status != -1)	again=true;						
					}
					else if (status < Dist(_rslt[*num_of_data-1],o))
					{							
						status = -1;
					}
					else
					{
						again = true;
					}
				}
				else
				{
					if(status>0)
					{
						if(status<he->key)
							status = -1;
						else
							son=he->son1;
					}
					else
					{
						son=he->son1;
					}
				}
			}
		}
		
		delete he;
	}
	delete heap;
}

void RTree::Point_BFN_NNQ(Point2D o, double *_rslt)
{
	
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		io_access++;
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float o1[2];
			o1[0]=(float)o[0];
			o1[1]=(float)o[1];
			float edist = MINDIST(o1, rtn->entries[i].bounces, dimension);
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;

		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->level == 0) //p is an object 
				{
					_rslt[0] = he->x1;
					_rslt[1] = he->y1;
					son=-1;
														
				}
				else
				{
					son=he->son1;
					
				}
			}
		}
		
		delete he;
	}
	delete heap;
}
//END




// The following code was copied from the implementation of TP-kNN queries from Tony
void RTree::TPNN_TP(float *_qline, int _k, Entry *_nn, Entry *_rslt, float _max_trvl)
{
	float key = _max_trvl; 
	  // the minimum distance that the query point must travel 

//we comment these lines to avoid initing the heap everytime--
//this function is called
//Heap *heap = new Heap();
//heap->init(dimension);
//------------------------------------------------------------
	if (tpheap==NULL)
		error("tpheap is not initialized\n", true);
	tpheap->used=0;
	
	int son = root;
	while (son != -1)
	{
		RTNode *rtn = new RTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			//first create an array m2 for e to cal dist----------
			float *m2 = new float [4 * dimension];
			memset(m2, 0, 4 * dimension * sizeof(float));
			memcpy(m2, rtn -> entries[i].bounces, 2 * dimension * sizeof(float));
			//----------------------------------------------------
//testing--------------------------
//if (rtn->entries[i].son==573673 && cnt==84-1)
//	printf("testing...\n");
//---------------------------------

			float edist = (float) MAXREAL;
			if (rtn -> level == 0)
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^leaf node case^^^^^^^^^^^^^^^^^^^^^
			{	
				int eNNsub=-1;
				//if e (i.e., m2) is one of the current NN points-
				//eNNsub stores its subsript in _nn; otherwise, 
				//eNNsub=-1
				for (int j=0; j<_k; j++)
					if (rtn->entries[i].son==_nn[j].son)
					{ eNNsub=j; j=_k;}
				//------------------------------------------------
		
				if (eNNsub==-1 || SEQUENCE_SENSITIVE)
					//namely, if sequence insensitive and e is indeed a NN
					//then we do not handle this entry
				{

					float *m1 = new float [4 * dimension];
					//find the NN that leads to the minimum-------
					//influence time
					int nn_num=-1; //the subsript of the NN to be found
					for (int j = 0; j < _k; j ++)
						//for each of the NN found in the 1st step
					{
						bool yesdo=true; //whether to compute
						  //the inflence time of nn[j] and e
						if (j==eNNsub) yesdo=false;
						
						//check if this pair has produced --------
						//influence time before
						for (int l=0; l<PAST_PAIR; l++)
							if (min(_nn[j].son, rtn->entries[i].son)==min(last_pair[l][0], last_pair[l][1]) &&
								max(_nn[j].son, rtn->entries[i].son)==max(last_pair[l][0], last_pair[l][1]))
							{	yesdo=false; l=PAST_PAIR; }
						//----------------------------------------

						if (yesdo)
						{
//these codes use NNinf===========================================
/*
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							//get the influence time of m2-------- 
							//(for entry e) with respect to(nn[j])
							float this_inf = NNinf(m1, m2, _qline, dimension); 
							//------------------------------------
*/
//================================================================

//these codes use NNinf2==========================================
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							m1[1]=m1[2]; m2[1]=m2[2];
							//create an arry m3 for _qline----------------
							float *m3=new float[2*dimension];
							m3[0]=_qline[0]; m3[1]=_qline[2];
							m3[2]=_qline[4]; m3[3]=_qline[6];
							//--------------------------------------------
							//get the influence time of m2-------- 
							//(for entry e) with respect to(nn[j])
							float this_inf = NNinf2(m1, m2, m3); 
							//------------------------------------
							delete []m3;
//================================================================

							if (this_inf>0 && this_inf<edist)
								//this_inf=0 means that there is another point that has the same distance
								//to the current query position as the current NN. In this implementation,
								//we choose to ignore handling such special cases, which, however, may cause
								//problems for datasets with big cardinality
//							if (this_inf>=0 && this_inf<edist)
							{
								edist=this_inf; nn_num=j;
							}
						}  //END if (yesdo)
					}//END checking all neighbors
					//-------------------------------------------------
					//if (edist<key && edist!=0)
					if (edist<key)
					{
						update_rslt(&(rtn->entries[i]), edist, _rslt, &key, 1);
						_rslt->nn_num=nn_num;
					}
					delete []m1;
				}
			}
//^^^^^^^^^^^^^^^^^^^^^^^^^non-leaf node case^^^^^^^^^^^^^^^^^^^^^
			else
				//Next handle non-leaf node case
			{
				float *m1 = new float [4 * dimension];
				for (int j = 0; j < _k; j ++)
				{
					//first create an array m1 to cal dist--------
					memset(m1, 0, 4 * dimension * sizeof(float));
					memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
					//--------------------------------------------
					float this_mininf = NNmininf(m1, m2, _qline, dimension);
					if (this_mininf < edist)
						edist = this_mininf;
				}
				delete [] m1;

				if (edist < key)
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					tpheap -> insert(he);
					delete he;
				}
			}
			delete [] m2;

		}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	
		delete rtn;

		//get next entry from the heap
		bool again = true;
		while (again)
		{
			again = false;
			HeapEntry *he = new HeapEntry();
			if (!tpheap -> remove(he))  //heap is empty
				son = -1;
			else
				if (he -> key > key)
					//the algorithm can terminate
					son = -1;
				else
					son = he -> son1;
			delete he;
		}
	}
//delete heap;
}



//RTree::private_kGNN --- Created by Tanzima for kGNN



void RTree::private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{	
	//variable initialization
	int end=0;

	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;
	
	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;
	
	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	
	while (son != -1 && end==0)
	{		
		RTNode *rtn = new RTNode(this, son);
		
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0;
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);	
			if((g_size*edist1)>maxdist[k-1])	continue;
			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{				
				edist1 += MINRECTDIST(R[j], rtn->entries[i].bounces);
			}

			if(edist1>maxdist[k-1])	continue;
			
			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{	
				edist2 += MAXRECTDIST(R[j], rtn->entries[i].bounces);
			}

			//update maxdistk		
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];						
					}
					maxdist[l]=edist2;
					break;
				}
			}
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

			}
			else
			{	
				if (he->level == 0) //p is an object 
				{
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;


						*num_of_data=*num_of_data+1;
						
						//get next data  from heap
						again =true;						
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}
		
		delete he;
	}
	delete heap;
}


//END


void RTree::private_kGNN_max(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data)
{
	//open debug file
	FILE * dFile, * dFile1, * dFile2;
	dFile = fopen( "C:/GRPLBSQ/Results/debug.txt", "w");
	if (dFile == NULL)
	{
		printf("Error writing debugfile\n");		
	}
	dFile1 = fopen( "C:/GRPLBSQ/Results/debug1.txt", "w");
	if (dFile1 == NULL)
	{
		printf("Error writing debug1file\n");		
	}
	dFile2 = fopen( "C:/GRPLBSQ/Results/debug2.txt", "w");
	if (dFile2 == NULL)
	{
		printf("Error writing debug1file\n");		
	}
	//

	//variable initialization
	int end=0;


	//Variables
	int i, j;

	double maxdist[kMAX];
	for(i=0; i<k; i++)
		maxdist[i] = MAXDOUBLE;
	


	//Find the MBR from provided rectangles
	Rectangle1 mbr;
	mbr.x1=R[0].x1;
	mbr.x2=R[0].x2;
	mbr.y1=R[0].y1;
	mbr.y2=R[0].y2;
	
	for(i=1; i<g_size; i++)
	{
		if(R[i].x1 < mbr.x1)	mbr.x1=R[i].x1;
		if(R[i].x2 > mbr.x2)	mbr.x2=R[i].x2;
		if(R[i].y1 < mbr.y1)	mbr.y1=R[i].y1;
		if(R[i].y2 > mbr.y2)	mbr.y2=R[i].y2;
	}
	//end

	// debug
	fprintf(dFile,"%lf\t%lf\t%lf\t%lf\n",mbr.x1,mbr.x2,mbr.y1,mbr.y2);
	//

	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next

	//debug

		fprintf(dFile,".............");
	//...
	while (son != -1 && end==0)
	{
		
		
		RTNode *rtn = new RTNode(this, son);
		//Experiment
		//io_access++;
		for (int i = 0; i < rtn -> num_entries; i++)
		{
			double edist1=0, edist2=0, tdist=0;
			edist1 = MINRECTDIST(mbr, rtn->entries[i].bounces);	
			if((edist1)>maxdist[k-1])	continue;

			//debug

			if(rtn -> level==0)
			{
			//fprintf(dFile,"\nI am in inserting heap mode\n");
			fprintf(dFile,"\nmaxdistk:%lf\tedist1:%lf\tedist2:%lf\tlevel:%d", maxdist[k-1],edist1,edist2,rtn -> level);
			//fprintf(dFile,"%lf\t%lf\t%lf\t%lf\n",rtn->entries[i].bounces[0],rtn->entries[i].bounces[1],rtn->entries[i].bounces[2],rtn->entries[i].bounces[3]);
			//}
			//...

			//debug
		
			fprintf(dFile,"\n\n%lf\t%lf\t%lf\t%lf",rtn->entries[i].bounces[0],rtn->entries[i].bounces[1],rtn->entries[i].bounces[2],rtn->entries[i].bounces[3]);
			}
			edist1 = MINRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				
				//debug
				fprintf(dFile, "\n%lf\t%lf\t%lf\t%lf", R[j-1].x1, R[j-1].x2,R[j-1].y1,R[j-1].y2);
				fprintf(dFile,"\nmaxdistk:%lf\ttdist:%lf\tedist1:%lf\tlevel:%d", maxdist[k-1],tdist,edist1,rtn -> level);
				//

				tdist = MINRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist1)
					edist1=tdist;

				
			}

			//debug
				fprintf(dFile, "\n%lf\t%lf\t%lf\t%lf", R[j-1].x1, R[j-1].x2,R[j-1].y1,R[j-1].y2);
				fprintf(dFile,"\nmaxdistk:%lf\ttdist:%lf\tedist1:%lf\tlevel:%d", maxdist[k-1],tdist,edist1,rtn -> level);
			//debug

			//if(rtn -> level==0)
			//{
			//fprintf(dFile,"\nI am in inserting heap mode\n");
			
			
			//}
			//...

			if(edist1>maxdist[k-1])	continue;
			
			edist2 = MAXRECTDIST(R[0], rtn->entries[i].bounces);
			for (j=1; j < g_size; j++)
			{
				//debug
				fprintf(dFile, "\n%lf\t%lf\t%lf\t%lf", R[j-1].x1, R[j-1].x2,R[j-1].y1,R[j-1].y2);
				fprintf(dFile,"\nmaxdistk:%lf\ttdist:%lf\tedist2:%lf\tlevel:%d", maxdist[k-1],tdist,edist2,rtn -> level);
				//
				
				tdist= MAXRECTDIST(R[j], rtn->entries[i].bounces);
				if(tdist>edist2)
					edist2=tdist;

				
			}

			//debug
				fprintf(dFile, "\n%lf\t%lf\t%lf\t%lf", R[j-1].x1, R[j-1].x2,R[j-1].y1,R[j-1].y2);
				fprintf(dFile,"\nmaxdistk:%lf\ttdist:%lf\tedist2:%lf\tlevel:%d", maxdist[k-1],tdist,edist2,rtn -> level);
			//debug
			//debug

			//if(rtn -> level==0)
			//{
			//fprintf(dFile,"\nI am in inserting heap mode\n");
			fprintf(dFile,"\nmaxdistk:%lf\tedist1:%lf\tedist2:%lf\tlevel:%d\n", maxdist[k-1],edist1,edist2,rtn -> level);
			//fprintf(dFile,"%lf\t%lf\t%lf\t%lf\n",rtn->entries[i].bounces[0],rtn->entries[i].bounces[1],rtn->entries[i].bounces[2],rtn->entries[i].bounces[3]);
			//}
			//...


			//update maxdistk		
			for(int l=0; l<k; l++)
			{
				if (edist2 < maxdist[l])
				{
					for(int j=k-1; j>l; j--)
					{
						maxdist[j]=maxdist[j-1];						
					}
					maxdist[l]=edist2;
					break;
				}
			}

			//debug
			for(int l=0; l<k; l++)
				fprintf(dFile2,"\n%lf", maxdist[l]);
			fprintf(dFile2,"\n\n\n");
			//...
			
			
			
			HeapEntry *he = new HeapEntry();
			he -> key = edist1;
			he -> key1 = edist2;
			he -> level = rtn -> level;
			he -> son1 = rtn->entries[i].son;
			he-> x1 = rtn->entries[i].bounces[0];
			//he-> x2 = rtn->entries[i].bounces[1];
			he-> y1 = rtn->entries[i].bounces[2];
			//he-> y2 = rtn->entries[i].bounces[3];
			heap -> insert(he);
			delete he;			

		}
	
		delete rtn;

		//debug

		fprintf(dFile,".............");
		//...
		
		//get next entry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{			
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else if (he->key>maxdist[k-1])
			{
				end = 1;

				//debug
				fprintf(dFile,"\nI am in end mode\n");
				fprintf(dFile,"\nmaxdistk:%lf\tedist:%lf\tlevel:%d\n", maxdist[k-1],he->key,he->level);
				//...
			}
			else
			{	

				//debug
				fprintf(dFile,"\nmaxdistk:%lf\tedist1:%lf\tedist2:%lf\tlevel:%d\n", maxdist[k-1],he->key,he->key1,he->level);
				fprintf(dFile1,"%lf\t%lf\n",he->x1,he->y1);
				//...

				if (he->level == 0) //p is an object 
				{
						//enter into result set
						if(*num_of_data==MAXDATALIMIT)
							//printf("\nGreater than 10000\n");
							error("Rect_kNNQ:DataPoint maximum Limit exceeded\n",TRUE);
						
						rslt[*num_of_data].x=he->x1;
						rslt[*num_of_data].y=he->y1;
						rslt[*num_of_data].dmin=he->key;
						rslt[*num_of_data].dmax=he->key1;

						
					

						

						//debug
						
						//fprintf(dFile,"%lf\t%lf\t%lf\t%lf\n",he->x1,he->y1,he->key,edist);
						
						//...

						*num_of_data=*num_of_data+1;
						//if(*num_of_data%1000==0) printf("\n%d",*num_of_data);
						
						//get next data  from heap
						again =true;						
				}
				else //not leaf node
				{
					son=he->son1;
				}
			}
		}
		
		delete he;
	}
	delete heap;

	//debug
	fclose(dFile);
	fclose(dFile1);
	//
}
