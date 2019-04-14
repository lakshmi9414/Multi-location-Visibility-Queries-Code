/*rtree.cpp
  this file implements the RTree class*/
//#undef min
//#undef max

#include <math.h>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "rtree.h"
#include "entry.h"
#include "rtnode.h"
#include "distance.h"
#include "../global.h"
#include "../blockfile/cache.h"
#include "../blockfile/blk_file.h"
#include "../linlist/linlist.h"
#include <array>
#include <list>
#include <math.h>
#include <queue>
#include <string.h>
#include <vector>
#include<iostream>
#include <map>
#include <utility>  
using namespace std;


//Added by Tanzima
//extern int io_access;
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
   
  
   
   inpname = "british_data.txt";
   fname = "mytree.tree"; 	
   
    Entry *d;
    FILE *fp;
    file = new BlockFile(fname, _b_length);
   	
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);
	delete [] header;

    
    //changing for sending numof obstacles	
    //dimension = _dimension;
	dimension = 2;
	int numobs=_dimension;	
	

    
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
      perror("Error");
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
	
      //while (!feof(fp)) //by lakshmi bcz some waste charcters it is taking
      while(record_count<=numobs)//should be changed
      {
		
		//cout << "record finserting -- " << record_count << endl;
		record_count ++;

		d = new Entry(dimension, NULL);

    	        fscanf(fp, "%d", &(d -> son));
        	
		//printf("ID=%d ",d->son);
		for (int i = 0; i < 2 * dimension; i ++)
		{
			fscanf(fp, " %f", &(d -> bounces[i]));
			
		    
		}
		for (int i = 0; i < 2 * dimension; i=i+2)
		{
			d->bounces[i+1] = d->bounces[i]+100;  // 100 is obstacle size here.. i set it like that
		}
		
    	       insert(d);
		
		
	       //cout << "inserted into r tree of record count " << record_count << endl;


		/*if (record_count % 100 == 0)
		{
			for (int i = 0; i < 79; i ++)  //clear a line
				printf("\b");

			printf("inserting object %d", record_count);
		}*/
      }
      //cout << "YAYYY " << record_count << endl;
      	
    }
	//TANZIMA
	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
	//cout << "BUILDING R-TREE IS OVER " << endl;
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

	//cout << "completed inserted function for one object " << endl; 
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
	//io_access++;
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
		//io_access++;
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
		//io_access++;
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





// method functions - lakshmi

bool ccw(Pointlocation A,Pointlocation B,Pointlocation C)
{
	return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x);
}

bool intersect(Pointlocation A,Pointlocation B,Pointlocation C,Pointlocation D)//if line segments intersect
{
        return ccw(A,C,D) != ccw(B,C,D) && ccw(A,B,C) != ccw(A,B,D);
}

double slope(double x1,double y1,double x2,double y2)//2 points of a line is given
{
       double s= (double)(y1-y2)/(x1-x2);
       return s;
}
bool intersect_ext(Pointlocation a,Pointlocation b, Pointlocation c, Pointlocation d)
{
	double m1,m2;
	if(a.x!=b.x)
	{
		m1 = slope(a.x,a.y,b.x,b.y);
	}

	if(c.x!=d.x)
	{
		m2 = slope(c.x,c.y,d.x,d.y);
	}

	if( (a.x==b.x) && (c.x==d.x) )
	{
		return false;//both ar perpendicular on x axis 
	}
	else if((a.x==b.x) && (c.x!=d.x) )
	{
		return true;
	}
	else if((a.x!=b.x) && (c.x==d.x) )
	{
		return true;
	}
	else if(m1 == m2)
	{
		return false;//parallal lines
	}
	else return true;
}



double angle_between_2_lines(double x1,double y1,double x2,double y2, double x11,double y11,double x21,double y21)
{
	double m1,m2;
	int flag1=0; int flag2=0;
    if(x1 != x2)
	{
		flag1=1;
		m1=slope(x1,y1,x2,y2);
		//cout << "SLOPE " << m1 << endl;
	}
	if(x11!=x21)
	{
		flag2=1;
		m2=slope(x11,y11,x21,y21);
		///cout << "SLOPE " << m2 << endl;
	}
	if((x1 == x2) && (x11 == x21))
	{
		return 0.0;
	}

	else if(x1 == x2)
	{	
		double ans=90 - (atan((float)fabs(m2)) * (180/M_PI)); //changed to fabs by lakshmi
		return ans;
	}
	else if(x11 == x21)
	{
		double ans = 90 - (atan((float)fabs(m1)) * (180/M_PI)); //changed to fabs by lakshmi
		//cout << "ANGLE .. " << ans << " " << M_PI << " " << (atan((float)m1) * (180/M_PI)) << endl;
		return ans;
	}
	else
	{
		double tan_theta=((double)(m1-m2)/(1+m1*m2));
		if(tan_theta<0.0)tan_theta=(tan_theta*-1);//acute angle only
		return (atan((float)tan_theta) * (180/M_PI));
	}
       
}


Pointlocation mid_of_line_segment(Pointlocation a,Pointlocation b)
{
	Pointlocation mid;
	mid.x=(a.x+b.x)/2.0;
	mid.y=(a.y+b.y)/2.0;
	return mid;
}

/*Pointlocation intersection_point(Pointlocation p1,Pointlocation p2,Pointlocation p3,Pointlocation p4) //farah function
{
	Pointlocation intersect;
	intersect.x = (((p1.x*p2.y- p2.x*p1.y)*(p3.x-p4.x)) - ((p1.x-p2.x)*(p3.x*p4.y - p3.y*p4.x))) / (((p1.x-p2.x)*(p3.y-p4.y)) - ((p1.y-p2.y)*(p3.x-p4.x)));

	intersect.y = (((p1.x*p2.y- p2.x*p1.y)*(p3.y-p4.y)) - ((p1.y-p2.y)*(p3.x*p4.y - p3.y*p4.x))) / (((p1.x-p2.x)*(p3.y-p4.y)) - ((p1.y-p2.y)*(p3.x-p4.x)));
	return intersect;
}*/
Pointlocation intersection_point(Pointlocation p1,Pointlocation p2,Pointlocation p3,Pointlocation p4)
{
	//p1,p2 is line segment
	double x1=p1.x,y1=p1.y,x2=p2.x,y2=p2.y,x3=p3.x,y3=p3.y,x4=p4.x,y4=p4.y;
	
	//cout << x1 << " " << y1 << endl;
	//cout << x2 << " " << y2 << endl;
	//cout << x3 << " " << y3 << endl;
	//cout << x4 << " " << y4 << endl;
	Pointlocation ans;
	double ml,ms,il,is;
	//ml,il;
	ml= (double)(y3-y4)/(x3-x4);
	il= y3- (ml*x3);

	
	if(x2==x1)
	{
		ans.x=x1;
		ans.y = (ml*x1)+il;
	}
	else if(y1==y2)
	{
		ans.y=y1;
		ans.x=(double)(ans.y - il)/(ml);
	}	
	else
	{	
	ms=(double)(y2-y1)/(x2-x1);
	is= y2- (ms*x2);
	
	ans.x=  is-il/ml-ms;
	ans.y= (ms*ans.x)+il;
	}
	return ans;
}

double reduced_visibility_factor4angle(double angle)
{
	return angle/90.0;
}
double Lakarea(double x1, double y1, double x2, double y2, double x3, double y3)
{
   return fabs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
}
 

bool LakIsPointInsideTriangle(Pointlocation a,Pointlocation b,Pointlocation c,Pointlocation x)
{   
	double x1=a.x, y1=a.y, x2=b.x,  y2=b.y,  x3=c.x,  y3=c.y,  x4=x.x,  y4=x.y;
   /* Calculate area of triangle ABC */
  double A = Lakarea (x1, y1, x2, y2, x3, y3);
 
   /* Calculate area of triangle PBC */  
   double A1 = Lakarea (x4, y4, x2, y2, x3, y3);
 
   /* Calculate area of triangle PAC */  
   double A2 = Lakarea (x1, y1, x4, y4, x3, y3);
 
   /* Calculate area of triangle PAB */   
   double A3 = Lakarea (x1, y1, x2, y2, x4, y4);
   
   //cout << "AREAAA ..."  << A1 << " " << A2 << " " << A3 << " " <<  A << endl;	
   /* Check if sum of A1, A2 and A3 is same as A */
   return (A == A1 + A2 + A3);
}




bool IsPointInsideTriangle(Pointlocation a,Pointlocation b,Pointlocation c,Pointlocation x)
{	
	
	double angle1=angle_between_2_lines(x.x,x.y,a.x,a.y,x.x,x.y,c.x,c.y);
	double angle2=angle_between_2_lines(x.x,x.y,a.x,a.y,x.x,x.y,b.x,b.y);
	double angle3=angle_between_2_lines(x.x,x.y,b.x,b.y,x.x,x.y,c.x,c.y);
	
	//cout << "indvi angle ... " << angle1 << " " << angle2 << " " << angle3 << endl;
	//cout << "Angle .. " << angle1+angle2+angle3 << endl;
	if(fabs(angle1+angle2+angle3-180)<=0.0001) // lakshmi changed
	{
		return true;
	}
	return false;
}

double Distance(Pointlocation p1, Pointlocation p2)
{
	return sqrt( ((p1.x-p2.x)*(p1.x-p2.x)) + ((p1.y-p2.y)*(p1.y-p2.y)) );
}
/////////////////////////////
vector<Portion>  Calc_portions_on_side(pair<Pointlocation,Pointlocation> cur_side_ends,int side_id,int no_of_portions)
{

         
	double num_portions=no_of_portions,t_num_portions;  // no.of portions per side , for now number of portions per side is fixed.
	int cnt=side_id*num_portions;	
	vector<Portion> cur_side_portions;

	Pointlocation side_start,side_end; // instead of pa,pb
	side_start=cur_side_ends.first;
	side_end=cur_side_ends.second;
	

	double len_side = Distance(side_start,side_end);
	double len_portion = len_side/num_portions; /// rewrite for checking double divison.
	//len_portion = floorf(len_portion * 100) / 100; // rewrite
		 
	 Portion cur_portion; 
	 Pointlocation old_p;
	 t_num_portions= num_portions;
	 
	 //cout << "num_portions= " << t_num_portions << endl;
	// std::cout << std::fixed;
	 //cout << "side ka start, end points= " << side_start.x << " " << side_start.y << "," << side_end.x << " " << side_end.y << endl;

	 if(side_start.x==side_end.x and side_start.y>side_end.y)
	 {	
	       // cout << "vertical line down " << endl;
	 	old_p=side_start;
	 	while(t_num_portions)
	 	{
	 		t_num_portions--;
	 				
	 	 Pointlocation new_p;
	 	 new_p.x=side_start.x;
	 	 new_p.y=old_p.y - len_portion;
	 	 if(new_p.y < side_end.y)  //check this again
	 	 	new_p.y = side_end.y;
		
		cur_portion.por_id=cnt ;// as of now;
		cur_portion.p1_overT=old_p;
		cur_portion.p2_overT=new_p;
		cur_side_portions.push_back(cur_portion);
		//std::cout << std::fixed;
		//cout << "portion: " << cur_portion.por_id << "," << cur_portion.p1_overT.x << " " << cur_portion.p1_overT.y << "," << cur_portion.p2_overT.x << " " << cur_portion.p2_overT.y << endl;
		cnt++;
		old_p=new_p;
	 	}
	 	

	 }
          else if(side_start.x==side_end.x and side_start.y<side_end.y)
	 {	
	       // cout << "vertical line down " << endl;
	 	old_p=side_start;
	 
	 	while(t_num_portions)
	 	{
	 		t_num_portions--;
	 				
	 	 Pointlocation new_p;
	 	 new_p.x=side_start.x;
	 	 new_p.y=old_p.y + len_portion;
	 	 if(new_p.y > side_end.y)
	 	 	new_p.y = side_end.y;
		
		cur_portion.por_id=cnt ;// as of now;
		cur_portion.p1_overT=old_p;
		cur_portion.p2_overT=new_p;
		cur_side_portions.push_back(cur_portion);
		//std::cout << std::fixed;
		//cout << "portion: " << cur_portion.por_id << "," << cur_portion.p1_overT.x << " " << cur_portion.p1_overT.y << "," << cur_portion.p2_overT.x << " " << cur_portion.p2_overT.y << endl;
		cnt++;
	 	 
		old_p=new_p;
	 	}
	 	

	 }
	 
		
	 else if(side_start.y==side_end.y and side_start.x>side_end.x)
	 {	
	       // cout << "horizontal line left " << endl;
		t_num_portions= num_portions;
      		old_p=side_start;
	 	while(t_num_portions)
	 	{
	 	 t_num_portions--;	
	 	 Pointlocation new_p;
	 	 new_p.y=side_start.y;
	 	 new_p.x=side_start.x;
	 	 
	 	 new_p.x = old_p.x - len_portion;
	 	 if(new_p.x<side_end.x) //as of now
	 	 	new_p.x=side_end.x;
		cur_portion.por_id=cnt ;// as of now;
		cur_portion.p1_overT=old_p;
		cur_portion.p2_overT=new_p;
		cur_side_portions.push_back(cur_portion);
		//std::cout << std::fixed;
		//cout << "portion: " << cur_portion.por_id << "," << cur_portion.p1_overT.x << " " << cur_portion.p1_overT.y << "," << cur_portion.p2_overT.x << " " << cur_portion.p2_overT.y << endl;
		cnt++;
		old_p=new_p;
	 	}
	 	
	 }
	else if(side_start.y==side_end.y and side_start.x<side_end.x)
	 {	
	        //cout << "horizontal line right " << endl;
		t_num_portions= num_portions;
      		old_p=side_start;
	 	while(t_num_portions)
	 	{
	 	 t_num_portions--;	
	 	 Pointlocation new_p;
	 	 new_p.y=side_start.y;
	 	 new_p.x=side_start.x;
	 	 
	 	 new_p.x = old_p.x + len_portion;
	 	 if(new_p.x>side_end.x) //as of now
	 	 	new_p.x=side_end.x;
		cur_portion.por_id=cnt ;// as of now;
		cur_portion.p1_overT=old_p;
		cur_portion.p2_overT=new_p;
		cur_side_portions.push_back(cur_portion);
		//std::cout << std::fixed;
		//cout << "portion: " << cur_portion.por_id << "," << cur_portion.p1_overT.x << " " << cur_portion.p1_overT.y << "," << cur_portion.p2_overT.x << " " << cur_portion.p2_overT.y << endl;
		cnt++;
		old_p=new_p;
	 	}
	 	
	 }
			
   
    return cur_side_portions;

}
//////////////////////////////


map< int ,vector<Portion> > RTree :: Calc_portions_on_target(Rectangle Target, map< int ,vector<Portion> > portion_data,int no_of_portions)
{
	
	
	Pointlocation p1=Target.upper_left;
	Pointlocation p2(Target.upper_left.x,Target.lower_right.y);
	Pointlocation p3(Target.lower_right.x,Target.upper_left.y);
	Pointlocation p4 = Target.lower_right;
	
	
	
	map< int, pair<Pointlocation,Pointlocation> > side_ends;
	pair<Pointlocation,Pointlocation> side_ends_1,side_ends_2,side_ends_3,side_ends_4;
	
	side_ends_1=make_pair(p1,p2);
	side_ends[1]=side_ends_1;
	
	side_ends_2=make_pair(p2,p4);
	side_ends[2]=side_ends_2;
	
	side_ends_3=make_pair(p4,p3);
	side_ends[3]=side_ends_3;
	
	side_ends_4=make_pair(p3,p1);
	side_ends[4]=side_ends_4;
	

	vector<Portion> cur_side_portions;
	portion_data[1] = cur_side_portions;
	

	
	
	for(int i=0;i<side_ends.size();i++)
	{
	cur_side_portions=Calc_portions_on_side(side_ends[i+1],i,no_of_portions);
	/*for(int k=0;k<cur_side_portions.size();k++)
	{
	cout << cur_side_portions[k].por_id << " " << cur_side_portions[k].p1_overT.x << "," << cur_side_portions[k].p1_overT.y << " " << cur_side_portions[k].p2_overT.x << "," << cur_side_portions[k].p2_overT.y << endl;
	}*/
        portion_data[i+1]=cur_side_portions;
	}
	return portion_data;	
}
////////////////////////////////////


//considering no obstacle - to be replaced my own
/*double QueryPoint :: init_visibility(Rectangle Target,map<int, vector<Portion> > portion_data)
{
	Pointlocation p1=Target.upper_left;
	Pointlocation p2(Target.upper_left.x,Target.lower_right.y);
	Pointlocation p3(Target.lower_right.x,Target.upper_left.y);
	Pointlocation p4 = Target.lower_right;

	bool p1_seen,p2_seen,p3_seen,p4_seen;
	p1_seen=p2_seen=p3_seen=p4_seen=false;

	total_visibility = 0.0;
	//initially kon kon side dekhte pay seta check korte hobe

	//check for line SEGMENT intersect 
	if(intersect(position, p1 , p2,p4)== false &&  intersect(position, p1 , p3,p4)== false) p1_seen=true;
	if(intersect(position, p4 , p1,p2)== false &&  intersect(position, p4 , p1,p3)== false) p4_seen=true;

	if(intersect(position, p2 , p3,p4)== false &&  intersect(position, p2 , p1,p3)== false) p2_seen=true;
	if(intersect(position, p3 , p1,p2)== false &&  intersect(position, p3 , p2,p4)== false) p3_seen=true;
	
	Pointlocation mid1, point1, point2;

	double angle_with_rect=0.0;
	double visibility_reduce_factor=1.0;
	double partial_visibility=0.0;

	if((p1_seen==true && p2_seen==true) || (p3_seen==true && p4_seen==true))
	{
		if(p1_seen==true && p2_seen==true){
			point1=p1; point2=p2;
		}
		else{		// if p1,p2 is seen, p3,p4 cannot be seen
			point1=p3; point2=p4;
		}
		
		mid1=mid_of_line_segment(point1,point2); //rectangle target, so p1,p2 length == p3,p4 length 
		angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, point1.x,point1.y,point2.x,point2.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
		visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

		partial_visibility = Distance(point1, point2) * visibility_reduce_factor;
		total_visibility +=  partial_visibility;

		//construct the visible region list
		Pointlocation a;
		struct VisibleRegionOverT segment;
		segment.p1_overT=point1;
		segment.p2_overT=point2;
		segment.partial_visibility = partial_visibility;
		VisibleRegion.push_back(segment);
	}

	
	if((p1_seen==true && p3_seen==true) || (p2_seen==true && p4_seen==true))
	{
		if(p1_seen==true && p3_seen==true){
			point1=p1; point2=p3;
		}
		else{		// if p1,p3 is seen, p2,p4 cannot be seen
			point1=p2; point2=p4;
		}
		
		mid1=mid_of_line_segment(point1,point2); //rectangle target, so p1,p2 length == p3,p4 length 
		angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, point1.x,point1.y,point2.x,point2.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
		visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

		partial_visibility = Distance(point1, point2) * visibility_reduce_factor;
		total_visibility +=  partial_visibility;

		//construct the visible region list
		Pointlocation a;
		struct VisibleRegionOverT segment;
		segment.p1_overT=point1;
		segment.p2_overT=point2;
		segment.partial_visibility = partial_visibility;
		VisibleRegion.push_back(segment);
	}

	return total_visibility;
}*/
//*************************
//considering no obstacle - my function
double QueryPoint :: init_visibility(Rectangle Target, map<int, vector<Portion> > portion_data,int dis_threshold,double vis_threshold)
{
        
        //cout << endl;
       // std::cout << std::fixed;
        //cout << "query points = " << position.x << " " << position.y << endl;
	//cout << "CENTROID .. " << dis_to_tar << " " << dis_threshold << endl;
	total_visibility = 0.0;
	//if((int)dis_to_tar>dis_threshold+500 or (int)dis_to_tar<dis_threshold-500)
	if((int)dis_to_tar>dis_threshold)
	return total_visibility;
	//side 1 - p1,p2,side -2 - p2,p4,side -3 - p4,p3, side -4 - p3,p1 (anti clkwise)
	Pointlocation p1=Target.upper_left;
	Pointlocation p2(Target.upper_left.x,Target.lower_right.y);
	Pointlocation p3(Target.lower_right.x,Target.upper_left.y);
	Pointlocation p4 = Target.lower_right;
	int side_seen=-1,num_portions=portion_data[1].size();
	vector<Portion> cur_side_portions;

	bool p1_seen,p2_seen,p3_seen,p4_seen;
	p1_seen=p2_seen=p3_seen=p4_seen=false;

	
	//initially kon kon side dekhte pay seta check korte hobe

	//check for line SEGMENT intersect 
	if(intersect(position, p1 , p2,p4)== false &&  intersect(position, p1 , p3,p4)== false) p1_seen=true;
	if(intersect(position, p4 , p1,p2)== false &&  intersect(position, p4 , p1,p3)== false) p4_seen=true;

	if(intersect(position, p2 , p3,p4)== false &&  intersect(position, p2 , p1,p3)== false) p2_seen=true;
	if(intersect(position, p3 , p1,p2)== false &&  intersect(position, p3 , p2,p4)== false) p3_seen=true;
	
	Pointlocation mid1, point1, point2;

	double angle_with_rect=0.0;
	double visibility_reduce_factor=1.0;
	double partial_visibility=0.0,seg_visibility=0.0;

	if((p1_seen==true && p2_seen==true) || (p3_seen==true && p4_seen==true))
	{
		if(p1_seen==true && p2_seen==true){
			point1=p1; point2=p2;
			side_seen=0;
		}
		else
		{		
			point1=p3; point2=p4;
			side_seen=2;
		}
		cur_side_portions=portion_data[side_seen+1];
		
		for(int i=0;i<num_portions;i++)
		{
		
			
		point1=cur_side_portions[i].p1_overT;
		point2=cur_side_portions[i].p2_overT;
		
		mid1=mid_of_line_segment(point1,point2);  
		angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, point1.x,point1.y,point2.x,point2.y);
		visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);
		seg_visibility=Distance(point1, point2) * visibility_reduce_factor;
		
		//if(visibility_reduce_factor > dis_to_tar/dis_threshold)
		if(visibility_reduce_factor>vis_threshold)	
		{
                //cout << "***1 " << visibility_reduce_factor << endl;
		//construct the visible region list
		struct VisibleRegionOverT segment;
		segment.vr_id = (side_seen*num_portions)+i;
		//cout << "IDDDD " << i << " " << segment.vr_id << endl;
		segment.p1_overT=point1;
		segment.p2_overT=point2;
		segment.partial_visibility = seg_visibility;
		segment.actual_dis = Distance(point1, point2);
		VisibleRegion.push_back(segment);
		total_visibility +=  seg_visibility;
		}
		}

	//cout << "sides seen aaa " << side_seen << endl;
	}
	

	
	if((p1_seen==true && p3_seen==true) || (p2_seen==true && p4_seen==true))
	{
		if(p1_seen==true && p3_seen==true){
			point1=p1; point2=p3;
			side_seen=3;
		}
		else
		{
			point1=p2; point2=p4;
			side_seen=1;
		}
		cur_side_portions=portion_data[side_seen+1];
		for(int i=0;i<num_portions;i++)
		{	
		point1=cur_side_portions[i].p1_overT;
		point2=cur_side_portions[i].p2_overT;
		//std::cout << std::fixed;
		//cout << "portion points = " << point1.x << " " << point1.y << " " << point2.x << " " << point2.y << endl;
		//cout << point1.x << "," << point1.y << " " << point2.x << "," << point2.y << endl;
		mid1=mid_of_line_segment(point1,point2); 
		angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, point1.x,point1.y,point2.x,point2.y);
		visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);
		seg_visibility=Distance(point1, point2) * visibility_reduce_factor;
		//total_visibility +=  seg_visibility; // chaging this line from here to inside if condition
		//construct the visible region list
		//cout << "** " << visibility_reduce_factor << endl;
		//std::cout << std::fixed;
		//cout << "reduce factors =  "   << visibility_reduce_factor << endl;
		if(visibility_reduce_factor>vis_threshold)
		{
                 //cout << "***2 " << visibility_reduce_factor << endl;
		struct VisibleRegionOverT segment;
		segment.vr_id = (side_seen*num_portions)+i;
		segment.p1_overT=point1;
		segment.p2_overT=point2;
		segment.partial_visibility = seg_visibility;
		segment.actual_dis = Distance(point1, point2);
		VisibleRegion.push_back(segment);
		total_visibility +=  seg_visibility;
		}
		}
	//cout << "sides seen bbb" << side_seen << endl;
	}
	
	//std::cout << std::fixed;
        //cout << "VR = " << position.x << " " << position.y << " - " << VisibleRegion.size() << " " << total_visibility << endl; 	
	return total_visibility;
}


/////////////////////////////



//visible region is the list of triangles where a triangle is formed of the query point and 2 points over the target
//check if obstacle is intersecting visibile regions of a query point.
//need not change as of now
bool QueryPoint::IsInVisibleRegion(Rectangle rect)
{
	//check if rect intersects with any of the items in VisibleRegion list
	Pointlocation p1=rect.upper_left;
	Pointlocation p2;
	p2.x=rect.upper_left.x;
	p2.y=rect.lower_right.y;
	Pointlocation p3;
	p3.x=rect.lower_right.x;
	p3.y=rect.upper_left.y;
	Pointlocation p4= rect.lower_right;

	//cout << "checking pehla .... " << p1.x << " " << p1.y << endl;
	int itr;

	for ( itr = 0; itr < VisibleRegion.size(); itr++ )
	{
		//if rect is in the region (fully or partially) of triangle formed by querypoint, VisibleRegion.p1,VisibleRegion.p2 
		//return true
		//cout << "vis reg traiangle... " <<  VisibleRegion[itr].p1_overT.x << " " << VisibleRegion[itr].p1_overT.y << " " << VisibleRegion[itr].p2_overT.x << " " << VisibleRegion[itr].p2_overT.y << endl;
		
		if(intersect(position, VisibleRegion[itr].p1_overT, p1,p2)) return true;
		else if(intersect(position, VisibleRegion[itr].p1_overT, p1,p3)) return true;
		else if(intersect(position, VisibleRegion[itr].p1_overT, p2,p4)) return true;
		else if(intersect(position, VisibleRegion[itr].p1_overT, p3,p4)) return true;

		else if(intersect(position, VisibleRegion[itr].p2_overT, p1,p2)) return true;
		else if(intersect(position, VisibleRegion[itr].p2_overT, p1,p3)) return true;
		else if(intersect(position, VisibleRegion[itr].p2_overT, p2,p4)) return true;
		else if(intersect(position, VisibleRegion[itr].p2_overT, p3,p4)) return true;

		else if(intersect(VisibleRegion[itr].p2_overT, VisibleRegion[itr].p1_overT, p1,p2)) return true;
		else if(intersect(VisibleRegion[itr].p2_overT, VisibleRegion[itr].p1_overT, p1,p3)) return true;
		else if(intersect(VisibleRegion[itr].p2_overT, VisibleRegion[itr].p1_overT, p2,p4)) return true;
		else if(intersect(VisibleRegion[itr].p2_overT, VisibleRegion[itr].p1_overT, p3,p4)) return true;

		//is rect fully inside the visible region
		//else if 
		else if(LakIsPointInsideTriangle(position,VisibleRegion[itr].p1_overT,VisibleRegion[itr].p2_overT, p1)==true)//rect er jeno ekta point (p1 here) check korlei hobe, jehetu ager else if gula check kora ache
		{
			
			return true;
		}
		//else return false; // may b this shouldnot be there

	}
	
     return false;
                                        
}


// dont know wat it does for now
bool if_intersects_target(Pointlocation position, Pointlocation obstacle_p, Rectangle target, Pointlocation& intersectPoint,Pointlocation p1,Pointlocation p2)
{
	
	Pointlocation target_p1 = target.upper_left;
	Pointlocation target_p2;
	target_p2.x=target.upper_left.x;
	target_p2.y=target.lower_right.y;
	Pointlocation target_p3;
	target_p3.x=target.lower_right.x;
	target_p3.y=target.upper_left.y;
	Pointlocation target_p4 = target.lower_right;

	float lowx,lowy,highx,highy;
	lowx =min(p1.x, p2.x);
	lowy = min(p1.y,p2.y);
	highx = max(p1.x,p2.x);
	highy = max(p1.y,p2.y);

	Pointlocation Tpoints[2]={p1,p2};	
	//Pointlocation Tpoints[8] = {target_p1,target_p2, target_p1, target_p3, target_p3,target_p4, target_p2,target_p4};
	for(int i=0;i<2;i+=2)
	{
		if(intersect_ext(position,obstacle_p, Tpoints[i],Tpoints[i+1])==true)
		{	
			//cout << "CHECK INSIDE " << "YAYYY " << endl;
			
			intersectPoint = intersection_point(Tpoints[i],Tpoints[i+1],position,obstacle_p); // check once gaina
			//cout << "ulllll 2222 akks " << Tpoints[i].x << " " << Tpoints[i].y << " " << Tpoints[i+1].x << " " << Tpoints[i+1].y << endl; 
			//cout << "ultiii 22 laks.. " << obstacle_p.x << " " << obstacle_p.y  << " " << position.x << " " << position.y << endl;
			//cout << "ulttt answer.. " << intersectPoint.x << " " << intersectPoint.y << endl;
			if(intersectPoint.x<=highx && intersectPoint.x>=lowx && intersectPoint.y<=highy && intersectPoint.y>=lowy)//intersect point is in that region
			{

				return true;
			}
		
		}
	}
	return false;
	
	
}

// dont know wat it oes for now
bool change_vRegion_one_side(Pointlocation position, Pointlocation obstacle_p1, Rectangle target, 
	Pointlocation unchanged_side_val, Pointlocation changed_side_val,Pointlocation& intersectPoint, 
	VisibleRegionOverT itr)
{
	bool check = if_intersects_target(position, obstacle_p1 , target, intersectPoint,unchanged_side_val,changed_side_val);
	//cout << "ONE SIDE... " << check << endl;
	//cout << "ultiii , lak .. " << intersectPoint.x << " " << intersectPoint.y << endl;
	if(check == true)
	{
		float a=Distance(unchanged_side_val, (intersectPoint));
		float b=Distance(unchanged_side_val,changed_side_val); 
		if(a<= b)
		{
						
			VisibleRegionOverT v;
			v.p1_overT = unchanged_side_val;
			v.p2_overT = (intersectPoint);
			
			return true;
		}
	
	}
	return false;
}


//dont know wat it does for now
bool isPointOf_VR_changed(Pointlocation p1_overT, Pointlocation position, Pointlocation obstacle_p1,
	Pointlocation obstacle_p2,Pointlocation obstacle_p3,Pointlocation obstacle_p4)
	
{
	//cout << "one point checkkk .. " << endl; 
	Pointlocation Opoints[8] = {obstacle_p1,obstacle_p2,obstacle_p1,obstacle_p3,obstacle_p2,obstacle_p4,obstacle_p3,obstacle_p4};
	
	for(int i=0;i<8;i+=2)
	{
		bool z = intersect(position, p1_overT, Opoints[i],Opoints[i+1]);
		//cout << p1_overT.x << " " << p1_overT.y << " " << Opoints[i].x << " " << Opoints[i].y << " " << Opoints[i+1].x << " " << Opoints[i+1].y << endl;
		//cout << "ANS .. " << z << endl;
		if(z)
		{
			return true;
		}
	}
	
	return false;
}



//update visibility for one obstacle - lakshmi
//VisibleRegion[itr].partial_visibility - when ever this is getting updated - it means if still not crossing threshold.
//if after considering ful obstacle, if croeses the thresold, delete this segment form visibleregion.
//no need of total visibility updation, simply delete those lines. 
void QueryPoint :: update_visibliliyRegion(Rectangle obstacle, Rectangle target,double vis_threshold)
{
	
	
	Pointlocation obstacle_p1 = obstacle.upper_left;
	Pointlocation obstacle_p2(obstacle.upper_left.x,obstacle.lower_right.y);
	Pointlocation obstacle_p3(obstacle.lower_right.x,obstacle.upper_left.y);
	Pointlocation obstacle_p4 = obstacle.lower_right;
	//cout << "START UPDATING............... " << "query point - " <<  position.x << " " << position.y << " " << endl;
	//cout << "initial visibile segments size - " << VisibleRegion.size() << " " << total_visibility << ",obstacle - " << obstacle_p1.x << " " << obstacle_p1.y << endl;

	Pointlocation target_p1 = target.upper_left;
	Pointlocation target_p2(target.upper_left.x,target.lower_right.y);
	Pointlocation target_p3(target.lower_right.x,target.upper_left.y);
	Pointlocation target_p4 = target.lower_right;

	bool p1_seen,p2_seen,p3_seen,p4_seen,p1_overT_changed,p2_overT_changed;
	p1_seen=p2_seen=p3_seen=p4_seen=p1_overT_changed=p2_overT_changed=false;

	int itr;
	Pointlocation temp_p;
	Pointlocation& intersectPoint=temp_p;
	
	Pointlocation& intersectPoint1 =temp_p;Pointlocation& intersectPoint2=temp_p;
	

	//construct temp with the changed v. regions. finally make visibleRegion=temp
	vector<VisibleRegionOverT> temp;
	itr=0;
	double temp_visibility=0.0,por_actual_dis;
	while(itr<VisibleRegion.size())
	{ 
		//cout << "POINTS on target " << VisibleRegion[itr].p1_overT.x << ","  << VisibleRegion[itr].p1_overT.y << " " << VisibleRegion[itr].p2_overT.x << "," <<  VisibleRegion[itr].p2_overT.y << endl;
		//cout << "total visibility before VR " << VisibleRegion[itr].partial_visibility << " " <<  total_visibility << endl;
		//cout << "ITR " << itr << " " << VisibleRegion.size() << endl;
		
		
		//actual_dis = Distance(VisibleRegion[itr].p1_overT, VisibleRegion[itr].p2_overT);
		por_actual_dis = VisibleRegion[itr].actual_dis;
		//cout << "DONT CHANGE " << VisibleRegion[itr].vr_id << " " << VisibleRegion[itr].actual_dis << endl;
		temp_visibility=0.0;
		bool entirely_inside_region= false;
		bool one_side_changed=false;

		p1_overT_changed = isPointOf_VR_changed((VisibleRegion[itr]).p1_overT,position, obstacle_p1,obstacle_p2, obstacle_p3,obstacle_p4);
		p2_overT_changed = isPointOf_VR_changed((VisibleRegion[itr]).p2_overT,position, obstacle_p1,obstacle_p2, obstacle_p3,obstacle_p4);		
				
		
		Pointlocation unchanged_side_val;
		Pointlocation changed_side_val;

		//condition for case 2 --- check if obstacle in inside triangle visible region
		if(p2_overT_changed == false && p1_overT_changed ==false && LakIsPointInsideTriangle(position,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT, obstacle_p1)==true)//rect er jeno ekta point (p1 here) check korlei hobe, jehetu ager else if gula check kora ache
		{
			
			entirely_inside_region= true;
		}
		
		//condition for case 3
		if(p1_overT_changed == true && p2_overT_changed == false)
		{
			one_side_changed=true;
			unchanged_side_val = (VisibleRegion[itr]).p2_overT;
			changed_side_val = (VisibleRegion[itr]).p1_overT;
		}
		//condition for case 3
		if(p2_overT_changed == true && p1_overT_changed == false)
		{
			one_side_changed= true;
			unchanged_side_val = (VisibleRegion[itr]).p1_overT;
			changed_side_val = (VisibleRegion[itr]).p2_overT;
		}
		//cout << "conditions..**** " << p1_overT_changed << " " << p2_overT_changed << " " << entirely_inside_region << " " << one_side_changed << endl;
		
		////**************************///// case 1
		if(p1_overT_changed==true && p2_overT_changed==true)//totally blocked this visibilitysegment
		{
			
			//total_visibility -= VisibleRegion[itr].actual_dis;
			total_visibility -= VisibleRegion[itr].partial_visibility;
			//cout << "temp vis after check " << total_visibility << endl;
			VisibleRegion.erase(VisibleRegion.begin()+itr);
			
			//cout << "CASE 1-output " << VisibleRegion[itr].partial_visibility << " " <<  total_visibility << endl;
		}


		///////*****************************///// case 2		
		else if( entirely_inside_region == true)
		{//case 2 start bracket
			//cout << "AMULYA - " << "init vis " << total_visibility << endl;
			//cout << "1111WWWW " << total_visibility << " ";
			bool intersectPoint1_found,intersectPoint2_found;
			intersectPoint1_found=intersectPoint2_found=false;
			
			//split the region
			if(intersect(position, obstacle_p1, obstacle_p3, obstacle_p4)==false &&
			intersect(position, obstacle_p1, obstacle_p2, obstacle_p4)==false)//doesn't intersect with own arm, for point obstacle_p1
			{
				//cout << "amull - check -1 " << endl;
				bool check;
				if(intersectPoint1_found==false)
				{
					check = if_intersects_target(position, obstacle_p1 , target, intersectPoint1,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					if(check == true)
					{
						intersectPoint1_found=true;
					}
				}

				if(intersectPoint2_found==false)
				{
					check = if_intersects_target(position, obstacle_p1 , target, intersectPoint2,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					if(check == true)
					{
						intersectPoint2_found=true;
					}
				}	
				
	

			}

			if(intersect(position, obstacle_p2, obstacle_p1, obstacle_p3)==false &&
				intersect(position, obstacle_p2, obstacle_p3, obstacle_p4)==false)//doesn't intersect with own arm, for point obstacle_p2
			{
				//cout << "amull - check -2 " << endl;
				bool check;
				if(intersectPoint1_found==false)
				{
					check = if_intersects_target(position, obstacle_p2 , target, intersectPoint1,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					//cout << "amul - check -2.inter1 "  << check << endl; 
					if(check == true)
					{
						intersectPoint1_found=true;
					}
				//cout << "OMGGG before if " << intersectPoint2_found << endl;
				}
				
				if(intersectPoint2_found==false)
				{
					//cout << "OMGGG after if " << intersectPoint2_found << endl;
					check = if_intersects_target(position, obstacle_p2 , target, intersectPoint2,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					//cout << "amul - check -2.inter2 "  << check << endl;
					if(check == true)
					{
						intersectPoint2_found=true;
					}
				}
				//cout << "amul yayy -2 " << intersectPoint1_found << " " << intersectPoint2_found << endl;

				
			}

			if( intersect(position, obstacle_p3, obstacle_p2, obstacle_p4)==false &&
				intersect(position, obstacle_p3, obstacle_p1, obstacle_p2)==false)//doesn't intersect with own arm, for point obstacle_p3
			{
				//cout << "amull - check -3 " << endl;				
				bool check;
				if(intersectPoint1_found==false)
				{
					check = if_intersects_target(position, obstacle_p3 , target, intersectPoint1,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					if(check == true)
					{
						intersectPoint1_found=true;
					}
				}

				if(intersectPoint2_found==false)
				{
					check = if_intersects_target(position, obstacle_p3 , target, intersectPoint2,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					if(check == true)
					{
						intersectPoint2_found=true;
					}
				}


			}

			if(intersect(position, obstacle_p4, obstacle_p1, obstacle_p2)==false &&
				intersect(position, obstacle_p4, obstacle_p1, obstacle_p3)==false)//doesn't intersect with own arm, for point obstacle_p4
			{
				//cout << "amull - check -4 " << endl;
				bool check;
				if(intersectPoint1_found==false)
				{
					check = if_intersects_target(position, obstacle_p4 , target, intersectPoint1,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					//cout << "amul - check -4.inter1 "  << check << endl;
					if(check == true)
					{
						intersectPoint1_found=true;
					}
				}

				if(intersectPoint2_found==false)
				{
					check = if_intersects_target(position, obstacle_p4 , target, intersectPoint2,(VisibleRegion[itr]).p1_overT,(VisibleRegion[itr]).p2_overT );
					//cout << "amul - check -4.inter2 "  << check << endl;
					if(check == true)
					{
						intersectPoint2_found=true;
					}
				}
				//cout << "amul yayy -4 " << intersectPoint1_found << " " << intersectPoint2_found << endl;
			}

			VisibleRegionOverT v;
			if(Distance(VisibleRegion[itr].p1_overT, intersectPoint1) < Distance(VisibleRegion[itr].p1_overT, intersectPoint2))
			{
				//cout << "AMUL - aaaa" << endl;
				//edit visibility and construct new segments
				total_visibility -= VisibleRegion[itr].partial_visibility;
				Pointlocation mid1;double angle_with_rect;double visibility_reduce_factor;
				if(intersectPoint1_found==true)
				{
					v.p1_overT=VisibleRegion[itr].p1_overT;
					v.p2_overT=intersectPoint1;

					mid1=mid_of_line_segment(v.p1_overT,v.p2_overT); //rectangle target, so p1,p2 length == p3,p4 length 
					angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, v.p1_overT.x,v.p1_overT.y,v.p2_overT.x,v.p2_overT.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
					visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					v.partial_visibility = Distance(VisibleRegion[itr].p1_overT, intersectPoint1)*visibility_reduce_factor;
					total_visibility +=v.partial_visibility;

					temp_visibility +=v.partial_visibility; // lakshmi changing
					//cout << "AMUL -1 " << v.partial_visibility << " " << temp_visibility << endl;
					//temp.push_back(v);
				}
				if(intersectPoint2_found==true)
				{
					v.p1_overT= VisibleRegion[itr].p2_overT;
					v.p2_overT=intersectPoint2;

					mid1=mid_of_line_segment(v.p1_overT,v.p2_overT); 
					angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, v.p1_overT.x,v.p1_overT.y,v.p2_overT.x,v.p2_overT.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
					visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					v.partial_visibility = Distance(VisibleRegion[itr].p2_overT, intersectPoint2)*visibility_reduce_factor;
 
					total_visibility +=v.partial_visibility;
					temp_visibility +=v.partial_visibility;
					//cout << "AMUL -2 " << v.partial_visibility << " " << temp_visibility << endl;
					//temp.push_back(v);
				}
			}

			else
			{	
				//cout << "AMUL - bbb " << intersectPoint1_found << " " << intersectPoint2_found << endl;
				// edit visibility and construct new segments
				total_visibility -= VisibleRegion[itr].partial_visibility;
				Pointlocation mid1;double angle_with_rect;double visibility_reduce_factor;
				if(intersectPoint1_found==true)
				{
					v.p1_overT=VisibleRegion[itr].p1_overT;
					v.p2_overT=intersectPoint2;

					mid1=mid_of_line_segment(v.p1_overT,v.p2_overT); //rectangle target, so p1,p2 length == p3,p4 length 
					angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, v.p1_overT.x,v.p1_overT.y,v.p2_overT.x,v.p2_overT.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
					visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					v.partial_visibility = Distance(VisibleRegion[itr].p1_overT, intersectPoint2)*visibility_reduce_factor;

					total_visibility +=v.partial_visibility;
					temp_visibility +=v.partial_visibility;
					//cout << "AMUL -3 " << v.partial_visibility << " " << temp_visibility << endl;
					//temp.push_back(v);
				}
				if(intersectPoint2_found==true) //changed 1 to 2 by lakshmi.. dont know if its correct, just checking
				{
					v.p1_overT= VisibleRegion[itr].p2_overT;
					v.p2_overT=intersectPoint1;
					mid1=mid_of_line_segment(v.p1_overT,v.p2_overT); //rectangle target, so p1,p2 length == p3,p4 length 
					angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, v.p1_overT.x,v.p1_overT.y,v.p2_overT.x,v.p2_overT.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
					visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					v.partial_visibility = Distance(VisibleRegion[itr].p2_overT, intersectPoint1)*visibility_reduce_factor;

					total_visibility +=v.partial_visibility;
					temp_visibility +=v.partial_visibility;
					//cout << "AMUL -4" << v.partial_visibility << " " << temp_visibility << endl;
					//temp.push_back(v);
				}
			}
			//cout << "CASE 2 output " << temp_visibility << " " <<  total_visibility << endl;
			if(temp_visibility!= 0.0 and temp_visibility/por_actual_dis >vis_threshold)
			{
			VisibleRegion[itr].partial_visibility= temp_visibility;
			itr++;
			}
			else
			{
			//cout << "CASE - 2.1 " << endl;
			VisibleRegion.erase(VisibleRegion.begin()+itr);
			total_visibility -= temp_visibility;
			}
			
				
			
		}//case 2 end bracket 

		/////******************************///// case 3
		else if( one_side_changed == true)
		{ 
			if(intersect(position, obstacle_p1, obstacle_p3, obstacle_p4)==false &&
			intersect(position, obstacle_p1, obstacle_p2, obstacle_p4)==false)//doesn't intersect with own arm, for point obstacle_p1
			{
				//cout << "lakkk case 1 " << endl;
				// is position-obstacle_p1 line intersects target? 
				//if so, is the point nearer the unchanged_side_val than the changed_side_val?
				//cout << "AMULYA 11aa" <<  " " << total_visibility << endl;
				p1_seen = change_vRegion_one_side(position, obstacle_p1, target, unchanged_side_val, changed_side_val, intersectPoint, 
					 VisibleRegion[itr]);
				//if so, then change the value of the v.region
				if(p1_seen)
				{
					VisibleRegion[itr].p1_overT= unchanged_side_val;
					VisibleRegion[itr].p2_overT= intersectPoint;
					
					
					//total_visibility -= VisibleRegion[itr].partial_visibility;
					
					Pointlocation mid1=mid_of_line_segment(unchanged_side_val,intersectPoint); 
					double angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, unchanged_side_val.x,unchanged_side_val.y,intersectPoint.x,intersectPoint.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
					double visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);
                                        
					//VisibleRegion[itr].partial_visibility = Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;
					//total_visibility += VisibleRegion[itr].partial_visibility;
					double x_vis= Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;					
					//total_visibility += x_vis;
					temp_visibility +=x_vis;
					//cout << "AMULYA 11" <<  " " << total_visibility << endl;
					
				}

			}

			if(p1_seen == false && intersect(position, obstacle_p2, obstacle_p1, obstacle_p3)==false &&
				intersect(position, obstacle_p2, obstacle_p3, obstacle_p4)==false)//doesn't intersect with own arm, for point obstacle_p2
			{
				//cout << "lakkk case 2 " << endl;
				p2_seen = change_vRegion_one_side(position, obstacle_p2, target, unchanged_side_val, changed_side_val, intersectPoint, 
					 VisibleRegion[itr]);
				//cout << "AMULYA 22aa" <<  " " << total_visibility << endl;
				if(p2_seen)
				{
					VisibleRegion[itr].p1_overT= unchanged_side_val;
					VisibleRegion[itr].p2_overT= intersectPoint;
					
					//total_visibility -= VisibleRegion[itr].partial_visibility;
					Pointlocation mid1=mid_of_line_segment(unchanged_side_val,intersectPoint); 
					double angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, unchanged_side_val.x,unchanged_side_val.y,intersectPoint.x,intersectPoint.y);//angle between the line connecting query point and mid of the p1,p2 and the line p1,p2
					double visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					//VisibleRegion[itr].partial_visibility = Distance(unchanged_side_val,intersectPoint)* visibility_reduce_factor;
					//total_visibility += VisibleRegion[itr].partial_visibility;
					 
					double x_vis= Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;
					//total_visibility += x_vis;
					temp_visibility +=x_vis;
					//cout << "AMULYA 22" <<  " " << total_visibility << endl;					
	
				}

			}

			if(p1_seen == false && p2_seen == false && intersect(position, obstacle_p3, obstacle_p2, obstacle_p4)==false &&
				intersect(position, obstacle_p3, obstacle_p1, obstacle_p2)==false)//doesn't intersect with own arm, for point obstacle_p3
			{
				//cout << "lakkk case 3 " << endl;
				p3_seen = change_vRegion_one_side(position, obstacle_p3, target, unchanged_side_val, changed_side_val, intersectPoint, 
					 VisibleRegion[itr]);
				//cout << "AMULYA 33aa" <<  " " << total_visibility << endl;
				if(p3_seen)
				{
					VisibleRegion[itr].p1_overT= unchanged_side_val;
					VisibleRegion[itr].p2_overT= intersectPoint;
					

					//total_visibility -= VisibleRegion[itr].partial_visibility;
					Pointlocation mid1=mid_of_line_segment(unchanged_side_val,intersectPoint); 
					double angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, unchanged_side_val.x,unchanged_side_val.y,intersectPoint.x,intersectPoint.y);
					double visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					//VisibleRegion[itr].partial_visibility = Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;
					//total_visibility += VisibleRegion[itr].partial_visibility;
					double x_vis= Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;
					//total_visibility += x_vis;
					temp_visibility +=x_vis;
					//cout << "AMULYA 33" <<  " " << total_visibility << endl;
					
				}
			}

			if(p1_seen == false && p2_seen == false && p3_seen == false && intersect(position, obstacle_p4, obstacle_p1, obstacle_p2)==false &&
				intersect(position, obstacle_p4, obstacle_p1, obstacle_p3)==false)//doesn't intersect with own arm, for point obstacle_p4
			{
				//cout << "lakkk case 4 " << endl;
				p4_seen = change_vRegion_one_side(position, obstacle_p4, target, unchanged_side_val, changed_side_val, intersectPoint, 
					 VisibleRegion[itr]);
				//cout << "AMULYA 44aa" <<  " " << total_visibility << endl;
				if(p4_seen)
				{
					VisibleRegion[itr].p1_overT= unchanged_side_val;
					VisibleRegion[itr].p2_overT= intersectPoint;
					
					//total_visibility -= VisibleRegion[itr].partial_visibility;
					Pointlocation mid1=mid_of_line_segment(unchanged_side_val,VisibleRegion[itr].p2_overT); 
					double angle_with_rect = angle_between_2_lines(mid1.x,mid1.y,position.x,position.y, unchanged_side_val.x,unchanged_side_val.y,intersectPoint.x,intersectPoint.y);
					double visibility_reduce_factor = reduced_visibility_factor4angle(angle_with_rect);

					//VisibleRegion[itr].partial_visibility = Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;
					//total_visibility += VisibleRegion[itr].partial_visibility;
					double x_vis= Distance(unchanged_side_val,intersectPoint)*visibility_reduce_factor;
					//total_visibility += x_vis;
					temp_visibility +=x_vis;
					//cout << "AMULYA 44" <<  " " << total_visibility << endl;
				}
			}
			//cout << "CASE 3 Before " << VisibleRegion[itr].partial_visibility << endl;
			if(temp_visibility!= 0.0 and temp_visibility/por_actual_dis > vis_threshold)
			{	
			VisibleRegion[itr].partial_visibility= temp_visibility;
			total_visibility -= VisibleRegion[itr].partial_visibility;
			total_visibility += temp_visibility;
			itr++;
			}
			else
			{
			VisibleRegion.erase(VisibleRegion.begin()+itr);
			//total_visibility -= VisibleRegion[itr].actual_dis;
			total_visibility -= temp_visibility;
			
			}
			//cout << "CASE 3 output " << VisibleRegion[itr].vr_id << " " << position.x << " " << position.y << " " << por_actual_dis << " " << temp_visibility << " " <<  total_visibility << endl;
			//cout << "CASE 3 threshold.. " << temp_visibility/por_actual_dis << endl;
			

		}//case 3 end bracket
		else
		{
		itr++;
		}
	}
	
	//cout << "END UPDATING............... " << "after visible segmenst size - " <<  VisibleRegion.size() << ",total visibility now -  " << total_visibility << endl;
}

class CompareRect 
{
    public:
    bool operator()(const Rectangle& t1, const Rectangle& t2) 
    {
		if(t1.upper_left.x<t2.upper_left.x) return true;
		return false;
    }
};

class CompareVisibility 
{
    public:
    bool operator()(QueryPoint& t1, QueryPoint& t2) // descending order
    {
		if(t2.total_visibility <= t1.total_visibility) return false;
		return true;
    }
};

Rectangle HeapEntry_to_rectangle(HeapEntry *he)
{
	Rectangle obj;
	float leftx,rightx,topy,lowy;	
	leftx = min(he->x1,he->x2);
	rightx = max(he->x1,he->x2);
	topy = max(he->y1,he->y2);
	lowy = min(he->y1,he->y2);
	
	obj.upper_left.x=leftx;
	obj.upper_left.y=topy;
	obj.lower_right.x=rightx;
	obj.lower_right.y=lowy;
	return obj;

}

bool Rectangle::operator==(Rectangle a)
{
	return (upper_left.x==a.upper_left.x) && (upper_left.y==a.upper_left.y) && (lower_right.x==a.lower_right.x) && (lower_right.y==a.lower_right.y);
}




/*bool line::operator==(line l)
{
	return (l.a == a && l.b == b);
}*/

int relative_oct_rectangles(Rectangle a, Rectangle b)
{
	int oct=0;
	if(a.lower_right.x<=b.upper_left.x)
	{
		if(a.lower_right.y>=b.upper_left.y)
		{
			oct=0;
		}
		else if(a.upper_left.y<=b.lower_right.y)
		{
			oct=6;
		}
		else 
		{
			oct=7;
		}
	}
	else if(a.upper_left.x >= b.lower_right.x)
	{
		if(a.lower_right.y>=b.upper_left.y)
		{
			oct=2;
		}
		else if(a.upper_left.y<=b.lower_right.y)
		{
			oct=4;
		}
		else
		{
			oct=3;
		}
	}
	else if(a.lower_right.y>=b.upper_left.y)
	{
		oct=1;
	}
	else if(a.upper_left.y<=b.lower_right.y)
	{
		oct=5;
	}
	else return -1;//overlapping rect
	return oct;
	
}

double MinDistBetweenRect(Rectangle a, Rectangle b)
{
	Rectangle1 a1;
	a1.x1 = a.upper_left.x;
	a1.x2 = a.lower_right.x;
	a1.y1 = a.lower_right.y;
	a1.y2 = a.upper_left.y;

	float *bounces =(float *)(malloc(sizeof(float)*4));
	*bounces = b.upper_left.x; bounces++;
	*bounces = b.lower_right.x; bounces++;
	*bounces = b.lower_right.y; bounces++;
	*bounces = b.upper_left.y; 
	bounces--;bounces--;bounces--;
	
	float dist = MINRECTDIST(a1,bounces);
	free(bounces);
	return dist;
	
}


void RTree::MOV(QueryPoint q[],int  NumOfObstacles, double vis_threshold,int dis_threshold, Rectangle T,int k, vector<QueryPoint> &k_answers,int& page, int& obstacle_considered,map<int, vector<Portion> > portion_data)
{

	int cost =0;
       //cout << "YESSS " << endl;
	int num_of_query_points = k;
	int k_position=0;
	int k_original=k;
	if(k<0) return;
	int end=0;
	int i,j;
	map<Rectangle,bool,CompareRect> obstacle_checked;
	priority_queue<QueryPoint, vector<QueryPoint>, CompareVisibility> k_ans; 

	map<int,bool> rnodeTraversed;
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);

	//priority queue of query points according to their visibility metric val
	priority_queue<QueryPoint, vector<QueryPoint>, CompareVisibility> qp_priority_queue; 

	//initially, the visibility is calculated without considering the obstacles
	//init the priority queue
	//cout << "calculating initial visibilty....... - printing query point and visible region size..." << endl;
	for(i=0;i<num_of_query_points;i++)
	{
		q[i].init_visibility(T,portion_data,dis_threshold,vis_threshold);
		//cout << q[i].position.x << " " << q[i].position.y << " - " << q[i].VisibleRegion.size() << " " << q[i].total_visibility << endl;
		qp_priority_queue.push(q[i]);	
		
	}
	//return;

	//------------------------------------------------------------
	
	int son = root; //this entry is to be visited next
	float leftx, topy, rightx, lowy;
	//son == -1 means empty && end = 0 means exit

	//while(false)
	while (son != -1 && end==0)
	{	
		
		map<int,bool>::iterator rnodeMapItr;
		rnodeMapItr=rnodeTraversed.find(son);
		if(rnodeMapItr!=rnodeTraversed.end() && rnodeTraversed.find(son)->second==true)
		{
			//traversed before,so do nothing
		}
		else
		{
			
			RTNode *rtn = new RTNode(this, son);
			rnodeTraversed[son]=true;
			page++;	
			cost++;
			for (i = 0; i < rtn -> num_entries; i++)
			{
				Rectangle nd;
			
				leftx = min(rtn->entries[i].bounces[0] , rtn->entries[i].bounces[1]);
				lowy = min(rtn->entries[i].bounces[2] , rtn->entries[i].bounces[3]);
				rightx = max(rtn->entries[i].bounces[0] , rtn->entries[i].bounces[1]);
				topy=max(rtn->entries[i].bounces[2] , rtn->entries[i].bounces[3]);
				//cout << "----------------------------------------" << endl;
				//cout << "MBR..." << leftx << " " << lowy << " " << rightx << " " << topy << endl;
				
				Rectangle1 nd_1; 

				nd_1.x1=nd.upper_left.x=leftx;
				nd_1.x2=nd.upper_left.y=topy;
				nd_1.y1=nd.lower_right.x=rightx;
				nd_1.y2=nd.lower_right.y=lowy;

				double edist1;
				//check for all query points, if nd is in the visibility region of any q, insert than into corresponding qp_heap
				for(j=0;j<num_of_query_points;j++)
				{
					//cout << "id of query point " << j << endl;
					Rectangle temp_rect;
					temp_rect.upper_left=temp_rect.lower_right=q[j].position;
					edist1 = MinDistBetweenRect(temp_rect, nd); // mindist from querypoint to node
					float dist2 = MinDistBetweenRect(T,nd);
					bool lak = q[j].IsInVisibleRegion(nd);
					//cout << "pehla check " << lak << endl;
					//if((edist1 <0.000001 && dist2 <0.000001) || (q[j].IsInVisibleRegion(nd)==true))
					if(q[j].IsInVisibleRegion(nd)==true)
					{
						//q and target are inside that mbr, or mbr is inside the visible region
						HeapEntry *h = new HeapEntry();
						h -> key = edist1;//sort the obs_heap according to mindist of obstacle from that qp
						h -> level = rtn -> level;
						h -> son1 = rtn->entries[i].son;
						h-> x1 = nd.upper_left.x;
						h-> x2 = nd.lower_right.x;
						h-> y1 = nd.upper_left.y;
						h-> y2 = nd.lower_right.y;
						(q[j]).obstacleList->insert(h);

						delete h;			

					}
				}
			
			}// end of num_entries's for 
			delete rtn;
		}	
		
		//get next entry from the heap of the top query point--------------------------

		if(qp_priority_queue.empty()==true)
		{
			end=1;
			break;
		}

		QueryPoint current_best_point;// = qp_priority_queue.top();
		//qp_priority_queue.pop();//visibility value change kore abar pore push korte hobe 

		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{	
                        //cout << "************************" << endl;		
			again = false;
			current_best_point = qp_priority_queue.top();
			//cout << "visss " << current_best_point.total_visibility << endl;
			//cout << "number of segments " << current_best_point.VisibleRegion.size() << endl;
			qp_priority_queue.pop();
			bool isempty = !current_best_point.obstacleList->remove(he);
			Rectangle obj;
			if (isempty)  //heap is empty, current_best_point is the answer
			{
				//cout << "KARMAAA 11 " << endl;
				k_ans.push(current_best_point);
				k--;

				if(k==0) end=1; //end;

				else
				{
					//qp_priority_queue theke best point already popped, so go for the next best point
					again=true;
					//break;
				}


			}
			//else if, check end condition, if reached the target
			////////////////////////////////////////////
			else if((obj=HeapEntry_to_rectangle(he)).equals_to(T))
			{
				//cout << "KARMAAA 222 " << endl;
				//end = 1; 
				k_ans.push(current_best_point);
				k--;

				if(k==0)end = 1; 

				else
				{
					//qp_priority_queue theke best point already popped, so go for the next best point
					again=true;
					//break;
				}
			}
			else
			{	
				//cout << "KARMAA " << endl;
				Rectangle1 p;
				p.x1 = he-> x1;
				p.y1 = he-> y1;
				p.x2 = he-> x2;
				p.y2 = he-> y2;

				obj= HeapEntry_to_rectangle(he);
				//check if obj is already been considered, if yes, don't consider it again;
			
		

				if (he->level == 0) //p is an object ; code for data object
				{
					//cout << "FOUND OBJECTTT " << endl;
					
					map<Rectangle,bool,CompareRect>::iterator mapItr;
					mapItr=obstacle_checked.find(obj);
					if(mapItr!=obstacle_checked.end())
					{
						if(obstacle_checked.find(obj)->second==true)
						{
							again=true;
							qp_priority_queue.push(current_best_point);
							continue;
						}
						
					}
					else obstacle_considered++;
					//calculate visibility for current_best_point

					//change visibility region of current_best_point

					if(current_best_point.IsInVisibleRegion(obj)==true)
						{
						//cout << "HAMMAYAAA " << endl;
						current_best_point.update_visibliliyRegion(obj,T,vis_threshold);
						//cout << current_best_point.total_visibility << endl;	
						}
					//check for all query points, if obj is in visibility region of q[p], 
					//pop q[p],calculate visibility for q[p] and then insert into the priority queue
					
					//change visibility region of q[p]

					int size1=qp_priority_queue.size();
					priority_queue<QueryPoint, vector<QueryPoint>, CompareVisibility> qp_temp; 
					for(int x=0; x<size1; x++)
					{
						//cout << "NOOO BABY " << endl;
						QueryPoint qp_check = qp_priority_queue.top();
						qp_priority_queue.pop();
						if(qp_check.IsInVisibleRegion(obj) == true)
						{
							//cout << "HAMMAYAAA " << endl;
							qp_check.update_visibliliyRegion(obj,T,vis_threshold);
						}
						qp_temp.push(qp_check);
					}
					qp_priority_queue = qp_temp; 
					
					//mark obj has been checked
					obstacle_checked[obj]=true;

					//get next data  from heap
					again =true;	
					qp_priority_queue.push(current_best_point);
				}
				else //not leaf node
				{
					qp_priority_queue.push(current_best_point);
					son=he->son1;
				}
			}
		}
		
		delete he;
	}
	// final ones after all obstacles and all q points are considered
	//cout << "everything done........ " << endl;
	QueryPoint final_qp;
	int size = k_ans.size();
	//QueryPoint k_answers[]
	vector<double> peace;
	for(i=0;i<size;i++)
	
		{	
		final_qp=k_ans.top();
		k_answers.push_back(final_qp);
		//cout << "&&&&&&" << final_qp.position.x << " " << final_qp.position.y << " - " << final_qp.VisibleRegion.size() << " " << final_qp.total_visibility << endl; 	
		
		if(final_qp.total_visibility!=0)
		peace.push_back(final_qp.total_visibility);
		
		for(int k =0;k<final_qp.VisibleRegion.size();k++)
		{
		 //cout << final_qp.VisibleRegion[k].vr_id << " ";	
		}
		//cout << endl;
		k_ans.pop();
			
		
	}
	//cout << "yayyy " << peace[peace.size()-1] << endl; 
        cout << "cost " << cost << endl; 
	delete heap;
	return;
	
}

