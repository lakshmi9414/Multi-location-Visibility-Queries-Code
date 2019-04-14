/* rtree.h
   this file defines the class RTree*/

//#undef min
//#undef max
#ifndef __RTREE
#define __RTREE
#define PI 3.14159265
//------------------------------------------------------------
#include "../func/gendef.h"
#include "../heap/heap.h"
#include <vector>
#include <map>
#include <list>
#include <utility>  

//Added by Tanzima
//#include "../global.h"
//...............

//Added by Tanzima
#include <vector>
using namespace std;


//added for TP KNN------------------------------------------------
#define SEQUENCE_SENSITIVE false
  //set it true if you want a new partition point every time the order of the NNs change
#define PAST_PAIR 1000

//Added for kGNN by Tanzima


const double MAXDOUBLE = numeric_limits<double>::infinity();
const int kMAX=100; // changed by Eunus : it was set to 500
const int MAXDATALIMIT=20000;
//const int MAXGROUP = 128;
//const int MINSUBGROUP = MAXGROUP/2;
//------------------------------------------------------------
class LinList;
class SortedLinList;
class Cache;
class RTNode;
class Entry;
//Added by Tanzima
struct Rectangle1;
struct DistfromPOI;
//------------------------------------------------------------


//added by fa
struct VisibleRegionOverT //
{
     int vr_id;	
     Pointlocation p1_overT;
     Pointlocation p2_overT;
     double actual_dis;	

     double partial_visibility;
	// int operator==(VisibleRegionOverT b);
};
struct Portion
{
       int por_id;
       Pointlocation p1_overT;
      Pointlocation p2_overT;
	
};


class Rectangle
{
	public:
      Pointlocation upper_left;
	  Pointlocation lower_right;
	  bool operator==(Rectangle a);
	  bool operator<(Rectangle a);
	  bool operator<=(Rectangle a);
	  bool operator>(Rectangle a);
	  bool equals_to(Rectangle a)
	  {
		  return (upper_left.x==a.upper_left.x) && (upper_left.y==a.upper_left.y) && (lower_right.x==a.lower_right.x) && (lower_right.y==a.lower_right.y);
	  }
};

class Rectangle3D
{
public:
	Point3D p1;
	Point3D p2;
	Point3D p3;
	Point3D p4;
	Rectangle3D(Point3D a, Point3D b,Point3D c,Point3D d){p1=a;p2=b;p3=c;p4=d;}
	Rectangle3D()
	{
	}
};

struct line
{
	Point3D a;
	Point3D b;
	line(){}
	line(Point3D _a,Point3D _b){a=_a;b=_b;}

	bool operator==(line a);
};
struct polygon
{
	vector <line> sides;
	float Polyarea;
	polygon(){};
	polygon(vector<line> lines)
	{
		sides = lines;
	}
	bool operator==(polygon poly);
};

struct VisibleRegion3D
{
	Rectangle3D boundary;
	float partial_visibility;
	list<polygon> invisible_parts;
	VisibleRegion3D()
	{
	}
};

struct inVR
{
	Point3D obstacle_corner;
	VisibleRegion3D vRegion;
	Point3D intersection_p;
};

/*int VisibleRegionOverT::operator==(VisibleRegionOverT b)
{
	return (p1_overT.x==b.p1_overT.x) && (p1_overT.y==b.p1_overT.y) && (p2_overT.x==b.p2_overT.x) && (p2_overT.y==b.p2_overT.y); 
}*/

class Box
{
public:
	Point3D a;
	Point3D b;
	Box(Point3D _a,Point3D _b){a=_a;b=_b;}
	Box(){}
	bool equals_to(Box T)
	{
		if(T.a.x == a.x && T.a.y == a.y && T.a.z == a.z && T.b.x == b.x && T.b.y == b.y && T.b.z ==b.z)
			return true;
		return false;
	}
};


class BoxPlanes
{
public:
	int a[4];//corner point number
	bool visible;
	BoxPlanes(int _a,int _b,int _c,int _d)
	{
		a[0]=_a;a[1]=_b;a[2]=_c;a[3]=_d;
	}
	BoxPlanes()
	{
		a[0]=a[1]=a[2]=a[3]=0;
	}
	
	bool isCornerInBoxPlanes(int p)
	{
		if(p == a[0] || p == a[1] || p == a[2] || p == a[3]) return true;
		else return false;
	}
};

class QueryPoint3D
{
public:
	Point3D position;
	float total_visibility;
	Heap *obstacleList; //obstacles that affects visiblilty, according to their mindist
	vector <VisibleRegion3D> visiblePlanes;
	float init_visibility(Box Target);//considering no obstacle
	bool IsInsideVisibleRegion(Box T,Box obj); 
	void update_visibilityRegion(Box obstacle, Box target);
	//float update_naiveVisibility(Box obj,Box T);
};

class QueryPoint
{
    public:  
		
	  Pointlocation position;
	  double dis_to_tar;	
	  
	  double total_visibility; //for this query_point only
	  
	  Heap *obstacleList; //obstacles that affects visiblilty, according to their mindist
         vector <VisibleRegionOverT>  VisibleRegion;//segments over Target object T, that are visible from this query point
	 
	  
         // list <Rectangle> ObstacleList; 
          bool IsInVisibleRegion(Rectangle rect); 

	  bool  del_from_ObstacleList(Rectangle rect);

	  void update_visibliliyRegion(Rectangle obstacle, Rectangle target,double vis_threshold);//IsInVisibleRegion true hole oi obstacle er jonno update korte hobe
	  double visibiliy(Rectangle Target);
      
	  double init_visibility(Rectangle Target,map<int, vector<Portion> > portion_data,int dis_threshold,double vis_threshold);//jokhon kono obstacle nai consider kore
      
};


//-------------------------------------------------------------
class RTree : public Cacheable
{
public:
//--===on disk===--
	int dimension;                       
	int num_of_data;	                 
    int num_of_dnodes;	                 
    int num_of_inodes;	                 
	int root;                            
	bool root_is_data;                   
//--===others===--
	RTNode *root_ptr;
    bool *re_level;  
    LinList *re_data_cands; 
	LinList *deletelist;

//--===added for TP KNN===--
	int last_pair[PAST_PAIR][2]; //records the pairs of points that produce the minimum inflence time
	int lastcnt; //next last pair to be replaced
	Heap *tpheap;

//--===functions===--
    RTree(char *fname, int _b_length, Cache* c, int _dimension);
    RTree(char *fname, Cache* c);
    RTree(char *inpname, char *fname, int _blength, Cache* c, int _dimension);

    ~RTree();
	void del_root();
	bool delete_entry(Entry *d);
	bool FindLeaf(Entry *e);
    int get_num() { return num_of_data; }
	void insert(Entry *d);
	void load_root();  
	void rangeQuery(float *mbr, SortedLinList *res);
	void read_header(char *buffer);      
	void write_header(char *buffer);
	int update_rslt(Entry *_e, float _dist, Entry *_rslt, 
					 float *_key, int _k);

	// This function was added to perform TP-kNN queries by Bobby
	void TPNN_TP(float *_qline, int _k, Entry *_nn, Entry *_rslt, float _max_trvl);
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
	void BFNN(float *_qpt, int _k, Entry *_rslt);
	void BFNNCont(float *qmbr, float *qmbr2, Heap *heap, HeapEntry *e, int k);
	//Added by Tanzima
	float KMin(Point2D m, int c, float cl, int k, Pointlocation _rslt[], int *num_of_data,float d_safe_corner,DistfromPOI _cornertoPOI[]);
	void UpdateCount(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, float r, int *count, DistfromPOI _cornertoPOI[]);
	int UpdateStatus(Rectangle1 R, float cl, int k, Pointlocation _rslt[], int *num_of_data, int *count, DistfromPOI _cornertoPOI[]);
	void Rect_kNNQ(Rectangle1 R, float cl, int k, Pointlocation *_rslt, int *num_of_data);
	void Point_BFN_NNQ(Point2D o, double *_rslt);

	//Added by Tanzima for kGNN
	//float RTree::KMax(int k, Pointlocation _rslt[], int *num_of_data);
	//void private_kGNN_sum(Rectangle1 R[], int g_size, int k, Pointlocation _rslt[], int *num_of_data);
	//void private_kGNN_max(Rectangle1 R[], int g_size, int k, Pointlocation rslt[], int *num_of_data);

	//Added by Eunus for Consensus
	//void RTree::consensus_kGNN(Rectangle1 R[], int g_size, int sg_size, int k, int f, SG L[][32], Pointlocation rslt[], int *num_of_data);


	//fa
	map< int ,vector<Portion> > Calc_portions_on_target(Rectangle T, map<int, vector<Portion> > portion_data,int no_of_portions);
	void MOV(QueryPoint q[],int  NumOfObstacles, double vis_threshold,int dis_threshold, Rectangle T,int k, vector<QueryPoint> &k_answers,int& page, int& obstacle_considered,map<int, vector<Portion> > portion_data);
	float VCM_visibility(QueryPoint q,Rectangle T);

};
bool ifIntersect3D_lineSegment(Point3D a, Point3D b, Point3D c, Point3D d, Point3D &intersect_Point);
#endif // __RTREE
