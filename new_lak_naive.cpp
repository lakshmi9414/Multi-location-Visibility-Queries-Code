//***************************************************
//This is the implementation of R*-tree v0.1

//Last revised July 4.
//***************************************************

#include <ctime>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <iostream>
#include <map>
#include <utility>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <vector>


#include<iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include<map>
#include <utility>  
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <iostream>
#include <map>
#include <utility>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <climits>


#include "./rtree/rtree.h"
#include "./rtree/rtnode.h"
#include "./rtree/entry.h"
#include "./blockfile/blk_file.h"
#include "./blockfile/cache.h"
#include "./linlist/linlist.h"
#include "./rtree/rtree_cmd.h"
#include "./rtree/distance.h"
#include "./his/histogram.h"
#include "./global.h"
#include "./func/gendef.h"

//Added by Tanzima
//#include "global.h"
//...............
using namespace std;



//----------------------------------------------------------------------------------------------------------
// SN 21/10/2007 001 <Start>

#define	DIMENSION 2
#define DEFAULT_C 0
const double FLOAT_INFTY =numeric_limits<double>::infinity();
//const double MAXDOUBLE = numeric_limits<double>::infinity();

const double MINX=0;
const double MINY=0;
const double MAXX=10000;
const double MAXY=10000;
const double THRESHOLD = 0.000000001;

//Experiments Parameters

const int SAMPLE=100; //was 1000
//DEFAULT
const int DEFAULT_GRPSIZE=64;
const int DEFAULT_K=8;
const double DEFAULT_M_AREAPART=0.08;
const double DEFAULT_M_RATIO=1;
const double DEFAULT_R_AREAPART=0.00005;
const double DEFAULT_R_RATIO=1;
//MIN
const int MIN_GRPSIZE=4;
const int MIN_K=2;
const double MIN_M_AREAPART=0.02;
const double MIN_M_RATIO=1;
const double MIN_R_AREAPART=0.00001;
const double MIN_R_RATIO=1;
//MAX
const int MAX_GRPSIZE=64; //WAS 1024
const int MAX_K=32;
const double MAX_M_AREAPART=0.32;
const double MAX_M_RATIO=16;
const double MAX_R_AREAPART=0.0001;
const double MAX_R_RATIO=16;
//INTERVAL
const int INTERVAL_GRPSIZE=4;
const int INTERVAL_K=2;
const double INTERVAL_M_AREAPART=2;
const double INTERVAL_M_RATIO=2;
const double INTERVAL_R_AREAPART=0.00001;
const double INTERVAL_R_RATIO=2;

//


//Tanzima
//char *TREEFILE = "F:/CSGLBSQ/Datasets/ca.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip20000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni20000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni15000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni10000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/uni5000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip15000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip10000.tree";
//char *TREEFILE = "C:/CSGLBSQ/Datasets/zip5000.tree";
//char *DATAFILE = "C:/CSGLBSQ/Datasets/zip5000.txt";
//char *RESFILE  = "F:/CSGLBSQ/Results/result.txt";
//char *RECTFILE  = "F:/CSGLBSQ/Results/rect.txt";
//char *RECTRATIOFILE  = "F:/CSGLBSQ/Results/rectratio.txt";
int io_access;
int disktime;
int updatetime;
int counttime;
int kmintime;
//--------------------------------------------------
//char *TREEFILE = "../Datasets/ca.24.tree";
//char *DATAFILE = "../Datasets/ca.txt";
//char *RESFILE  = "../Results/result.txt";

const int MAXPOINTS=10000;
const int MAXTRAJECTORYPOINTS=5000;
const int MAXTRAJECTORY=5;
//double PI = (double) 3.1415926535;
//typedef double Point2D[DIMENSION];

vector< vector<int> > total_transac,cov_locations;
//float minRF=0.2;
//float minCS = 0.3;
//float maxOR=0.05;

map<int, vector<int> > q_por_dict={};
map< vector<int>, vector<int> > cs_dict={}; 
vector< vector<int> > cov_patterns;

//-----------------------MOV-----------
void build_transac(vector<QueryPoint> q,int total_portions_cnt,int num_queries)
{
	
	int i,j;
	vector <VisibleRegionOverT>  cur_vis_segments;
	for(i=0;i<q.size();i++)
	{
	  
	      cur_vis_segments =  q[i].VisibleRegion;
	       // q_por_dict[i]=[];
		for(j=0;j<cur_vis_segments.size();j++)
		{
		q_por_dict[i].push_back(cur_vis_segments[j].vr_id);
		}
	}
	
	

}
double cal_distance(Pointlocation p1, Pointlocation p2)
{
	return sqrt( ((p1.x-p2.x)*(p1.x-p2.x)) + ((p1.y-p2.y)*(p1.y-p2.y)) );
}

void Calc_query_to_target(QueryPoint q[],Rectangle T,int num_queries,int *min_dis,int *max_dis)
{
	Pointlocation centroid;
	//T.upper_left.x=442148.29088;T.upper_left.y=111650.921465; // a random target from small british ordnance
	//T.lower_right.x=442174.190996;T.lower_right.y=111642.190447
	centroid.x = (T.upper_left.x+T.lower_right.x)/2;
	centroid.y = (T.upper_left.y+T.lower_right.y)/2;
	int i;
	double cur_dis;
	for(i=0;i<num_queries;i++)
	{
	  cur_dis = cal_distance(q[i].position,centroid);
	 *min_dis=min(*min_dis,(int)cur_dis);
	 *max_dis=max(*max_dis,(int)cur_dis);
	  q[i].dis_to_tar=cur_dis;	
	  //cout << "Cennnn... " << cur_dis << " " << *min_dis << " " << *max_dis << endl;	
	}
	return;
	
	

}
void check_cov_conditions(vector<int> pattern,int num_portions,float minRF, float minCS, float maxOR)
{
    int i;
    float rf, cs, orr;
    vector< pair<int,int> > pattern_or;
    pair <int,int> a;
    vector<int> cs_vec,or_vec;
    int len=pattern.size();
    
    bool rf_flag=false,cs_flag=false;
   
   for(i=0;i<pattern.size();i++)
    {
      
      a.first=q_por_dict[pattern[i]].size();
      a.second=pattern[i];
      pattern_or.push_back(a);
      cs_dict[{a.second}] = q_por_dict[pattern[i]];
      rf=float(a.first)/float(num_portions);
      //cout << "&& " << a.second << "," << a.first << " " << rf << endl;
      
    }
    if(len==1)
    {
     if(rf<minRF or rf< minCS)
      {
       return;
      }
    }
    else
    {
    int new_item = pattern[pattern.size()-1];
    
    //cs calculation
    sort(pattern.begin(), pattern.end(), std::greater<int>());
    
    cs_vec=pattern;
    cs_vec.pop_back();
    vector<int> A = cs_dict[cs_vec];
    vector<int> B = q_por_dict[new_item];
    
    
    vector<int> cs_un(A.size()+B.size());
    vector<int>::iterator it= set_union(A.begin(),A.end(),B.begin(),B.end(),cs_un.begin());
    cs_un.resize(it-cs_un.begin());
    cs_dict[pattern] = cs_un;
    
    cs=float(float(cs_un.size())/float(num_portions));
    
    //cout << "cs " << cs << endl;
    
    if(cs < minCS)
    cs_flag=true;
    
    //ol calculation
    sort(pattern_or.begin(), pattern_or.end());
    reverse(pattern_or.begin(), pattern_or.end());
    B=q_por_dict[pattern_or[pattern_or.size()-1].second];
    
    //cout << "&&&& " << B.size() << endl;
    
    for(i=0;i<pattern_or.size()-1;i++)
    {
      if(pattern_or[i].second!=new_item)
      or_vec.push_back(pattern_or[i].second);
    }
    if (or_vec.size()==pattern_or.size()-1)
    {
       A = cs_dict[or_vec];     
    }
    else
    {
     vector<int> tA = cs_dict[or_vec], tB=q_por_dict[new_item];
     vector<int> temp_un(tA.size()+tB.size());
    vector<int>::iterator it= set_union(tA.begin(),tA.end(),tB.begin(),tB.end(),temp_un.begin());
     temp_un.resize(it-temp_un.begin());
     A=temp_un;
     vector<int> store=or_vec;
     store.push_back(new_item);
     sort(store.begin(), store.end());
     cs_dict[store]=A;
    }
   
    vector<int> or_in(A.size()+B.size());
    std::vector<int>::iterator it1;
    it1=std::set_intersection (A.begin(),A.end(),B.begin(),B.end(), or_in.begin());
    or_in.resize(it1-or_in.begin());  
    
    if(B.size()==0)
    orr=0;
    else
    orr = float(float(or_in.size())/float(B.size()));
    
    /*for(int d=0;d<pattern.size();d++)
    {
     cout << pattern[d] << " ";
    }
    cout << ": " << cs << "," << orr << endl;*/
    
 
	    
    
    if(rf_flag == true or cs_flag == true or orr>maxOR)
    return;
    }
    
    cov_patterns.push_back(pattern);
    return;
}




void powerset(vector<int> q_ids,int num_queries,int num_portions,float minRF, float minCS, float maxOR)
{
   unsigned int pow_set_size = pow(2, num_queries);
    int counter, j,i;
    vector<int> pattern;
    for(counter = 0; counter < pow_set_size; counter++)
    { 
      pattern.clear();
      for(j = 0; j < num_queries; j++)
       {
          if(counter & (1<<j))
          {
            pattern.push_back(q_ids[j]);
          }
       }
       if(pattern.size()!=0)
       check_cov_conditions(pattern,num_portions,minRF,minCS,maxOR);
       
    }
    cout << "YAYYYY " << pow_set_size  << endl;
    /*for(i=0;i<cov_patterns.size(); i++)
    {
     vector<int> x = cov_patterns[i];
     //cout << "** " << x.size() << endl;
     for(j=0;j<x.size();j++)
     {
      cout << x[j] << " ";
     } 
    cout << endl;
    }*/


}

    





bool Mov_exp(Rectangle T, int NumOfObstacles,double vis_threshold,int no_of_portions,int num_queries,int dis_threshold,float minRF, float minCS, float maxOR)//just name, no extension needed
{
	// my own k deciding.
	int k=num_queries;
	
	//tree file
	char TREEFILE[200]="mytree.tree";
	
	//input file
	char DATAFILE[200]="/home/lakshmi/Desktop/obstacles/rtree/british_data.txt";
	
	//queries file
	char QUERYFILE[200]="british_queries.txt";
	
	//output file
	char OUTFILE[200]="test_out.txt";
	
	int blocksize = 1024; //4096;//1024;//4096;
	int b_length = 1024;
	//Cache *cache = new Cache(0, blocksize);
	
	//changing for num of obstacles	
	//int dimension=2;
	int dimension = NumOfObstacles;
	
		
	RTreeCmdIntrpr *r=new RTreeCmdIntrpr();
	r->build_tree(DATAFILE,TREEFILE,b_length,dimension,0);

	QueryPoint q[num_queries];
	vector<QueryPoint> k_ans;
	
	float temp;
	int i=0;
	FILE *fp; 
	FILE *fpDataSetSize;
	FILE *fp2;
	if((fp = fopen(QUERYFILE,"r")) == NULL)
	{
	//cout << "error opening file" << endl;
	//error("Cannot open query text file", TRUE);
        }
	else
	{
		i=0;
		while (i<=num_queries)//condition for restricting number of query points
		{
			fscanf(fp, "%f %f %f", &temp, &(q[i].position.x), &(q[i].position.y));
			i++; 
		}
		//cout << "over loading qurey points " << endl;
	}
	
	//cout << "NUmber of query points = " << (sizeof(q)/sizeof(*q)) << endl;
		
		int page,obstacle_considered;  page=obstacle_considered=0;
		fp2 = fopen(OUTFILE,"w+");
		float avg=0.0;
		float overallAvgTime=0.0;
		float overallAvgPage=0.0;
		float overallAvgObstacle=0.0;
		int count=0;
		
		
		for(int l=0;l<num_queries;l++)
		{
		q[l].obstacleList = new Heap();
		q[l].obstacleList->init(2);
		q[l].VisibleRegion.clear();
		}

		page=0;obstacle_considered=0;
                //cout << "before entering MOV funcion " << endl;
		map<int, vector<Portion> > portion_data;
		
		portion_data=r->tree->Calc_portions_on_target(T,portion_data,no_of_portions);
		
		
		r->tree->MOV(q, NumOfObstacles,vis_threshold,dis_threshold,T,k,k_ans,page, obstacle_considered,portion_data);
		
		
		QueryPoint final_qp;
		
		build_transac(k_ans,4*no_of_portions,i-1);
		int s;
		vector<int> q_ids;
		for(s=0;s<num_queries;s++)
                q_ids.push_back(s);
		powerset(q_ids,num_queries,4*no_of_portions,minRF,minCS,maxOR);

		for(int l=0;l<num_queries;l++)
		{
		delete q[l].obstacleList;
		}
		fclose(fp2);
		fclose(fp);
	return true;
	
}
void gen_input(int num_queries)
{
	double x,y;

  	random_device rand_dev;
  	mt19937 generator(rand_dev());

  	double minx=441834.0, maxx=442377.0,miny=110917.0,maxy=111719.0; // minx,maxx,miny,maxy

  	std::uniform_real_distribution<double> distributionx (minx-50.0, maxx+50.0);
  	std::uniform_real_distribution<double> distributiony (miny-50.0,maxy+50.0);
   	ofstream outputfile("british_queries.txt");
  	for (int i=0; i<num_queries; ++i)
  	{
  	x=distributionx(generator);
  	y=distributiony(generator);

    	std::cout << std::fixed;
    	outputfile << i+1 << " " << x << " " << y << "\n";
   	} 
	return;
}

//----------------------------------- main -----------------------------------
int main(int argc, char* argv[])
{
	Rectangle T;
	
	//2500--> 443.75 537.1875 3140.125 3175.281 46.52697 169.3154
	//T.upper_left.x=443.75;T.upper_left.y=3175.281; // a random target from british ordnance
	//T.lower_right.x=537.1875;T.lower_right.y=3140.125;

	//T.upper_left.x=500.0;T.upper_left.y=600.0; // my own for testing purpose
	//T.lower_right.x=600.0;T.lower_right.y=400.0;
	

	//T.upper_left.x=(442148.29088)*100;T.upper_left.y=(111650.921465)*100; // a random target from small british ordnance
	//T.lower_right.x=(442174.190996)*100;T.lower_right.y=(111642.190447)*100;

	//important query info:
	//441834 442377 110917 111719 - minx,maxx,miny,maxy

	T.upper_left.x=44214800.29088;T.upper_left.y=11165000.921465; // a random target from small british ordnance
	//T.lower_right.x=442174.190996;T.lower_right.y=111642.190447;
	T.lower_right.x=44257400.190996;T.lower_right.y=11144200.190447;

	//int NumOfObstacles = 550;
	//double vis_threshold = 0.4;
	//int no_of_portions = 10;
	//int num_queries = 500;

	int NumOfObstacles = atoi(argv[1]);
	double vis_threshold = atof(argv[2]);
	int no_of_portions = atoi(argv[3]);
	int num_queries = atoi(argv[4]);
	int dis_threshold = atoi(argv[5]);
	float minrf_n = atoi(argv[6]);
	float minrf_d = atoi(argv[7]);
	float mincs_n = atoi(argv[8]);
	float mincs_d = atoi(argv[9]);
	float maxor_n = atoi(argv[10]);
	float maxor_d = atoi(argv[11]);
	
	float minrf= minrf_n/minrf_d, mincs= mincs_n/mincs_d, maxor=maxor_n/maxor_d;
	


	
	
	//gen_input(num_queries);
		
	
	Mov_exp(T, NumOfObstacles, vis_threshold,no_of_portions,num_queries,dis_threshold,minrf, mincs, maxor);
        return 0;
	
}
