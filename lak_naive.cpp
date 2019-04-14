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
float minRF=0.2;
//float minCS = 0.3;
float maxOR=0.05;


//-----------------------MOV-----------
void build_transac(vector<QueryPoint> q,int total_portions_cnt,int num_queries)
{
	//cout << "BUILDDDing.. transac... " << total_portions_cnt << endl;
	//cout << "Q SIZE " <<  q.size() << endl;
 	 //vector< vector<int> > transac;
	 int i,j,por_id;
	
	for(i=0;i<total_portions_cnt;i++)
    	{
    	vector<int> temp_v;
    	total_transac.push_back(temp_v);
    	}	

	for(i=0;i<q.size();i++)
	{
		//cout << "&&&& " << i << endl;
		vector <VisibleRegionOverT>  cur_vis_segments =  q[i].VisibleRegion;
		for(j=0;j<cur_vis_segments.size();j++)
		{
		
		por_id = cur_vis_segments[j].vr_id;
		//cout << "JJJJ " << j << " " << por_id << endl;
           	total_transac[por_id].push_back(i);
		}	
	}
	//cout << "TRAAANSSS " << endl;
	//cout << "transactions size = " << total_transac.size() << endl;
	vector<int> qps;
	/*for(i=0;i<transac.size();i++)
	{
	qps=transac[i];
	for(j=0;j<qps.size();j++)
	{
	 cout << qps[j] << " ";	
	}
	cout << endl;
	}*/
	
	ofstream outputfile("transactions.txt");
	for(i=0;i<total_transac.size();i++)
	{
        //cout << "*** " << i << endl;
	qps=total_transac[i];
	for(j=0;j<qps.size();j++)
	{
	 outputfile << qps[j] << " ";	
	}
	outputfile << "\n";
	}
	//cout << "OVERRR " << endl;
	outputfile.close();
	//cout << "OVERRRRRRRRRRR " << endl;
	
	return;
  	
	

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
void get_smallparts_subset(vector<int> subset,int no_of_portions,float minCS)
{ 

	/*for(int d=0;d<subset.size();d++)
    {
  		cout << subset[d] << ",";
    }	
    cout << endl;*/

	//get smallest one
	vector<int> min_one=total_transac[subset[0]];
	int min_size=total_transac[subset[0]].size(),min_ind=0;
	for(int j=0;j<subset.size();j++)
	{
		if(total_transac[subset[j]].size()<min_size)
		{
			min_size=total_transac[subset[j]].size();
			min_ind=j;
			min_one=total_transac[subset[j]];

		}	
	}	


	float rel_freq;

	
    //return ;
	vector<int> A,B;
	int sel_ind;
	for(int j=0;j<subset.size();j++)
	{
		if(j!=min_ind)
		{	
			A=total_transac[subset[j]];
			sel_ind=j;
		    break;
	}

	}
	
	sort(A.begin(),A.end());
	rel_freq = float(float(A.size())/float(no_of_portions));
	if(rel_freq<minRF)
		return;    

    float cs_ratio,or_ratio;

	for(int i=0;i<subset.size();i++)
	{
		if(i!=min_ind and i!=sel_ind)
		{	
		B=total_transac[subset[i]];
		sort(B.begin(),B.end());
		rel_freq = float(float(B.size())/float(no_of_portions));
		if(rel_freq<minRF)
		return;

		vector<int> C(A.size()+B.size());
		vector<int>::iterator it= set_union(A.begin(),A.end(),B.begin(),B.end(),C.begin());
    	C.resize(it-C.begin());
    	A=C;
    	sort(A.begin(),A.end());
    }

    
	}
	
	B=total_transac[subset[min_ind]];	
	sort(B.begin(),B.end());
	rel_freq = float(float(B.size())/float(no_of_portions));
	if(rel_freq<minRF)
		return;


	//cout << "sizes " << A.size() << " " << B.size() << endl;
	
	vector<int> cs_un(A.size()+B.size());
	vector<int>::iterator it= set_union(A.begin(),A.end(),B.begin(),B.end(),cs_un.begin());
    cs_un.resize(it-cs_un.begin());
    


    vector<int> or_in(A.size()+B.size());
    std::vector<int>::iterator it1;
    it1=std::set_intersection (A.begin(),A.end(),B.begin(),B.end(), or_in.begin());
    or_in.resize(it1-or_in.begin());  
    
    //for(int b=0;b<A.size();b++)
    //cout << A[b] << " ";
    //cout << "___" << endl;
     //for(int b=0;b<A.size();b++)
    //cout << B[b] << " ";
    //cout << endl;
    
    //cout << "CHECKKKK " << or_in.size() << " " << B.size() << endl;    

    
    
    cs_ratio=float(float(cs_un.size())/float(no_of_portions)); //denominaror = no.of small portions.
    or_ratio = float(float(or_in.size())/float(B.size()));

    //cout << "ratios check - " << cs_ratio <<  " , " << or_ratio  << endl;

    if(float(cs_ratio)<float(minCS))
    {
         //cout << "ayyoo - " << cs_ratio << "-"<< minCS << " , " << or_ratio  << endl;
    	return;
    }
    if(float(or_ratio)>float(maxOR))
    {
       //cout << "yyyyyayyoo - " << cs_ratio << "-"<< minCS << " , " << or_ratio  << endl;
    	return;
    }

    cov_locations.push_back(subset);
    //cout << "yayyy" << subset.size() << endl;
    //cout << "yayyy" << cov_locations.size() << "ratios - " << cs_ratio << "," << or_ratio << endl;
    return ;


}


void getAllSubsets(vector<int> o_set,int no_of_portions,float minCS)
{
	vector<int> temp,set;
	vector< vector<int> > old_set,new_set;
	int s,p,item;
	for(s=0;s<o_set.size();s++)
	{
		temp.push_back(o_set[s]);
		old_set.push_back(temp);
		temp.clear();
	}
     /*for(s=0;s<old_set.size();s++)
     {
     	set=old_set[s];
     	//cout << set.size() << endl;
     	item=set[set.size()-1];
     	cout << item << endl;

     }*/	


	int len=1;	
	while(len<o_set.size()-1)
	{	
    for(s=0;s<old_set.size();s++)
    {
        
    	set = old_set[s];
    	item=set[set.size()-1];
    	for(p=0;p<o_set.size();p++)
    	{
    		temp.clear();
    		if(item < o_set[p])
    		{
    			temp=set;
    			temp.push_back(o_set[p]);
    			new_set.push_back(temp);
    			/*for(int d=0;d<temp.size();d++)
    			{
    				cout << temp[d] << ",";
    			}	
    			cout << endl;*/

    			get_smallparts_subset(temp,no_of_portions,minCS);
    		}	
    	}	
    }
    len++;
    old_set=new_set;
    new_set.clear();
}

    return;
}




bool Mov_exp(Rectangle T, int NumOfObstacles,double vis_threshold,int no_of_portions,int num_queries,int dis_threshold,float minCS)//just name, no extension needed
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
		
		//int min_dis=INT_MAX, max_dis=INT_MIN,avg_dis,dis_threshold;
		//Calc_query_to_target(q,T,i-1,&min_dis,&max_dis);
		//avg_dis=(min_dis+max_dis)/2;
		//cout << "MIN DIS " << min_dis << " " << "MAX DIS " << max_dis << endl;
		//dis_threshold = max_dis/2;
		//dis_threshold=50000;
		r->tree->MOV(q, NumOfObstacles,vis_threshold,dis_threshold,T,k,k_ans,page, obstacle_considered,portion_data);
		
		//cout << "finalllyyy :D :D ... " << endl;
		QueryPoint final_qp;
		
		/*for(int l=0;l<k_ans.size();l++)	
		{	
		final_qp=k_ans[l];
		cout << "-----------" << endl;
		cout << final_qp.position.x << " " << final_qp.position.y << " - " << final_qp.VisibleRegion.size() << " " << final_qp.total_visibility << endl;
		
		for(int k =0;k<final_qp.VisibleRegion.size();k++)
		{
		 cout << final_qp.VisibleRegion[k].vr_id << " ";	
		}
		cout << endl;
		}*/
		
		build_transac(k_ans,4*no_of_portions,i-1);

		vector<int> set,temp_one;
      		int s;
      		float rf;


      for(s=0;s<num_queries;s++)
      set.push_back(s);	

  	  for(s=0;s<set.size();s++)	
  	  {
          rf=float(float(total_transac[s].size())/float(no_of_portions));
          //cout << "&& " << total_transac[s].size() << "&& " << float(total_smallparts.size()) << " " << rf << endl;

	  //cout << "&& " << rf << endl; 
          if(rf>=minRF and rf>=minCS)
          {
          	//cout << "*** " << endl;
          	temp_one.clear();
          	temp_one.push_back(s);
          	cov_locations.push_back(temp_one);
          }	
  	  }	
      getAllSubsets(set,no_of_portions,minCS);
     //cout << "no.of subsets - " << pow(2,t_cnt) << endl;
      //cout << "no.of coverage locations - " << cov_locations.size() << endl;
      vector<int> cov_loc;
     		
     		/*for(int g=0;g<cov_locations.size();g++)
      		{
      	cov_loc=cov_locations[g];
      	for(int h=0;h<cov_loc.size()-1;h++)
      	{
      		cout << cov_loc[h] << ",";
      	}
        cout << cov_loc[cov_loc.size()-1];
      	cout << endl;
     	 }*/
      
			


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
	int dis_threshold = atoi(argv[4]);
	float mincs = atoi(argv[6]);
	mincs=mincs/100;
	


	
	
	//gen_input(num_queries);
		
	
	Mov_exp(T, NumOfObstacles, vis_threshold,no_of_portions,num_queries,dis_threshold,mincs);
        return 0;
	
}
