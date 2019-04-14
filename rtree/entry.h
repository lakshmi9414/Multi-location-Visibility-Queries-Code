/*entry.h
  this file defines class Entry*/
#ifndef __Entry 
#define __Entry
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class RTNode;
class RTree;
struct Linkable;
//------------------------------------------------------------
class Entry 
{
public:
//--===on disk===--
	int son;
	float *bounces;                     
//--===others===--
	int dimension;                      
	int level;
    RTree *my_tree;                     
    RTNode *son_ptr;                    
//for TP queries
	float dist; // the threshold of the entry when output
	int nn_num; // the nn that should be replaced 
   
//--===functions===--
	Entry();
	Entry(int dimension, RTree *rt);
    ~Entry();

	void del_son();
	Linkable *gen_Linkable();
	int get_size(); 
	RTNode *get_son();
	void init_entry(int _dimension, RTree *_rt);
	void read_from_buffer(char *buffer);// reads data from buffer
    SECTION section(float *mbr);        // tests, if mbr intersects the box
	bool section_circle(float *center, float radius);
	void set_from_Linkable(Linkable *link);
    void write_to_buffer(char *buffer); // writes data to buffer

    virtual Entry & operator = (Entry &_d);
	bool operator == (Entry &_d);
};

#endif