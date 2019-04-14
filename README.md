# musical-goggles
Commands to run:
A. Creating Portion Transactions.
  1.g++ -std=gnu++11 lak_main.cpp ./rtree/rtree.cpp ./blockfile/cache.cpp ./blockfile/blk_file.cpp ./linlist/linlist.cpp ./heap/heap.cpp ./rtree/distance.cpp ./rtree/entry.cpp ./rtree/rtnode.cpp ./his/histogram.cpp ./rtree/rtree_cmd.cpp ./func/gendef.cpp
  2./a.out $NumOfObstacles $visthreshold $numportions $numqueries $disthreshold
B. Extract multi-location visibility queries
	1.python run_cppg.py $minrf_n $minrf_d $mincs_n $mincs_d $maxor_n $maxor_d
	
