import cppg
import sys
f = open('out.txt','w')
f.close()

#minrf=0.1
minrf=float(sys.argv[1])/float(sys.argv[2])

#mincs=0.1
mincs=float(sys.argv[3])/float(sys.argv[4])

#maxor=0.025
maxor=float(sys.argv[5])/float(sys.argv[6]) # 1000 for max visibility in 3d files


p = cppg.cppg('transactions.txt','out.txt',float(minrf),float(mincs),float(maxor))
num_p, cand_count = p.getcoveragepatterns()
#print minrf,",",mincs,",",maxor,",",cand_count,",",num_p

#print num_p,
#p.flushcovpatterns() ##for getting coverapettens in a file no fluah for time
#no_p = sum(1 for line in open('out.txt'))
#print no_p,minrf,mincs,maxor
#f.close()
del p

