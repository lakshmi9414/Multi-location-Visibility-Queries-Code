import cppgcontinue
import sys
f = open('out.txt','w')
f.close()

#minrf=0.1
minrf=float(sys.argv[1])/float(sys.argv[2])

#mincs=0.1
mincs=float(sys.argv[3])/float(sys.argv[4])

#maxor=0.025
maxor=float(sys.argv[5])/float(sys.argv[6]) # 1000 for max visibility in 3d files


p = cppgcontinue.cppgcontinue('transactions.txt','out.txt',float(minrf),float(mincs),float(maxor))
num_p = p.getcoveragepatterns()
#print minrf,",",mincs,",",maxor,",",num_p
#p.flushcovpatterns()
del p

