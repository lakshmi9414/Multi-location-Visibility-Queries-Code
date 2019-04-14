import re
import sys
import os
from operator import itemgetter
import numpy as np

class cppgcontinue:
    def __init__(self, inputfile, outputfile, minRF, minCS, maxOR):
        self.minRF = minRF
        self.minCS = minCS
        self.maxOR = maxOR
        
        self.dataset = dict()
        self.nots = self.getlines(inputfile)
        self.covpatterns = {}
        self.candidatecount=0
        [self.items, self.bitpattern] = self.dbscan(inputfile)
       

        sorteditems = sorted(self.items.items(), key = lambda a: (-a[1],a[0]))
        mintracs = self.minRF * 1.0 * self.nots
        freqitems = filter(lambda x: (x[1] >= mintracs), sorteditems)
        #print map(lambda x: list([x[0],x[1]*1.0/self.nots]), freqitems)
        self.freqitems = map(lambda x: x[0], freqitems)
        
        
        f = open(inputfile, 'r')
        counter = 1
        for i in f:
            i = i.rstrip('\n')
            l = i.split(" ")
            if len(l[-1]) == 0:
                l.pop()
            l = filter(lambda x: (self.items[x] >= mintracs), l)
            l.sort(key = lambda x: (-self.items[x],x))
            if(len(l)) > 0:
                self.dataset[counter] = l
                counter += 1
        f.close()
        #print 'Projected database formed'
	self.out = outputfile
    # self.getcoveragepatterns() call if not used in function



    def dbscan(self, db):
        f = open(db, 'r')
        a = {}
        bitpattern = {}
        counter = 0
        for i in f:
            counter += 1
            pattern = np.array([0] * (counter-1)  + [1] + [0] * (self.nots - counter))
            i = i.rstrip('\n')
            x = i.split(' ')
            if len(x[-1]) == 0:
                x.pop()
            for j in x:
                if a.has_key(j):
                    a[j] += 1
                    bitpattern[j] = np.bitwise_or(bitpattern[j], pattern)
                else:
                    a[j] = 1 
                    bitpattern[j] = pattern
                   
        #bitpattern[-1]=bitpattern[self.nots-1]
        #bitpattern[self.nots]=bitpattern[0]
        return [a, bitpattern]

    def getlines(self,filename):
        with open(filename,"r") as f:
            return sum(1 for _ in f)
    
    def flushcovpatterns(self):
        f = open(self.out, 'a')
	#print "*** " , len(self.covpatterns)
        for key in self.covpatterns:
            l =  str(self.covpatterns[key]) + "," + str(key[:-1]) +'\n'
            #print l,
            f.write(l)	
            #f.write(",")
            #f.write(key)
        
        #f.writelines(str(self.covpatterns)) 
        self.covpatterns = {}
        f.close()

    def getcoveragepatterns(self):
        count=0
        dataset = {}
        for key in self.dataset.iterkeys():
            dataset[key] = 0
        self.candidatecount = self.candidatecount + len(self.freqitems)        
        for item in self.freqitems:
            cs = (self.items[item]*1.0)/self.nots
            if (cs)>=self.minCS:
		#print item,(self.items[item]*1.0)/self.nots;
                #self.covpatterns.append(item+'\n')
                #self.covpatterns[item+'\n']=(self.items[item]*1.0)

                self.covpatterns[item+'\n']=cs
            projds = self.noProjection(dataset, [item])
            count += 1
            #print "-----"
            #print count,item,self.dataset.__len__(),projds.__len__(),self.coveragepatterns.__len__()
            self.cppgrec([item], 1, projds)
            del projds
	return len(self.covpatterns),self.candidatecount

    def noProjection(self, dataset, pattern):
        projdataset = dict()
        item = pattern[-1]
        #outfile = inpfile + item

        for key, value in dataset.iteritems():
            patternset = set(pattern)
            inputpattern = self.dataset[key][value:] 
            if set(inputpattern).isdisjoint(patternset):
                index = value
                itemindex = self.freqitems.index(item)
                for j in inputpattern:
                    if self.freqitems.index(j) < itemindex:
                        index += 1
                if index < len(self.dataset[key]):
                    projdataset[key] = index
	return projdataset

    def getfqs(self, ds):
        items={}
        for key, value in ds.iteritems():
            pattern = self.dataset[key][value:]
            for w in pattern:
                if not(items.has_key(w)):
                    items[w] = 1
        fqs = []
        for key,value in items.iteritems():
            fqs.append(key)
        return fqs 

    def patterntolist(self, pattern):
        init = np.array([0]*self.nots)
        #c1_init=np.array([0]*self.nots)
        #c2_init=np.array([0]*self.nots)
        for i in pattern:
            init = np.bitwise_or(init, self.bitpattern[i])
            #j1=i+1
            #j2=i-1
        return init

    def continuous(self,strpattern):
        c1_strpattern = np.array([0]*self.nots)
        c2_strpattern = np.array([0]*self.nots)
        c=np.where(strpattern==1)
        c1=c[0]+1
        c2=c[0]-1
        for x in c1:
            if x==-1:
                x=self.nots-1
            if x==self.nots:
                x=0    
            c1_strpattern[x]=1
        for y in c2:
            if y==-1:
                y=self.nots-1
            if y==self.nots:
                y=0    
            c2_strpattern[y]=1    

        return c1_strpattern,c2_strpattern
    def cppgrec(self, pattern, length, projds):

	
        fqs = self.getfqs(projds)
        self.candidatecount = self.candidatecount + len(fqs)
        #print "eneterred ", length, len(fqs) 
        strpattern = self.patterntolist(pattern)
        c1_strpattern,c2_strpattern=self.continuous(strpattern)
        #print "^^",strpattern,c1_strpattern,c2_strpattern
        
        for i in fqs:
            newpattern = pattern + [i]
            
            newstrpattern = np.bitwise_or(strpattern, self.bitpattern[i])
            c1_newstrpattern =np.bitwise_and(c1_strpattern, self.bitpattern[i])
            c2_newstrpattern =np.bitwise_and(c2_strpattern, self.bitpattern[i])
            
            newostrpattern = np.bitwise_and(strpattern, self.bitpattern[i])
            
            # overlap ratio
            ovratio = (np.sum(newostrpattern) * 1.0)/(np.sum(self.bitpattern[i]))
            
            #if(ovratio <= self.maxOR and ovratio >= self.minOR):
            #if(ovratio <= self.maxOR):
            if not(ovratio <= self.maxOR):
                continue
            
            #continus ratio
            c1_r=np.sum(c1_newstrpattern)
            c2_r=np.sum(c2_newstrpattern)
              
            if not(c1_r>=1 or c2_r>=1):
                continue

            # coverage support
            cs = np.sum(newstrpattern)/(1.0 * self.nots)
            #print newspattern, np.sum(newstrpattern)
            if cs >= self.minCS:
                #self.covpatterns[reduce(lambda x,y: str(x) +',' + str(y), newpattern)+'\n']=np.sum(newstrpattern)
                self.covpatterns[reduce(lambda x,y: str(x) +',' + str(y), newpattern)+'\n']=cs
                #print  reduce(lambda x,y: str(x) +',' + str(y), newpattern)+'\n'
                #print  np.where( strpattern==1 )
                #print np.where( self.bitpattern[i]==1)
                #print np.where( newstrpattern==1 )
                #print c1_r,c2_r
                #print "****"
		#print "for extra information", ovratio : keyword for overlap ratio;
                #print ",".join(newpattern),cs, ovratio

            newprojds = self.noProjection(projds, newpattern)
            #print "pattern:",newpattern,"proj dataset",projds,newprojds
            
            self.cppgrec(newpattern, length+1, newprojds)
            #print length
            #if sys.getsizeof(self.covpatterns)>102400:
                #self.flushcovpatterns()
            # deleting new projected dataset
            del newprojds
            

