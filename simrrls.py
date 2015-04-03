#!/usr/bin/env python2

"""
##################################################
##  Script to simulate RADseq-like data         ##
##  author:  Deren Eaton                        ##
##  contact: deren.eaton@yale.edu               ##
##  date:    4/2/15                             ##
##  version: 1.04                               ##
##################################################

##################################################
##  change log
##  1.04: - added option to change sequencing depth
##        - dropout is now controlled by the data type
##            and sequence length with insert
##  1.03: - paired data file names have _R1_ & _R2_
##  1.02: - dropout is only affected by arg 2 not data type
##  1.01: - barcodes differ by 2 bp by default
##################################################
"""

## load modules                                        
import egglib
import sys
import numpy as np
import gzip


""" load/parse arguments """
indel   = float(sys.argv[1])
sitemut = int(sys.argv[2])
nLoci   = int(sys.argv[3])     
Ninds   = int(sys.argv[4])     ## per tiptax
depth   = str(sys.argv[5])
frags   = str(sys.argv[6])
datatype = str(sys.argv[7])
outhandle = str(sys.argv[8])



## check args
if datatype in ['rad','gbs','pairgbs','ddrad','pairddrad']:
    print '\tsimulating '+datatype+" data"
else:
    print 'datatype not recognized'
    sys.exit()

if "," in depth:
    stringdepth = "N("+depth+")"
else: 
    stringdepth = depth

print "\t"+str(nLoci)+" loci at "+stringdepth+"X coverage"
print "\t"+str(Ninds)+" samples per taxon across 12 tip taxa"
print "\tindels arise at frequency of %.4f per mutation" % indel
print "\tlocus dropout =",bool(sitemut)
print "\tsize selection window =",frags
print "\tsequencing error rate = 0.0005"

## fixed args    
locuslength = 1000
N = 50000
u = (7e-9)*locuslength
theta = 4*N*u
print '\ttheta=4Nu=',theta/locuslength

if int(frags.split(",")[0]) < 200:
    print "\tmin fragment length allows read overlaps/adapter sequences "


""" cuts at Pst1 (and) at EcoR1, attaches barcode to Pst1 """
CUT1 = "CTGCAG"    ## C|TGCAG
CUT2 = "GAATT"     ## complement of G|AATT
cut1 = "TGCAG"
cut2 = "AATT"

tiptax = list("ABCDEFGHIJKLX")
tiptax = [i+j for i,j in zip(list("1111222233334"),list("ABCDEFGHIJKLX"))]
outgroup = "4X"

"random number seeds used in test example"
np.random.seed(998877)

"simulate species with divergence times"
A=B=C=D=E=F=G=H=I=J=K=L=X=Ninds   ## number of diploid individuals to sample
"clade 1 = (((A,B),C), D)"
"clade 2 = (((E,F),G), H)"
"clade 3 = (((I,J),K), L)"
"tree = (((((A,B),C),D),(((E,F),G),H)),(((I,J),K),L))"
newick = "(((((A:2,B:2):2,C:4):4,D:8):4,(((E:2,F:2):2,G:4):4,H:8):4):4,(((I:2,J:2):2,K:4):4,L:8):8,X:16):16;"

divscale=1.0

nodeABCDEFGHIJKLX  = 4.*divscale
nodeABCDEFGHIJKL   = 4.*divscale
nodeABCDEFGH       = 3.*divscale
nodeIJKL           = 2.*divscale
nodeEFGH           = 2.*divscale
nodeABCD           = 2.*divscale
nodeIJK            = 1.*divscale
nodeEFG            = 1.*divscale
nodeABC            = 1.*divscale
nodeIJ             = 0.5*divscale
nodeEF             = 0.5*divscale
nodeAB             = 0.5*divscale

# sets the two parameter classes
paramSet = egglib.simul.CoalesceParamSet(singleSamples=None, doubleSamples=[A,B,C,D,E,F,G,H,I,J,K,L,X],M=0.0)
"clade 1"
paramSet.populationFusion(nodeABCD, 0,3)     ## D into A
paramSet.populationFusion(nodeABC, 0,2)      ## C into A
paramSet.populationFusion(nodeAB, 0,1)       ## B into A
"clade 2"
paramSet.populationFusion(nodeEFGH, 4,7)     ## H into E
paramSet.populationFusion(nodeEFG, 4,6)      ## G into E
paramSet.populationFusion(nodeEF, 4,5)       ## F into E
"clade 3"
paramSet.populationFusion(nodeIJKL, 8,11)    ## L into I
paramSet.populationFusion(nodeIJK, 8,10)     ## K into I
paramSet.populationFusion(nodeIJ, 8,9)       ## J into I
"together and outgroup"
paramSet.populationFusion(nodeABCDEFGHIJKLX, 0,12) ## X into A
paramSet.populationFusion(nodeABCDEFGHIJKL, 0,8)   ## I into A
paramSet.populationFusion(nodeABCDEFGH, 0,4)       ## E into A


mutator = egglib.simul.CoalesceFiniteAlleleMutator(theta=theta,
                                                   alleles= 4,
                                                   randomAncestralState=True)
mutator.setSites(locuslength)


def twodiffs(a,b):
    "requires two base differences between barcodes"
    if len(a) == len(b):
        t = [a[i]==b[i] for i in range(len(a))]
        if t.count(False) > 1:
            return True


def barcoder(aligns, Ninds, tiptax, outgroup, outhandle, Barcodes):
    """ takes an align object and puts names to sequences and
    makes a barcode map """

    " append names"
    names = []
    for name in tiptax:
        for i in range(Ninds):
            names.append(name+str(i))
            names.append(name+str(i))

    " appends names for allele 1 vs. allele 2 for each diploid sample "
    l = 1
    for i in aligns:
        for s,j in zip(i,names):
            s.name = j
    
    " make dictionary with list of loci for each sample "
    D = {}
    for i in set(names):
        D[i] = []
    for samp in D:
        for i in aligns:
            l = []
            for j in i:
                if j.name == samp:
                    l.append(j.sequence)
            D[samp].append(l) 

    " generate barcodes to append to start of loci "
    bases = list("ATGC")
    if not Barcodes:
        print "creating new barcode map"
        Barcodes[names[0]] = "CATCAT"
        for name in names:
            while name not in Barcodes.keys():
                bbs = np.random.randint(0,3,6)
                bar = "".join([bases[i] for i in bbs])
                if all([twodiffs(bar,bb) for bb in Barcodes.values()]):
                    if not any([i in bar for i in [cut1,cut2]]):
                        Barcodes[name] = bar
        " creates random barcodes and writes map to file "
        with open(outhandle+".barcodes",'w') as barout:
            bnames = list(Barcodes.keys())
            bnames.sort()
            for bar in bnames:
                if outgroup not in bar:
                    print >>barout, "\t".join([bar,Barcodes[bar]])
    return D,Barcodes



def mutate(base):
    "introduce point sequencing errors"
    nbase = list("ATGC")[np.random.randint(0,4)]
    while nbase == base:
        nbase = list("ATGC")[np.random.randint(0,4)]
    return nbase


def revcomp(seq):
    S = seq.replace('A','t').replace("T",'a').replace("C",'g').replace("G",'c')
    return S[::-1].upper()


def locus_dropout(datatype,seq,cut1,cut2):
    if datatype=='rad':
        return bool( any( [i in seq for i in cut1,revcomp(cut1)]) )
    if 'gbs' in datatype:
        return bool( any( [i in seq for i in cut1,revcomp(cut1)]) )
    if 'ddrad' in datatype:
        return bool( any( [i in seq for i in cut1, revcomp(cut1),
                                             cut2, revcomp(cut2)]) )


def dressitup(D, Barcodes, sitemut, indel, locuslength, frags, counter, Nloci):
    """ return a list of sequences all dressed up """
    SEQ1s     = []
    SEQ2s     = []
    Illumina_P1_adapter = "ACGACGCTCTTCCGATCT"
    Illumina_P2_adapter = "AGATCGGAAGAGCTCGTATG"

    for loc in range(100):
        ## sample fragment in size selection window for this locus
        frag1,frag2 = frags.split(',')
        insert = np.random.randint(int(frag1),int(frag2))
        frag = 200+insert

        ## do we keep this locus
        keptoutg = 0
        if not locus_dropout(datatype,D['4X0'][loc][0][:frag],cut1,cut2):
            counter += 1
            keptoutg = 1

        ## sample depth of sequencing at this locus
        if ',' in depth:
            meandepth,stdevdepth = depth.split(",")
            copy1 = int(np.random.normal(float(meandepth)/2.,float(stdevdepth)/2.,1))
            copy2 = int(np.random.normal(float(meandepth)/2.,float(stdevdepth)/2.,1))
        else:
            ## return even number copies
            copies = round(int(depth)/2.)*2
            copy1 = copy2 = int(copies/2)

        if indel:
            skip = 5
            "indel = probability site is checked for indel-mutation, skip first 'skip' bp"
            inds = np.random.binomial(locuslength-skip,indel)
            "if indel, for each indel site randomly choose location on locus"
            wh = [np.random.randint(skip,locuslength) for i in range(inds)]
            "randomly select size of indel (default size 1)"
            le = [np.random.randint(1,2) for i in wh]
            where = [wh,le]

        for samp in D:
            """ exclude the outgroup taxon used for detecting restriction site
            mutations and indels use both chromosomes """
            if sitemut:
                if keptoutg:
                    " remove locus if mutation arises in restriction site"
                    L1 = [list(D[samp][loc][0]) for copy in range(copy1) if not \
                          locus_dropout(datatype,D[samp][loc][0][:frag],cut1,cut2)]

                    L2 = [list(D[samp][loc][1]) for copy in range(copy2) if not \
                          locus_dropout(datatype,D[samp][loc][0][:frag],cut1,cut2)]
                else:
                    L1 = L2 = []
            else:
                L1 = [list(D[samp][loc][0]) for copy in range(copy1)] 
                L2 = [list(D[samp][loc][1]) for copy in range(copy2)]
            L = L1+L2

            for copy in range(len(L)):
                ll = locuslength
                """ introduce indel if mutation present at site
                relative to the outgroup """
                if indel:
                    for ww,ee in zip(where[0],where[1]):
                        df = locuslength-ll
                        ww -= df
                        if L[copy][ww] != D['4X0'][loc][0][ww]:
                            L[copy] = L[copy][:ww] + L[copy][ww+ee:]
                            ll -= ee

                " introduce sequencing errors"
                Es = np.random.binomial(ll,0.0005)
                Elocs = np.random.randint(0,ll,Es)
                eL = L[copy]
                for err in Elocs:
                    eL[err] = mutate(eL[err])

                " shorten fragment to potentially overlapping length"
                eL = eL[:frag]

                " make fastQ formatted"
                if datatype in ['ddrad','pairddrad']:
                    eL1 = ("A"*30)+Illumina_P1_adapter+Barcodes[samp]+cut1+\
                          "".join(eL)+revcomp(cut2)+Illumina_P2_adapter+"A"*30
                else:
                    eL1 = ("A"*30)+Illumina_P1_adapter+Barcodes[samp]+cut1+\
                          "".join(eL)+revcomp(cut1)+Illumina_P2_adapter+"A"*30

                if datatype in ['gbs']:
                    eL2 = ("A"*30)+Illumina_P1_adapter+Barcodes[samp]+cut1+\
                          revcomp("".join(eL))+revcomp(cut1)+Illumina_P2_adapter+"A"*30
                    
                " sequence read1 from 5' end "
                startL = eL1.index("CGATCT")
                "---->"
                ii = startL+len("CGATCT")
                sss = eL1[ii:ii+100]

                if "X" not in samp:
                    if keptoutg:
                        if counter <= Nloci:
                            SEQ1s.append("@lane1_locus"+str(counter)+"_"+str(D.keys().index(samp))+"_R1_"+str(copy)+\
                                         " 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

                if datatype in ['gbs']:
                    " sequence read1 from 3' end for GBS"
                    startL = eL2.index("CGATCT")
                    "---->"
                    ii = startL+len("CGATCT")
                    sss = eL2[ii:ii+100]
                    if "X" not in samp:
                        if keptoutg:
                            if counter <= Nloci:
                                SEQ1s.append("@lane1_fakedata"+str(counter)+"_"+str(D.keys().index(samp))+\
                                             "_R1_"+str(copy)+" 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

                " sequence read2 from 3' end for paired read "
                startR = eL1.rindex("AGATCG")
                "<----"
                ii = len(eL1)-startR
                sss = eL1[-1*(ii+100):-ii]
                if "X" not in samp:
                    if keptoutg:
                        if counter <= Nloci:
                            SEQ2s.append("@lane1_locus"+str(counter)+"_"+str(D.keys().index(samp))+\
                                         "_R2_"+str(copy)+" 1:N:0:"+"\n"+revcomp(sss)+"\n+\n"+("B"*len(sss))+"\n")

                if datatype in ['pairgbs']:
                    ' double check this '
                    " sequence read2 from 5' end for GBS"
                    startL = eL2.index("CGATCT")
                    "---->"
                    ii = startL+len("CGATCT")
                    sss = eL2[ii:ii+100]
                    if "X" not in samp:
                        if keptoutg:
                            if counter <= Nloci:
                                SEQ2s.append("@lane1_locus"+str(counter)+"_"+str(D.keys().index(samp))+\
                                             "_R1_"+str(copy)+" 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

    return SEQ1s,SEQ2s,counter




print '\t',
D = []
Barcodes = {}

out1 = gzip.open(outhandle+"_R1_.fastq.gz",'wb')
if datatype in ['pairddrad','pairgbs']:
    out2 = gzip.open(outhandle+"_R2_.fastq.gz",'wb')

if datatype in ['gbs']:
    nrepets = nLoci/2  # because forward and reverse sequenced

counter = 0
while counter < nLoci: #for i in range(0, nrepets, divis):
    " performs simulation of 1K loci at a time"
    aligns = egglib.simul.coalesce(paramSet, mutator, 100)

    " creates barcodes if not made yet and fill D dict."
    D,Barcodes = barcoder(aligns, Ninds, tiptax, outgroup, outhandle, Barcodes)

    " dresses up data to be fastq and puts in errors, indels, etc."
    SEQ1s,SEQ2s,counter = dressitup(D, Barcodes, sitemut, indel,
                                    locuslength, frags, counter, nLoci)
    out1.write("".join(SEQ1s))
    if datatype in ['pairddrad','pairgbs']:
        out2.write("".join(SEQ2s))

out1.close()
if datatype in ['pairddrad','pairgbs']:
    out2.close()
