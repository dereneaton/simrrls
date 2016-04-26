#!/usr/bin/env python2

"""
##################################################
##  Script to simulate RADseq-like data         ##
##  author:  Deren Eaton                        ##
##  contact: deren.eaton@yale.edu               ##
##  date:    4/16/15                            ##
##  version: 1.05                               ##
##################################################

###############################################################
##  change log                                               ##
##  1.05: - major re-write/cleanup, and opt-parsing added    ##
##            and user-supplied topology                     ##
##  1.04: - added option to change sequencing depth          ##
##        - dropout is now controlled by the data type       ##
##            and sequence length with insert                ##
##  1.03: - paired data file names have _R1_ & _R2_          ##
##  1.02: - dropout is only affected by arg 2 not data type  ##
##  1.01: - barcodes differ by 2 bp by default               ##
###############################################################
"""

## load modules                                        
import egglib
import sys
import numpy as np
import gzip
import getopt


def parseopts(opts):
    """ parses command line arguments """

    ## defaults
    params = {'indels': 0,
              'dropout': 0,
              'nLoci': 1000,
              'Ninds': 1,
              'error': 0.0005,
              'N': 1e5,
              'mu': 5e-9,
              'length': 100,
              'insert': "300,800",
              'datatype': 'rad',
              'tree': "default 12 tip tree",
              'outname': "out",
              'depth': "10,1e-9",
              'cut1': 'CTGCAG',
              'cut2': 'GAATTC',
              'seed1': 123456,
              'seed2': 987654,
              'verbose': '',
              'logfile': ''
              }


    for opt, arg in opts:
        if opt in ["-I"]:
            params["indels"] = float(arg)
        if opt in ["-D"]:
            params["dropout"] = int(arg)
        if opt in ["-L"]:
            params["nLoci"] = int(arg)
        if opt in ["-l"]:
            params["length"] = int(arg)
        if opt in ["-i"]:
            params["Ninds"] = int(arg)
        if opt in ["-N"]:
            params["N"] = int(float(arg))
        if opt in ["-u"]:
            params["mu"] = float(arg)
        if opt in ["-s"]:
            params["insert"] = str(arg)
        if opt in ["-d"]:
            params["depth"] = str(arg)
        if opt in ["-f"]:
            params["datatype"] = str(arg)
        if opt in ["-e"]:
            params["error"] = float(arg)
        if opt in ["-t"]:
            params["tree"] = str(arg)
        if opt in ["-o"]:
            params["outname"] = str(arg)
        if opt in ["-s1"]:
            params["seed1"] = int(arg)
        if opt in ["-s2"]:
            params["seed2"] = int(arg)
        if opt in ["-c1"]:
            params["cut1"] = str(arg)
        if opt in ["-c2"]:
            params["cut2"] = str(arg)
        if opt in ["-v"]:
            params["verbose"] = bool(arg)
        if opt in ["-log"]:
            params["logfile"] = bool(arg)

    return params


def checkopts(params):
    """ print errors if options are incorrect"""

    ## should move tree checking functions here...
    if params["datatype"] in ['rad', 'gbs', 'pairgbs',
                              'ddrad', 'pairddrad']:
        print '\n\tsimulating '+params["datatype"]+" data"
    else:
        print 'datatype not recognized'
        sys.exit()

    ## warning if overlap occurs
    #if int(params["insert"].split(",")[0]) < 200:
    #    if params['datatype'] not in ['rad', 'ddrad']:
    #       print "\tmin fragment length allows reads to overlap "
    if int(params["insert"].split(",")[0]) < 0:
        print "\tmin insert size allows read overlaps/adapter sequences "

    ## which tree to use
    if params["tree"] == "default 12 tip tree": 
        tiptax, paramset = defaulttree(params)
    else:
        ## check the tree string or file
        is_ultrametric(egglib.Tree(fname=params["tree"]))
        tiptax, paramset = usertree(params)
    return tiptax, paramset


def usage():
    """
    brief description of various flags and options for this script
    """
    print "\nHere is how you can use this script\n"
    print "Usage: python %s"%sys.argv[0]
    print "\t -o  <str>   output file prefix name (def 'out')"
    print "\t -I  <float> probability mutation is an indel (default 0.0) "
    print "\t -D  <bool>  allow locus dropout (def 0) "
    print "\t -L  <int>   Number of loci to simulate (def 1000) "
    print "\t -l  <int>   locus length (def 100) "    
    print "\t -u  <float> per-site mutation rate (def 1e-9) "
    print "\t -N  <int>   effective population size (def 1e5) "
    print "\t -i  <int>   individuals sampled per tip taxon (def 1) "
    print "\t -s  <str>   insert size per locus : U(min,max) def=300,800 "
    print "\t -e  <int>   sequencing error rate (def 0.0005) "    
    print "\t -d  <str>   depth of sequencing (def=10,0). N(mean,std) per copy"
    print "\t -t  <str>   tree w br lens (file or text) or use default"
    print "\t -f  <str>   format (def=rad), opts: ddrad,gbs,pairddrad,pairgbs)"
    print "\t -s1 <int>   random seed 1 (def 123456)"
    print "\t -s2 <int>   random seed 2 (def 987654)"
    print "\t -c1 <int>   restriction cut site 1 (def CTGCAG) ->TGCAG"
    print "\t -c2 <int>   restriction cut site 2 (def CAATTG) ->AATTG"
    print "\t -v  <int>   verbose: 1=screen, 2=log, 3=both"
    

def lengthtotip(node):
    """ returns node height """
    dec = 0.
    while node.descendants():
        dec += node.branch_to(node.descendants()[0])
        node = node.descendants()[0]
    return dec


def is_ultrametric(tree):
    """ check if tree is ultrametric """
    tips = tree.get_terminal_nodes()
    dist_to_root = []
    for node in tips:
        leaf = node.get_label()
        edgesums = 0
        while node.ascendants():
            edgesums += node.branch_from(node.ascendants()[0])
            node = node.ascendants()[0]
        dist_to_root.append((leaf, edgesums))
    if len(set([round(i[1], 3) for i in dist_to_root])) == 1:
        return 1
    else:
        print "\nError: tree is not ultrametric"
        for node in dist_to_root:
            print node[0], "dist_to_root:", node[1]
        sys.exit(2)



def defaulttree(params):
    """ supply default tree """
    tiptax = [i+j for i, j in zip(list("1111222233334"),
                                  list("ABCDEFGHIJKLX"))]

    ## simulate species with divergence times"
    # newick = """(((((A:2,B:2):2,C:4):4,D:8):4,
    #             (((E:2,F:2):2,G:4):4,H:8):4):4,
    #             (((I:2,J:2):2,K:4):4,L:8):8,X:16):16;"""

    divscale = 1.0
    nodeABCDEFGHIJKLX  = 2.5*divscale
    nodeABCDEFGHIJKL   = 2.5*divscale
    nodeABCDEFGH       = 2.*divscale
    #nodeIJKL           = 2.*divscale
    #nodeEFGH           = 2.*divscale
    nodeABCD           = 1.5*divscale
    #nodeIJK            = 1.*divscale
    #nodeEFG            = 1.*divscale
    nodeABC            = 1.*divscale
    #nodeIJ             = 0.5*divscale
    #nodeEF             = 0.5*divscale
    nodeAB             = 0.5*divscale

    # sets the two parameter classes
    paramSet = egglib.simul.CoalesceParamSet(singleSamples=None, 
                            doubleSamples=[params["Ninds"]]*len(tiptax),
                            M=0.0)
    ## clade 1"
    paramSet.populationFusion(nodeABCD, 0, 3)     ## D into A
    paramSet.populationFusion(nodeABC, 0, 2)      ## C into A
    paramSet.populationFusion(nodeAB, 0, 1)       ## B into A
    ## clade 2"
    paramSet.populationFusion(nodeABCD, 4, 7)     ## H into E
    paramSet.populationFusion(nodeABC, 4, 6)      ## G into E
    paramSet.populationFusion(nodeAB, 4, 5)       ## F into E
    ## clade 3"
    paramSet.populationFusion(nodeABCD, 8, 11)    ## L into I
    paramSet.populationFusion(nodeABC, 8, 10)     ## K into I
    paramSet.populationFusion(nodeAB, 8, 9)       ## J into I
    ## together and outgroup"
    paramSet.populationFusion(nodeABCDEFGHIJKLX, 0, 12)   ## X into A
    paramSet.populationFusion(nodeABCDEFGHIJKL, 0, 8)     ## I into A
    paramSet.populationFusion(nodeABCDEFGH, 0, 4)         ## E into A
    return tiptax, paramSet


def usertree(params):
    """ simulate data on user input topology """
    ## get tip names
    tree = egglib.Tree(fname=params["tree"])
    tiptax = tree.all_leaves()+["4X"]

    ## sets the two parameter classes
    paramset = egglib.simul.CoalesceParamSet(
                                singleSamples=None,
                                doubleSamples=[params["Ninds"]]*len(tiptax),
                                M=0.0)

    ## add outgroup tip for polarization
    treelength = lengthtotip(tree.root_node())
    tree.add_node(tree.root_node(), label="4X", brlen=treelength)

    ## traverse tree fusing populations at each node
    for node in tree:
        if node.descendants():
            date = lengthtotip(node)
            tip1 = min([tiptax.index(l) for l in \
                        node.descendants()[0].leaves_down()])
            tip2 = min([tiptax.index(l) for l in \
                        node.descendants()[1].leaves_down()])
            paramset.populationFusion(date, tip1, tip2)
    ## fuse outgroup "4X0"
    paramset.populationFusion(treelength, 0, len(tiptax)-1)
    return tiptax, paramset


def twodiffs(bara, barb):
    "requires two base differences between barcodes"
    if len(bara) == len(barb):
        sames = [bara[i] == barb[i] for i in range(len(bara))]
        if sames.count(False) > 1:
            return True


def barcoder(aligns, params, tiptax, outgroup, Barcodes):
    """ takes an align object and puts names to sequences and
    makes a barcode map """

    ## append names"
    names = []
    for name in tiptax:
        for i in range(params["Ninds"]):
            names.append(name+str(i))
            names.append(name+str(i))

    ## appends names for allele 1 vs. allele 2 for each diploid sample "
    for align in aligns:
        for alobj, name in zip(align, names):
            alobj.name = name
    
    ## make dictionary with list of loci for each sample "
    loldic = {sname: [] for sname in set(names)}
    for samp in loldic:
        for align in aligns:
            loc = []
            for seq in align:
                if seq.name == samp:
                    loc.append(seq.sequence)
            loldic[samp].append(loc) 

    ## generate barcodes to append to start of loci "
    over1 = params['cut1'][1:]
    over2 = params['cut2'][1:]

    bases = list("ATGC")
    if not Barcodes:
        print "\tcreating new barcode map"
        Barcodes[names[0]] = "CATCAT"
        for name in names:
            while name not in Barcodes.keys():
                bbs = np.random.randint(0, 3, 6)
                bcd = "".join([bases[i] for i in bbs])
                if all([twodiffs(bcd, bb) for bb in Barcodes.values()]):
                    if not any([i in bcd for i in [over1, over2]]):
                        Barcodes[name] = bcd
        ## creates random barcodes and writes map to file "
        with open(params["outname"]+".barcodes", 'w') as barout:
            bnames = list(Barcodes.keys())
            bnames.sort()
            for bcd in bnames:
                if outgroup not in bcd:
                    print >>barout, "\t".join([bcd, Barcodes[bcd]])
    return loldic, Barcodes



def mutate(base):
    "introduce point sequencing errors"
    nbase = list("ATGC")[np.random.randint(0, 4)]
    while nbase == base:
        nbase = list("ATGC")[np.random.randint(0, 4)]
    return nbase


def revcomp(seq):
    """ retrurn the reverse complement of a seq """
    rcseq = seq.replace('A', 't').replace('T', 'a').\
                replace('C', 'g').replace('G', 'c')
    return rcseq[::-1].upper()


def locus_dropout(datatype, ingroupseq, outgroupseq, cut1, cut2):
    """ check if cut site occurs in the sequence fragment with the 
    cut sites determined by data type """
    muts1 = [i for i, j in zip(ingroupseq[:len(cut1)], 
                              outgroupseq[:len(cut1)]) if i != j]
    muts2 = [i for i, j in zip(ingroupseq[-len(cut2):], 
                              outgroupseq[-len(cut2):]) if i != j]
    muts3 = [i for i, j in zip(ingroupseq[-len(cut1):], 
                              outgroupseq[-len(cut1):]) if i != j]
    if datatype == 'rad':
        drop1 = bool(any([i in ingroupseq for i in cut1, revcomp(cut1)]))
        drop2 = bool(muts1) 
        return any([drop1, drop2])
    if 'gbs' in datatype:
        drop1 = bool(any([i in ingroupseq for i in cut1, revcomp(cut1)]))
        drop2 = any([muts1, muts3])
        return any([drop1, drop2])
    if 'ddrad' in datatype:
        drop1 = bool(any([i in ingroupseq for i in cut1, revcomp(cut1),
                                             cut2, revcomp(cut2)])) 
        drop2 = any([muts1, muts2]) 
        return any([drop1, drop2])

# def locus_dropout2(datatype, ingroupseq, outgroupseq, cut1, cut2):
#     """ check if mutation occurs in cut site """
#     if datatype == 'rad':
#         return 
#     if 'gbs' in datatype:
#         return 
#     if 'ddrad' in datatype:
#         return 


def dressitup(loldic, Barcodes, params,
              locuslength, counter):
    """ return a list of sequences all dressed up """
    seqs1 = []
    seqs2 = []

    ## parse allele depth sampling
    meandepth, stdevdepth = params["depth"].split(",")
    if stdevdepth == "0":
        stdevdepth = 1e-9

    for loc in range(100):
        ## sample fragment in size selection window for this locus
        frag1, frag2 = params["insert"].split(',')
        insert = np.random.randint(int(frag1), int(frag2))
        frag = (2*params["length"])+insert

        ## do we keep this locus in the outgroup sample?
        if not locus_dropout(params["datatype"],
                             loldic['4X0'][loc][0][:frag],
                             loldic['4X0'][loc][0][:frag],
                             params['cut1'], params['cut2']):
            counter += 1

            if params["indels"]:
                skip = 0
                ## indel = prob. site is checked for indel-mut, skip 'skip' bp"
                inds = np.random.binomial(locuslength-skip, params["indels"])
                ## for each indel site randomly choose location on locus"
                wh = np.random.randint(skip, locuslength, inds)
                ## randomly select size of indel (default size 1)"
                le = np.random.randint(1, 2, len(wh))
                where = [wh, le]

            for samp in loldic:
                ## sample depth of sequencing at this locus    
                copy1 = int(round(np.random.normal(float(meandepth),
                                                  float(stdevdepth))))
                copy2 = int(round(np.random.normal(float(meandepth), 
                                                  float(stdevdepth))))
                ## exclude outgroup taxon used for detecting restriction site
                ## mutations and indels use both chromosomes 
                if params["dropout"]:
                    ## remove locus if mutation arises in restriction site"
                    passed1 = [list(loldic[samp][loc][0]) for \
                               copy in range(copy1) if not \
                               locus_dropout(params["datatype"], 
                                        loldic[samp][loc][0][:frag],
                                        loldic['4X0'][loc][0][:frag],                                        
                                        params['cut1'], params['cut2'])]

                    passed2 = [list(loldic[samp][loc][1]) for \
                               copy in range(copy2) if not \
                               locus_dropout(params["datatype"],
                                        loldic[samp][loc][1][:frag],
                                        loldic['4X0'][loc][0][:frag],                                        
                                        params['cut1'], params['cut2'])]
                else:
                    passed1 = [list(loldic[samp][loc][0]) for \
                                    copy in range(copy1)] 
                    passed2 = [list(loldic[samp][loc][1]) for \
                                    copy in range(copy2)]

                allpass = passed1+passed2

                for copy in range(len(allpass)):
                    nowlength = locuslength
                    ## introduce indel if mutation present at site
                    ## relative to the outgroup
                    if params["indels"]:
                        for iloc, isize in zip(where[0], where[1]):
                            difflength = locuslength-nowlength
                            iloc -= difflength
                            if allpass[copy][iloc] != loldic['4X0'][loc][0][iloc]:
                                allpass[copy] = allpass[copy][:iloc] +\
                                                allpass[copy][iloc+isize:]
                                nowlength -= isize

                    ## introduce sequencing errors"
                    errors = np.random.binomial(nowlength, params["error"])
                    errlocs = np.random.randint(0, nowlength, errors)
                    thiscopy = allpass[copy]
                    for error in errlocs:
                        thiscopy[error] = mutate(thiscopy[error])

                    ## shorten fragment to potentially overlapping length"
                    ## 0------->insert<-------[frag]
                    thiscopy = thiscopy[:frag]

                    seqs1, seqs2, counter = stacklist(thiscopy, params, 
                                                Barcodes, loldic, copy, 
                                                counter, samp, seqs1, seqs2)
    return seqs1, seqs2, counter



def stacklist(thiscopy, params, Barcodes,
               loldic, copy, counter,
               samp, seqs1, seqs2):
    """ returns seqs formatted for the given data type """

    illumina_p1 = "ACGACGCTCTTCCGATCT"
    illumina_p2 = "AGATCGGAAGAGCTCGTATG"
    over1 = params['cut1'][1:]
    over2 = params['cut2'][1:]

    if params["datatype"] in ['ddrad', 'pairddrad']:
        thiscopy1 = ("A"*30)+illumina_p1+Barcodes[samp]+over1+\
                    "".join(thiscopy)+revcomp(over2)+illumina_p2+"A"*30
    else:
        thiscopy1 = ("A"*30)+illumina_p1+Barcodes[samp]+over1+\
                     "".join(thiscopy)+revcomp(over1)+illumina_p2+"A"*30

    if params["datatype"] == 'gbs':
        thiscopy2 = ("A"*30)+illumina_p1+Barcodes[samp]+over1+\
                    revcomp("".join(thiscopy))+revcomp(over1)+\
                    illumina_p2+"A"*30
                    
    ## sequence read1 from 5' end "
    startleft = thiscopy1.index("CGATCT")
    ## ---->
    here = startleft+len("CGATCT")
    sss = thiscopy1[here:here+100]

    if "X" not in samp:
        if counter <= params["nLoci"]:
            seqs1.append("@lane1_locus"+str(counter)+"_"+\
                         str(loldic.keys().index(samp))+"_R1_"+str(copy)+\
                        " 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

            ## add the second orientation of the read...
            if params["datatype"] == 'gbs':
                ## sequence fragment from 3' end for GBS
                startleft = thiscopy2.index("CGATCT")
                ## ---->
                here = startleft+len("CGATCT")
                sss = thiscopy2[here:here+100]
                if "X" not in samp:
                    #if counter <= params["nLoci"]:
                    seqs1.append("@lane1_fakedata"+str(counter)+\
                                 "_"+str(loldic.keys().index(samp))+\
                                 "_R1_"+str(copy)+" 1:N:0:"+"\n"+\
                                 sss+"\n+\n"+("B"*len(sss))+"\n")

            ## sequence read2 from 3' end for paired read "
            if 'pair' in params['datatype']:
                startright = thiscopy1.rindex("AGATCG")
                ## <----
                here = len(thiscopy1)-startright
                sss = thiscopy1[-1*(here+100):-here]
                if "X" not in samp:
                    #if counter <= params["nLoci"]:
                    seqs2.append("@lane1_locus"+str(counter)+"_"+\
                                str(loldic.keys().index(samp))+\
                               "_R2_"+str(copy)+" 1:N:0:"+"\n"+\
                                revcomp(sss)+"\n+\n"+("B"*len(sss))+"\n")

            if params["datatype"] == 'pairgbs':
                ##  double check this 
                ##  sequence read2 from 5' end for GBS"
                startleft = thiscopy2.index("CGATCT")
                ## ---->
                here = startleft+len("CGATCT")
                sss = thiscopy2[here:here+100]
                if "X" not in samp:
                    #if counter <= params["nLoci"]:
                    seqs2.append("@lane1_locus"+str(counter)+"_"+\
                                str(loldic.keys().index(samp))+\
                                "_R1_"+str(copy)+" 1:N:0:"+"\n"+\
                                sss+"\n+\n"+("B"*len(sss))+"\n")
    return seqs1, seqs2, counter



def createlog(params, tiptaxa):
    """ prints log file to log (and to screen) """

    if params["verbose"]:
        print "\tSequencing error rate = 0.0005 "
        print "\tNloci =", params["nLoci"]
        print "\tDepth per allele = N("+params["depth"]+")"
        print "\t"+str(len(tiptaxa)-1), "taxa with", params["Ninds"], "samples/taxon"
        print "\tTopology =", params["tree"]
        print "\tIndels arise at frequency of", params["indels"]
        print "\tLocus dropout =", bool(params["dropout"])
        maxfrag = int(params["insert"].split(',')[1])+(2*params["length"])
        print "\tInsert range=("+params["insert"]+"),", "maxfrag=", maxfrag


if __name__ == "__main__":

    ## parse command line options
    ARGS = sys.argv[1:]
    SMALLFLAGS = "I:D:L:u:N:i:e:s:d:t:f:o:s1:s2:c1:c2:v:l:"
    BIGFLAGS = []
    try:
        OPTS, ARGS = getopt.getopt(ARGS, SMALLFLAGS, BIGFLAGS)
        if not OPTS:
            usage()
            sys.exit(2)
    except getopt.GetoptError:
        print "Incorrect options passed"
        usage()
        sys.exit()

    ## get params
    PARAMS = parseopts(OPTS)
    TIPTAX, PARAMSET = checkopts(PARAMS)

    ## overhang from cuts at Pst1 (and) at EcoR1
    OVER1 = PARAMS['cut1'][1:]
    OVER2 = PARAMS['cut2'][1:]
    OUTGROUP = "4X0"

    ## get theta
    LOCUSLENGTH = int(PARAMS["insert"].split(',')[1])+(2*PARAMS["length"])
    THETA = 4.*PARAMS["N"]*PARAMS["mu"]*LOCUSLENGTH
    print "\tTHETA=", THETA/LOCUSLENGTH

    ## setup mutator
    MUTATOR = egglib.simul.CoalesceFiniteAlleleMutator(theta=THETA,
                                                      alleles=4,
                                                      randomAncestralState=True)
    MUTATOR.setSites(LOCUSLENGTH)

    ## create log file and print to screen
    createlog(PARAMS, TIPTAX)
 
    ## simulate data
    OUT1 = gzip.open(PARAMS["outname"]+"_R1_.fastq.gz", 'wb')
    if PARAMS["datatype"] in ['pairddrad', 'pairgbs']:
        OUT2 = gzip.open(PARAMS["outname"]+"_R2_.fastq.gz", 'wb')

    ## b/c single cutter will sequence fragment either way
    if PARAMS["datatype"] in ['gbs']:
        PARAMS["nLoci"] = PARAMS["nLoci"]/2   

    BARCODES = {}
    COUNTER = 0
    while COUNTER < PARAMS["nLoci"]: 
        SEED1 = PARAMS["seed1"] * (COUNTER+1)
        SEED2 = PARAMS["seed2"] * (COUNTER+1)
        ## performs simulation of 100 loci at a time"
        ALIGNS = egglib.simul.coalesce(PARAMSET, MUTATOR, 100,
                                       random=(SEED1, SEED2))

        ## creates barcodes if not made yet and fill D dict."
        LOLDIC, BARCODES = barcoder(ALIGNS, PARAMS, TIPTAX, 
                                    OUTGROUP, BARCODES)

        ## dresses up data to be fastq and puts in errors, indels, etc."
        SEQS1, SEQS2, COUNTER = dressitup(LOLDIC,
                                          BARCODES, 
                                          PARAMS,
                                          LOCUSLENGTH, 
                                          COUNTER)

        OUT1.write("".join(SEQS1))
        if PARAMS["datatype"] in ['pairddrad', 'pairgbs']:
            OUT2.write("".join(SEQS2))

    OUT1.close()
    if PARAMS["datatype"] in ['pairddrad', 'pairgbs']:
        OUT2.close()
