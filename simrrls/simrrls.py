#!/usr/bin/env python2

"""
Script to simulate RADseq-like data 
"""

import sys
import numpy as np
import egglib
import gzip
# pylint: disable=E1101



def checkopts(params):
    """ print errors if options are incorrect"""

    ## should move tree checking functions here...
    if params.datatype not in ['rad', 'gbs', 'pairgbs',
                               'ddrad', 'pairddrad']:
        sys.exit('\n\tdatatype not recognized')

    ## warning
    if params.min_insert < 0:
        print "\tmin insert size allows read overlaps/adapter sequences "

    ## which tree to use
    if not params.tree:
        tiptax, paramset = defaulttree(params)
    else:
        ## check the tree string or file
        is_ultrametric(egglib.Tree(fname=params["tree"]))
        tiptax, paramset = usertree(params)
    return tiptax, paramset



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

    # """ simulate species with divergence times
    # newick = (((((A:2,B:2):2,C:4):4,D:8):4,
    #          (((E:2,F:2):2,G:4):4,H:8):4):4,
    #          (((I:2,J:2):2,K:4):4,L:8):8,X:16):16;"""

    divscale = 1.0
    nodeABCDEFGHIJKLX  = 2.5*divscale
    nodeABCDEFGHIJKL   = 2.5*divscale
    nodeABCDEFGH       = 2.*divscale
    nodeABCD           = 1.5*divscale
    nodeABC            = 1.*divscale
    nodeAB             = 0.5*divscale

    # sets the two parameter classes
    paramset = egglib.simul.CoalesceParamSet(
                            singleSamples=None, 
                            doubleSamples=params.Ninds*len(tiptax),
                            M=0.0)
                            #N=[int(params["N"])]*(len(tiptax)-1)+[1e6])
    ## clade 1"
    paramset.populationFusion(nodeABCD, 0, 3)     ## D into A
    paramset.populationFusion(nodeABC, 0, 2)      ## C into A
    paramset.populationFusion(nodeAB, 0, 1)       ## B into A
    ## clade 2"
    paramset.populationFusion(nodeABCD, 4, 7)     ## H into E
    paramset.populationFusion(nodeABC, 4, 6)      ## G into E
    paramset.populationFusion(nodeAB, 4, 5)       ## F into E
    ## clade 3"
    paramset.populationFusion(nodeABCD, 8, 11)    ## L into I
    paramset.populationFusion(nodeABC, 8, 10)     ## K into I
    paramset.populationFusion(nodeAB, 8, 9)       ## J into I
    ## together and outgroup"
    paramset.populationFusion(nodeABCDEFGHIJKLX, 0, 12)   ## I into A    
    paramset.populationFusion(nodeABCDEFGHIJKL, 0, 8)     ## I into A
    paramset.populationFusion(nodeABCDEFGH, 0, 4)         ## E into A
    return tiptax, paramset



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
            while name not in Barcodes:
                bbs = np.random.randint(0, 3, 6)
                bcd = "".join([bases[i] for i in bbs])
                if all([twodiffs(bcd, bb) for bb in Barcodes.values()]):
                    if not any([i in bcd for i in [over1, over2]]):
                        Barcodes[name] = bcd
        ## creates random barcodes and writes map to file "
        with open(params["outname"]+"_barcodes.txt", 'w') as barout:
            bnames = list(Barcodes)
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



def dressitup(loldic, barcodes, params,
              locuslength, counter, stepsize):

    """ return a list of sequences all dressed up """
    seqs1 = []
    seqs2 = []

    ## parse allele depth sampling
    meandepth, stdevdepth = params["depth"].split(",")
    if stdevdepth == "0":
        stdevdepth = 1e-9

    for loc in range(stepsize):
        ## sample fragment in size selection window for this locus
        frag1, frag2 = params["insert"].split(',')
        insert = np.random.randint(int(frag1), int(frag2))
        frag = (2*params["length"])+insert    
        #print 'fragment length', frag

        ## make sure locus passes filter in the ancestor/outgroup sample?
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
                #print loldic[samp][loc][0][:20], 'what'
                ## sample depth of sequencing at this locus    
                # pylint: disable=E1101
                copy1 = int(round(np.random.normal(
                            float(meandepth), float(stdevdepth))))
                copy2 = int(round(np.random.normal(
                            float(meandepth), float(stdevdepth))))
                ## exclude outgroup taxon used for detecting restriction site
                ## mutations and indels use both chromosomes 
                if params["dropout"]:
                    #print loldic[samp][loc][0][:10]
                    #print loldic['4X0'][loc][0][:10] 

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
                #print len(allpass)

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
                    #print 'frag', "".join(thiscopy[:20]), '...', "".join(thiscopy[frag-20:frag])
                    thiscopy = thiscopy[:frag]

                    ## formats reads for the appropriate data type
                    seqs1, seqs2, counter = stacklist(thiscopy, params, 
                                                Barcodes, loldic, copy, 
                                                counter, samp, seqs1, seqs2)
    return seqs1, seqs2, counter



def stacklist(thiscopy, params, Barcodes,
               loldic, copy, counter,
               samp, seqs1, seqs2):
    """ returns seqs formatted for the given data type """

    ## start of adapters/primers
    illumina_p1 = "ACGACGCTCTTCCGATCT"
    illumina_p2 = "AGATCGGAAGAGCTCGTATG"

    ## assume cutter have one base cut
    over1 = params['cut1'][1:]
    over2 = params['cut2'][1:]

    ## embed seq data beween adapters and with barcodes
    ## and restriction sites to look like raw data
    if params["datatype"] in ['ddrad', 'pairddrad']:
        thiscopy1 = ("A"*30)+illumina_p1+Barcodes[samp]+over1+\
                    "".join(thiscopy)+revcomp(over2)+illumina_p2+"A"*30
    else:
        thiscopy1 = ("A"*30)+illumina_p1+Barcodes[samp]+over1+\
                     "".join(thiscopy)+revcomp(over1)+illumina_p2+"A"*30

    if 'gbs' in params["datatype"]: #  == 'gbs':
        thiscopy2 = ("A"*30)+illumina_p1+Barcodes[samp]+over1+\
                    revcomp("".join(thiscopy))+revcomp(over1)+\
                    illumina_p2+"A"*30

    ## print out reminder of the read structure                    
    # print " ".join(['1', thiscopy1[30:48],
    #                      thiscopy1[48:54],
    #                      thiscopy1[54:59],
    #                      thiscopy1[59:66],                         
    #                      '...',
    #                      thiscopy1[-62:-55],                         
    #                      thiscopy1[-55:-50],
    #                      thiscopy1[-50:-28]])

    # print " ".join(['2', thiscopy2[30:48],
    #                      thiscopy2[48:54],
    #                      thiscopy2[54:59],
    #                      thiscopy2[59:66],                         
    #                      '...',
    #                      thiscopy2[-62:-55],                         
    #                      thiscopy2[-55:-50],
    #                      thiscopy2[-50:-28]])

    if "4X" not in samp:
        if counter <= params["nLoci"]:
            ## always add to seqs1
            ## sequence read1 from 5' end "
            ## ----> start from primer
            start = thiscopy1.index("CGATCT")+len("CGATCT")
            sss = thiscopy1[start:start+params['length']]
            seqs1.append("@lane1_locus"+str(counter)+"_"+\
                         str(loldic.keys().index(samp))+"_R1_"+str(copy)+\
                        " 1:N:0:"+"\n"+sss+"\n+\n"+("B"*len(sss))+"\n")

            ## if gbs add fragment2 to the list of first reads
            if 'gbs' in params["datatype"]:
                ## sequence fragment from 3' end for GBS
                start = thiscopy2.index("CGATCT")+len("CGATCT")
                sss = thiscopy2[start:start+params["length"]]
                seqs1.append("@lane1_clocus"+str(counter)+\
                             "_"+str(loldic.keys().index(samp))+\
                             "_R1_"+str(copy)+" 1:N:0:"+"\n"+\
                             sss+"\n+\n"+("B"*len(sss))+"\n")

            ## sequence read2 from 3' end for paired read "
            if 'pair' in params['datatype']:
                start = thiscopy1.rindex("AGATCG")
                ## <----
                sss = thiscopy1[start-params['length']:start]
                #sss = thiscopy1[-1*(start+params['length']):-start]
                seqs2.append("@lane1_locus"+str(counter)+"_"+\
                            str(loldic.keys().index(samp))+\
                            "_R2_"+str(copy)+" 1:N:0:"+"\n"+\
                            revcomp(sss)+"\n+\n"+("B"*len(sss))+"\n")

            ## add both R1 and R2 to second reads
            ## r2 was just added above, so here just add R1
            if params["datatype"] == 'pairgbs':
                ##  double check this 
                ##  sequence read2 from 5' end for GBS"
                #startleft = thiscopy2.index("CGATCT")
                ## ---->
                #here = startleft+len("CGATCT")
                #sss = thiscopy2
                #sss = thiscopy2[-1*(start+params['length']):-start]                
                sss = thiscopy2[start-params['length']:start]                
                #sss = thiscopy2[here:here+params['length']]
                #if counter <= params["nLoci"]:
                seqs2.append("@lane1_clocus"+str(counter)+"_"+\
                            str(loldic.keys().index(samp))+\
                            "_R2_"+str(copy)+" 1:N:0:"+"\n"+\
                            revcomp(sss)+"\n+\n"+("B"*len(sss))+"\n")
    return seqs1, seqs2, counter



def run(params):
    ## get params
    tiptax, paramset = checkopts(params)

    ## overhang from cuts at Pst1 (and) at EcoR1
    outgroup = "4X0"

    ## get theta
    ## all data are simulated at the max fragment size of window
    locuslength = params.min_insert+(2*params.length)
    theta = 4.*params.N*params.mu*locuslength
    print("\tTHETA=", theta/locuslength)

    ## initialize random number sampling for numpy
    np.random.seed(params.seed1)  # pylint: disable=E1101

    ## setup mutator
    mutator = egglib.simul.CoalesceFiniteAlleleMutator(
                    theta=theta,
                    alleles=4,
                    randomAncestralState=True)
    mutator.setSites(locuslength)

    ## create log file and print to screen
    #createlog(params, tiptax)
 
    ## open files to write seq data to
    out1 = gzip.open(params.outfile+"_R1_.fastq.gz", 'wb')
    if 'pair' in params.datatype:
        out2 = gzip.open(params.outfile+"_R2_.fastq.gz", 'wb')

    ## b/c single cutter will sequence fragment either way
    if params.datatype == 'gbs':
        params.nloci = params.nloci/2

    ## set stepsize for simulation depending on Nloci
    if params.nloci >= 100:
        stepsize = 100
    elif params.nloci >= 1000:
        stepsize = 500
    elif params.nloci >= 2000:
        stepsize = 1000
    else:
        stepsize = 20

    ## simulate the data
    barcodes = {}
    counter = 0
    while counter < params.nloci: 
        ##  this is a problem if cut is too frequent..
        ## seed has to change each iteration in a known way
        localseed1 = params.seed1 * (counter+1)
        localseed2 = params.seed2 * (counter+1)        

        ## performs simulation of setpsize loci at a time"
        aligns = egglib.simul.coalesce(
                    paramset, mutator, stepsize,
                    random=(localseed1, localseed2))

        ## creates barcodes if not made yet and fill D dict."
        loldic, barcodes = barcoder(aligns, params, tiptax, outgroup, barcodes)

        ## dresses up data to be fastq and puts in errors, indels, etc."
        seqs1, seqs2, counter = dressitup(loldic, barcodes, params,
                                          locuslength, counter, stepsize)

        out1.write("".join(seqs1))
        if params["datatype"] in ['pairddrad', 'pairgbs']:
            out2.write("".join(seqs2))

    out1.close()
    if params["datatype"] in ['pairddrad', 'pairgbs']:
        out2.close()

