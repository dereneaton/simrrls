#!/usr/bin/env python2

"""
Script to simulate RADseq-like data 
"""

import sys
import numpy as np
import egglib
import gzip
# pylint: disable=E1101
# pylint: disable=W0142


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
        is_ultrametric(egglib.Tree(fname=params.tree))
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
    tiptax = [i+j for i, j in zip(list("111122223333")+["OUT"],
                                  list("ABCDEFGHIJKL")+[""])]

    # """ simulate species with divergence times
    # newick = (((((A:2,B:2):2,C:4):4,D:8):4,
    #          (((E:2,F:2):2,G:4):4,H:8):4):4,
    #          (((I:2,J:2):2,K:4):4,L:8):8,X:16):16;"""

    divscale = 1.0
    node_abcdefghijkl = 2.5*divscale
    node_abcdefgh = 2.*divscale
    node_abcd = 1.5*divscale
    node_abc = 1.*divscale
    node_ab = 0.5*divscale

    # sets the two parameter classes
    paramset = egglib.simul.CoalesceParamSet(
                            singleSamples=None, 
                            doubleSamples=[params.Ninds]*len(tiptax),
                            M=0.0)
    ## clade 1"
    paramset.populationFusion(node_abcd, 0, 3)     ## D into A
    paramset.populationFusion(node_abc, 0, 2)      ## C into A
    paramset.populationFusion(node_ab, 0, 1)       ## B into A
    ## clade 2"
    paramset.populationFusion(node_abcd, 4, 7)     ## H into E
    paramset.populationFusion(node_abc, 4, 6)      ## G into E
    paramset.populationFusion(node_ab, 4, 5)       ## F into E
    ## clade 3"
    paramset.populationFusion(node_abcd, 8, 11)    ## L into I
    paramset.populationFusion(node_abc, 8, 10)     ## K into I
    paramset.populationFusion(node_ab, 8, 9)       ## J into I
    ## together and outgroup"
    paramset.populationFusion(node_abcdefghijkl, 0, 12)     ## X into A    
    paramset.populationFusion(node_abcdefghijkl, 0, 8)     ## I into A
    paramset.populationFusion(node_abcdefghijkl, 0, 8)     ## I into A
    paramset.populationFusion(node_abcdefgh, 0, 4)         ## E into A
    return tiptax, paramset



def usertree(params):
    """ simulate data on user input topology """
    ## get tip names
    tree = egglib.Tree(fname=params.tree)
    tiptax = tree.all_leaves()+["OUT"]

    ## sets the two parameter classes
    paramset = egglib.simul.CoalesceParamSet(
                                singleSamples=None,
                                doubleSamples=[params.Ninds]*len(tiptax),
                                M=0.0)

    ## add outgroup tip for polarization
    treelength = lengthtotip(tree.root_node())
    tree.add_node(tree.root_node(), label="OUT", brlen=treelength)

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
        if sames.count(False) > 2:
            return True



def namer(aligns, params, tiptax):
    """ takes an align object and puts names to sequences and
    makes a barcode map if not already made"""

    ## append names for N individuals if >1
    names = []
    for name in tiptax:
        for i in range(params.Ninds):
            names.append(name+"_"+str(i))
            names.append(name+"_"+str(i))

    ## appends names for allele 1 vs. allele 2 for each diploid sample "
    for align in aligns:
        for alobj, name in zip(align, names):
            alobj.name = name
    return aligns, names



def barcoder(names, params, barcodes):
    """ dlll"""
    ## generate barcodes to append to start of loci "
    over1 = params.cut1[1:]
    over2 = params.cut2[1:]

    bases = list("ATGC")
    barcodes[names[0]] = "CATCATCAT"
    for name in names:
        while name not in barcodes:
            bbs = np.random.randint(0, 4, 9)
            bcd = "".join([bases[i] for i in bbs])
            if all([twodiffs(bcd, bb) for bb in barcodes.itervalues()]):
                if not any([i in bcd for i in [over1, over2]]):
                    barcodes[name] = bcd
    ## creates random barcodes and writes map to file "
    if params.outname:
        outname = params.outname+"_barcodes.txt"
    else:
        outname = "out_barcodes.txt"
    with open(outname, 'w') as barout:
        bnames = list(barcodes)
        bnames.sort()
        for bcd in bnames:
            if "OUT_" not in bcd:
                print >>barout, "\t".join([bcd, barcodes[bcd]])
    return barcodes



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



def mutation_new_cut(params, aligns, dropcheck):
    """ check if cut site occurs in the sequence fragment with the 
    cut sites determined by data type """
    ## check if restriction site is in the sequence
    cutlist1 = [params.cut1, revcomp(params.cut1)]
    cutlist2 = [params.cut2, revcomp(params.cut2)]                    

    tt1 = 0
    tt2 = 0
    keepgrp = []
    for locus in range(len(aligns)):
        ## get size for this locus
        #insert = np.random.randint(params.min_insert, 
        #                           params.max_insert)
        #frag = (2*params.length)+insert

        if aligns[locus]:
            ## get the outgroup seq at this locus
            #outseq = aligns[locus].sequenceByName("OUT_0")
            keeps = []
            ## iterate over ingroup samples
            for haplo in range(len(aligns[locus])):
                if "OUT_" not in aligns[locus][haplo][0]:
                    ## get this ingroup seq
                    inseq = aligns[locus][haplo][1]
                    drop = dropcheck[locus][haplo][1]
                    hits1 = 0
                    hits2 = 0
                    ## find occurrence of cut site
                    where = np.array([inseq.find(cut) for cut \
                                                      in cutlist1])
                    if any(where > 0):
                        #print where, '1'                        
                        check1 = drop[20:20+len(params.cut1)]
                        ## if mutation occurred to give rise to this cut site
                        if len(set(check1)) > 1:
                            ## if new frag is less than min frag length
                            if any(where > 0):
                                hits1 = 1

                    if 'ddrad' in params.datatype:
                        where = np.array([inseq.find(cut) for cut \
                                                      in cutlist2])
                        if any(where > 0):
                            #print where, '2'
                            check2 = drop[-20:-20+len(params.cut2)]
                            if len(set(check2)) > 1:
                                #print where, check2                                
                                hits2 += 1

                    if not hits1:
                        if not hits2:
                            keeps.append(aligns[locus][haplo])
                        else:
                            tt2 += 1
                    else:
                        tt1 += 1
                #else:
                #    keeps.append(aligns[locus][haplo])
            if keeps:
                keepgrp.append(egglib.Align.create(keeps))
            else:
                keepgrp.append([])
    #print tt1, tt2, 'dropped'
    return keepgrp



def mutation_in_cut(params, aligns, dropcheck):
    """ uses drop mutations to remove allelic dropout """
    keepgrp = []

    ## get length of area to check for mutation
    lentocheck = len(params.cut1)
    if 'gbs' in params.datatype:
        lentocheck += len(params.cut1)
    if 'ddrad' in params.datatype:
        lentocheck += len(params.cut2)

    ## iterate over loci
    for locus in range(len(dropcheck)):
        keeps = []
        outg = []
        for haplo in range(len(dropcheck[locus])):
            if "OUT_" not in aligns[locus][haplo][0]:
                check = dropcheck[locus][haplo][1][:lentocheck]
                if len(set(check)) == 1:
                    keeps.append(aligns[locus][haplo])
            else:
                ## always keep the outgroup b/c needed for other mut
                outg.append(aligns[locus][haplo])
        if keeps:
            keepgrp.append(egglib.Align.create(keeps+outg))
        else:
            keepgrp.append([])
    return keepgrp



def seq_copies(aligns, barcodes, params, counter, stepsize):
    """ return a list of sequences all dressed up """
    seqs1 = []
    seqs2 = []

    ## parse allele depth sampling
    if params.depthstd == 0:
        params.depthstd = 1e-9

    ## iterate over each locus 
    for loc in aligns:
        ## sample fragment in uniform size window for this locus
        insert = np.random.randint(params.min_insert, params.max_insert)
        #insert = np.random.randint(0, params.max_insert - params.min_insert)
        frag = (2*params.length)+insert

        ## make sure all reads did not get disrupted
        if loc:
            #print aligns[loc], "ALIGNS"
            
            ## make indels, sample copies and introduce seq errors
            if params.indels:
                skip = 0
                ## indel = prob. site is checked for indel-mut, skip 'skip' bp"
                inds = np.random.binomial(frag-skip, params.indels)
                ## for each indel site randomly choose location on locus"
                iwhere = np.random.randint(skip, frag, inds)
                ## randomly select size of indel (default size 1)"
                ilength = np.random.randint(1, 2, len(iwhere))
                where = [iwhere, ilength]

            reads = {}
            for samp in loc:
                tempseq = samp.sequence

                ## introduce indel if mutation at site relative to outgroup
                if params.indels:
                    outseq = loc.sequenceByName("OUT_0")
                    for iloc, isize in zip(where[0], where[1]):
                        ## both tempseq and outseq must have iloc index
                        if len(tempseq) > iloc:
                            try:
                                if tempseq[iloc] != outseq[iloc]:
                                    tempseq = tempseq[:iloc] + tempseq[iloc+isize:] 
                            except IndexError:
                                print iloc, isize, where
                                print len(tempseq), tempseq
                                print len(outseq), outseq
                                print ""

                ## get copies
                if params.depthfunc == 'norm':
                    allele1copies = int(round(np.random.normal(
                                        float(params.depthmean),
                                        float(params.depthstd))))
                    allele2copies = int(round(np.random.normal(
                                        float(params.depthmean),
                                        float(params.depthstd))))
                else:
                    ## TODO: exponential function 
                    pass

                ## build dict
                if samp.name in reads:
                    reads[samp.name] += [tempseq]*allele2copies
                else:                    
                    reads[samp.name] = [tempseq]*allele1copies

                ## introduce sequencing errors into copies
                for copy in range(len(reads[samp.name])):
                    coplen = len(reads[samp.name][copy])
                    errors = np.random.binomial(coplen, params.error)
                    errlocs = np.random.randint(0, coplen, errors)
                    for err in errlocs:
                        ## also do final read size trim after mutate
                        tochange = list(reads[samp.name][copy])
                        tochange[err] = mutate(tochange[err])
                        reads[samp.name][copy] = "".join(tochange)
                    ## shorten fragment to potentially overlapping length"
                    ## 0------->insert<-------[frag]
                    reads[samp.name][copy] = reads[samp.name][copy][:frag]
                    #print len(reads[samp.name][copy][:frag]), 'gotcha', frag

            if counter < params.nLoci:
                ## formats reads for the appropriate data type
                seqs1, seqs2, counter = stacklist(params, reads, barcodes, 
                                                  counter, seqs1, seqs2)
                counter += 1                    
        ## skipped
    return seqs1, seqs2, counter



def stacklist(params, reads, barcodes, counter, seqs1, seqs2):
    """ returns seqs formatted for the given data type """

    ## start of adapters/primers
    #illumina_p1 = "ACGACGCTCTTCCGATCT"
    illumina_p1 = "ACGACGCTCTTCCGATCT"
    illumina_p2 = "AGATCGGAAGAGCTCGTATG"

    ## assume cutter has one base overhang
    over1 = params.cut1[1:]
    over2 = params.cut2[1:]

    for key, copies in reads.iteritems():
        if "OUT_" not in key:
            for copy in range(len(copies)):
                ## embed seq data beween adapters and with barcodes
                ## and restriction sites to look like raw data
                if 'ddrad' in params.datatype:
                    reads[key][copy] = ("A"*30)+illumina_p1+\
                    barcodes[key]+over1+\
                    copies[copy]+revcomp(over2)+illumina_p2+("A"*30)
                else:
                    reads[key][copy] = ("A"*30)+illumina_p1+\
                    barcodes[key]+over1+\
                    copies[copy]+revcomp(over1)+illumina_p2+("A"*30)
                ## add reverse read direction for gbs
                if 'gbs' in params.datatype:
                    reads[key].append(("A"*30)+illumina_p1+\
                    barcodes[key]+over1+\
                    revcomp(copies[copy])+revcomp(over1)+illumina_p2+"A"*30)

                ## sequence read1 from 5' end "
                ## ----> start from primer
                start = reads[key][copy].index("CGATCT")+len("CGATCT")
                sss = reads[key][copy][start:start+params.length]
                seqs1.append("@lane1_locus{}_{}_{} 1:N:0\n{}\n+\n{}\n"\
                             .format(str(counter), key, str(copy),
                                     sss,
                                     "B"*len(sss)))

                ## sequence read2 from 3' end for paired read "
                if 'pair' in params.datatype:
                    start = reads[key][copy].rindex("AGATCG")
                    ## <----
                    sss = reads[key][copy][start-params.length:start]
                    seqs2.append("@lane1_locus{}_{}_{} 2:N:0\n{}\n+\n{}\n"\
                                 .format(str(counter), key, str(copy),
                                         revcomp(sss),
                                         "B"*len(sss)))
                                         
    return seqs1, seqs2, counter



def run(params):
    """ the main function for calling egglib, samplers and filters """
    ## get params
    tiptax, paramset = checkopts(params)

    ## data are initially simulated at the min fragment size of window
    locuslength = params.max_insert+(2*params.length)
    theta = 4.*int(params.N)*params.mu*locuslength
    #print "\n\tTHETA={}".format(theta/locuslength)

    ## initialize random number sampling for numpy
    np.random.seed(params.seed1)

    ## setup seq mutator
    seqmutator = egglib.simul.CoalesceFiniteAlleleMutator(
                        theta=theta,
                        alleles=4,
                        randomAncestralState=True)
    seqmutator.setSites(locuslength)

    ## setup drop mutator
    dropmutator = egglib.simul.CoalesceFiniteAlleleMutator(
                        theta=theta,
                        randomAncestralState=False)
    dropmutator.setSites(locuslength)   #(len(params.cut1)+len(params.cut2))

    ## open files to write seq data to
    if params.outname:
        out1 = gzip.open(params.outname+"_R1_.fastq.gz", 'wb')
        if 'pair' in params.datatype:
            out2 = gzip.open(params.outname+"_R2_.fastq.gz", 'wb')
    else:
        out1 = out2 = sys.stderr

    ## b/c single cutter will sequence twice as many copies
    #if 'gbs' in params.datatype:
    #    params.nLoci = params.nLoci/2

    ## set stepsize for simulation depending on Nloci to make faster
    if params.nLoci >= 100:
        stepsize = 100
    elif params.nLoci >= 1000:
        stepsize = 500
    elif params.nLoci >= 2000:
        stepsize = 1000
    else:
        stepsize = 20

    ## simulate the data
    barcodes = {}
    counter = 0
    while counter < params.nLoci: 
        ##  this is a problem if cut is too frequent..
        ## seed has to change each iteration in a known way
        localseed1 = params.seed1 * (counter+1)
        localseed2 = params.seed2 * (counter+1)        

        ## performs simulation of setpsize loci at a time"
        aligns = egglib.simul.coalesce(
                    paramset, seqmutator, stepsize,
                    random=(localseed1, localseed2))
        ## put data in a dict
        aligns, names = namer(aligns, params, tiptax)

        ## make barcodes file if not already
        if not barcodes:
            barcodes = barcoder(names, params, barcodes)

        if params.dropout_cut or params.dropout_seq:
            dropcheck = egglib.simul.coalesce(
                            paramset, dropmutator, stepsize,
                            random=(localseed1, localseed2))

        ## mutation-disruption (allelic dropout)
        if params.dropout_cut:
            ## sim cut site length of data for dropout
            ## filter data dict based on dropcheck
            aligns = mutation_in_cut(params, aligns, dropcheck)

        if params.dropout_seq:
            ## filter for new cut sites generated by mutation
            aligns = mutation_new_cut(params, aligns, dropcheck)

        ## dresses up data to be fastq and puts in errors, indels, etc
        seqs1, seqs2, counter = seq_copies(aligns, barcodes, params, 
                                           counter, stepsize)
        #print counter-1, "Count"
        out1.write("".join(seqs1))
        if 'pair' in params.datatype:
            out2.write("".join(seqs2))

    out1.close()
    if "pair" in params.datatype:
        out2.close()

