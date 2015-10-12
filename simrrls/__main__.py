#!/usr/bin/env python2

""" the main CLI for simrrls """

from __future__ import print_function, division ## requires python 2.7

from .simrrls import run
import argparse 
import pkg_resources


def main():
    """ parse args and run simrrls """
    params = parse()
    #print(params)
    ## run simrrls
    run(params)


def parse():
    """ parse the CLI arguments """
    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  Example usage: \n
  simrrls -o test1 -N 1e6 -n 4 -L 1000       
  simrrls -o test2 -f ddrad -c1 CCTGCAGG -c2 CCGG 
  simrrls -o test3 -f pairgbs -c1 CTGCAG -D 1 -i1 100 -i2 400 
  echo "((a:1,b:1):1,c:2);" > treefile 
  simrrls -o test4 -t treefile -df norm -dm 5 -ds 1
        """)

    ## add arguments 
    parser.add_argument('--version', action='version', 
        version=str(pkg_resources.get_distribution('simrrls')))

    parser.add_argument('-o', metavar="outname", dest='outname', 
        type=str, default=None, 
        help="[str] output file name prefix (default 'out')")

    parser.add_argument('-mc', metavar="dropout", dest="dropout_cut",
        type=int, default=0, 
        help="[0/1] allelic dropout from mutation to cut sites (default 0)")

    parser.add_argument('-ms', metavar="dropout", dest="dropout_seq",
        type=int, default=0, 
        help="[0/1] allelic dropout from new cut sites in seq (default 0)")

    parser.add_argument('-e', metavar="error", dest="error", 
        type=float, default=0.0005, 
        help="[float] sequencing error rate (default 0.0005)")

    parser.add_argument('-f', metavar="datatype", dest="datatype", 
        type=str, default='rad', 
        help="[str] datatype (default rad)                          \
        (options: rad, gbs, ddrad, pairddrad, pairgbs)")

    parser.add_argument('-I', metavar="indels", dest="indels", 
        type=float, default=0, 
        help="[float] rate of indel mutations (default 0) ex: 0.001")

    parser.add_argument('-l', metavar="length", dest='length', 
        type=int, default=100,
        help="[int] length of simulated sequences (default 100)")

    parser.add_argument('-L', metavar="nLoci", dest="nLoci", 
        type=int, default=100, 
        help="[int] number of loci to simulate (default 100)")

    parser.add_argument('-n', metavar="Ninds", dest="Ninds",
        type=int, default=1, 
        help="[int] N individuals from each taxon (default 1)")

    parser.add_argument('-N', metavar="Ne", dest="N", 
        type=float, default=int(5e5), 
        help="[int] pop size (Ne for all lineages; default 5e5)")

    parser.add_argument('-t', metavar="tree", dest='tree', 
        type=str, default="",
        help="[str] file name or newick string of ultrametric tree \
             (default 12 taxon balanced tree w/ bls=1)")

    parser.add_argument('-u', metavar="mu", dest='mu', 
        type=float, default=1e-9, 
        help="[float] per site mutation rate (default 1e-9)")

    parser.add_argument('-df', metavar="depthfunc", dest="depthfunc",
        type=str, default='norm', 
        help="[str] model for sampling copies (default norm, other=exp)")

    parser.add_argument('-dm', metavar="depthmean", dest="depthmean",
        type=int, default=10, 
        help="[int] mean sampled copies in norm, 1/m for exp (default 10)")

    parser.add_argument('-ds', metavar="depthstd", dest="depthstd",
        type=int, default=0, 
        help="[int] stdev sampled copies, used with norm model (default 0)")

    parser.add_argument('-c1', metavar="cut_1", dest='cut1', 
        type=str, default='CTGCAG',
        help="[str] restriction site 1 (default CTGCAG)")

    parser.add_argument('-c2', metavar="cut_2", dest='cut2', 
        type=str, default='CCGG',
        help="[str] restriction site 1 (default CCGG)")

    parser.add_argument('-i1', metavar="min_insert", dest='min_insert', 
        type=int, default=100,
        help="[int] total frag len = (2*l)+insert (default 100)")

    parser.add_argument('-i2', metavar="max_insert", dest='max_insert', 
        type=int, default=400,
        help="[int] total frag len = (2*l)+insert (default 400)")

    parser.add_argument('-r1', metavar="seed_1", dest='seed1', 
        type=int, default=1234567, 
        help="[int] random seed 1 (default 1234567)")

    parser.add_argument('-r2', metavar="seed_2", dest='seed2', 
        type=int, default=7654321,
        help="[int] random seed 2 (default 7654321)")

    ## parse sys.argv 
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()

