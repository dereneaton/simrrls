#!/usr/bin/env python2.7
""" the main CLI for simrrls """

from __future__ import print_function, division ## requires python 2.7

from .simrrls import run
import argparse 
import pkg_resources


def main():
    """ parse args and run simrrls """
    params = parse()
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
  simrrls -o test4 -t treefile -d 10,0
        """)

    ## add arguments 
    parser.add_argument('--version', action='version', 
        version=str(pkg_resources.get_distribution('simrrls')))

    parser.add_argument('-D', metavar="dropout", dest="dropout",
        type=bool, default=0, 
        help="[0/1] allow allelic dropout (default 0)")

    parser.add_argument('-e', metavar="error", dest="error", 
        type=float, default=0.0005, 
        help="[float] sequencing error rate (default 0.0005)")

    parser.add_argument('-f', metavar="datatype", dest="datatype", 
        type=str, default='rad', 
        help="[str] datatype (default rad)                          \
        (options: rad, gbs, ddrad, pairddrad, pairgbs)")

    parser.add_argument('-I', metavar="indels", dest="indels", 
        type=bool, default=0, 
        help="[0/1] allow indels (default 0)")

    parser.add_argument('-l', metavar="length", dest='length', 
        type=int, default=100,
        help="[int] length of simulated sequences (default 100)")

    parser.add_argument('-L', metavar="nLoci", dest="nloci", 
        type=int, default=100, 
        help="[int] number of loci to simulate (default 100)")

    parser.add_argument('-n', metavar="Ninds", dest="Ninds",
        type=int, default=1, 
        help="[int] N individuals from each taxon (default 1)")

    parser.add_argument('-N', metavar="Ne", dest="N", 
        type=int, default=5e5, 
        help="[int] pop size (Ne for all lineages; default 5e5)")

    parser.add_argument('-o', metavar="outfile", dest='outfile', 
        type=str, default="out.fastq.gz", 
        help="[str] output file name prefix (default out)")

    parser.add_argument('-t', metavar="tree", dest='tree', 
        type=str, default="",
        help="[str] file name or newick string of ultrametric tree \
             (default 12 taxon balanced tree w/ bls=1)")

    parser.add_argument('-u', metavar="mu", dest='mu', 
        type=float, default=1e-9, 
        help="[float] per site mutation rate (default 1e-9)")

    parser.add_argument('-c1', metavar="cut_1", dest='cut1', 
        type=str, default='CTGCAG',
        help="[str] restriction site 1 (default CTGCAG)")

    parser.add_argument('-c2', metavar="cut_2", dest='cut2', 
        type=str, default='GAATTC',
        help="[str] restriction site 1 (default GAATTC)")

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

