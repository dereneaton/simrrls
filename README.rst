**simrrls**
===========

A program to simulate raw RADseq-like data with options to modify tree, population, sequencing, and formatting parameters. 

Requirements
------------
+ Python 2.7+
+ Numpy Python module
+ Egglib Python module .. _http://egglib.sourceforge.net/

Example usage: 
---------------

Using default parameters::

    $ simrrls -o test1

Modified population parameters::

    $ simrrls -o test2 -N 1e6 -u 2e-8 

Modifications to sequencing parameters::

    $ simrrls -o test3 -L 5000 -l 200 -e 0.001 -d 10,2 

More modifications related to library prep methods  
(In this case making paired-end reads and allowing read overlap)::

    $ simrrls -o test4 -d 10,2 -f pairddrad -s -50,200 

Examples of using a non-default input topology::
    $ echo "((a:1,b:1):1,c:2);" > treefile  
    $ simrrls -o test4 -d 10,2 -f pairddrad -s -50,200 -t treefile
