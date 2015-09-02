## _simrrls_

A program to simulate raw radseq-like sequence data with options to modify tree, population, sequencing, and formatting parameters. More detailed documentation coming soon. 

#### Requirements
+ Python 2.7+
+ Numpy python module
+ Egglib python module

### Example usage: 

#### Using default parameters
```bash
simrrls -o test1
```

#### Modified population parameters
```bash
simrrls -o test2 -N 1e6 -u 2e-8 
```

#### Modifications to sequencing parameters
```bash
simrrls -o test3 -L 5000 -l 200 -e 0.001 -d 10,2 
```

#### More modifications related to library prep methods  
(In this case making paired-end reads and allowing read overlap)
```bash
simrrls -o test4 -d 10,2 -f pairddrad -s -50,200 
```

#### Examples of using a non-default input topology
```bash
echo "((a:1,b:1):1,c:2);" > treefile  
simrrls -o test4 -d 10,2 -f pairddrad -s -50,200 -t treefile
```