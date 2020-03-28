# Graph Pattern Sampling (Output Space Sampling)

The MCMC approach was extended in [2009-graphsampling] to mine a sample of all frequent patterns, to mine support-biased patterns and also to mine a sample of discriminative patterns.

See also https://github.com/zakimjz/Origami which can mine a sample of maximal graph patterns, but does not provide any uniformity guarantee. 

See also https://github.com/zakimjz/MUSK that proposes a Markov Chain Monte Carlo based approach to guarantee a uniform sample of all maximal patterns. 

**Relevant Publications**

* [2009-graphsampling] Mohammad Al Hasan and Mohammed J. Zaki. Output space sampling for graph patterns. Proceedings of the VLDB Endowment (35th International Conference on Very Large Data Bases), 2(1):730–741, 2009.


## TO COMPILE:

 make

## TO RUN:

./uniform_sampling -d ./dataset/GRAPH_large.dat -c 1000 -s 30

        -d  = input file  (see the dataset directory for input graph format)
        -c  = maximum iteration
        -s  = minimum support


