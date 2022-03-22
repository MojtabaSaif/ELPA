========================================================================
    ELPA: Evolutionary Label Propagation Algorithm
========================================================================

The example implements evolutionary label propagation algorithm (ELPA) that is 
an efficient algorithm based on LPA and the evolutionary algorithm framework for overlapping 
community detection. ELPA intelligently searches the different node-processing order 
and proposes a tradeoff between local information and global information. 

Fitting procedure and the community-Affiliation Graph Model are described in the following paper:
DETECTING OVERLAPPING COMMUNITIES IN COMPLEX NETWORKS: AN EVOLUTIONARY LABEL PROPAGATION APPROACH.

The code works under Windows with Visual Studio. Make sure that a
C++ compiler is installed on the system. Visual Studio project files
is provided. 

/////////////////////////////////////////////////////////////////////////////
Parameters:
   -f, --graph arg  Edge list of graph file name (default: ..\Data\karate.txt)
   -t, --gt arg     partitions file name: (default: ..\Data\gt\gt_karate.txt)
   -o, --out arg    partitions file name: (default: None)
      --gao arg    GA setting file name (default: None)
      --gaf arg    GA result file name (default: None)
      --run arg    The number of runs (default: 1)
  -h, --help       Print usage
   
   
