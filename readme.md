### Bioinformatics project
=========

This repo contains our implementation Navaro's algorithm for finding the minimum edit distance of strings (gene sequences), this implementation uses only one thread and achieves dubious results. This project was done for the bioinformatics course at fer [course page](https://www.fer.unizg.hr/predmet/bio/).
The algorithm is presented in: [this paper](https://www.sciencedirect.com/science/article/pii/S0304397599003333)
We have generated the test data from this repo: [Here](https://github.com/maickrau/GraphAligner/tree/PaperExperiments/WabiExperimentSnake) and compare our results to their implementation of the Navaro algorithm, we have also used their functions for reading from fastq files as well as their gfagraph class, which we heavily modified for reading the snp and tangle topologies since our approach was different from theirs.
Test files are: ref10000_simulatedreads.fastq, ref10000_snp.gfa, ref10000_onechar.gfa, ref10000_twopath.gfa, ref10000_tangle.gfa.

Running our project:
--------------
* 1. Compile the code using this command: g++ -std=c++11 main.cpp gfagraph.cpp fastqloader.cpp -o test
* 2. Run the code using this command: ./test 0 1 1 ref10000_simulatedreads.fastq ref10000_onechar.gfa 1
(you can comment/uncomment line 145 depending if you want output after every sequence (there are 74 sequences in our test file))

Arguments explained:
- 1. match cost
- 2. miss cost
- 3. indel cost
- 4. sequences file (.fastq format)
- 5. the graph you want to run the algorithm on (.gfa format)
- 6. display boolean - passed as integer; 1:display results in console after every sequence, 0:write all sequence results in the _results.txt

Results:
-------------
Our algorithm works ~6 times slower than the on in maickrau's repo for every graph topology except for the tangled graph. we were not able to run the algortihm on tangle graph because it has cyclic regions.
You can find out results in the results_ours folder and the results from maickrau's repo in results_theirs. Both have been run on my PC which has an Iintel i5-4670k, 16GB of 1333MHz RAM (the project had a restriction of 16GB RAM but to me it seems that it is redundamt since the whole point of Navaro's algorithm is that it only stores 2 rows in memory which is SIGNIFICANTLY less than 16GB) on ubuntu 18.04


| Runtime [μs]   | onechar graph | twopath graph | snp graph |
| ------------- | ------------- | ------------- | ------- |
| Our runtime  | 970988279  | 2119507841 | 1170734973 |
| Their runtime  | 140008518  | 321483028 | 127856765 |
| Runtimes ratio | 6.93 | 6.59 | 9.15 |
| Memory consumption [MiB] | 6.5 | 12.6 | 7.5 |

Future work:
-------------
Expand the algorithm for cyclic graphs.

Contact:
--------
jeronim.matijevic@fer.hr or
ema.puljak@fer.hr or
luka.cagalj@fer.hr
