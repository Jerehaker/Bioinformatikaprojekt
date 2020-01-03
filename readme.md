### Bioinformatics project
---------------------
This repo contains our implementation Navaro's algorithm for finding the minimum edit distance of strings (gene sequences), this implementation uses only one thread and achieves dubious results.
the algorithm is presented in: https://www.sciencedirect.com/science/article/pii/S0304397599003333
We have generated the test data from this repo: https://github.com/maickrau/GraphAligner/tree/PaperExperiments/WabiExperimentSnake and compare our results to their implementation of the Navaro algorithm, we have also used their functions for reading from fastq files as well as their gfagraph class.
Running our project:
--------------
* 1. Compile the code using this command: g++ -std=c++11 main.cpp gfagraph.cpp fastqloader.cpp -o test
* 2. Run the code using this command: ./test 0 1 1 ref10000_simulatedreads.fastq ref10000_onechar.gfa 
(you can comment/uncomment line 145 depending if you want output after every sequence (there are 74 sequences in our test file))

Arguments explained:
- match cost
- miss cost
- indel cost
- sequences file (.fastq format)
- the graph you want to run the algorithm on (.gfa format)

Results:
-------------
Our algorithm works 6 times slower than the on in maickrau's repo for every graph topology except for the tangled graph, we are not sure why that is happening
You can find out results in the results_nasi folder and the results from maickrau's repo in results_njihovi. Both have been run on my PC which has an Iintel i5-4670k, 16GB of 1333MHz RAM(the project had a restriction of 16GB RAM but to me it seems that it is redundat since the whole point of Navaro's algorithm is taht it onnly stores 2 rows in memory which is SIGNIFICANTLY less than 16GB) on ubuntu 18.04

Contact:
--------
jeronim96@fer.hr or here
