# paralel_relief_algorithm
This repository contains source code for parallel relief algorithm implementation in cpp. There is 1 master 
processor and (P-1) slave processors. (P is the total number 
processors). Master procesor handles I/O operations and send the 
given data to slave processor (same number of instance to each of 
them, number of iterations and number of feature need to be detected 
by each slave ). Slave processors applies relief alghorithm. Relief is an 
algorithm that enables us to detect most important and varying 
features in our data set. I will explain relief alghorithm in detail 
incoming sections. At end slaves prints the feature ids they find and 
sends this ids to master. Master takes this ids, removes dublicates and 
print result. 

# Program Execution:
my mpi version (2.1.1)
gcc (Ubuntu 7.4.0-1ubuntu1~18.04.1) 
7.4.0)
There should be a cpp compiler in your computer to run this program. 
Following comments will be explain the setup needed in lubuntu 
environment. Then download libopenmpi-dev library with "apt install" 
command (version 2.1.1). After this you are ready to run this program. 
Open command line and go to file that includes source file. After that 
run the following command to create to object file.
1. mpic++ -o hellopp ./cmpe300_mpi_2017400234.cpp
Let's say name of our input file is "input.txt" and we have 11 
processors. Then to run this object file 
2. mpirun -oversubscribe -np 11 ./hellopp input.txt
If program does not stop, you can use Ctrl^c to stop.
cmpe300_mpi_2017400234.cpp

