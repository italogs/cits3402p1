cits3402p1
==========

Programming project 1


General comments:

To avoid confusion in the code, I removed some useless comments in the code, to make it more readable and cleaner.

All the parallelizations made here were implemented in small code blocks whereas I could not parallelize the main loop because of its complexity and its dependency of pior information.

As stated on the project description, the code is suitably commented in order to the reader can have a general understanding about the code and its OpenMP directives.

One approach that I have made use of was using OpenMP sections deal with writing files. So, each section was assigned to write an output file.
This approach really speed up the code.
Also, a good approach was using reduction, what had a significant speed up.

Results:

The obtained results from the modified code(parallelized code) were really promising.
In this document, I am going to provide a considerable amount of results for three set of variables that I considered suitable for this project.
I tried to combine different parameters for the variables: chainlngth and dt.




First set 
chainlngth 10
dt  0.001

Second set
chainlngth 2000
dt  0.001

Third set
chainlngth 100
dt 	0.00000001


Fourth set
chainlngth 10
dt  0.00000001


