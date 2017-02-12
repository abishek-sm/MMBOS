# MMBOS - MinMin Based on Order Statistics 
This is a novel algorithm designed and implemented in Java which forms a part of my thesis project. 
# Introduction
The paper for this algorithm is currently being prepared. This is a novel algorithm that uses Distance Correlation to build networks from gene expression data. I am presenting a sample of the code here written by me in src/MinMinAlg.java. This code depends on classes and methods which are open source and some written by a post doc in the lab. 

For convenience I am also providing the entire JAR file which will be made available on the lab website along with the paper on the algorithm details. 

# Usage of the JAR file
java -jar MMBOSg.jar [-path to data-] [-path to distance correlation matrix-] ["extension of the result file"] [number of edges to be recovered in thousands]

for example: 
java -jar ELMM_MinMinAlg2.jar "/home/usr/data.txt" "/home/usr/dc_matrix.txt" "bsdc300" 300
will recover 300,000 edges in the network





