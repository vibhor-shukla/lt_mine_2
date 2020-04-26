# Deployment Strategies in Underground Coal Mines



## Significance of various files 

#### _Deploy.m_ : 
For deployment of sensor nodes\
Inputs are length of interest and communication range\
Returns the total number of nodes and their positions and position of sink


#### _dist.m_ : 
returns the distance between 1 pair of points\
Input is the three dimensional coordinates of the points

#### _bfs\_connectivity.m_ : 
Implementation of Energy Model\
Assigns the amount of data transmitted and received


#### _MeanEnergy.m_ :
Calculates the mean energy of live nodes

#### _tpEnergy\_fun.m_ :
Assigns Energy level in each iteration according to their receiving and transmitting levels

#### _bfs.m_ :
Checks whether the live nodes are connected or not

#### _Destroy.m_ :
Input is the percentage of live nodes to be destroyed\
Changes the status from live to dead of certain percentage

#### _cpavsimulate.m_ :
Simulates for Average Deployment

#### _cpsimulate.m_ : 
Simulates for Second Deployment

#### _cptlayeredsimulate.m_ :
Simulates for Layered Deployment
