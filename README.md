###An Atmospheric model implementation with MPI.

The atmospheric model simulates atmospheric processes (wind speed, clouds, rain, etc.) 
which influence the weather or the climate. Resolves a set of differential equations 
describing, in this case, the dynamics of the atmosphere. The behavior of the equations 
in continuous real space, approached from the behavior of a finite set of points in space.
These points form a three-dimensional (3D) mesh:

![Atmo cube](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/asd.jpg)

The grid is periodic in the x, y axes, whereby the point (0,0,0) is adjacent to (0, Ny-1, 0) 
and (Nx-1, 0, 0). A vector of variables (pressure, temperature, wind speed and humidity) 
contained in each grid point. 

In the present work three variables:
* Xi,j	 -> through nine point stencil
* Yz,	 -> via Y = pow (X, n)
* Mi,j,k -> by Total Mass / # of processes

![9 point Stencil](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/stencil.jpg)

Calculations and distinct communications:

**1**.Horizontal nine-point stencil
	Calculation of value in position xi, j of the grid in year t + 1 (update) as follows:
	(Jacobi finite difference method)
	
![9 point Stencil](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/stencil2.jpg)
	
![equation](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/x_i_j.jpg)

Vertically, in our problem, there is no communication, only calculation based on the position in the grid:
 
y (x) = a * pow(x, n), n = 10
	
**2**.Calculation magazines of the total mass of the atmosphere, to confirm that the simulation works correctly.
	
![equation 2](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/total%20mass.jpg)
	
Where the mass at grid point (i, j, k).
	
Details: We consider the mass of the atmosphere is divided into processes. It does not change during the execution, but every few steps-repetitions will become reduce to determine that there is no difference.
	
**3**.Calculations in Physics.
y(x) = a * pow(x, n), n = 10
If any operation takes from one grid point then an important communication.
	
The communication related to the stencil is distributed, so it may proceed simultaneously. 
The same applies for the operation of the universal communications (air mass calculation).

In the process of domain decomposition atmospheric model created a point per task, ie Nx x Ny x Nz 
~ 10^5 processes (rather demanding), which leads us to the process of agglomeration.



**Agglomeration**

* From one point of the matrix per process 4 (the surface x, y) per process reduce communication
 requirements as regards the nine-point stencil, from 8 to 4 messages per process per step.

 ![agglo](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/agglomeration.jpg)
 
 ![agglo2](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/agglomeration2.jpg)
 
* Similarly it is possible to reduce the demands on communication and the vertical axis using batteries
 - No need, since there is no communication along the vertical axis z.
 
**Mapping**

According to our data, there is not a workload imbalance issue, so mapping the following form is deemed sufficient.
 ![grid](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/grid.jpg)

Combining the foregoing:
Having as the picture 3D mesh from above:

* each processor will take a piece of the surface (eg. 100x50 grid points / # processors) and their respective columns.

* In the vertical z axis does not occur any form of communication, however, the corresponding size of partition 
multiplies the size of messages and ultimately plays a significant role in the choice of partition, as we will see below in the analysis of performance.


**Theoretical calculation**

In our program there is full support three-dimensional segmentation of data. In each dimension, the segmentation units are called partitions (and by extension the crowds xpc / ypc / zpc of x / y / z partition count). The condition to operate the partition is the amount of data in each dimension is a multiple of the corresponding plurality of partition unit so that each unit getting equal data with a simple division.

As mentioned above, the typical dimensions of the problem are:

![nxnynz](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/nxnynz.jpg)

Let's define and counterparts of segmentation units, namely those that exist in each MPI process:

![xyz](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/xyz.jpg)

Processing times are:

![Tc](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/Tc.jpg)

where **internal** are the elements that do not require previous communication to calculate and the rest busywait / external.

Considering the network processor has sufficient capacity in relation to our needs, we will assume that all communications may be done in parallel. Therefore, the communication time will not be equal to the set of messages, but with the maximum message:

![Tcomm](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/Tcomm.jpg)

where ts is the time constant of each message, tw the message time per byte, 3 is the number of doubles in our struct, and 8 is the sizeof(double).

The calculation of processing and communication time, in the non-blocking communication, starts immediately. The calculation of internal data follows, then there are
MPI_Wait to calculate the rest. Thus, the internal processing time spans a portion of the communication time:

![T](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/T.jpg)

At this point the parallel communication and computation makes a further theoretical calculation (speedup and efficiency) difficult. We shall confine ourselves to comparing the theoretical time with the sampled time.


The processing time is inevitable since the arithmetic operations must occur.
Therefore, we try to reduce communication. Let's rethink the equation of Tcomm.

This equation refers to a 3-dimensional partition. For a better communication parallelization the size of messages should be equal:

![equation](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/equation.jpg)

So we get a criterion for the division.
We also note the z partition indirectly affects the problem (as it forces a change in the other partitions).
In addition, removing only the x- or y- partition has a real effect in communication time. One message will be large, and the other would virtually not exist, eliminating the possibility of parallel communication.

However, the abolition of both is indeed very advantageous: no communication is required! This results from the complete absence of communication in the z-axis.

**Code implementation and Tiles**

An example of the elements of each rank (halo elements).

![tiles](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/tiles.jpg)

**Experimental calculations**

Estimations:

* The total time who ran each process have been calculated with MPI_Reduce and mode MPI_SUM. To find the average we have divided by the number of total processes.
* Calculate tc time by running the program with one process. We choose a value for the larger problem, which is more representative (orange fill)
* Again, we estimate the running time of a single process, if the memory fits the data (gray cells represent impotence of execution).
* Calculation of Speedup and Efficiency. (orange fill use aforementioned calculations).

The graph Efficiency(P) is for the scaling of processors. As we expected, we note the decline in performance by increasing processor which is also smaller for larger number of data.

The graph Efficiency(Data) is for the scaling data. We observe increased performance by increasing data. The increase is lower for most processes.

![E(data)](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/E%28data%29.jpg)
![E(p)](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/E%28p%29.jpg)

**Theoretical calculations**

It is important to note that the aforementioned theoretical calculations, use a fixed tc, tx, ts (in yellow cell A35) that will help us for theoretical calculations and conclusions.

The times internal/comm/busywait-external/theoretical(sum) correspond to the theoretical analysis.

Ttheoritical / Texperimental * 100%: Here we compare the theory with practice. We would have a successful match for values that are ​​close to 100. We note two ways in which the situation differs:

**1**.In areas of high performance (small number of processes and / or a large number of data) instead of values ​​are near 100, there are values near 200.

**2**.In areas of low performance, experimental time is much greater than the theoretical. This correlates with the constant ts and tw. 
In small problem sizes/multiple processors we expect a larger portion of time to be in communication.

*In the MPI measurements ts = 0.8ms, tw = 0.004ms / byte
*From colleague Measurements: ts = 82ms, tw = 0.2ms / byte

The difference is up to two orders of magnitude. However, this is not enough to "fill the gap". By testing at the entrance of ts / tw, to reach the values approaching theoretical practices should be 3.5 orders of magnitude larger. Of course, there is still variation in the values, but not so great! This is confirmed by the fact that running the program without calculations, waiting times (only for send / recv wait) is much larger than the values calculated from mpptest.
The last two tables are constructed to expose the real problem: For the given problem sizes used, although the theory predicts that the calculation of the internal elements will occupy 99% of the total time in all cases, in practice the results are much more different.

**Reduce**

Here are measurements and graphs relative to the Reduce in the total mass of the problem. We observe an expected linearity. Apparently this calculation step has a fixed time cost, which has linearity both in terms of frequency and in the size of the remaining problem. 
The more often it is called, the more time it requires in overall, and the bigger the data is, the calculations and communications of the remaining problem, the more negligible time it takes for reduce.

![reduce](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/reduce.jpg)

**Conclusions**

Nowadays that the speed of processors have ceased following the evolution of the law of Moore and with multiple cores even in mobile phones, the parallelism at a data level with threads and processes, is a very important field in software engineering.
In this project we experienced the base level processes of MPI. We used the Foster methodology for the development and provisioning of the speed of a parallel program. There were differences between the theoretical expectations and the experimental results. However, experimental results in terms of scalability and data processors were as expected.
