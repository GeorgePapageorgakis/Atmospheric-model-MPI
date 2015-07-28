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

1.Horizontal nine-point stencil
	Calculation of value in position xi, j of the grid in year t + 1 (update) as follows:
	(Jacobi finite difference method)
![9 point Stencil](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/stencil2.jpg)
	
![equation](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/x_i_j.jpg)

Vertically, in our problem, there is no communication, only calculation based on the position in the grid: 
y (x) = a * pow(x, n), n = 10
	
2.Calculation magazines of the total mass of the atmosphere, to confirm that the simulation works correctly.
	
![equation 2](https://github.com/GeorgePapageorgakis/Atmospheric-model-MPI/blob/master/figures/total%20mass.jpg)
	
	Where the mass at grid point (i, j, k).
	
	Details: We consider the mass of the atmosphere is divided into processes. It does not change during the 
	execution, but every few steps-repetitions will become reduce to determine that there is no difference.
	
3.Calculations in Physics.
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
multiplies the size of messages and ultimately plays his role in the choice of partition, as we will see below 
in the analysis of performance.
