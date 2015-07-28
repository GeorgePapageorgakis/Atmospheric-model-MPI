/**
    atmosphere.c
    Purpose: An MPI-openMP implementation of 3D atmospheric model
    
    @version 1.0 3/2014
*/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _OPENMP
	#include <omp.h>
#else
	int omp_get_thread_num	(void) 	{ return 0; }
	int omp_get_num_threads	(void) 	{ return 1; }
#endif //Failsafe for systems without openmp

int STEPS,NX,NY,NZ; 
// N(x/y/z) is the size of the problem in elements

//Note:they refer to the relative position between partitions, NOT the data position within the partitions
#define DTUTAG   2       //down-to-up tag. 
#define UTDTAG   3       //up-to-down tag
#define LTRTAG   4       //left-to-right tag
#define RTLTAG   5       //right-to-left tag

struct element{
	double mass;
	double xy_value;
	double z_value;
};

// Constants
const 	double total_mass = 1000.0;
const 	double a		  = 2.718281828; 
	double mass_reduced;
	double rank_mass;
		
void prtdat	(char *fnam);
void update	(int, int, int, struct element *, struct element *);

// c coordinates to simplify the 3-dimensional array handling
int 	c	(int x, int y, int z); 
char filename	[200];
/* array of 2 pointers, each pointer represents three-dimensional array */
struct element	*u[2];
/* task's unique id */     
int	rank,
/* number of tasks */                   
numtasks,             
/* loop variables */     
i, ix, iy, iz, iu, it,

//PARTITIONS
x_length,	x_partition_count,	u_x_length,	x_rank,	x_start, 	x_end,
y_length,	y_partition_count,	u_y_length,	y_rank,	y_start,	y_end,
z_length,	z_partition_count,	u_z_length,	z_rank,	z_start,	z_end,

// Vars to check if a section adjacent to a neighbour is received and calculated
done_count,	done[4],	done_index[4], 
// For subarray datatype
sizes[3], 	subsizes_right_left[3], subsizes_down_up[3], starts_right_left[3]; 

// The datatypes, the dimension in which their elements progress (the other dimension is constant)
MPI_Request 	req_recv[4],	req_send[4];
MPI_Datatype 	xmargindata,	ymargindata,	right_left_type,down_up_type; 
MPI_Status 		stat_recv[4],	stat_send[4];

double	time_start_internal,	time_start_busywait,	time_start_external,	time_start_mass_reduce,
	time_start_sendrecv,	time_start_waitsend,	duration_internal,  	duration_busywait,  
	duration_external,  	duration_mass_reduce,  	duration_sendrecv,  	duration_waitsend,
	total_internal,     	total_busywait,     	total_external,     	total_mass_reduce,     
	total_sendrecv,     	total_waitsend,		total_accounted_for;
FILE * fp;
int reduce_frequency = 5;


int main (int argc, char *argv[]){
	MPI_Init	(&argc,&argv);
    	MPI_Comm_size	(MPI_COMM_WORLD,&numtasks);
    	MPI_Comm_rank	(MPI_COMM_WORLD,&rank);
	STEPS	= 50;
	NX	= 800;
	NY	= 400;
	NZ	= 30;
	x_partition_count = 2;
	y_partition_count = 1;
	z_partition_count = 1;
	
	for(i=1; i<argc; ++i) {
		if(!strcmp(argv[i], "s"))
			STEPS = atoi(argv[i+1]);
		else if(!strcmp(argv[i], "x"))
			NX = atoi(argv[i+1]);
		else if(!strcmp(argv[i],"y"))
			NY = atoi(argv[i+1]);
		else if(!strcmp(argv[i],"z"))
			NZ = atoi(argv[i+1]);
		else if(!strcmp(argv[i],"px"))
			x_partition_count = atoi(argv[i+1]);
		else if(!strcmp(argv[i],"py"))
			y_partition_count = atoi(argv[i+1]);
		else if(!strcmp(argv[i],"pz"))
			z_partition_count = atoi(argv[i+1]);
	}
	
	rank_mass = total_mass/(double)numtasks;
	
	if(rank == 0) {
		printf("\n\n==================================\n==========SESSION START===========\n==================================\n"
	        "Program size (x/y/z/steps): %d x %d x %d x %d\n"
			"Partition grid (x/y/z): %d x %d x %d\n", NX, NY, NZ, STEPS, x_partition_count, y_partition_count, z_partition_count);
		fp = fopen("log.txt", "a");
		fprintf(fp,"%dx%dx%dx%d\n%dx%dx%d\n%d processes\n\n",NX, NY, NZ, STEPS, x_partition_count, y_partition_count, z_partition_count,
															x_partition_count * y_partition_count * z_partition_count);
		if(!(NX && NY && NZ && x_partition_count && y_partition_count && z_partition_count)){
			puts("Elements/Grid zero, cannot continue\n"\
				"Use -x  <number> to input x elements count\n"\
				"Use -y  <number> to input y elements count\n"\
				"Use -z  <number> to input z elements count\n"\
				"Use -s  <number> to input step count\n"\
				"Use -px <number> to input x-dimension partition count\n"\
				"Use -py <number> to input y-dimension partition count\n"\
				"Use -pz <number> to input z-dimension partition count\n"\
				);
			//getchar();
			return;
		}
	}
	
	#pragma omp parallel
	{
		if((rank==0) && (omp_get_thread_num() == 0)) 
			printf("Internal element processing threads (OpenMP parallelization): %d\n", omp_get_num_threads());
	}
	MPI_Barrier(MPI_COMM_WORLD); //for printf to apper in order
	// ================================START Data partition assignment================================
	
	x_length = NX/x_partition_count;		//Divide elements uniformly among partitions
	x_rank	 = rank % x_partition_count;		//rank, as is "partition rank"
	x_start	 = x_rank * x_length ;			//min
	x_end	 = (x_rank+1) * x_length - 1;

	y_length = NY/y_partition_count;
	y_rank	 = (rank  / x_partition_count ) % y_partition_count; //rank, as is "partition rank"
	y_start	 = y_rank * y_length;			//min
	y_end	 = (y_rank+1) * y_length - 1;
		
	z_length = NZ/z_partition_count;	
	z_rank	 = rank / (x_partition_count*y_partition_count);
	z_start	 = z_rank * z_length;			//min
	z_end	 = (z_rank+1) * z_length - 1;
	
	printf("Rank %d range: x(%d-%d) of %d, y(%d-%d) of %d, z(%d-%d) of %d\n",
					rank, x_start, x_end, NX, y_start, y_end, NY, z_start, z_end, NZ);
	//================================END Data partition assignment================================

	
	//=====================================START Initialization====================================
	duration_sendrecv   = 0.0;
	duration_internal   = 0.0;
	duration_busywait   = 0.0;
	duration_external   = 0.0;
	duration_mass_reduce= 0.0;
	duration_waitsend   = 0.0;
	//Each of the arrays needs to have a size of (x+4)*(y +4)* z elements.
	//The size must be identical for all so that MPI datatype column will work correctly.
	//The +4 is the Halo zone for receives.
	u[0]	= (struct element*) malloc((x_length+4)*(y_length+4)*z_length*sizeof(struct element));
	u[1]	= (struct element*) malloc((x_length+4)*(y_length+4)*z_length*sizeof(struct element));
	for (iz=0; iz<z_length; ++iz)
		for (iy=0; iy<y_length+4; ++iy)
			for (ix=0; ix<x_length+4; ++ix){
				(u[0]+c(ix,iy,iz))->mass 	= total_mass/numtasks/NX/NY/NZ;
				(u[0]+c(ix,iy,iz))->xy_value 	= (double)(rand() % 100);
				(u[0]+c(ix,iy,iz))->z_value 	=  a*pow((u[0]+c(ix,iy,iz))->xy_value,10.0);
				(u[1]+c(ix,iy,iz))->mass 	= total_mass/numtasks/NX/NY/NZ;
				(u[1]+c(ix,iy,iz))->xy_value 	= 0.0;
			}
	iu = 0; //iz: Track which of the u arrays is the "old"
	//sprintf(filename,"atm%ds%d.txt", rank, it);
	//prtdat(filename);
	//printf("Rank %d saving in %s\n", rank, filename);
	// for printf to apper in order
	MPI_Barrier(MPI_COMM_WORLD);
	
	// DATATYPE: Notice how column size(1st arg) depends on partition size
	// an element consists of 3 doubles
	MPI_Type_vector (2,  3*x_length, 3*(x_length+4), MPI_DOUBLE, &xmargindata); 
	MPI_Type_commit (&xmargindata);
	MPI_Type_vector (y_length, 6, 3*(x_length+4), MPI_DOUBLE, &ymargindata);
	MPI_Type_commit (&ymargindata);
	sizes[2]		= 3*(x_length+4);
	sizes[1]		= y_length+4;
	sizes[0]		= z_length;
	subsizes_right_left[2]	= 3*2;
	subsizes_right_left[1]	= y_length;
	subsizes_right_left[0]	= z_length;
	starts_right_left[0]	= 0;
	starts_right_left[1]	= 0;
	starts_right_left[2]	= 0;
	MPI_Type_create_subarray(3, sizes, subsizes_right_left, starts_right_left, MPI_ORDER_C, MPI_DOUBLE, &right_left_type);
	MPI_Type_commit (&right_left_type);
	subsizes_down_up[2]	= 3*x_length;
	subsizes_down_up[1]	= 2;
	subsizes_down_up[0]	= z_length;
	MPI_Type_create_subarray(3, sizes, subsizes_down_up, starts_right_left, MPI_ORDER_C, MPI_DOUBLE, &down_up_type);
	MPI_Type_commit (&down_up_type);
	
	printf("Rank %d has finished initialisation\n", rank);
	//==============================================END Initialization==============================================

	//Main Computation
	for (it = 1; it <= STEPS; ++it)	{
		if(rank == 0) printf("Step %d\n", it);
		time_start_sendrecv = MPI_Wtime();
		
		//printf("Rank %d starts iteration %d\n",rank,it);
	/*	if(STEPS==1) printf("Rank %d neighbours: U %d  D %d  L %d  R %d\n",rank,
		((rank+x_partition_count) % (x_partition_count*y_partition_count) 
			+ (x_partition_count*y_partition_count)*(rank/x_partition_count/y_partition_count)),
				((rank-x_partition_count+x_partition_count*y_partition_count) % (x_partition_count*y_partition_count) 
					+ (x_partition_count*y_partition_count)*(rank/x_partition_count/y_partition_count)),
						(rank+(rank % x_partition_count ? 0 : x_partition_count)-1),
							(rank+((rank+1) % x_partition_count ? 0 : -x_partition_count)+1));*/
		
		MPI_Isend(u[iu]+c(2, y_length, 0), 1, down_up_type,\
			(rank+x_partition_count) % (x_partition_count*y_partition_count) + 
				(x_partition_count*y_partition_count)*(rank/x_partition_count/y_partition_count), 
					DTUTAG, MPI_COMM_WORLD, req_send + 0);
		MPI_Irecv(u[iu]+c(2, y_length+2, 0), 1, down_up_type,\
			(rank+x_partition_count) % (x_partition_count*y_partition_count) +
				(x_partition_count*y_partition_count)*(rank/x_partition_count/y_partition_count), 
					UTDTAG, MPI_COMM_WORLD, req_recv + 0);

		MPI_Isend(u[iu]+c(2,2,0), 1, down_up_type,\
			(rank-x_partition_count+x_partition_count*y_partition_count) % (x_partition_count*y_partition_count) + 
				(x_partition_count*y_partition_count)*(rank/x_partition_count/y_partition_count), 
					UTDTAG, MPI_COMM_WORLD, req_send + 1);
		MPI_Irecv(u[iu]+c(2,0,0), 1, down_up_type,\
			(rank-x_partition_count+x_partition_count*y_partition_count) % (x_partition_count*y_partition_count) +
				(x_partition_count*y_partition_count)*(rank/x_partition_count/y_partition_count), 
					DTUTAG, MPI_COMM_WORLD, req_recv+1);

		// use % due to spatial wraparound
		MPI_Isend(u[iu]+c(2,2,0), 1, right_left_type,\
			rank+(rank % x_partition_count ? 0 : x_partition_count)-1, RTLTAG, MPI_COMM_WORLD, req_send+2); 
		MPI_Irecv(u[iu]+c(0,2,0), 1, right_left_type,\
			rank+(rank % x_partition_count ? 0 : x_partition_count)-1, LTRTAG, MPI_COMM_WORLD, req_recv+2);

		MPI_Isend(u[iu]+c(x_length,2,0), 1, right_left_type,\
			rank+((rank+1) % x_partition_count ? 0 : -x_partition_count)+1, LTRTAG, MPI_COMM_WORLD, req_send+3);
		MPI_Irecv(u[iu]+c(x_length+2,2,0), 1, right_left_type,\
			rank+((rank+1) % x_partition_count ? 0 : -x_partition_count)+1, RTLTAG, MPI_COMM_WORLD, req_recv+3);
		//printf("Rank %d has finished nonblocking sendrecvs\n", rank);
		duration_sendrecv += MPI_Wtime() - time_start_sendrecv;
			
		//begin update of internal elements
		time_start_internal = MPI_Wtime();
		#pragma omp parallel
		{
			#pragma omp for //collapse (2)
				for(iz=0; iz<z_length; ++iz){ //full z range
					//printf("Iteration %d is assigned to thread %d\n", iz, omp_get_thread_num());
					//disregard both the data waiting to be received (width 2 perimeter) and the ones 
					//who need them to be calculated (another 2 width perimeter)(central elements)
					for(iy=4; iy<y_length; ++iy) 
						for(ix=4; ix<x_length; ++ix)
							update(ix, iy, iz, u[iu], u[1-iu]);
				}
		}
		duration_internal += MPI_Wtime() - time_start_internal;
		// printf("Rank %d has finished internal elements\n", rank);
		// finished update of internal elements
	
		time_start_busywait = MPI_Wtime();
		done_count = 0;                  
		memset(done, 0, 4*sizeof(int)); 
		while(done_count<4){
			for(i=0; i<4; ++i)
				if(!done[i]){
					MPI_Test(req_recv+i, done+i, MPI_STATUS_IGNORE);
					if(done[i]){
						switch(i){
						case 0:
							for(iz=0; iz<z_length; ++iz) //full z range
								for(iy=y_length; iy<y_length+2; ++iy)
									for(ix=2; ix<x_length+2; ++ix)
										update(ix,iy,iz,u[iu],u[1-iu]);//update top row except corners
							break;
						case 1:
							for(iz=0; iz<z_length; ++iz) //full z range
								for(iy=2; iy<4; ++iy)
									for(ix=2; ix<x_length+2; ++ix)
										update(ix,iy,iz,u[iu],u[1-iu]);//update bottom row except corners
							break;
						case 2:
							for(iz=0; iz<z_length; ++iz) //full z range
								for(ix=2; ix<4; ++ix)
									for(iy=2; iy<y_length+2; ++iy)
										update(ix,iy,iz,u[iu],u[1-iu]);//update left column except corners
							break;
						case 3:
							for(iz=0;iz<z_length;iz++) //full z range
								for(ix=x_length;ix<x_length+2;ix++)
									for(iy=2;iy<y_length+2;iy++)
										update(ix,iy,iz,u[iu],u[1-iu]);//update right column except corners
						}
						++done_count;
					}//end if(done[i])
				}//end if(!done[i]).
		}//end while(done_count<4)
		//printf("Rank %d has finished busywait phase\n", rank);
		duration_busywait += MPI_Wtime() - time_start_busywait;

		time_start_external = MPI_Wtime();
		for(iz=0; iz<z_length; ++iz) //full z range
			for(iy=2*y_length-2; iy<2*y_length+2; ++iy)
				for(ix=2*x_length-2; ix<2*x_length+2; ++ix)
					update(ix%x_length+2,iy%y_length+2,iz,u[iu],u[1-iu]);//update the four corners
		duration_external += MPI_Wtime() - time_start_external;

		time_start_mass_reduce = MPI_Wtime();
		if(it % reduce_frequency == 0){
			MPI_Reduce(&rank_mass, &mass_reduced, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank == 0)
				printf("Step %d: Rank %d reduced total mass of %f\n", it, rank, mass_reduced);
		}
		duration_mass_reduce += MPI_Wtime() - time_start_mass_reduce;
			
		time_start_waitsend = MPI_Wtime();
		for(i=0; i<4; ++i)
			MPI_Wait(req_send+i, MPI_STATUS_IGNORE);//Wait for the sends
		//MPI_Barrier(MPI_COMM_WORLD);
		//printf("rank %d finished MPI_Waits at step %d\n", rank, it);	
		//Revert arrays
		iu = 1-iu;
		duration_waitsend += MPI_Wtime() - time_start_waitsend;
		
		//sprintf(filename,"atm%ds%d.txt", rank, it);
		//prtdat(filename);
				
	}//end STEPS iteration
	
	MPI_Reduce(&duration_sendrecv   ,&total_sendrecv   ,1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&duration_busywait   ,&total_busywait   ,1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&duration_internal   ,&total_internal   ,1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&duration_external   ,&total_external   ,1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&duration_mass_reduce,&total_mass_reduce,1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&duration_waitsend   ,&total_waitsend   ,1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	total_accounted_for = total_sendrecv + total_internal + total_external + total_busywait + total_mass_reduce + total_waitsend;
	
	if(!rank) printf("Time elapsed: %f seconds\n", total_accounted_for);
	if(!rank) printf("Durations:\nSend/Recv    = %f\nInternal    = %f\nBusywait    = %f\nExternal    = %f\nReduce mass = %f\nWait sends    = %f\n\n"
		"Respective percentages (based on accounted for):\nSend/Recv    = %f\nInternal    = %f\nBusywait    = %f\nExternal    = %f\nReduce mass = %f\nWait sends    = %f\n\n",
		total_sendrecv,total_internal,total_busywait,total_external,total_mass_reduce,total_waitsend,
		100.0 * total_sendrecv   /total_accounted_for,
		100.0 * total_internal   /total_accounted_for,
		100.0 * total_busywait   /total_accounted_for,
		100.0 * total_external   /total_accounted_for,
		100.0 * total_mass_reduce/total_accounted_for,
		100.0 * total_waitsend   /total_accounted_for);
	if(rank==0)	{
		fprintf(fp,"Total/Sendrecv/internal/busywait/external/mass reduce/waitsend durations:\n%f\t%f\t%f\t%f\t%f\t%f\t%f\n\n\n",
				total_accounted_for, total_sendrecv, total_internal, total_busywait, total_external, total_mass_reduce, total_waitsend);
		fclose(fp);
	}

	if(!rank) printf("\n\n==================================\n===========SESSION END============\n==================================\n\n\n");
    MPI_Finalize();
	return 0;
}

// c coordinates to simplify the 3-dimensional array handling
int c(int x,int y,int z){
	return z*(y_length+4)*(x_length+4)+y*(x_length+4)+x;
}
	
/**************************************************************************
 *  update elements
 ****************************************************************************/
void update(int x, int y, int z, struct element *u1, struct element *u2){
	(u2+c(x,y,z))->xy_value = (	(8*(u1+c(x,y,z))->xy_value)+
					(u1+c(x-2,y,z))->xy_value+
					(u1+c(x-1,y,z))->xy_value+
					(u1+c(x+1,y,z))->xy_value+
					(u1+c(x+2,y,z))->xy_value+
					(u1+c(x,y-2,z))->xy_value+
					(u1+c(x,y-1,z))->xy_value+
					(u1+c(x,y+1,z))->xy_value+
					(u1+c(x,y+2,z))->xy_value)/16.0;	
	(u2+c(x,y,z))->z_value = a*pow((u2+c(x,y,z))->xy_value,10.0);
	return;
}

/**************************************************************************
 * Print data
 **************************************************************************/
void prtdat(char *fnam){
	int ix, iy, iu;
	FILE *fp;

	fp = fopen(fnam, "w");
	fprintf(fp,"This file contains the data of partition %d, sized %d x %d x %d:\n", rank, x_length, y_length, z_length);
	for(iu=0; iu<2; ++iu)
		for (iz = 0; iz < z_length; ++iz){
			fprintf(fp,"u[%d], z = %d:\n\n   x|        0",iu,iz);
			for (ix = 1; ix < x_length+4; ++ix) fprintf(fp,"%7d",ix);
			fprintf(fp,"\n-----------------------------------------------------------------------------------\n");
			for (iy = y_length + 3; iy >=0; --iy){
				fprintf(fp, "y %3d|    ",iy);
				for (ix = 0; ix < x_length+4; ++ix){
					fprintf(fp, "%6.0f", (u[iu]+c(ix,iy,iz))->xy_value);
					if (ix != (x_length+4)-1) 
						fprintf(fp, " ");
					else
						fprintf(fp, "\n");
				}
			}
			fprintf(fp,"\n\n");
		}
	fclose(fp);
	return;
}

/*
void calculate_grid(int* x_partitions, int* y_partitions, int x_data, int y_data, int partitions){
	int best,i,j;
	double target_ratio,current_ratio,best_ratio;

	if(partitions>x_data*y_data) partitions=x_data*y_data;//Drop cpus to data count. As if.
	target_ratio=(double)x_data/(double)y_data;
	best_ratio=0;
	for(i=1;i<=(int)sqrt((double)partitions);i++)
		if(!(partitions%i)){ //valid divisor
			current_ratio=(double)i/(double)(partitions/i);
			if(fabs(target_ratio-best_ratio)>fabs(target_ratio-current_ratio)){
				best_ratio=current_ratio;
				best=i;
			}
			current_ratio=(double)(partitions/i)/(double)i;//invert
			if(fabs(target_ratio-best_ratio)>fabs(target_ratio-current_ratio)){
				best_ratio=current_ratio;
				best=partitions/i;
			}
		}
	*x_partitions=best;
	*y_partitions=partitions/best;
	return;
}
*/
