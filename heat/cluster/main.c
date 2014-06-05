/*
 *	HOST program - Cluster 
 *	Heat transfer calculation based on Bondeli's algorithm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>	//fabs
#include <sys/time.h>

#include <mpi.h>

#include "e-hal.h"
#include "common.h"

shared_buf_t 	Mailbox;
e_mem_t		DRAM;

float d[2][_M];		// global D vector
float v[_M];
int   prev=0, next=1;

float	ld[_M/_NBOARDS];	
float	ldtmp[_M/_NBOARDS];

float	hyy[_M/_NBOARDS];	
float	hza[_M/_NBOARDS];	
float	hzb[_M/_NBOARDS];

float	p_ld, f_ld;	// preceding and following temperature value

float gy_alti[_P];	
float gy_bassi[_P];	
float gz_alti[(_P*2)];	// first and last element to be ignored
float gz_bassi[(_P*2)];	// first and last element to be ignored

float lgy_alti[_P/_NBOARDS];	
float lgy_bassi[_P/_NBOARDS];	
float lgz_alti[(_P*2)/_NBOARDS];
float lgz_bassi[(_P*2)/_NBOARDS];

float ss[(_P*2)-2];
float rr[(_P*2)-2];
float tt[(_P*2)-2];
float alphas[(_P*2)-2];

char nodename[6];
int rank, numprocs;

struct timeval	timer[4], start, stop, brr;	
double		tdiff[2], speedup;
double 		diff, tot_pd=0, tot_s1=0, tot_s2=0, tot_lc=0, tot_gv=0, tot_al=0, tot_ba=0, tot_cu=0, tot_gu=0, tot_bc=0, tot_rd=0, tot_mu=0, tot_gb=0;
double		tot_l1=0, tot_l2=0, tot_l3=0, tot_l4=0, tot_ed=0;
float perc;



//TEST
float a[_M], b[_M], c[_M];


void initialize_d_vector();
void initialize_local_d_vector();
void rearrange_d_vector();
void send_rank_info();
void local_computation();
void solve_tridiagonal( float* a,  float* b,  float* c,  float* d, unsigned p);
void validate_results();
int  compare_res(float *d,float *v);
void print_matrix(float *d, int s, int r);


int main(int argc, char *argv[]) 
{
	const int MASTER = 0;
	int TAG = 1;

	int t;
	int i,j,k,pos;
	int rc;
	int sweep;
	int row,col;
	int muld=1;
	int dst,src;

	// openMPI
	MPI_Request 	dreq[_MSIDE];
	MPI_Request 	sreq[2];
	MPI_Request 	rreq[2];
	
	// Epiphany
	e_platform_t 	platform;
	e_epiphany_t	dev;

	// get nodename from file
	FILE * fp;	
	fp = fopen( "./bin/id.txt" , "r");
	rc = fscanf(fp, "%s", nodename);
	if(rc==0) return 0;
	fclose(fp);

	// Init openMPI
	MPI_Init(&argc,&argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	//vector type
	MPI_Datatype gcolumn,lcolumn; 
	MPI_Type_vector(_MSIDE, 1, _MSIDE, MPI_FLOAT, &gcolumn);
	MPI_Type_commit(&gcolumn);
	MPI_Type_vector(_BSIDE, 1, _MSIDE, MPI_FLOAT, &lcolumn);
	MPI_Type_commit(&lcolumn); 


	// Initalize Epiphany device
	e_init(NULL);                      
	e_reset_system();                                      
	e_get_platform_info(&platform);
	e_alloc(&DRAM, 0x00000000 , sizeof(Mailbox));
	send_rank_info();
	e_open(&dev, 0, 0, _NROWS, _NCOLS);
	if (rank == MASTER) fprintf(stdout, "\n[%d BOARDS] Loading programs to [%d] Epiphany core(s)... [%dx%d] [%d eq.]\n\n",_NBOARDS,_NCORES,_MSIDE,_MSIDE,_M);
	e_load_group("./bin/core.srec", &dev, 0, 0, _NROWS, _NCOLS, E_TRUE);

	// set timer
	gettimeofday(&timer[0], NULL);

	// Prepare initial Temperature array LOCALLY
	initialize_local_d_vector();


for (t=0; t<_ITERATIONS; t++) {

	// 2 sweeps: Implicit over X-direction and Explicit over Y-direction
	for (sweep=0; sweep<2; sweep++) {

	
		// local computation
		//gettimeofday(&start, NULL);
		local_computation();
		//gettimeofday(&stop, NULL);
		//tot_lc += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

		//gettimeofday(&start, NULL);
		MPI_Barrier(MPI_COMM_WORLD);
		//gettimeofday(&stop, NULL);
		//tot_s1 += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

		// collect Yalti/bassi, Zalti/bassi
		//gettimeofday(&start, NULL);
		MPI_Gather(&lgy_alti, _P/_NBOARDS, MPI_FLOAT, &gy_alti, _P/_NBOARDS, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
		MPI_Gather(&lgy_bassi, _P/_NBOARDS, MPI_FLOAT, &gy_bassi, _P/_NBOARDS, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
		MPI_Gather(&lgz_alti, (_P*2)/_NBOARDS, MPI_FLOAT, &gz_alti, (_P*2)/_NBOARDS, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
		MPI_Gather(&lgz_bassi, (_P*2)/_NBOARDS, MPI_FLOAT, &gz_bassi, (_P*2)/_NBOARDS, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
		//gettimeofday(&stop, NULL);
		//tot_gv += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

		
		//gettimeofday(&start, NULL);
		if (rank == MASTER) {
		
			// define tridiagonal system
			for (i=0; i<(2*_P)-2; i++) {

				if (i%2 == 0) {

					ss[i]=gz_bassi[i+1]; 
					if(i!=0) rr[i]=gz_bassi[i];
					if(i!=(2*_P)-3) tt[i]= 1/(-1*_GAMMA);
					alphas[i]=-1*(gy_bassi[i/2]);

				} else {

					ss[i]=gz_alti[i+1]; 
					if(i!=0) rr[i]= 1/(-1*_GAMMA); 
					if(i!=(2*_P)-3) tt[i]=gz_alti[i+2]; 
					alphas[i]=-1*(gy_alti[(i+1)/2]);

				}
			}

			// calculate alphas
			solve_tridiagonal(rr,ss,tt,alphas,(2*_P)-2);
		}
		
		// broadcast alphas
		MPI_Bcast(&alphas, (2*_P)-2, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
		//gettimeofday(&stop, NULL);
		//tot_al += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);


		// calculate unknowns
		//gettimeofday(&start, NULL);	
		j=0;
		for (i=0; i<(_M/_NBOARDS); i+=_MCORE) {

			//position in alphas array
			pos=(((_P/_NBOARDS)*rank)+j)*2-1;	
		
			if((j==0)&&(rank==0)) 
			{
				// absolute first core
				for (k=0; k<_MCORE; k++) ld[k]=hyy[k] + alphas[0]*hza[k];
			}
			else if ((j==(_P/_NBOARDS)-1) && (rank==numprocs-1)) 
			{
				// absolute last core
				for (k=(j*_MCORE); k<((j+1)*_MCORE); k++) ld[k]=hyy[k] + alphas[pos]*hza[k];
			}
			else {
				// all others
				for (k=(j*_MCORE); k<((j+1)*_MCORE); k++) ld[k]=hyy[k] + alphas[pos]*hza[k] + alphas[pos+1]*hzb[k];
			}

			j++;
		}
		//gettimeofday(&stop, NULL);
		//tot_cu += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);



		// Apply Boundary conditions and Heat source LOCALLY
		if (sweep==1) {
			//gettimeofday(&start, NULL);

			// Apply Boundary conditions
			if(rank==MASTER)     for (i=0; i<_MSIDE; i++) ld[i]=0; // top
			if(rank==numprocs-1) for (i=0; i<_MSIDE; i++) ld[(_M/_NBOARDS)-1-i]=0; // bottom
			for (i=0; i<_BSIDE; i++) {
				ld[i*_MSIDE]=0; // left
				ld[(i*_MSIDE)+_MSIDE-1]=0; // right
			}

			// Apply Heat source
			pos=_MSIDE/16;
			for (i=0; i<_M; i++) {
				row=i/_MSIDE;
				col=i-(row*_MSIDE);
				if ((row>(7*pos))&&(row<(9*pos))&&(col>(3*pos))&&(col<(13*pos)) ) {
					if( (i>=rank*(_M/_NBOARDS)) && (i<(rank+1)*(_M/_NBOARDS)) ) ld[i-(rank*_M/_NBOARDS)] = _HEATSOURCE; 
				}
			}
		
			//gettimeofday(&stop, NULL);
			//tot_bc += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);
		}


	
		if((t==_ITERATIONS-1) && (sweep==1)) muld=0;
		if(muld) {

			//gettimeofday(&start, NULL);

			// exchange border values between boards
			p_ld=0; f_ld=0;
/*
			if(rank!=0) MPI_Isend(&ld[0], 1, MPI_FLOAT, rank-1, TAG, MPI_COMM_WORLD, &req);
			if(rank!=(numprocs-1)) MPI_Isend(&ld[(_M/_NBOARDS)-1], 1, MPI_FLOAT, rank+1, TAG, MPI_COMM_WORLD, &req);
			if(rank!=0) MPI_Recv(&p_ld, 1, MPI_FLOAT, rank-1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(rank!=(numprocs-1)) MPI_Recv(&f_ld, 1, MPI_FLOAT, rank+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
*/
/*
			// SEND-RECEIVE Switching
			if(rank!=0) MPI_Sendrecv(&ld[0], 1, MPI_FLOAT, rank-1, TAG,
						 &p_ld, 1, MPI_FLOAT, rank-1, TAG,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(rank!=(numprocs-1)) MPI_Sendrecv(&ld[(_M/_NBOARDS)-1], 1, MPI_FLOAT, rank+1, TAG,
						 &f_ld, 1, MPI_FLOAT, rank+1, TAG,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
*/
			// ASYNCHRONOUS Send-Receive
			if(rank!=0) MPI_Isend(&ld[0], 1, MPI_FLOAT, rank-1, TAG, MPI_COMM_WORLD, &sreq[0]);
			if(rank!=(numprocs-1)) MPI_Isend(&ld[(_M/_NBOARDS)-1], 1, MPI_FLOAT, rank+1, TAG, MPI_COMM_WORLD, &sreq[1]);
			if(rank!=0) MPI_Irecv(&p_ld, 1, MPI_FLOAT, rank-1, TAG, MPI_COMM_WORLD, &rreq[0]);
			if(rank!=(numprocs-1)) MPI_Irecv(&f_ld, 1, MPI_FLOAT, rank+1, TAG, MPI_COMM_WORLD, &rreq[1]);

			//gettimeofday(&stop, NULL);
			//tot_gb += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

			// Multiply Unknowns by right hand matrix LOCALLY
			//gettimeofday(&start, NULL);
			for (i=1; i<(_M/_NBOARDS)-1; i++) {
				ldtmp[i]= _GAMMA*ld[i-1]+(1-(2*_GAMMA))*ld[i]+_GAMMA*ld[i+1];
			}
			//gettimeofday(&stop, NULL);
			//tot_mu += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

			//gettimeofday(&start, NULL);
			// Wait for border values to be received
			if(rank!=0) MPI_Wait(&rreq[0], MPI_STATUSES_IGNORE);
			if(rank!=(numprocs-1)) MPI_Wait(&rreq[1], MPI_STATUSES_IGNORE);
			// Complete right hand matrix multiplication
			ldtmp[0]=_GAMMA*p_ld+(1-(2*_GAMMA))*ld[0]+_GAMMA*ld[1];
			ldtmp[(_M/_NBOARDS)-1]=_GAMMA*ld[(_M/_NBOARDS)-2]+(1-(2*_GAMMA))*ld[(_M/_NBOARDS)-1]+_GAMMA*f_ld;
			//gettimeofday(&stop, NULL);
			//tot_gb += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);


			// Reordering: Redistridute D vector slices through boards
			//gettimeofday(&start, NULL);	
			src=-1;
			for (i=0; i<_MSIDE; i++) {
				dst=i/(_MSIDE/_NBOARDS);
				MPI_Isend(&ldtmp[i], 1, lcolumn, dst, TAG, MPI_COMM_WORLD, &dreq[i]);
			}
			for (i=0; i<_MSIDE; i++) {
				src++; if(src==numprocs) src=0;
				MPI_Recv(&ld[i*_BSIDE], _BSIDE, MPI_FLOAT, src, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			//gettimeofday(&stop, NULL);
			

			//TEST
			MPI_Barrier(MPI_COMM_WORLD);
			//gettimeofday(&brr, NULL);

			//tot_rd += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);
			//tot_s2 += (brr.tv_sec - stop.tv_sec) * 1000 + ((double) (brr.tv_usec - stop.tv_usec) / 1000.0);


		}
		else {
			// COLLECT RESULT
			//gettimeofday(&start, NULL);	
			MPI_Gather(&ld, _M/_NBOARDS, MPI_FLOAT, &d[next], _M/_NBOARDS, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
			//gettimeofday(&stop, NULL);
			//tot_gu += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);
		}


	} // end sweep

}// end iterations

	// stop timer
	gettimeofday(&timer[1], NULL);

if (rank == MASTER) {
	// Calculate Elapsed time
	tdiff[0] = (timer[1].tv_sec - timer[0].tv_sec) * 1000 + ((double) (timer[1].tv_usec - timer[0].tv_usec) / 1000.0);
	fprintf(stdout, "CLUSTER -  time: %9.3f msec \n", tdiff[0]);


	perc=(tot_lc/tdiff[0])*100;
	printf("\n           local calc: %9.3f msec  (%.1f)\n", tot_lc,perc);

	perc=(tot_l1/tdiff[0])*100;
	printf("\n                ->  write D: %9.3f msec  (%.1f)\n", tot_l1,perc);
	perc=(tot_l2/tdiff[0])*100;
	printf("                -> epiphany: %9.3f msec  (%.1f)\n", tot_l2,perc);
	perc=(tot_l3/tdiff[0])*100;
	printf("                -> read YZZ: %9.3f msec  (%.1f)\n\n", tot_l3,perc);

	perc=(tot_s1/tdiff[0])*100;
	printf("               sync_1: %9.3f msec  (%.1f)\n", tot_s1,perc);
	perc=(tot_gv/tdiff[0])*100;
	printf("           gather y,z: %9.3f msec  (%.1f)\n", tot_gv,perc);
	perc=(tot_al/tdiff[0])*100;
	printf("     calc/send alphas: %9.3f msec  (%.1f)\n", tot_al,perc);
	perc=(tot_cu/tdiff[0])*100;
	printf("         calc unknown: %9.3f msec  (%.1f)\n", tot_cu,perc);
	perc=(tot_bc/tdiff[0])*100;
	printf("           boundaries: %9.3f msec  (%.1f)\n", tot_bc,perc);
	perc=(tot_gb/tdiff[0])*100;
	printf("         exc. borders: %9.3f msec  (%.1f)\n", tot_gb,perc);
	perc=(tot_mu/tdiff[0])*100;
	printf("          mul D right: %9.3f msec  (%.1f)\n", tot_mu,perc);
	perc=(tot_rd/tdiff[0])*100;
	printf("       redistribute D: %9.3f msec  (%.1f)\n", tot_rd,perc);
	perc=(tot_s2/tdiff[0])*100;
	printf("               sync_1: %9.3f msec  (%.1f)\n", tot_s2,perc);
	perc=(tot_gu/tdiff[0])*100;
	printf("        collect final: %9.3f msec  (%.1f)\n", tot_gu,perc);
	
	printf("\n");


	// save previous result for later comparison
	for (i=0; i<_M; i++) v[i]=d[next][i];

}
	// close epiphany device
	e_close(&dev);
	e_finalize();

 	
	if (rank == MASTER) {
		
		validate_results();
		
		// Write to file
		fp = fopen("heat.csv", "w");
		if (fp == NULL) { fprintf(stderr, "Can't open input file in.list!\n"); exit(1); }
		for (i=0; i<_M-1; i++) fprintf(fp, "%.1f;", v[i]);
		fprintf(fp, "%.1f", v[_M-1]);
		fclose(fp);
	}


	MPI_Finalize();
	return 0;
}




void initialize_d_vector()
{
	int i,p,row,col;

	p=_MSIDE/16;
	
	// Vector d
	for (i=0; i<_M; i++) {
		row=i/_MSIDE;
		col=i-(row*_MSIDE);
		if ((row>(7*p))&&(row<(9*p))&&(col>(3*p))&&(col<(13*p)) )
		{ d[next][i] = _HEATSOURCE; }
		else { d[next][i] = 0; }
	}

	//for (i=0; i<_M; i++) d[next][i]=i+1;
	for (i=0; i<_M; i++) d[next][i]=0;
}

void initialize_local_d_vector()
{
	int i;

	for (i=0; i<_M/_NBOARDS; i++) ld[i]=0;	
	//for (i=0; i<_M/_NBOARDS; i++) ld[i]=1+i+rank*(_M/_NBOARDS);
}


/*
 *	Multiply D vector with right hand matrix and rearrange it
 */
void rearrange_d_vector() 
{
	int i,row,col,pos;

	for (i=1; i<_M-1; i++) {
		row=i/_MSIDE;
		col=i-(row*_MSIDE);
		pos=(col*_MSIDE)+row;
		d[next][pos]= _GAMMA*d[prev][i-1]+(1-(2*_GAMMA))*d[prev][i]+_GAMMA*d[prev][i+1];
	}
	d[next][0]=(1-(2*_GAMMA))*d[prev][0]+_GAMMA*d[prev][1];
	d[next][_M-1]=_GAMMA*d[prev][_M-2]+(1-(2*_GAMMA))*d[prev][_M-1];
}

/*
 *	Send rank Info to Epiphany chip
 */
void send_rank_info()
{
	unsigned int addr;

	Mailbox.myrank=rank;
	addr = offsetof(shared_buf_t, myrank);
	e_write(&DRAM, 0, 0, addr, (void*)&Mailbox.myrank, sizeof(int));

	Mailbox.islastrank=0;
	if(rank==(numprocs-1)) Mailbox.islastrank=1;
	addr = offsetof(shared_buf_t, islastrank);
	e_write(&DRAM, 0, 0, addr, (void*)&Mailbox.islastrank, sizeof(int));
}


/*
 *	LOCAL COMPUTATION:
 *	Each board sends a portion of the local D vector to the ephipany chip.
 *	When epihany computation is completed, gZalti and Gbassi are calulated
 */
void local_computation()
{
	int s,i,j,h;
	unsigned int addr;

	// int nchunks=(_M/_NBOARDS)/_MSHR;

	//split in 1024*1024bytes chunks
	for (s=0; s<(_M/_NBOARDS)/_MSHR;s++) { 

		// Write d vector to Mailbox
		addr = offsetof(shared_buf_t, d[0]);
		e_write(&DRAM, 0, 0, addr, (void *)&ld[s*_MSHR], _size_MSHR); 



		// send go signal
		Mailbox.completed=0;
		addr = offsetof(shared_buf_t, completed);
		e_write(&DRAM, 0, 0, addr, (void*)&Mailbox.completed, sizeof(int));
		Mailbox.go=1;
		addr = offsetof(shared_buf_t, go);
		e_write(&DRAM, 0, 0, addr, (void*)&Mailbox.go, sizeof(int));

		// wait for completed signal
		addr = offsetof(shared_buf_t, completed);
		while (Mailbox.completed == 0)
			e_read(&DRAM, 0, 0, addr, (void*)&Mailbox.completed, sizeof(int));


		// get Y, ZA, ZB from Mailbox
		addr = offsetof(shared_buf_t, tyy);
		e_read(&DRAM, 0, 0, addr, (void*)&hyy[s*_MSHR], _size_MSHR);
		addr = offsetof(shared_buf_t, tza);
		e_read(&DRAM, 0, 0, addr, (void*)&hza[s*_MSHR], _size_MSHR);
		addr = offsetof(shared_buf_t, tzb);
		e_read(&DRAM, 0, 0, addr, (void*)&hzb[s*_MSHR], _size_MSHR);
	}

	
	// extract GY values to be sent to the coordinator
	j=0;
	for (i=0; i<(_M/_NBOARDS); i+=_MCORE) {
		lgy_alti[j]=hyy[i];
		lgy_bassi[j]=hyy[i+_MCORE-1];
		j++;
	}
	
	// extract GZ values to be sent to the coordinator
	j=0; h=0;
	for (i=0; i<(_M/_NBOARDS); i+=_MCORE) {

		if((j==0)&&(rank==0)) {					// absolute fisrt core
			lgz_alti[0]=0;
			lgz_alti[1]=hza[i];
			lgz_bassi[0]=0;
			lgz_bassi[1]=hza[i+_MCORE-1];
		}
		else if ( (j==(_P/_NBOARDS)-1)&&(rank==(numprocs-1)) ) { // absolute last core
			lgz_alti[h]=hza[i];
			lgz_alti[h+1]=0;
			lgz_bassi[h]=hza[i+_MCORE-1];
			lgz_bassi[h+1]=0;
		}
		else {							// All other cases
			lgz_alti[h]=hza[i];
			lgz_alti[h+1]=hzb[i];
			lgz_bassi[h]=hza[i+_MCORE-1];
			lgz_bassi[h+1]=hzb[i+_MCORE-1];
		}

		j++; h=(j*2);
	}

}




/*
 *	SOLVE LOCAL TRIDIAGONAL SYSTEM with THOMAS alghorithm
 *	Vector C and D are overwritten, local result is stored in vector D
 */
void solve_tridiagonal( float* a,  float* b,  float* c,  float* d, unsigned p) {

	int i;

	c[0] /= b[0];
	d[0] /= b[0];

	for (i = 1; i < p-1; i++) {
		c[i] /= b[i] - a[i]*c[i-1];
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}
	
	d[p-1] = (d[p-1] - a[p-1]*d[p-2]) / (b[p-1] - a[p-1]*c[p-2]);

	for (i = p-2; i >= 0; i--) {
		d[i] -= c[i]*d[i+1];
	}

}




void validate_results()
{
	int i,t,sweep,row,col,pos;

	prev=0; next=1;

	

	// ARM CALCULATION ------------------------------------------------------------------------------

	// set timer
	gettimeofday(&timer[2], NULL);

	// Prepare Matrix [T] and initial temperature vector [d]
	for (i=0; i<_M; i++) a[i] = -1*_GAMMA;
	for (i=0; i<_M; i++) b[i] =  1+(2*_GAMMA);
	for (i=0; i<_M; i++) c[i] = -1*_GAMMA;

	initialize_d_vector();

	// Iterations 
	for (t=0; t<_ITERATIONS; t++) {

		for (sweep=0; sweep<2; sweep++) {

			// switch temperature vectors
			next=1-next;
			prev=1-prev;

			// multiply vector d by right hand triband matrix and
			// re-arrange vector d as [d]= [d(1,1), d(2,1), d(3,1), ..d(1,2),d(2,2),d(3,2)..];
			for (i=1; i<_M-1; i++) {
				row=i/_MSIDE;
				col=i-(row*_MSIDE);
				pos=(col*_MSIDE)+row;
				d[next][pos]= _GAMMA*d[prev][i-1]+(1-(2*_GAMMA))*d[prev][i]+_GAMMA*d[prev][i+1];
			}
			d[next][0]=(1-(2*_GAMMA))*d[prev][0]+_GAMMA*d[prev][1];
			d[next][_M-1]=_GAMMA*d[prev][_M-2]+(1-(2*_GAMMA))*d[prev][_M-1];

			// solve tridiagonal system;
			solve_tridiagonal(a, b, c, d[next], _M); 
			for (i=0; i<_M; i++) c[i] = -1*_GAMMA;

		}

		// Apply Boundary conditions
		for (i=0; i<_MSIDE; i++) {
			d[next][i]=0;			  // top
			d[next][_M-1-i]=0;		  // bottom
			d[next][i*_MSIDE]=0;		  // left
			d[next][(i*_MSIDE)+_MSIDE-1]=0;   // right
		}

		// Apply Heat source
		pos=_MSIDE/16;
		for (i=0; i<_M; i++) {
			row=i/_MSIDE;
			col=i-(row*_MSIDE);
			if ((row>(7*pos))&&(row<(9*pos))&&(col>(3*pos))&&(col<(13*pos)) ) d[next][i] = 50; 
		}

	} // END ITERATIONS


	// stop timer
	gettimeofday(&timer[3], NULL);

	// Calculate Elapsed time and speedup
	tdiff[1] = (timer[3].tv_sec - timer[2].tv_sec) * 1000 + ((double) (timer[3].tv_usec - timer[2].tv_usec) / 1000.0);
	fprintf(stdout, "Thomas on arm    -  time: %9.3f msec \n", tdiff[1]);
	speedup=tdiff[1]/tdiff[0];
	fprintf(stdout, "\nSpeedup: %5.2fx \n", speedup);


	// COMPARE ----------------------------------------------------------------------------------------
	compare_res(d[next],v);

}

/*
 * Compare Two Vectors, returns 1 if equal
 */
int compare_res(float *d,float *v)
{
	int i;
	
	for (i=0; i<_M; i++) {
		
		if (fabs( d[i] - v[i] ) > 0.001) { 
		//if (d[i] != x[i]) {
			fprintf(stdout, "\n\nD[%d]: %.5f     V[%d]: %.5f\n", i, d[i], i, v[i]);
			fprintf(stdout, "\nCheck: MISTMATCHING RESULTS!\n\n");
			return 0;
		}
	}
		
	fprintf(stdout, "Check  : OK\n\n");	
	return 1;
}


void print_matrix(float *d, int s, int r)
{
	int i;
	for (i=0; i<s; i++) {
		if (i%r==0) printf("\n");
		printf("[%.2f]",d[i]);
	}
	printf("\n\n");
}

