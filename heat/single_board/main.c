/*
 *	HOST program - Single Parallella board
 *	Heat transfer calculation based on Bondeli's algorithm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>	//fabs
#include <sys/time.h>

#include "e-hal.h"
#include "common.h"

shared_buf_t 	Mailbox;
e_mem_t		DRAM;

e_platform_t 	platform;
e_epiphany_t	dev;
struct timeval	timer[4], start, stop;	
double		tdiff[2], speedup;
double 		diff, tot=0, tot_ep=0, tot_r=0, tot_pd0=0, tot_b=0, tot_gy=0, tot_al=0, tot_un=0;
unsigned int	addr;
FILE *fp;


float	a[_M];		// | b1 c1          |   |x1|   |d1|
float	b[_M];		// | a2 b2 c2       |   |x2|   |d2|
float	c[_M];		// |    a3 b3 c3    | * |x3| = |d3|
float	d[2][_M];	// |       a4 b4 c4 |   |x4|   |d4|

int 	prev=0, next=1;
float	v[_M];		// store result for later validation

float	hyy[_M];	
float	hza[_M];	
float	hzb[_M];


float gy_alti[_P];	
float gy_bassi[_P];	
float gz_alti[(_P*2)-2];
float gz_bassi[(_P*2)-2];

float ss[(_P*2)-2];
float rr[(_P*2)-2];
float tt[(_P*2)-2];
float alphas[(_P*2)-2];


void generate_matrix_T();
void generate_vector_d();
void prevtore_vector_c();
void send_go_signal();
void solve(float* a, float* b, float* c, float* d, unsigned p);
int  validate_prevult(float *a,float *b);


int main(int argc, char *argv[])
{
	int t,i;
	int j,h;
	int sweep, s;
	int row,col,pos;

float perc;
int nchunks;


	// Initalize Epiphany device
	e_init(NULL);                      
	e_reset_system();                                      
	e_get_platform_info(&platform);
	
	// Allocate memory for the Mailbox on the epiphany DRAM (max 32MB)
	e_alloc(&DRAM, 0x00000000 , sizeof(Mailbox));

	// set flag 
	Mailbox.go=0;
	addr = offsetof(shared_buf_t, go);
	e_write(&DRAM, 0, 0, addr, &Mailbox.go, sizeof(Mailbox.go));
	
	// Open Coprev on epiphany chip
	e_open(&dev, 0, 0, _NROWS, _NCOLS);

	// Load program to epiphany coprev and start execution
	fprintf(stdout, "\nLoading programs to [%d] Epiphany core(s)... [%dx%d] [%d eq.]\n",_NCORES,_MSIDE,_MSIDE,_M);
	e_load_group("./bin/core.srec", &dev, 0, 0, _NROWS, _NCOLS, E_TRUE);
	fprintf(stdout, "running.. \n\n");


/* -----------------------------------------------------------------------------------------*/
/* HOST+EPIPHANY CALCULATION								    */
/* -----------------------------------------------------------------------------------------*/

	// set timer
	gettimeofday(&timer[0], NULL);

	// Prepare Matrix [T] and initial temperature vector [d]
	generate_matrix_T();
	generate_vector_d();

	// Iterations 
	for (t=0; t<_ITERATIONS; t++) {

	// 2 sweeps: Implicit over X-direction and Explicit over Y-direction
	for (sweep=0; sweep<2; sweep++) {

		// Switch temperature arrays
		next=1-next;
		prev=1-prev;

		// multiply vector d by right hand triband matrix and
		// re-arrange vector d as [d]= [d(1,1), d(2,1), d(3,1), ..d(1,2),d(2,2),d(3,2)..];
		gettimeofday(&start, NULL);
		for (i=1; i<_M-1; i++) {
			row=i/_MSIDE;
			col=i-(row*_MSIDE);
			pos=(col*_MSIDE)+row;
			d[next][pos]= _GAMMA*d[prev][i-1]+(1-(2*_GAMMA))*d[prev][i]+_GAMMA*d[prev][i+1];
		}
		d[next][0]=(1-(2*_GAMMA))*d[prev][0]+_GAMMA*d[prev][1];
		d[next][_M-1]=_GAMMA*d[prev][_M-2]+(1-(2*_GAMMA))*d[prev][_M-1];
		gettimeofday(&stop, NULL);
		tot_pd0 += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

		//split in 1024*1024bytes chunks
		for (s=0; s<(_M/_MSHR);s++) {

			// Write d vector to Mailbox
			gettimeofday(&start, NULL);
			addr = offsetof(shared_buf_t, d[0]);
			e_write(&DRAM, 0, 0, addr, (void *)&d[next][s*_MSHR], _size_MSHR); 
			gettimeofday(&stop, NULL);
			tot += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

			// send go signal
			send_go_signal();

			// wait for completed signal
			gettimeofday(&start, NULL);
			addr = offsetof(shared_buf_t, completed);
			while (Mailbox.completed == 0)
				e_read(&DRAM, 0, 0, addr, &Mailbox.completed, sizeof(int));	
			gettimeofday(&stop, NULL);
			tot_ep += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);

			// get Y, ZA, ZB from Mailbox
			gettimeofday(&start, NULL);
			addr = offsetof(shared_buf_t, tyy);
			e_read(&DRAM, 0, 0, addr, (void*)&hyy[s*_MSHR], _size_MSHR);
			addr = offsetof(shared_buf_t, tza);
			e_read(&DRAM, 0, 0, addr, (void*)&hza[s*_MSHR], _size_MSHR);
			addr = offsetof(shared_buf_t, tzb);
			e_read(&DRAM, 0, 0, addr, (void*)&hzb[s*_MSHR], _size_MSHR);
			gettimeofday(&stop, NULL);
			tot_r += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);
		}


		// calculate GY and GZ
		gettimeofday(&start, NULL);
		j=0;
		for (i=0; i<_M; i+=_MCORE) {
			gy_alti[j]=hyy[i];
			gy_bassi[j]=hyy[i+_MCORE-1];
			j++;
		}
		j=0; h=0;
		for (i=0; i<_M; i+=_MCORE) {
			gz_alti[h]=hza[i];
			gz_bassi[h]=hza[i+_MCORE-1];
			
			if( (j!=0) && (j!=(_P-1)) ) {
				gz_alti[h+1]=hzb[i];
				gz_bassi[h+1]=hzb[i+_MCORE-1];
			}
			j++; h=(j*2)-1;
		}
		gettimeofday(&stop, NULL);
		tot_gy += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);	


		// calculate alphas
		gettimeofday(&start, NULL);
		for (i=0; i<(2*_P)-2; i++) {

			if (i%2 == 0) {
				ss[i]=gz_bassi[i];
				if(i!=0) rr[i]=gz_bassi[i-1];
				if(i!=(2*_P)-3) tt[i]= 1/(-1*_GAMMA);
				alphas[i]=-1*(gy_bassi[i/2]);
			} else {
				ss[i]=(gz_alti[i]);
				if(i!=0) rr[i]= 1/(-1*_GAMMA);
				if(i!=(2*_P)-3) tt[i]=gz_alti[i+1];
				alphas[i]=-1*(gy_alti[(i+1)/2]);
			}
		}
		solve(rr,ss,tt,alphas,(2*_P)-2);
		gettimeofday(&stop, NULL);
		tot_al += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);


		// calculate unknowns
		gettimeofday(&start, NULL);
		for (i=1; i<_P-1; i++) {
			for (h=0; h<_MCORE; h++) {
				j=(i*_MCORE)+h;
				d[next][j]=hyy[j] + alphas[(2*i)-1]*hza[j] + alphas[(2*i)]*hzb[j];
			}
		}
		i=0;
		for (h=0; h<_MCORE; h++) {
			j=(i*_MCORE)+h;
			d[next][j]=hyy[j] + alphas[0]*hza[j]; 
		}
		i=_P-1;
		for (h=0; h<_MCORE; h++) {
			j=(i*_MCORE)+h;
			d[next][j]=hyy[j] + alphas[(2*i)-1]*hza[j];
		}
		gettimeofday(&stop, NULL);
		tot_un += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);	

		
	}


	
		// -------------------------------------------------------------------------------------
		// APPLY BOUNDARY CONDITIONS AND HEAT SOURCE
		// -------------------------------------------------------------------------------------
		gettimeofday(&start, NULL);

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

		gettimeofday(&stop, NULL);
		tot_b += (stop.tv_sec - start.tv_sec) * 1000 + ((double) (stop.tv_usec - start.tv_usec) / 1000.0);


	} // END ITERATIONS



	// stop timer
	gettimeofday(&timer[1], NULL);

	// Calculate Elapsed time
	tdiff[0] = (timer[1].tv_sec - timer[0].tv_sec) * 1000 + ((double) (timer[1].tv_usec - timer[0].tv_usec) / 1000.0);
	fprintf(stdout, "Host+Epiphany -  time: %9.3f msec \n", tdiff[0]);

	perc=(tot_pd0/tdiff[0])*100;
	fprintf(stdout, "               prep D: %9.3f msec  (%.1f)\n", tot_pd0,perc);
	perc=(tot/tdiff[0])*100;
	fprintf(stdout, "              write D: %9.3f msec  (%.1f)\n", tot,perc);
	perc=(tot_ep/tdiff[0])*100;
	fprintf(stdout, "             epiphany: %9.3f msec  (%.1f)\n", tot_ep,perc);
	perc=(tot_r/tdiff[0])*100;
	fprintf(stdout, "             read Y,Z: %9.3f msec  (%.1f)\n", tot_r,perc);
	perc=(tot_gy/tdiff[0])*100;
	fprintf(stdout, "           calc GY,GZ: %9.3f msec  (%.1f)\n", tot_gy,perc);
	perc=(tot_al/tdiff[0])*100;
	fprintf(stdout, "          calc alphas: %9.3f msec  (%.1f)\n", tot_al,perc);
	perc=(tot_un/tdiff[0])*100;
	fprintf(stdout, "             unknowns: %9.3f msec  (%.1f)\n", tot_un,perc);
	perc=(tot_b/tdiff[0])*100;
	fprintf(stdout, "             boundary: %9.3f msec  (%.1f)\n", tot_b,perc);
	


	fprintf(stdout, "\n");


	// STORE prevULT FOR LATER VALIDATION
	for (i=0; i<_M; i++) v[i]=d[next][i];
 

	
/* -----------------------------------------------------------------------------------------*/
/* HOST ONLY CALCULATION								    */
/* -----------------------------------------------------------------------------------------*/
	
	// set timer
	gettimeofday(&timer[2], NULL);

	// Prepare Matrix [T] and initial temperature vector [d]
	generate_matrix_T();
	generate_vector_d();

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
			solve(a, b, c, d[next], _M); 
			prevtore_vector_c();

		}

		// Apply Boundary conditions
		for (i=0; i<_MSIDE; i++) {
			d[prev][i]=0;			  // top
			d[prev][_M-1-i]=0;		  // bottom
			d[prev][i*_MSIDE]=0;		  // left
			d[prev][(i*_MSIDE)+_MSIDE-1]=0; // right
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
	fprintf(stdout, "Host only     -  time: %9.3f msec \n", tdiff[1]);
	speedup=tdiff[1]/tdiff[0];
	fprintf(stdout, "\nSpeedup: %5.2fx \n", speedup);


/* -----------------------------------------------------------------------------------------*/
/* COMPARE prevULTS									    */
/* -----------------------------------------------------------------------------------------*/


	// Check the prevults
	validate_prevult(d[next],v);

	
	// Write to file
	fp = fopen("heat.csv", "w");
	if (fp == NULL) { fprintf(stderr, "Can't open input file in.list!\n"); exit(1); }
	for (i=0; i<_M-1; i++) fprintf(fp, "%.1f;", v[i]);
	fprintf(fp, "%.1f", v[_M-1]);
	fclose(fp);

	// Close down Epiphany device
	e_close(&dev);
	e_finalize();
	return 0;
}







void generate_matrix_T()
{
	int i;

	// Matrix T ( diagonals a,b,c )
	for (i=0; i<_M; i++) a[i] = -1*_GAMMA;
	for (i=0; i<_M; i++) b[i] =  1+(2*_GAMMA);
	for (i=0; i<_M; i++) c[i] = -1*_GAMMA;
}

void prevtore_vector_c()
{
	int i;
	for (i=0; i<_M; i++) c[i] = -1*_GAMMA;
}

void generate_vector_d()
{
	int i,p,row,col;

	p=_MSIDE/16;
	
	// Vector d
	for (i=0; i<_M; i++) {
		row=i/_MSIDE;
		col=i-(row*_MSIDE);
		if ((row>(7*p))&&(row<(9*p))&&(col>(3*p))&&(col<(13*p)) )
		{ d[next][i] = 50; }
		else { d[next][i] = 0; }
	}

	//for (i=0; i<_M; i++) d[i]=(i/4)+1;
	//for (i=0; i<_M; i++) d[i]=i+1;
}

void send_go_signal()
{
	// Completed [0]
	Mailbox.completed=0;
	addr = offsetof(shared_buf_t, completed);
	e_write(&DRAM, 0, 0, addr, &Mailbox.completed, sizeof(int));

	// Go [1]
	Mailbox.go=1;
	addr = offsetof(shared_buf_t, go);
	e_write(&DRAM, 0, 0, addr, &Mailbox.go, sizeof(int));
}



/*
 * THOMAS algorithm
 */
void solve(float* a, float* b, float* c, float* d, unsigned p) {

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


/*
 * Compare Two Vectors, returns 1 if equal
 */
int validate_prevult(float *d,float *v)
{
	int i;
	
	for (i=0; i<_M; i++) {
		
		if (fabs( d[i] - v[i] ) > 0.01) { 
		//if (d[i] != x[i]) {
			fprintf(stdout, "\n\nD[%d]: %.5f     V[%d]: %.5f\n", i, d[i], i, v[i]);
			fprintf(stdout, "\nCheck: MISTMATCHING RESULTS!\n\n");
			return 0;
		}
	}
		
	fprintf(stdout, "Check  : OK\n\n");	
	return 1;
}


