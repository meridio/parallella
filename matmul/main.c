/*
 *	HOST
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/time.h>

#include "e-hal.h"
#include "common.h"

shared_buf_t 	Mailbox;
e_mem_t		DRAM;

void matrix_init();
void copy_data_to_epiphany();
void calculate_local_result(float *a, float *b, float *c);
void print_matrix(float *m);
int  validate_reslt(float *a,float *b);


int main(int argc, char *argv[])
{
	e_platform_t 	platform;
	e_epiphany_t	dev;
	float		R[_M * _M];	// R matrix to store the result calculated by the host
	struct timeval	timer[4];	
	double		tdiff[2], speedup;
	unsigned int	addr;

	// max buffer size 32MB	
	// fprintf(stdout, "\nShared Buffer size: %u Bytes\n",sizeof(Mailbox));
	
	// Initialize operand matrices
	matrix_init();

	//Initalize Epiphany device
	e_init(NULL);                      
	e_reset_system();                                      
	e_get_platform_info(&platform);
	
	// Allocate memory for the Mailbox on the epiphany DRAM (max 32MB)
	e_alloc(&DRAM, 0x00000000 , sizeof(Mailbox));

	// Copy operand matrices and flags to DRAM
	copy_data_to_epiphany();
	
	// Open all Cores on epiphany chip
	e_open(&dev, 0, 0, _NCORES_SIDE, _NCORES_SIDE);

	// Load program to cores and run
	fprintf(stdout, "\nLoading programs to [%d] Epiphany core(s)... [%d]x[%d]\n\n",_NCORES,_M,_M);
	e_load_group("./bin/core.srec", &dev, 0, 0, _NCORES_SIDE, _NCORES_SIDE, E_TRUE);
	

	// Wait for dev_ready flag
	addr = offsetof(shared_buf_t, dev_ready);
	while (Mailbox.dev_ready == 0)
		e_read(&DRAM, 0, 0, addr, &Mailbox.dev_ready, sizeof(Mailbox.dev_ready));	


	// Start calculation on Epiphany system by sending the GO signal
	Mailbox.hst_ready = 1;
	addr = offsetof(shared_buf_t, hst_ready);
	e_write(&DRAM, 0, 0, addr, &Mailbox.hst_ready, sizeof(Mailbox.hst_ready));

	// Save time: START Epiphany
	gettimeofday(&timer[0], NULL);

	// Busy wait: Poll the Mailbox DONE signal every second
	addr = offsetof(shared_buf_t, completed);
	while (Mailbox.completed == 0)
		e_read(&DRAM, 0, 0, addr, &Mailbox.completed, sizeof(Mailbox.completed));

	// Save time: STOP Epiphany
	gettimeofday(&timer[1], NULL);

	// Calculate Elapsed time
	tdiff[0] = (timer[1].tv_sec - timer[0].tv_sec) * 1000 + ((double) (timer[1].tv_usec - timer[0].tv_usec) / 1000.0);
	fprintf(stdout, "Epiphany -  time: %9.3f msec \n", tdiff[0]);


	// Read the result from DRAM
	addr = offsetof(shared_buf_t, C[0]);
	e_read(&DRAM, 0, 0, addr, &Mailbox.C, sizeof(Mailbox.C));


	// local calculation
	gettimeofday(&timer[2], NULL);
	calculate_local_result(Mailbox.A, Mailbox.B, R);
	gettimeofday(&timer[3], NULL);


	// Calculate Elapsed time and speedup
	tdiff[1] = (timer[3].tv_sec - timer[2].tv_sec) * 1000 + ((double) (timer[3].tv_usec - timer[2].tv_usec) / 1000.0);
	fprintf(stdout, "Host     -  time: %9.3f msec \n", tdiff[1]);
	speedup=tdiff[1]/tdiff[0];
	fprintf(stdout, "\nSpeedup: %5.2fx \n", speedup);

	
	// Check the results
	if (validate_result(Mailbox.C,R)) { fprintf(stdout, "  check: OK\n\n"); }
	else { fprintf(stdout, "  check: MISTMATCHING RESULTS!\n\n");	}


	// Close down Epiphany device
	e_close(&dev);
	e_finalize();
	return 0;
}




//
// Initialize operand matrices
//
void matrix_init()
{
	int i, j;

	// Matrix initialization
	for (i=0; i<_M; i++)
		for (j=0; j<_M; j++)
			Mailbox.A[i*_M+j] = (i + j + 1) % 32;

	for (i=0; i<_M; i++)
		for (j=0; j<_M; j++)
			Mailbox.B[i*_M+j] = ((i + j) * 2) % 32;

	for (i=0; i<_M; i++)
		for (j=0; j<_M; j++)
			Mailbox.C[i*_M+j] = 0;

	// Flags initialization
	Mailbox.dev_ready=0;
	Mailbox.hst_ready=0;
	Mailbox.completed = 0;

	return;
}

//
// Copy operand matrices to Epiphany system (DRAM)
//
void copy_data_to_epiphany()
{
	int m=_M;
	size_t sz;
	unsigned int addr;

	// Copy Matrices to DRAM
	addr = offsetof(shared_buf_t, A[0]);
	sz = sizeof(Mailbox.A);
	//fprintf(stdout, "Copying Matrix A [%d x %d] (%u Bytes) to DRAM address %08x... \n", m, m, sz, addr);
	e_write(&DRAM, 0, 0, addr, (void *) Mailbox.A, sz);

	addr = offsetof(shared_buf_t, B[0]);
	sz = sizeof(Mailbox.B);
	//fprintf(stdout, "Copying Matrix B [%d x %d] (%u Bytes) to DRAM address %08x... \n", m, m, sz, addr);
	e_write(&DRAM, 0, 0, addr, (void *) Mailbox.B, sz);

	addr = offsetof(shared_buf_t, C[0]);
	sz = sizeof(Mailbox.C);
	//fprintf(stdout, "Copying Matrix C [%d x %d] (%u Bytes) to DRAM address %08x... \n", m, m, sz, addr);
	e_write(&DRAM, 0, 0, addr, (void *) Mailbox.C, sz);


	// Copy FLAGS to DRAM
	addr = offsetof(shared_buf_t, dev_ready);
	e_write(&DRAM, 0, 0, addr, &Mailbox.dev_ready, sizeof(Mailbox.dev_ready));
	
	addr = offsetof(shared_buf_t, hst_ready);
	e_write(&DRAM, 0, 0, addr, &Mailbox.hst_ready, sizeof(Mailbox.hst_ready));
	
	addr = offsetof(shared_buf_t, completed);
	e_write(&DRAM, 0, 0, addr, &Mailbox.completed, sizeof(Mailbox.completed));
}

//
// Calculare local result [R]=[A]*[B]
//
void calculate_local_result(float *a, float *b, float *c)
{
	int i, j, k;

	for (i = 0; i < _M; i++)
		for (j = 0; j < _M; j++) {
			c[i*_M + j] = 0;
			for (k = 0; k < _M; k++)
				c[i*_M + j] += a[i*_M + k] * b[k*_M + j];
		}
}

//
// Print Matrix
//
void print_matrix(float *m)
{
	int i, j;

	for (i=0; i<_M; i++)
	{
		for (j=0; j<_M; j++) { fprintf(stdout, "%2.0f ", m[i*_M+j]); }
		fprintf(stdout, "\n");
	}
}

//
// Compare Two Matrices, returns 1 if equal
//
int validate_result(float *a,float *b)
{
	int i, j;

	for (i=0; i<_M; i++)
		for (j=0; j<_M; j++)
			// if (fabs(a[j+j*i] - b[j+j*i]) > 0.00001) { return 0; }
			if (a[i*_M+j] != b[i*_M+j]) { return 0; }

	return 1;
}


