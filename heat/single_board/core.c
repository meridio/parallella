/*
 *	CORE
 */

#include "e-lib.h"
#include "common.h"


typedef struct {
	unsigned corenum;
	unsigned master;
} core_t;


volatile shared_buf_t *Mailbox;
volatile e_barrier_t  barriers[_NCORES];
	 e_barrier_t *tgt_bars[_NCORES];
	 core_t me;

// local
float aa[_MCORE];
float bb[_MCORE];
float cc[_MCORE];

// input/output
float dd[_MCORE];


// function prototypes
void init();
void generate_local_matrix_T();
void reset_cc();
void copy_vector(float* a, float* b);
void unit_vector_1(float* a);
void unit_vector_k(float* a);
void solve_tridiagonal(float* a, float* b, float* c, float* d, unsigned msize);


int main(int argc, char *argv[])
{
	int nsteps;
	int offset;
	int n;
	unsigned row,col;
	

	// Initialization (whoami and pointers)
	init();

	// calculate local input
	generate_local_matrix_T();

	nsteps=_MSHR/(_NCORES*_MCORE);
	
	while(1) {

		row = e_group_config.core_row;
		col = e_group_config.core_col;
		me.corenum = row * _NCOLS + col;

		// Wait for GO signal from Host
		if (me.master) {
			while (Mailbox->go == 0) {};
			Mailbox->go = 0;
		}

		// Sync with all other cores
		e_barrier(barriers, tgt_bars);


		// Iterate
		for(n=0; n < nsteps; n++) {

			// read a portion of vector [d] from Mailbox
			offset=(me.corenum*_MCORE);
			e_dma_copy(&dd, (void*)&Mailbox->d[offset], _size_MCORE);

			// calculate YY: solve Ty=d
			solve_tridiagonal(aa,bb,cc,dd,_MCORE);

			// send current portion of vector YY to Host
			offset=(me.corenum*_MCORE)*sizeof(float);
			e_dma_copy((void *)&(Mailbox->tyy)+offset, (void*)&dd[0], _size_MCORE);

			// calculate ZA and ZB
			offset=(me.corenum*_MCORE)*sizeof(float);
			if (me.corenum == 0) {					
				unit_vector_k(dd);
				solve_tridiagonal(aa,bb,cc,dd,_MCORE); //T1ZA=ek
				e_dma_copy((void *)&(Mailbox->tza)+offset, (void*)&dd[0], _size_MCORE);

			} else if (me.corenum == _P-1) {			
				unit_vector_1(dd);
				solve_tridiagonal(aa,bb,cc,dd,_MCORE); //TpZAp-2=e1 	
				e_dma_copy((void *)&(Mailbox->tza)+offset, (void*)&dd[0], _size_MCORE);

			} else {						
				unit_vector_1(dd);
				solve_tridiagonal(aa,bb,cc,dd,_MCORE); // TiZA=e1
				e_dma_copy((void *)&(Mailbox->tza)+offset, (void*)&dd[0], _size_MCORE);
				
				unit_vector_k(dd);
				solve_tridiagonal(aa,bb,cc,dd,_MCORE); // TiZB=ek
				e_dma_copy((void *)&(Mailbox->tzb)+offset, (void*)&dd[0], _size_MCORE);		
			}


			// set next step coreid
			me.corenum += _NCORES;

		} // End step iterations

		// Sync with all other cores
		e_barrier(barriers, tgt_bars);

		// Signal End-Of-Calculation to the host.
		if (me.master) Mailbox->debug = nsteps;
		if (me.master) Mailbox->completed = 1;
	}
	
	return 0;
}




/*
 *	INIT: detect core number and define pointers to core_0
 */
void init()
{
	unsigned row,col;

	// Init pointer to mailbox
	Mailbox = (shared_buf_t *) e_emem_config.base;
	//Mailbox = (shared_buf_t *) BASE_ADDRESS;

	// Core info
	row = e_group_config.core_row;
	col = e_group_config.core_col;
	me.corenum = row * _NCOLS + col;

	me.master=0;
	if(me.corenum==0) me.master=1;

	// Initialize the barriers
	e_barrier_init(barriers, tgt_bars);
}


/*
 *	Calculate local initial values
 */
void generate_local_matrix_T()
{
	int i;
	// Matrix T
	for (i=0; i<_MCORE; i++) aa[i] = -1*_GAMMA;
	for (i=0; i<_MCORE; i++) bb[i] =  1+(2*_GAMMA);
	for (i=0; i<_MCORE; i++) cc[i] = -1*_GAMMA;
}

void reset_cc()
{
	int i;
	for (i=0; i<_MCORE; i++) cc[i] = -1*_GAMMA;
}

void copy_vector(float* a, float* b) 
{
	int i;
	for (i=0; i<_MCORE; i++) b[i] = a[i];
}

void unit_vector_1(float* a)
{
	int i;
	for (i=1; i<_MCORE; i++) a[i] = 0;
	a[0] = 1;
}

void unit_vector_k(float* a)
{
	int i;
	for (i=0; i<_MCORE-1; i++) a[i] = 0;
	a[_MCORE-1] = 1;
}



/*
 *	SOLVE LOCAL TRIDIAGONAL SYSTEM with THOMAS alghorithm
 *	Vector C and D are overwritten, local result is stored in vector D
 */
//void solve_tridiagonal(volatile float* a, volatile float* b, volatile float* c, volatile float* d, unsigned p) {
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

	reset_cc();
}














