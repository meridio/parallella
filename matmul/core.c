/*
 *	CORE
 */

#include "e-lib.h"
#include "common.h"

#define _size_M		sizeof(float)*_M
#define _size_MCORE	sizeof(float)*_MCORE
#define _PING   0
#define _PONG   1

typedef struct {
	unsigned coreID;
	unsigned corenum;
	unsigned row, col;

	void *local_Bank_A[2]; 
	void *local_Bank_B[2];
	void *target_Bank_A[2];
	void *target_Bank_B[2];

	int    pingpong;  
} core_t;


volatile float  AA[2][_MCORE*_MCORE]	SECTION(".data_bank1");  // local A submatrix
volatile float  BB[2][_MCORE*_MCORE]	SECTION(".data_bank2");  // local B submatrix
volatile float  CC   [_MCORE*_MCORE]	SECTION(".data_bank3");  // local C submatrix

// Explicit placement of static objects in memory. The three matrices
// the linker may relocate the objects wherever within the bank. The core structure "me" 
// is specifically located at an explicit address - 0x7000.
// To do that, a custom linker file (LDF) was defined adding "section_core" data section

volatile shared_buf_t *Mailbox		SECTION("section_core"); // Pointers to Mailbox
volatile e_barrier_t  barriers[_Ncores] SECTION("section_core"); // barriers array
         e_barrier_t *tgt_bars[_Ncores] SECTION("section_core"); // barriers array
         core_t me                      SECTION("section_core"); // core data structure

// function prototypes
void init();
void do_matrix_multiplication();
void calculate_local_result(volatile float *a, volatile float *b, volatile float *c);
void data_copy(e_dma_desc_t *dma_desc, void *dst, void *src);


int main(int argc, char *argv[])
{
	// Initialization (whoami and addresses)
	init();

	// Initialize the barriers
	e_barrier_init(barriers, tgt_bars);
	
	// Wait for start signal (mailbox)
	if (me.corenum == 0)
	{
		Mailbox->dev_ready = 1;			// ALL cores initialized
		while (Mailbox->hst_ready == 0) {};	// Wait for GO signal
		Mailbox->hst_ready = 0;
	}

	// Sync with all other cores
	e_barrier(barriers, tgt_bars);

	// local execution
	do_matrix_multiplication();

	// Sync with all other cores
	e_barrier(barriers, tgt_bars);

	// Signal End-Of-Calculation to the host.
	if (me.corenum == 0) Mailbox->completed = 1;
	return 0;
}


void init()
{
	unsigned rowh,colh;	// coordinates of neighbor on the left
	unsigned rowv,colv;	// coordinates of neighbor on top

	// Init pointer to mailbox
	Mailbox = (shared_buf_t *) e_emem_config.base;  

	// Core info
	me.coreID  = e_get_coreid();
	me.row     = e_group_config.core_row;
	me.col     = e_group_config.core_col;
	//me.corenum = me.row * e_group_config.group_cols + me.col;
	me.corenum = me.row * _NCORES_SIDE + me.col;

	// Initialize pointers to the operand matrices ping-pong arrays
	me.local_Bank_A[_PING] = (void *) &(AA[_PING][0]); //= &AA[_PING];
	me.local_Bank_A[_PONG] = (void *) &(AA[_PONG][0]);
	me.local_Bank_B[_PING] = (void *) &(BB[_PING][0]);
	me.local_Bank_B[_PONG] = (void *) &(BB[_PONG][0]);

	// Pointers to external banks
	e_neighbor_id(E_PREV_CORE, E_ROW_WRAP, &rowh, &colh); 	
	e_neighbor_id(E_PREV_CORE, E_COL_WRAP, &rowv, &colv);	
	me.target_Bank_A[_PING] = e_get_global_address(rowh, colh, me.local_Bank_A[_PONG]);
	me.target_Bank_A[_PONG] = e_get_global_address(rowh, colh, me.local_Bank_A[_PING]);
	me.target_Bank_B[_PING] = e_get_global_address(rowv, colv, me.local_Bank_B[_PONG]);
	me.target_Bank_B[_PONG] = e_get_global_address(rowv, colv, me.local_Bank_B[_PING]);

	me.pingpong = _PING;
}


/*
 *	MATRIX MULTIPLICATION:
 *	Multiply Matrix A and Matrix B e write result back in DRAM
 */
void do_matrix_multiplication()
{
	int  i,j;
	int  im, jm, km;		
	int  ip, jp;		// coordinates in the global matrixes
	void *src, *dst;	// source and destination addresses for DMA transfers


	// For each chip block
	for (im=0; im<_M; im+=_MCHIP)
	{
		for (jm=0; jm<_M; jm+=_MCHIP)
		{
			// reset local block C result
			for (i=0; i<_MCORE*_MCORE; i++) CC[i] = 0;

			// For as many chip blocks in a row/col
			for (km=0; km<_M; km+=_MCHIP) {
			

				// get block A from DRAM
				ip= im + (_MCORE * me.row);
				jp= km + (_MCORE * me.col);
				src = (void *)&(Mailbox->A[(ip * _M) + jp]);
				dst = (void *) AA[0];
				for (i = 0; i < _MCORE; i++) e_dma_copy(dst + i*_size_MCORE, src + i*_size_M, _size_MCORE); // Time: 30ms
				
				// get block B from DRAM
				ip= km + (_MCORE * me.row);
				jp= jm + (_MCORE * me.col);
				src = (void *)&(Mailbox->B[(ip * _M) + jp]);
				dst = (void *) BB[0];
				for (i = 0; i < _MCORE; i++) e_dma_copy(dst + i*_size_MCORE, src + i*_size_M, _size_MCORE); // Time: 30ms

		
				// local submatrices multiplications
				for (i=0; i<_NCORES_SIDE; i++)
				{	
					// local matrix multiplication
					calculate_local_result(AA[me.pingpong], BB[me.pingpong], CC); // 158ms

					// Swap A_banks horizontally
					src = me.local_Bank_A[me.pingpong];
					dst = me.target_Bank_A[me.pingpong];
					e_dma_copy(dst,src,(_MCORE*_MCORE*4));  //Time: 1ms
				
					// Swap B_banks vertically
					src = me.local_Bank_B[me.pingpong];
					dst = me.target_Bank_B[me.pingpong];
					e_dma_copy(dst,src,(_MCORE*_MCORE*4));	//Time: 1ms

					// Sync with all other cores
					e_barrier(barriers, tgt_bars);

					me.pingpong = 1 - me.pingpong;
				}

			}

			// Write the block C result to DRAM
			ip= im + (_MCORE * me.row);
			jp= jm + (_MCORE * me.col);		
			dst = (void *)&(Mailbox->C[(ip * _M) + jp]);
			src = (void *) CC;
			for (i = 0; i < _MCORE; i++) e_dma_copy(dst + i*_size_M, src + i*_size_MCORE, _size_MCORE); // Time 7ms

		}
	} 
	
} 


/*
 *	CALCULATE LOCAL RESULT [CC]=[AA]*[BB]:
 *	Intermediate C values are not reset for cumulative result
 */
void calculate_local_result(volatile float *a, volatile float *b, volatile float *c)
{
	int i, j, k;

	for (i = 0; i < _MCORE; i++)
		for (j = 0; j < _MCORE; j++) {
			//c[i*_MCORES + j] = 0;
			for (k = 0; k < _MCORE; k++)
				c[i*_MCORE + j] += a[i*_MCORE + k] * b[k*_MCORE + j];
		}
}


