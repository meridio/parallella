
#define _NCLUSTERS	1	// # of epiphany board connected in a cluster [1= single board]
#define _NCORES		16	// # of cores in epiphany chip GROUP [options: 1,4,16]
#define _NCORES_SIDE	4	// # of cores in epiphany chip GROUP side [sqrt of _NCORES]

#define _M		1024	// Side of operand Matrices [M x M]
				// Multiple of: _MCORE * _NCORE_SIDE (in this case [128])
				// Max value to have 3 matrices fitting in 32MB DRAM: [1664]

#define _MCORE		32	// side size of per-core sub-submatrix 
				// Max value to store 2 copies of the matrix in a core: [32]
				// 32 x 32 x 4 Byte = 4KB, half the size of a core block memory
				// Max value to store a single copy of the matrix in a core: [44]

#define _MCHIP	(_MCORE * _NCORES_SIDE) // side of per chip sub-matrix 



// A Matrix [_M x _M] is split in [_MCHIP x _MCHIP] submatrices to fit a 16-cores epiphany chip.
// Inside the chip, the submatrix is then split in [_MCORE x _MCORE] sub-submatrices to fit the
// memory banks of the single cores (8Kb per memory bank).

// NOTE: All matrices are defined as arrays of 4-Bytes float numbers.


#define _Nchips 4                  // # of chips in operand matrix side
#define _Nside  4                  // # of cores in chip side
#define _Ncores (_Nside * _Nside)  // Num of cores = 16
#define _Score  32                 // side size of per-core sub-submatrix (max 32)
#define _Schip  (_Score * _Nside)  // side size of per-chip submatrix
#define _Smtx   (_Schip * _Nchips) // side size of operand matrix


typedef struct {
	float	A[_M * _M];	// A matrix 
	float	B[_M * _M];	// B matrix 
	float	C[_M * _M];	// C matrix 
	int	dev_ready;	// FLAG: [0,1] Epiphany initialized, ready to receive data in DRAM and start execution
	int 	hst_ready;	// FLAG: [0,1] Host has copied data to DRAM and is ready to receive result (GO signal)
	int 	completed;	// FLAG: [0,1] Epiphany completed execution, result is available in DRAM for host

} shared_buf_t;

