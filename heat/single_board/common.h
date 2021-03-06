

#define _NROWS		4		// # of rows in core group
#define _NCOLS		4		// # of cols in core group 
#define _NCORES		(_NROWS*_NCOLS) // # of cores in use



#define _MSIDE		1024 // 1024 max


#define _M		(_MSIDE*_MSIDE) // # of equations/unknowns in tridiagonal system
#define _MCORE		1024		// # of equations/unknowns handled by single core, MAX:1024
#define _MSHR 		_M		// 1048576 // maximum size of Mailbox arrays


#define _P		(_M/_MCORE)	// number of parallel processes required

#define _GAMMA		1
#define _ITERATIONS	10


#define _size_MCORE	sizeof(float)*_MCORE
#define _size_M		sizeof(float)*_M
#define _size_MSHR	sizeof(float)*_MSHR
#define _size_P 	sizeof(float)*((2*_P)-2)
#define _size_F 	sizeof(float)


typedef struct {

	float	d[_MSHR];		// right hand temperature vector

	float	tyy[_MSHR];	
	float	tza[_MSHR];	
	float	tzb[_MSHR];	

	int	go;		
	int	completed;
	int	debug;

} shared_buf_t;




