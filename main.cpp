#include "sudoku_header.h"
#define comment 0

int main()
{
	int *solved;

#if comment
	int v[N*N] = { 
		
		0, 3, 0,  4, 7, 0,  8, 0, 0,
		0, 0, 1,  8, 0, 0,  0, 0, 0,
		0, 0, 0,  0, 0, 0,  0, 0, 0,

		0, 6, 0,  0, 0, 0,  7, 0, 0,
		0, 8, 5,  0, 0, 2,  0, 9, 6,
		0, 0, 0,  0, 8, 0,  0, 0, 0,

		0, 4, 0,  0, 9, 6,  0, 0, 1,
		3, 5, 0,  0, 0, 0,  0, 2, 4,
		1, 0, 0,  3, 2, 0,  0, 0, 0
		
		};
#else
	int v[N*N];
	getFunc(v);
#endif // comment

	printFunc(v);

	solved = solveFunc(v);

	printFunc(solved);

	return 0;
}
