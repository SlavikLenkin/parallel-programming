#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int main(int argc, char *argv[]) {
	int size = 1000;
	char *string = (char*)malloc(size * sizeof(char));
	int sum = 0;
	int ProcNum, ProcRank;
	double times;
	int k;
	int tmp = 0;
	int TotalSum;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	//times = MPI_Wtime();
	if (ProcRank == 0) {
		gets(string);
		int k = size / ProcNum;
		MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
		for (int i = 1; i < ProcNum; i++)
			MPI_Send(string + (k * i),k, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		for (int i = 0; i < k; i++)
			if (string[i] == '.' || string[i] == '?' || string[i] == '!')
				tmp += 1;
		
		int TotalSum;
		MPI_Reduce(&tmp, &TotalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		printf("Total Sum = %d\n", TotalSum);
		
		
	}
	else {
		MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
		MPI_Status status;
	    string = (char*)malloc(k * sizeof(char));
		
		MPI_Recv(string, k, MPI_DOUBLE, MPI_ANY_TAG, 0, MPI_COMM_WORLD, &status);
		int tmp=0;
		for (int i = 0; i < k; i++)
			if (string[i] == '.' || string[i] == '?' || string[i] == '!')
				tmp += 1;
		MPI_Reduce(&tmp, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


	}
	

	
	free(string);

	//printf("     %d  - %20.10f\n", ProcRank, MPI_Wtime() - times);
	MPI_Finalize();
	
    return 0;
}