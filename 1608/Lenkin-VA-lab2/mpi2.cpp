#include <iostream>
#include <mpi.h>
using namespace std;


int main(int argc, char* argv[])
{
	
	int size = 4, rank = -1, start = 0, end = 0, task_size = 0, count = 0,q=0;
	double **Matrix;
	double *Vector;
	double start_time = 0, end_time = 0, max_el = 0, final = 0;
	int reord = 1;

	//size = atoi(argv[1]);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &count);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	task_size = size / count;
	Matrix = new double*[task_size];
	for (int i = 0; i<task_size; i++)
		Matrix[i] = new double[size + 1];
	Vector = new double[size + 1];

	for (int i = 0; i<task_size; i++)
	{
		for (int j = 0; j<size; j++)
		{
			if (rank + count * i == j)
			{
				Matrix[i][j] = 2;
			}
			else
			{
				Matrix[i][j] = 1;
			}
				cout << Matrix[i][j] << " ";
			
		}
		Matrix[i][size] = 1 * size + 1;
		cout <<" = " << Matrix[i][size] << "  " ;
	

	}
	start_time = MPI_Wtime();
	//прямой ход
	for (int k = 0; k<task_size; k++)
	{
		for (int p = 0; p<count; p++)
		{
			if (rank == p)
			{
				max_el = 1.0 / Matrix[k][count*k + p];
				for (int j = size; j >= count * k + p; j--)
				{
					Matrix[k][j] = Matrix[k][j] * max_el;
				}
				for (int j = 0; j <= size; j++)
				{
					Vector[j] = Matrix[k][j];
				}
				MPI_Bcast(Vector, size + 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
				for (int i = k + 1; i<task_size; i++)
				{
					for (int j = size; j >= count * k + p; j--)
					{
						Matrix[i][j] = Matrix[i][j] - Matrix[i][count*k + p] * Matrix[k][j];
					}
				}
			}
			else if (rank<p)
			{
				MPI_Bcast(Vector, size + 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
				for (int i = k + 1; i<task_size; i++)
					for (int j = size; j >= count * k + p; j--)
						Matrix[i][j] = Matrix[i][j] - Matrix[i][count*k + p] * Vector[j];

			}
			else if (rank>p)
			{
				MPI_Bcast(Vector, size + 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
				for (int i = k; i<task_size; i++)
					for (int j = size; j >= count * k + p; j--)
						Matrix[i][j] = Matrix[i][j] - Matrix[i][count*k + p] * Vector[j];
			}
		}
		
	}

	//обратный ход

	for (int k = task_size - 1; k >= 0; k--)
	{
		for (int p = count - 1; p >= 0; p--)
		{
			if (rank == p)
			{
				final = Matrix[k][size];
				MPI_Bcast(&final, 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
				for (int i = k - 1; i >= 0; i--)
					Matrix[i][size] -= Matrix[k][size] * Matrix[i][count*k + p];
			}

			else if (rank < p)
			{
				MPI_Bcast(&final, 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
				for (int i = k; i >= 0; i--)
				{
					Matrix[i][size] -= final * Matrix[i][count*k + p];
				}
			}
			else if (rank > p)
			{
				MPI_Bcast(&final, 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
				for (int i = k - 1; i >= 0; i--)
				{
					Matrix[i][size] -= final * Matrix[i][count*k + p];
					
				}
			}

		
		}
		

	}
	end_time = MPI_Wtime();
	cout << "time on " << rank << "proc is " << end_time - start_time << endl;
	for (int i=0;i<task_size;i++)
	{       
            cout << Matrix[i][size] << " ";
	}
	MPI_Finalize();
	

	return 0;
}