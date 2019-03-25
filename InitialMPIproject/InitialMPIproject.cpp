
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define HEAVY 100000
#define MASTER_ID 0
#define N 20

double heavy(int x, int y) {
	int i, loop = 1;
	double sum = 0;

	// Super heavy tasks
	if (x < 5 || y < 5)
		loop = 10;
	// Heavy calculations
	for (i = 0; i < loop*HEAVY; i++)
		sum += sin(exp(sin((double)i / HEAVY)));

	return sum;
}

int find_proc_with_min_weight(int *weight_per_process, int size)
{
	int i, min_index = 0;
	for (i = 1; i < size; i++)
	{
		if (weight_per_process[i] < weight_per_process[min_index])
			min_index = i;
	}
	return min_index;
}

int main(int argc, char *argv[])
{

	int numprocs, myid;
	int i = 0;


	int *my_x, *my_y;
	int my_num_of_x_y;
	double answer = 0, t1, t2;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


	if (myid == MASTER_ID)
	{
		int id_of_min_weight, weight_to_add, actual_num_of_tasks_to_send, x_y_col;
		int heavy_task_weight = 10, easy_task_weight = 1;
		int *weight_per_process;
		int x, y;
		int *actual_number_of_x_y_per_proc;
		int **x_per_proc, **y_per_proc;

		weight_per_process = (int*)calloc(numprocs, sizeof(int));
		actual_number_of_x_y_per_proc = (int*)calloc(numprocs, sizeof(int));
		x_per_proc = (int **)malloc(numprocs * sizeof(int *));
		y_per_proc = (int **)malloc(numprocs * sizeof(int *));

		t1 = MPI_Wtime();

		for (i = 0; i < numprocs; i++)
		{
			x_per_proc[i] = (int *)malloc(N*N * sizeof(int));
			y_per_proc[i] = (int *)malloc(N*N * sizeof(int));
		}


		for (x = 0; x < N; x++)
		{
			for (y = 0; y < N; y++)
			{
				id_of_min_weight = find_proc_with_min_weight(weight_per_process, numprocs);
				x_y_col = actual_number_of_x_y_per_proc[id_of_min_weight];
				x_per_proc[id_of_min_weight][x_y_col] = x;
				y_per_proc[id_of_min_weight][x_y_col] = y;
				actual_number_of_x_y_per_proc[id_of_min_weight] += 1;

				if (x < 5 || y < 5)
					weight_to_add = heavy_task_weight;
				else
					weight_to_add = easy_task_weight;

				weight_per_process[id_of_min_weight] += weight_to_add;
			}
		}

		for (i = 1; i < numprocs; i++)
		{
			actual_num_of_tasks_to_send = actual_number_of_x_y_per_proc[i];
			MPI_Send(&(actual_num_of_tasks_to_send), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(x_per_proc[i], actual_num_of_tasks_to_send, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(y_per_proc[i], actual_num_of_tasks_to_send, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		//for (i = 0; i < numprocs; i++)
		//	printf("\nfor id = %d, weight of tasks =  %d, actual_number_of_x_y_per_proc = %d", i, weight_per_process[i], actual_number_of_x_y_per_proc[i]);

		my_num_of_x_y = actual_number_of_x_y_per_proc[MASTER_ID];
		my_x = (int*)malloc(my_num_of_x_y * sizeof(int));
		my_x = x_per_proc[MASTER_ID];
		my_y = (int*)malloc(my_num_of_x_y * sizeof(int));
		my_y = y_per_proc[MASTER_ID];
	}
	else
	{
		MPI_Recv(&my_num_of_x_y, 1, MPI_INT, MASTER_ID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		my_x = (int*)malloc(my_num_of_x_y * sizeof(int));
		MPI_Recv(my_x, my_num_of_x_y, MPI_INT, MASTER_ID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		my_y = (int*)malloc(my_num_of_x_y * sizeof(int));
		MPI_Recv(my_y, my_num_of_x_y, MPI_INT, MASTER_ID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	}

	//printf("\n\nmyid = %d", myid);
	//for (i = 0; i < my_num_of_x_y; i++)
	//{
	//	printf("\nx = %d, y = %d", my_x[i], my_y[i]);
	//}


	for (i = 0; i < my_num_of_x_y; i++)
	{
		answer += heavy(my_x[i], my_y[i]);
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("\nmy id = %d, my_num_of_x_y = %d, my answer = %e", myid, my_num_of_x_y, answer);



	if (myid == MASTER_ID)
	{
		int temp_ans;

		for (i = 1; i < numprocs; i++)
		{
			MPI_Recv(&temp_ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			answer += temp_ans;
		}
		t2 = MPI_Wtime();
		printf("\n\nanswer = %e\ntime = %1.9f", answer, t2 - t1);
	}
	else
	{
		MPI_Send(&answer, 1, MPI_DOUBLE, MASTER_ID, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}

