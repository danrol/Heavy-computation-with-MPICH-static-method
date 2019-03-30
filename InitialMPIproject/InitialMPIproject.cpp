#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define HEAVY 100000

#define X_ARR_TAG 0
#define Y_ARR_TAG 1
#define GET_ANSWERS_TAG 3
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

int main(int argc, char *argv[])
{

	int numprocs, myid;
	int N = 20, i, x, y;
	int num_of_tasks = N * N, tasks_counter = 0, *x_arr, *y_arr;
	int tasks_per_proc, remainder;
	int start_from_task, end_task;
	double answer = 0, t1, t2;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	tasks_per_proc = num_of_tasks / numprocs;

	x_arr = (int*)malloc(N*N * sizeof(int));
	y_arr = (int*)malloc(N*N * sizeof(int));

	t1 = MPI_Wtime();
	remainder = num_of_tasks % numprocs;

	if (myid == 0)
	{
		for (x = 0; x < N; x++)
		{
			for (y = 0; y < N; y++)
			{
				x_arr[tasks_counter] = x; //put all x tasks in array
				y_arr[tasks_counter] = y; //put all y tasks in array
				tasks_counter++;
			}
		}
		for (i = 1; i < numprocs; i++)
		{
			MPI_Send(x_arr, num_of_tasks, MPI_INT, i, X_ARR_TAG, MPI_COMM_WORLD);
			MPI_Send(y_arr, num_of_tasks, MPI_INT, i, Y_ARR_TAG, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(x_arr, num_of_tasks, MPI_INT, 0, X_ARR_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(y_arr, num_of_tasks, MPI_INT, 0, Y_ARR_TAG, MPI_COMM_WORLD, &status);
	}


	//Each process calculates his tasks - if more processes than tasks, each calculates one
	if (numprocs >= num_of_tasks)
	{
		if (myid <= num_of_tasks)
		{
			answer += heavy(x_arr[myid], y_arr[myid]);
			remainder = 0;
		}
	}
	else
	{
		start_from_task = myid * tasks_per_proc;
		end_task = myid * tasks_per_proc + tasks_per_proc;
		for (i = start_from_task; i < end_task; i++)
		{
			answer += heavy(x_arr[i], y_arr[i]);
		}
	}

	if (myid == 0)
	{
		double temp_ans;
		//calculate answers for all remaining tasks
		for (i = num_of_tasks - remainder; i < num_of_tasks; i++)
		{
			answer += heavy(x_arr[i], y_arr[i]);
		}

		//get answers from slaves
		for (i = 1; i < numprocs; i++)
		{
			MPI_Recv(&temp_ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, GET_ANSWERS_TAG, MPI_COMM_WORLD, &status);
			answer += temp_ans;
		}
		t2 = MPI_Wtime();
		printf("answer = %e\ntime = %1.9f", answer, t2 - t1);

	}
	else
	{
		// send answers to master
		MPI_Send(&answer, 1, MPI_DOUBLE, 0, GET_ANSWERS_TAG, MPI_COMM_WORLD);
	}

	free(x_arr);
	free(y_arr);
	MPI_Finalize();
	return 0;
}