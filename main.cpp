#include <time.h>
#include "math.h"
#include <stdio.h>
#include "LR_cell.h"
#include "LR_lattice.h"
#include "parallel.h"
#include <assert.h>
#include "stdlib.h"
#include "fcntl.h"

#include <stdio.h>

static void convertRst(short *src, short *dst, int Size, int GridSize);
short *convert_buf = NULL;

int main(int argc, char *argv[])
{
	double t1,t2,duration;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);//Number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);//Rank of process
	t1 = MPI_Wtime();
	MPI_Group WorldGroup, CalculatorGroup;
	MPI_Comm Calculators;
	int ranks[1];
	ranks[0] = ProcNum-1;
	MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
	MPI_Group_excl(WorldGroup, 1, ranks, &CalculatorGroup);
	MPI_Comm_create(MPI_COMM_WORLD,CalculatorGroup,&Calculators);

	
	int GridSize=sqrt((double)(ProcNum-1));//size of virtual topology(Grid)
	N=Size/GridSize+2;//N-size of subsystem; +2 - for boundary condition
    assert(GridSize*GridSize+1 == ProcNum);
	

	if (ProcRank!=ProcNum-1)
	{
		
		MPI_Datatype column;
		MPI_Type_vector(N,1,N,MPI_DOUBLE,&column);
		MPI_Type_commit(&column);

		//creating datatype for square
		MPI_Datatype square;
		MPI_Type_vector(N-2,N-2,N,MPI_DOUBLE,&square);
		MPI_Type_commit(&square);

		MPI_Comm GridComm;
		CreateGrid(GridSize, &GridComm, &Calculators);
		//Latt_Init1(N,"200last_state_gk_d0p1.bin");
		LattInit(N);
        SolveEquations(250, GridComm, GridSize, column, square);
		MPI_Type_free(&square);
		MPI_Type_free(&column);
		MPI_Group_free(&CalculatorGroup);
		MPI_Comm_free(&Calculators);
		delete[] V;
		delete[] Cell;
		delete[] I_ext;
	}
	else
	{
        convert_buf = new short[(N-2)*(N-2)*ProcNum];
		double* V_all = new double[(N-2)*(N-2)*ProcNum];
		double* V_temp = new double[(N-2)*(N-2)];
		short* V_save = new short[(N-2)*(N-2)*ProcNum];
		int ii;
		//printf("Total number of frames: %i\n", DrawNum);
#ifdef OS_WINDOWS
		int fd = open("200snapshots_gk_d0p1.bin",O_RDWR|O_CREAT | O_BINARY,S_IREAD|S_IWRITE);
#else
        int fd = open("snapshots.bin",O_RDWR|O_CREAT ,S_IREAD|S_IWRITE);
#endif
        for (int i=0;i<DrawNum*2;i++)
		{
			
			MPI_Gather(V_temp,(N-2)*(N-2),MPI_DOUBLE,V_all,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
			for (ii=0;ii<(N-2)*(N-2)*(ProcNum-1);ii++)
			{
				V_save[ii]=short(V_all[ii]*250.);
			}
            convertRst(V_save,convert_buf,Size,GridSize);
            save(convert_buf,Size*Size,fd);
        }
		close(fd);

		MPI_Gather(V_temp,(N-2)*(N-2),MPI_DOUBLE,V_all,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
		for (ii=0;ii<(N-2)*(N-2)*(ProcNum-1);ii++)
		{
			V_save[ii]=short(V_all[ii]*250.);
		}
#ifdef OS_WINDOWS
		fd = open("200last_V_gk_d0p1.bin",O_RDWR|O_CREAT | O_BINARY,S_IREAD|S_IWRITE);
#else
        fd = open("lastV.bin",O_RDWR|O_CREAT,S_IREAD|S_IWRITE);
#endif
		save(V_save,(N-2)*(N-2)*(ProcNum-1),fd);
		close(fd);

#ifdef OS_WINDOWS
		fd = open("200last_state_gk_d0p1.bin",O_RDWR|O_CREAT | O_BINARY,S_IREAD|S_IWRITE);
#else
        fd = open("state.bin",O_RDWR|O_CREAT ,S_IREAD|S_IWRITE);
#endif
		for (int i=0; i<9; i++)
		{
		MPI_Gather(V_temp,(N-2)*(N-2),MPI_DOUBLE,V_all,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
		save_double(V_all,(N-2)*(N-2)*(ProcNum-1),fd);
		}
		close(fd);



		delete[] V_all;
		delete[] V_temp;
		delete[] V_save;
        delete[] convert_buf;
	}



	
	
	t2 = MPI_Wtime();
	if (ProcRank==0)
	{
		printf("Experiment duration: %f seconds\n",t2-t1);
	}
	MPI_Finalize();
	
	return 0;
}

static void convertRst(short *src, short *dst, int Size, int GridSize){
    int N = Size/GridSize;
    int counter = 0;
    for (int i=0; i<GridSize; i++){
        for (int j=0; j<GridSize; j++){
            for (int ii=0; ii<N; ii++){
                for (int jj=0; jj<N; jj++){
                    dst[(i*N+ii)*Size+j*N+jj]=src[counter];
                    counter++;
                }
            }
        }
    }
}
