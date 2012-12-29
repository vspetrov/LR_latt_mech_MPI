#include "LR_lattice.h"
#include "LR_cell.h"
#include "parallel.h" 
#include "math.h"
#include <time.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>

int Size=600;

int DrawNum = 5;

#ifdef OS_LINUX
#define _lseek lseek
#endif

CellVariables *Cell;
double *I_ext;
double *V, *gk1, *type;
double *Gk, *Gsi;

int ln, rn, un, dn; //variables for defining coupling neighbours

const double dt=0.1;
//double D=0.1; //was 0.1
double Do = 0.;
double De = 0.5;
double Deo = 0.0;//0.05
double Doe = 0.;

void LattInit(const int &size)
{
	srand((ProcRank+1)*10);
	Cell = new CellVariables[size*size];
	I_ext = new double[size*size];
	V = new double[size*size];
	gk1 = new double[size*size];
	Gk = new double[size*size];
	Gsi = new double[size*size];
	type = new double[size*size];
	//define all lattice variables
	int i,j; //int counters
	for (i=0;i<size*size;i++)
	{
		V[i]=-75.;
		Cell[i].m=0.0143576;
		Cell[i].h=0.745286;
		Cell[i].j=0.76125;
		Cell[i].d=0.00934794;
		Cell[i].f=0.999772;
		Cell[i].X=0.018801218;
		Cell[i].Cai=0.00025;
        gk1[i] = 0.6047;// 0.04*double(rand())/double(RAND_MAX);
        type[i] = 1; // 1 - excitable, -1 - oscillatory
        Gk[i] = 0.705;//0.705
        Gsi[i] = 0.07;//0.07
        Cell[i].tension = 0;
        Cell[i].s_tension = 0;
	}
	//define external currents
	for (i=0;i<size*size;i++)
		I_ext[i]=0;//-2.7+(double)rand()/double(RAND_MAX)*(2.7-3.2);
		
	


    int ii,jj;
    int I,J;
    for (i=1; i<size-1; i++)
    {
        for (j=1; j<size-1; j++)
        {
            ii = i-1;
            jj = j-1;
            I = OwnCoords[0]*(size-2)+ii;
            J = OwnCoords[1]*(size-2)+jj;
            if (I < 10 && J < 10){
                V[i*size+j] = -30;
            }
        }
    }

}


void SolveEquations(double MaxTime, MPI_Comm &GridComm, int &GridSize, MPI_Datatype &column, MPI_Datatype &square)
{
	int i,j,k; //counting variables
	int time, MT;//time itarator
	int counter=1;
	MT=int(MaxTime/dt);
//	int wall_flag=0;
	int drawCounter=0;
	int step = int(floor(((MaxTime-1.)/(dt*(DrawNum-1)))));
	
    double *buf = new double[N*N];
	//integrating on the interval from time to time+dt/2
	Exchange(&GridComm,GridSize,V,N,&column);
	Exchange(&GridComm,GridSize,type,N,&column);

	for (i=1; i<N-1; i++)
		for (j=1; j<N-1; j++)
			V[i*N+j]+=dt/2.*Coupling(i,j);

	int GS=sqrt((double)(ProcNum-1));
    int percent_step = MT/100;
	

	for (time=0; time<MT; time++) 
	{
        if ((time+1)/percent_step*percent_step == (time+1) && (0 == ProcRank)){
            printf("%d%%..\r",(time+1)/percent_step);
            fflush(stdout);
        }
		//solve system with NO coupling on the interval from time to time+dt
		for (i=1; i<N-1; i++)
			for (j=1; j<N-1; j++)
				OdeSolve(i, j);

				
		Exchange(&GridComm,GridSize,V,N,&column);
		//save data to the file
		drawCounter++;

		
        if (time/step*step == time){
			MPI_Gather(&V[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
            for (i=0;i<N*N; i++)
            {
                buf[i] = Cell[i].tension;
            }
            MPI_Gather(&buf[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
        }


		//again calculate coupling
		for (i=1; i<N-1; i++)
			for (j=1; j<N-1; j++)
				V[i*N+j]+=dt*Coupling(i,j);

		
	}


	MPI_Gather(&V[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	
	//saving state - 9 gathers
	//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
	MPI_Gather(&V[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
	double *temp=new double[N*N];
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].m;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].h;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].j;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//ddddddddddddddddddddddddddddddddddddddddddddddddddd
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].d;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//fffffffffffffffffffffffffffffffffffffffffffffffffff
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].f;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].X;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//CaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCaiCai
	for (i=0;i<N*N; i++)
	{
		temp[i] = Cell[i].Cai;
	}
	MPI_Gather(&temp[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
	//gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1gk1
	MPI_Gather(&gk1[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);

    delete[] temp;
    delete[] buf;
}

void OdeSolve(int &ii, int &jj)
{
	double vd;
	int kstep=1;
	double delta_t=dt;
	for(int i=0; i<kstep; i++)
	{
		vd=VFunction(ii,jj);
		
		if(i==0) //decide on the time substep
		{
			kstep=Substeps(vd);
			delta_t=dt/(double)kstep;
		}

		V[ii*N+jj]+=delta_t*vd;
		Cell[ii*N+jj].m=mFunction(delta_t);
		Cell[ii*N+jj].h=hFunction(delta_t);
		Cell[ii*N+jj].j=jFunction(delta_t);
		Cell[ii*N+jj].d=dFunction(delta_t);
		Cell[ii*N+jj].f=fFunction(delta_t);
		Cell[ii*N+jj].X=XFunction(delta_t);
		Cell[ii*N+jj].Cai+=delta_t*CaiFunction();
        Cell[ii*N+jj].s_tension += delta_t * s_tensionFunction();
	}
}

inline int Substeps(double &vd)
{
	const int kmax=100;
	int k;

	const int k0=vd>0. ? 5 : 1;
 	k=k0+(int)fabs(vd);
	
	return k<kmax ? k : kmax;
}

double Coupling(int &ii, int &jj)
{
	double S = 0;

	


	if (type[ii*N+jj] > 0)
	{
		if (type[ii*N+jj-1] > 0)
			S += De*(V[ii*N+jj-1] - V[ii*N+jj]);
		else
			S += Doe*(V[ii*N+jj-1] - V[ii*N+jj]);

		if (type[ii*N+jj+1] > 0)
			S += De*(V[ii*N+jj+1] - V[ii*N+jj]);
		else
			S += Doe*(V[ii*N+jj+1] - V[ii*N+jj]);

		if (type[(ii+1)*N+jj] > 0)
			S += De*(V[(ii+1)*N+jj] - V[ii*N+jj]);
		else
			S += Doe*(V[(ii+1)*N+jj] - V[ii*N+jj]);

		if (type[(ii-1)*N+jj] > 0)
			S += De*(V[(ii-1)*N+jj] - V[ii*N+jj]);
		else
			S += Doe*(V[(ii-1)*N+jj] - V[ii*N+jj]);
	}
	else if (type[ii*N+jj] < 0)
	{
		if (type[ii*N+jj-1] < 0)
			S += Do*(V[ii*N+jj-1] - V[ii*N+jj]);
		else
			S += Deo*(V[ii*N+jj-1] - V[ii*N+jj]);

		if (type[ii*N+jj+1] < 0)
			S += Do*(V[ii*N+jj+1] - V[ii*N+jj]);
		else
			S += Deo*(V[ii*N+jj+1] - V[ii*N+jj]);

		if (type[(ii+1)*N+jj] < 0)
			S += Do*(V[(ii+1)*N+jj] - V[ii*N+jj]);
		else
			S += Deo*(V[(ii+1)*N+jj] - V[ii*N+jj]);

		if (type[(ii-1)*N+jj] < 0)
			S += Do*(V[(ii-1)*N+jj] - V[ii*N+jj]);
		else
			S += Deo*(V[(ii-1)*N+jj] - V[ii*N+jj]);
	}

	return S;
	//return D*(V[ii*N+jj+1]+V[ii*N+jj-1]+V[(ii+1)*N+jj]+V[(ii-1)*N+jj]-4.*V[ii*N+jj]);
}

void LattInit(int &size, char a[])
{
	Cell = new CellVariables[size*size];
	I_ext = new double[size*size];
	V = new double[size*size];
	gk1 = new double[size*size];	
    int rc;
#ifdef OS_WINDOWS
	int file_id=open(a,O_RDONLY | O_BINARY, S_IREAD );
#else
    int file_id=open(a,O_RDONLY , S_IREAD );
#endif
	if (ProcRank != ProcNum - 1)
	{
		//read V
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=0; i < N-2; i++)
		{
            rc = read(file_id,&(V[N+i*N+1]),(N-2)*sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read m __1 - skip
		_lseek(file_id,Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].m),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read h __2 - skip
		_lseek(file_id,2*Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].h),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read j __3 - skip
		_lseek(file_id,3*Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].j),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read d __4 - skip
		_lseek(file_id,4*Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].d),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read f __5 - skip
		_lseek(file_id,5*Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].f),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read X __6 - skip
		_lseek(file_id,6*Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].X),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}

        //rc = read Cai __7 - skip
		_lseek(file_id,7*Size*Size*sizeof(double),SEEK_SET);
		_lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
		for (int i=1; i < N-1; i++)
		{
			for (int j=1;j<N-1;j++)
                rc = read(file_id,&(Cell[i*N+j].Cai),sizeof(double));
			_lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
		}
	}
	close(file_id);
	for (int i=0;i<size*size;i++)
	{
		I_ext[i]=0.;
	

	}

//	printf("Process %i: ready..\n", ProcRank);
}

void Latt_Init1(int &size, char a[])
{
	Cell = new CellVariables[size*size];
	I_ext = new double[size*size];
	V = new double[size*size];
	gk1 = new double[size*size];	
	Gk = new double[size*size];	
	Gsi = new double[size*size];	
	type = new double[size*size];	
	int i,j, n = size-2;
    int rc;
#ifdef OS_WINDOWS
	int file_id=open(a,O_RDONLY | O_BINARY, S_IREAD );
#else
    int file_id=open(a,O_RDONLY, S_IREAD );
#endif
	if (ProcRank != ProcNum - 1)
	{
		_lseek(file_id,n*n*8*ProcRank,SEEK_CUR);
		
		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(V[i*size+j]),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].m),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].h),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].j),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].d),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].f),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].X),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(Cell[i*size+j].Cai),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);

		for (i=1; i<size-1; i++)
			for (j=1; j<size-1; j++)
                rc = read(file_id,&(gk1[i*size+j]),8);
		_lseek(file_id,(ProcNum-2)*n*n*8,SEEK_CUR);


	}
		close(file_id);
	for (int i=0;i<size*size;i++)
	{
		I_ext[i]=0.0;
	
		Gk[i] = 0.282;
		Gsi[i] = 0.09;
		type[i] = -1.;
	}
}
