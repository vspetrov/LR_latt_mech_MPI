#ifndef _PARALLEL_H
#define _PARALLEL_H

#include "mpi.h"

#ifdef OS_WINDOWS
#include "io.h"
#include "sys\stat.h"
#elif defined(OS_LINUX)
#include <sys/stat.h>
#include <unistd.h>
#endif

extern int ProcNum;
extern int ProcRank;
extern int N;
extern double *null_p;
extern int OwnCoords[2];
void CreateGrid(int GridSize,MPI_Comm *Comm, MPI_Comm* old_comm);
void Exchange(MPI_Comm *Comm, int GridSize, double* V,int N, MPI_Datatype *column);
extern void save(short* V,int N,int fd);
extern void save_double(double* Vv,int Nn,int fd);

#endif
