/*
 * frank_energy.cpp
 * created by: Eric Koudhi Omori <koudhi@gmail.com>
 * Creation date: 2019-06-24
 * Last Modification: 2019-06-24
 * version: 0.0
 */


#include <iostream>
using namespace std;
#include <math.h>
#include "../include/dir.h"
#include "../include/io.h"
#include "../include/energy_eval.h"

int main(int argc, char **argv)
{
	if (argc < 6)
    {
        cout<<"Entre com: "<<argv[0]<<"input_file output_file Nx Ny Nz"<<endl;
        exit (1);
    }
    int Nx=atoi(argv[3]), Ny=atoi(argv[4]) ,Nz=atoi(argv[5]), N=Nx*Ny*Nz;
    if (N<=0){ cout<<"ERRO!!\n Tamanho de rede inválido"<<endl; exit (1);}
    vector *n, lat_size;
    double *elastic_energy;
    elastic_energy=(double*)calloc(3*N,sizeof (double));
    n=(vector*)calloc(N, sizeof (vector));
    if (n==NULL)
    {
        cout<<"erro ao alocar o vetor n"<<endl;
        exit (2);
    }
    read_file(argv[1],n,lat_size,Nx,Ny,Nz);
    elastic_energy_eval(Nx, Ny, Nz, n, elastic_energy, lat_size);
    print_file(argv[2],n,elastic_energy,lat_size,Nx,Ny,Nz);
    energy_integrator(elastic_energy, Nx, Ny, Nz);
    
    
    
	return 0;
}

