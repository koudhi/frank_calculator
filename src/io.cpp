/*
 * io.cpp
 * created by: Eric Koudhi Omori <koudhi@gmail.com>
 * Creation date: 2019-06-24
 * Last Modification: 2019-06-24
 * version: 0.0
 */


#include <iostream>
using namespace std;
#include "../include/dir.h"
#include "../include/io.h"
int get_size(FILE *input, vector &lat_size, int Nx, int Ny, int Nz)
{
    char escaper[400];
    double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z;
    int tester;
    while (fscanf(input, "%lf,%lf,%lf", &x, &y, &z)!= EOF)
    {
        {
            if ( x> xmax) xmax = x;
            if ( y> ymax) ymax = y;
            if ( z> zmax) zmax = z;
            if ( x< xmin) xmin = x;
            if ( y< ymin) ymin = y;
            if ( z< zmin) zmin = z;
        }
        fgets(escaper, 400, input);
    }
    lat_size.x=xmax-xmin;
    lat_size.y=ymax-ymin;
    lat_size.z=zmax-zmin;
    return 0;
}
int read_file(char* fname, vector *n, vector &lat_size, int Nx, int Ny, int Nz)
{
    char escaper[400];
    FILE *input=fopen(fname,"r");
    int tester, i, j, k;
    double nx, ny, nz, ii, jj, kk;
    if (input ==0) { cout<<"Não foi possível abrir o arquivo "<< fname << endl; exit (2);}
    get_size(input,lat_size, Nx, Ny, Nz);
    rewind(input);
    while (tester!=EOF)
    {
        tester=fscanf(input, "%lf,%lf,%lf,%lf,%lf,%lf", &ii, &jj, &kk, &nx, &ny, &nz);
        if(tester==6)
        {
            i=(int)((Nx-1)*(ii+0.5*lat_size.x)/lat_size.x+0.1);
            j=(int)((Ny-1)*(jj+0.5*lat_size.y)/lat_size.y+0.1);
            k=(int)((Nz-1)*(kk+0.5*lat_size.z)/lat_size.z+0.1);
            n[i+Nx*(j+k*Ny)].x=nx;
            n[i+Nx*(j+k*Ny)].y=ny;
            n[i+Nx*(j+k*Ny)].z=nz;
        }
        fgets(escaper,400,input);
    }
     
    
    return 0;
}
int print_file(char *fname, vector *n, double *K, vector lat_size, int Nx, int Ny, int Nz)
{
   FILE *output=fopen(fname, "w");
   double norm_x=Nx>1?lat_size.x/(Nx-1):1;
   double norm_y=Ny>1?lat_size.y/(Ny-1):1;
   double norm_z=Nz>1?lat_size.z/(Nz-1):1;
   fprintf(output,"i,j,k,K1,K2,K3\n");
   for (int i = 0; i < Nx; i++)
   {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int nn=i+Nx*(j+Ny*k);
                fprintf(output,"%g,%g,%g,%g,%g,%g\n",(i-0.5*(Nx-1))*norm_x,(j-0.5*(Ny-1))*norm_y,(k-0.5*(Nz-1))*norm_z,K[3*nn+0],K[3*nn+1],K[3*nn+2]);
               // fprintf(output,"%g,%g,%g,%lf,%lf,%lf\n",(i-0.5*(Nx-1))*norm_x,(j-0.5*(Ny-1))*norm_y,(k-0.5*(Nz-1))*norm_z,n[nn].x,n[nn].y,n[nn].z);
            }
        }
   }
    return 0;
}
