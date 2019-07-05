/*
 * energy_eval.cpp
 * created by: Eric Koudhi Omori <koudhi@gmail.com>
 * Creation date: 2019-06-24
 * Last Modification: 2019-06-24
 * version: 0.0
 */
#define n2(pos,x) n[pos].x*n[pos]
#include <iostream>
#include <math.h>
#include "../include/dir.h"
#include "../include/energy_eval.h"

int elastic_energy_eval (int Nx, int Ny, int Nz, vector *n, double *energy, vector lat_size)
{
    double k11=16.7e-12;
    double k22= 7.0e-12;
    double k33=18.1e-12;
    double q0=-8*M_PI/lat_size.z;
    
    int i, ip, im, j, jp, jm, k, kp, km, nn, dir;
    vector np, nm;
    /*ddi: 1-ni/x, 2-ni/y, 3-ni/z */
    double dnx[3],dny[3],dnz[3];
    long double scalep,scalem;
    double splay1, splay, bend1[3], bend, twist1, twist;
    
    //Bulk evaluation
    
    for (k = 1; k < Nz-1 ; k++)
    {
        kp=k+1; km=k-1;    
        for (j = 0; j < Ny; j++)
        {
            jp=(j+1)%Ny;
            jm=(j-1+Ny)%Ny;
            for (i = 0; i < Nx; i++)
            {
                nn=i+Nx*(j+Ny*k);
                ip=(i+1)%Nx;
                im=(i-1+Nx)%Nx;
                scalep=n[nn].x*n[ip+Nx*(j+Ny*k)].x + n[nn].y*n[ip+Nx*(j+Ny*k)].y + n[nn].z*n[ip+Nx*(j+Ny*k)].z;
                scalem=n[nn].x*n[im+Nx*(j+Ny*k)].x + n[nn].y*n[im+Nx*(j+Ny*k)].y + n[nn].z*n[im+Nx*(j+Ny*k)].z;
                dnx[0]=scalep*n[ip+Nx*(j+Ny*k)].x-scalem*n[im+Nx*(j+Ny*k)].x;
                dny[0]=scalep*n[ip+Nx*(j+Ny*k)].y-scalem*n[im+Nx*(j+Ny*k)].y;
                dnz[0]=scalep*n[ip+Nx*(j+Ny*k)].z-scalem*n[im+Nx*(j+Ny*k)].z;
                scalep=n[nn].x*n[i+Nx*(jp+Ny*k)].x + n[nn].y*n[i+Nx*(jp+Ny*k)].y + n[nn].z*n[i+Nx*(jp+Ny*k)].z;
                scalem=n[nn].x*n[i+Nx*(jm+Ny*k)].x + n[nn].y*n[i+Nx*(jm+Ny*k)].y + n[nn].z*n[i+Nx*(jm+Ny*k)].z;
                dnx[1]=scalep*n[i+Nx*(jp+Ny*k)].x-scalem*n[i+Nx*(jm+Ny*k)].x;
                dny[1]=scalep*n[i+Nx*(jp+Ny*k)].y-scalem*n[i+Nx*(jm+Ny*k)].y;
                dnz[1]=scalep*n[i+Nx*(jp+Ny*k)].z-scalem*n[i+Nx*(jm+Ny*k)].z;
                scalep=n[nn].x*n[i+Nx*(j+Ny*kp)].x + n[nn].y*n[i+Nx*(j+Ny*kp)].y + n[nn].z*n[i+Nx*(j+Ny*kp)].z;
                scalem=n[nn].x*n[i+Nx*(j+Ny*km)].x + n[nn].y*n[i+Nx*(j+Ny*km)].y + n[nn].z*n[i+Nx*(j+Ny*km)].z;
                dnx[2]=scalep*n[i+Nx*(j+Ny*kp)].x-scalem*n[i+Nx*(j+Ny*km)].x;
                dny[2]=scalep*n[i+Nx*(j+Ny*kp)].y-scalem*n[i+Nx*(j+Ny*km)].y;
                dnz[2]=scalep*n[i+Nx*(j+Ny*kp)].z-scalem*n[i+Nx*(j+Ny*km)].z;
                dnx[0]*=0.5*Nx/lat_size.x;
                dnx[1]*=0.5*Ny/lat_size.y;
                dnx[2]*=0.5*Nz/lat_size.z;
                dny[0]*=0.5*Nx/lat_size.x;
                dny[1]*=0.5*Ny/lat_size.y;
                dny[2]*=0.5*Nz/lat_size.z;
                dnz[0]*=0.5*Nx/lat_size.x;
                dnz[1]*=0.5*Ny/lat_size.y;
                dnz[2]*=0.5*Nz/lat_size.z;
                splay1=dnx[0]+dny[1]+dnz[2];
                
                twist1=n[nn].x*(dnz[1]-dny[2])
                      +n[nn].y*(dnx[2]-dnz[0])
                      +n[nn].z*(dny[0]-dnx[1])-q0;
                bend1[0]=n[nn].y*(dny[0]-dnx[1])-n[nn].z*(dnx[2]-dnz[0]);
                bend1[1]=n[nn].z*(dnz[1]-dny[2])-n[nn].x*(dny[0]-dnx[1]);
                bend1[2]=n[nn].x*(dnx[2]-dnz[0])-n[nn].y*(dnz[1]-dny[2]);      
                bend=bend1[0]*bend1[0]+bend1[1]*bend1[1]+bend1[2]*bend1[2];
                energy[3*(nn)+0]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k11*splay1*splay1;
                energy[3*(nn)+1]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k22*twist1*twist1;
                energy[3*(nn)+2]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k33*bend;
            }
        }
    }
    
    //bottom evaluation k=0
    km=k=0;kp=1; 
        for (j = 0; j < Ny; j++)
        {
            jp=(j+1)%Ny;
            jm=(j-1+Ny)%Ny;
            for (i = 0; i < Nx; i++)
            {
                nn=i+Nx*(j+Ny*k);
                ip=(i+1)%Nx;
                im=(i-1+Nx)%Nx;
                scalep=n[nn].x*n[ip+Nx*(j+Ny*k)].x + n[nn].y*n[ip+Nx*(j+Ny*k)].y + n[nn].z*n[ip+Nx*(j+Ny*k)].z;
                scalem=n[nn].x*n[im+Nx*(j+Ny*k)].x + n[nn].y*n[im+Nx*(j+Ny*k)].y + n[nn].z*n[im+Nx*(j+Ny*k)].z;
                dnx[0]=scalep*n[ip+Nx*(j+Ny*k)].x-scalem*n[im+Nx*(j+Ny*k)].x;
                dny[0]=scalep*n[ip+Nx*(j+Ny*k)].y-scalem*n[im+Nx*(j+Ny*k)].y;
                dnz[0]=scalep*n[ip+Nx*(j+Ny*k)].z-scalem*n[im+Nx*(j+Ny*k)].z;
                scalep=n[nn].x*n[i+Nx*(jp+Ny*k)].x + n[nn].y*n[i+Nx*(jp+Ny*k)].y + n[nn].z*n[i+Nx*(jp+Ny*k)].z;
                scalem=n[nn].x*n[i+Nx*(jm+Ny*k)].x + n[nn].y*n[i+Nx*(jm+Ny*k)].y + n[nn].z*n[i+Nx*(jm+Ny*k)].z;
                dnx[1]=scalep*n[i+Nx*(jp+Ny*k)].x-scalem*n[i+Nx*(jm+Ny*k)].x;
                dny[1]=scalep*n[i+Nx*(jp+Ny*k)].y-scalem*n[i+Nx*(jm+Ny*k)].y;
                dnz[1]=scalep*n[i+Nx*(jp+Ny*k)].z-scalem*n[i+Nx*(jm+Ny*k)].z;
                scalep=n[nn].x*n[i+Nx*(j+Ny*kp)].x + n[nn].y*n[i+Nx*(j+Ny*kp)].y + n[nn].z*n[i+Nx*(j+Ny*kp)].z;
                scalem=n[nn].x*n[i+Nx*(j+Ny*km)].x + n[nn].y*n[i+Nx*(j+Ny*km)].y + n[nn].z*n[i+Nx*(j+Ny*km)].z;
                dnx[2]=scalep*n[i+Nx*(j+Ny*kp)].x-scalem*n[i+Nx*(j+Ny*km)].x;
                dny[2]=scalep*n[i+Nx*(j+Ny*kp)].y-scalem*n[i+Nx*(j+Ny*km)].y;
                dnz[2]=scalep*n[i+Nx*(j+Ny*kp)].z-scalem*n[i+Nx*(j+Ny*km)].z;
                dnx[0]*=0.5*Nx/lat_size.x;
                dnx[1]*=0.5*Ny/lat_size.y;
                dnx[2]*=	Nz/lat_size.z;
                dny[0]*=0.5*Nx/lat_size.x;
                dny[1]*=0.5*Ny/lat_size.y;
                dny[2]*=	Nz/lat_size.z;
                dnz[0]*=0.5*Nx/lat_size.x;
                dnz[1]*=0.5*Ny/lat_size.y;
                dnz[2]*=	Nz/lat_size.z;
                splay1=dnx[0]+dny[1]+dnz[2];
                
                twist1=n[nn].x*(dnz[1]-dny[2])
                      +n[nn].y*(dnx[2]-dnz[0])
                      +n[nn].z*(dny[0]-dnx[1])-q0;
                bend1[0]=n[nn].y*(dny[0]-dnx[1])-n[nn].z*(dnx[2]-dnz[0]);
                bend1[1]=n[nn].z*(dnz[1]-dny[2])-n[nn].x*(dny[0]-dnx[1]);
                bend1[2]=n[nn].x*(dnx[2]-dnz[0])-n[nn].y*(dnz[1]-dny[2]);      
                bend=bend1[0]*bend1[0]+bend1[1]*bend1[1]+bend1[2]*bend1[2];
                energy[3*(nn)+0]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k11*splay1*splay1;
                energy[3*(nn)+1]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k22*twist1*twist1;
                energy[3*(nn)+2]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k33*bend;
            }
        }
    kp=k=Nz-1;km=Nz-2;
        for (j = 0; j < Ny; j++)
        {
            jp=(j+1)%Ny;
            jm=(j-1+Ny)%Ny;
            for (i = 0; i < Nx; i++)
            {
                nn=i+Nx*(j+Ny*k);
                ip=(i+1)%Nx;
                im=(i-1+Nx)%Nx;
                scalep=n[nn].x*n[ip+Nx*(j+Ny*k)].x + n[nn].y*n[ip+Nx*(j+Ny*k)].y + n[nn].z*n[ip+Nx*(j+Ny*k)].z;
                scalem=n[nn].x*n[im+Nx*(j+Ny*k)].x + n[nn].y*n[im+Nx*(j+Ny*k)].y + n[nn].z*n[im+Nx*(j+Ny*k)].z;
                dnx[0]=scalep*n[ip+Nx*(j+Ny*k)].x-scalem*n[im+Nx*(j+Ny*k)].x;
                dny[0]=scalep*n[ip+Nx*(j+Ny*k)].y-scalem*n[im+Nx*(j+Ny*k)].y;
                dnz[0]=scalep*n[ip+Nx*(j+Ny*k)].z-scalem*n[im+Nx*(j+Ny*k)].z;
                scalep=n[nn].x*n[i+Nx*(jp+Ny*k)].x + n[nn].y*n[i+Nx*(jp+Ny*k)].y + n[nn].z*n[i+Nx*(jp+Ny*k)].z;
                scalem=n[nn].x*n[i+Nx*(jm+Ny*k)].x + n[nn].y*n[i+Nx*(jm+Ny*k)].y + n[nn].z*n[i+Nx*(jm+Ny*k)].z;
                dnx[1]=scalep*n[i+Nx*(jp+Ny*k)].x-scalem*n[i+Nx*(jm+Ny*k)].x;
                dny[1]=scalep*n[i+Nx*(jp+Ny*k)].y-scalem*n[i+Nx*(jm+Ny*k)].y;
                dnz[1]=scalep*n[i+Nx*(jp+Ny*k)].z-scalem*n[i+Nx*(jm+Ny*k)].z;
                scalep=n[nn].x*n[i+Nx*(j+Ny*kp)].x + n[nn].y*n[i+Nx*(j+Ny*kp)].y + n[nn].z*n[i+Nx*(j+Ny*kp)].z;
                scalem=n[nn].x*n[i+Nx*(j+Ny*km)].x + n[nn].y*n[i+Nx*(j+Ny*km)].y + n[nn].z*n[i+Nx*(j+Ny*km)].z;
                dnx[2]=scalep*n[i+Nx*(j+Ny*kp)].x-scalem*n[i+Nx*(j+Ny*km)].x;
                dny[2]=scalep*n[i+Nx*(j+Ny*kp)].y-scalem*n[i+Nx*(j+Ny*km)].y;
                dnz[2]=scalep*n[i+Nx*(j+Ny*kp)].z-scalem*n[i+Nx*(j+Ny*km)].z;
                dnx[0]*=0.5*Nx/lat_size.x;
                dnx[1]*=0.5*Ny/lat_size.y;
                dnx[2]*=	Nz/lat_size.z;
                dny[0]*=0.5*Nx/lat_size.x;
                dny[1]*=0.5*Ny/lat_size.y;
                dny[2]*=	Nz/lat_size.z;
                dnz[0]*=0.5*Nx/lat_size.x;
                dnz[1]*=0.5*Ny/lat_size.y;
                dnz[2]*=	Nz/lat_size.z;
                splay1=dnx[0]+dny[1]+dnz[2];
                
                twist1=n[nn].x*(dnz[1]-dny[2])
                      +n[nn].y*(dnx[2]-dnz[0])
                      +n[nn].z*(dny[0]-dnx[1])-q0;
                bend1[0]=n[nn].y*(dny[0]-dnx[1])-n[nn].z*(dnx[2]-dnz[0]);
                bend1[1]=n[nn].z*(dnz[1]-dny[2])-n[nn].x*(dny[0]-dnx[1]);
                bend1[2]=n[nn].x*(dnx[2]-dnz[0])-n[nn].y*(dnz[1]-dny[2]);      
                bend=bend1[0]*bend1[0]+bend1[1]*bend1[1]+bend1[2]*bend1[2];
                energy[3*(nn)+0]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k11*splay1*splay1;
                energy[3*(nn)+1]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k22*twist1*twist1;
                energy[3*(nn)+2]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k33*bend;
            }
        }
        /*
    for (k = 1; k < Nz-1 ; k++)
    {
        kp=k+1; km=k-1;    
        for (j = 0; j < Ny; j++)
        {
            jp=(j+1)%Ny;
            jm=(j-1+Ny)%Ny;
            for (i = 0; i < Nx; i++)
            {
                nn=i+Nx*(j+Ny*k);
                ip=(i+1)%Nx;
                im=(i-1+Nx)%Nx;
                dir=(n[ip+Nx*(j+Ny*k)].x*n[im+Nx*(j+Ny*k)].x+n[ip+Nx*(j+Ny*k)].y*n[im+Nx*(j+Ny*k)].y+n[ip+Nx*(j+Ny*k)].z*n[im+Nx*(j+Ny*k)].z)>0?-1:1;
                dnx[0]=n[ip+Nx*(j+Ny*k)].x+dir*n[im+Nx*(j+Ny*k)].x;
                dny[0]=n[ip+Nx*(j+Ny*k)].y+dir*n[im+Nx*(j+Ny*k)].y;
                dnz[0]=n[ip+Nx*(j+Ny*k)].z+dir*n[im+Nx*(j+Ny*k)].z;
                dir=(n[i+Nx*(jp+Ny*k)].x*n[i+Nx*(jm+Ny*k)].x+n[i+Nx*(jp+Ny*k)].y*n[i+Nx*(jm+Ny*k)].y+n[i+Nx*(jp+Ny*k)].z*n[i+Nx*(jm+Ny*k)].z)>0?-1:1;
                dnx[1]=n[i+Nx*(jp+Ny*k)].x+dir*n[i+Nx*(jm+Ny*k)].x;
                dny[1]=n[i+Nx*(jp+Ny*k)].y+dir*n[i+Nx*(jm+Ny*k)].y;
                dnz[1]=n[i+Nx*(jp+Ny*k)].z+dir*n[i+Nx*(jm+Ny*k)].z;
                dir=(n[i+Nx*(j+Ny*kp)].x*n[i+Nx*(j+Ny*km)].x+n[i+Nx*(j+Ny*kp)].y*n[i+Nx*(j+Ny*km)].y+n[i+Nx*(j+Ny*kp)].z*n[i+Nx*(j+Ny*km)].z)>0?-1:1;
                dnx[2]=n[i+Nx*(j+Ny*kp)].x+dir*n[i+Nx*(j+Ny*km)].x;
                dny[2]=n[i+Nx*(j+Ny*kp)].y+dir*n[i+Nx*(j+Ny*km)].y;
                dnz[2]=n[i+Nx*(j+Ny*kp)].z+dir*n[i+Nx*(j+Ny*km)].z;
                dnx[0]*=0.5*Nx/lat_size.x;
                dnx[1]*=0.5*Ny/lat_size.y;
                dnx[2]*=0.5*Nz/lat_size.z;
                dny[0]*=0.5*Nx/lat_size.x;
                dny[1]*=0.5*Ny/lat_size.y;
                dny[2]*=0.5*Nz/lat_size.z;
                dnz[0]*=0.5*Nx/lat_size.x;
                dnz[1]*=0.5*Ny/lat_size.y;
                dnz[2]*=0.5*Nz/lat_size.z;
                splay1=dnx[0]+dny[1]+dnz[2];
                
                twist1=n[nn].x*(dnz[1]-dny[2])
                      +n[nn].y*(dnx[2]-dnz[0])
                      +n[nn].z*(dny[0]-dnx[1])-q0;
                bend1[0]=n[nn].y*(dny[0]-dnx[1])-n[nn].z*(dnx[2]-dnz[0]);
                bend1[1]=n[nn].z*(dnz[1]-dny[2])-n[nn].x*(dny[0]-dnx[1]);
                bend1[2]=n[nn].x*(dnx[2]-dnz[0])-n[nn].y*(dnz[1]-dny[2]);      
                bend=bend1[0]*bend1[0]+bend1[1]*bend1[1]+bend1[2]*bend1[2];
                energy[3*(nn)+0]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k11*splay1*splay1;
                energy[3*(nn)+1]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k22*twist1*twist1;
                energy[3*(nn)+2]=(lat_size.x*lat_size.y*lat_size.z/(Nx*Ny*Nz))*k33*bend;
            }
        }
    }*/
        
    return 0;
}
int energy_integrator (double *K, int Nx, int Ny, int Nz)
{
	double total_en[3]={0,0,0};
	for(int ii=0; ii< Nx*Ny*Nz; ii++)
	{
		total_en[0]+=K[ii*3+0];
		total_en[1]+=K[ii*3+1];
		total_en[2]+=K[ii*3+2];
	}
	printf("%g,%g,%g\n",total_en[0],total_en[1],total_en[2]);	
	
	return 0;
}
	
	
