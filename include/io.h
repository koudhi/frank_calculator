/*
 * io.cpp
 * created by: Eric Koudhi Omori <koudhi@gmail.com>
 * Creation date: 2019-06-24
 * Last Modification: 2019-06-24
 * version: 0.0
 */

#include <iostream>
int get_size(FILE *input, vector &lat_size, int Nx, int Ny, int Nz);
int read_file(char* fname, vector *n, vector &lat_size, int Nx, int Ny, int Nz);
int print_file(char *fname, vector *n, double *K, vector lat_size, int Nx, int Ny, int Nz);
