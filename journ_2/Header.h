#pragma once



void init();
void discretize_space();
void coefficients();
void bound();
void solver();
void force();

void print();


double jr, bl, cb, mi, rho, omega, beta, epsilon, pps, pla, pa, pdh;

int nx, nz, L, i, j, ismax, ie, iw, in, is, iter, maxiter, posmaxp;

double lx, lz, dx, dz, hx, hz, error, tol, relax, res, area;

double U, par1, par2, par3, hmin, thetal, thetar, maxp, angle, hminn,Ftotalx,Ftotaly,ftemp , fang,Re,Re_min;

const double pi = 3.14159265358979323846;

double x[20000000], z[20000000], fi[2000000], h[20000000];

double A[20000000], B[20000000], C[20000000], D[20000000], E[20000000], dsol[20000000], sol[20000000];

int ii[3000][2500], inlet[20000000];