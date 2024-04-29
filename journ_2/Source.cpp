#include <iostream>
#include <tgmath.h>
#include <math.h>
#include<fstream>




void init();
void discretize_space();
void coefficients();
void bound();
void solver();
void force();

void print();


double jr, bl, cb, mi, rho, omega, beta, epsilon, pps, pla, pa, pdh;

int nx, nz, L, i, j, ismax, ie, iw, in, is, iter, maxiter, posmaxp;

double lx, lz, dx, dz, hx, hz, error, tol, relax, res, area, time_before, time_after;

double U, par1, par2, par3, hmin, thetal, thetar, maxp, angle, hminn, Ftotalx, Ftotaly, ftemp, fang, Re, Re_min;

const double pi = 3.14159265358979323846;

double x[500000], z[500000], fi[500000], h[500000], A[500000], B[500000], C[500000], D[500000], E[500000], dsol[500000], sol[500000];;

int ii[5000][5000], inlet[500000];




void init() {

	jr = 63.5E-3;		//m
	bl = 127E-3;		//m
	cb = 0.0635E-3;		//m
	mi = 0.0138;		//Pa*s
	rho = 887;			//kg/m^3
	omega = 5500;		//rpm
	beta = 330.07;		//deg
	epsilon = 0.46743;	//E/cb
	pps = 100;			//Pa
	pla = 85;			//fi
	pa = 10;			//dfi
	pdh = 0.1E-3;		//m

	//Rotational Speed (RPM) --> Linear Speed (m/s) 
	U = omega * 2 * pi * jr / 60;

	//Dimentions of plane
	nx = 100;
	nz = 100;
	lx = 2 * pi * jr ;
	lz = bl ;
	dx = lx / (nx - 1);
	dz = lz / (nz - 1);

	//area of each cell
	area = dx * dz;

	//parameters used in the calculation of the coefficients
	par1 = (1 / pow(dx, 2)) + (1 / pow(dz, 2));
	par2 = 1 / pow(dx, 2);
	par3 = 1 / pow(dz, 2);

	hmin = 100000;
	maxp = 0;

	//actual min film thickness
	hminn = cb - cb * epsilon;


	//Angle between the eccentricity axis and xx'
	angle = 2 * pi - beta * (pi / 180);

	//angles of left and right sides of inlet
	thetal = pi + angle + (pla + pa) * (pi / 180);
	thetar = pi + angle + (pla) * (pi / 180);

	maxiter = 1000000;
	error = 1;
	tol = 1E-5;
	relax = 1;

	Re = rho * U * cb / mi;		//Reynolds number

	L = 0;
	for (i = 0; i <= nx+1; i++) {
		for (j = 0; j <= nz+1; j++) {
			ii[i][j] = 0;
			x[L] = 0;
			z[L] = 0;
			fi[L] = 0;
			h[L] = 0;
			A[L] = 0;
			B[L] = 0;
			C[L] = 0;
			D[L] = 0;
			E[L] = 0;
			sol[L] = 0;
			dsol[L] = 0;
			inlet[L] = 0;
			
			L++;
		}
	}

	ismax = 0;
	for (i = 1; i <= nx; i++) {
		for (j = 1; j <= nz; j++) {
			ismax = ismax + 1;
			ii[i][j] = ismax;
		}
	}

}

void discretize_space() {

	for (i = 1; i <= nx; i++) {
		for (j = 1; j <= nz; j++) {
			L = ii[i][j];
			x[L] = (i - 1) * dx;  //circumfrential coordinates
			z[L] = (j - 1) * dz; //axial coordinates
			fi[L] = (x[L] / jr); //angle relative to the eccentricity axis
			h[L] = cb*(1 + epsilon * cos(fi[L])); //film thickness
			if ( fi[L] < thetal && fi[L] > thetar ) { //check if in abissica of inlet 
				h[L] = h[L] +pdh; 
				inlet[L] = 1;
			}
			if (h[L] <= hmin && h[L]!=0) {
				hmin = h[L];
			}

		}
	}

}

void coefficients() {

	for (i = 2; i <= nx-1; i++) {
		for (j = 2; j <= nz-1; j++) {
			L = ii[i][j];
			in = ii[i + 1][j];
			is = ii[i - 1][j];
			ie = ii[i][j + 1];
			iw = ii[i][j - 1];
			
			hx = (h[in] - h[is]) / (2 * dx);		//derivative dh/dx
			hz = (h[ie] - h[iw]) / (2 * dz);		//derivative dh/dz


			A[L] = (par2 - 3 * hx / (2 * dx * h[L])) / (2 * par1);
			B[L] = (par3 - 3 * hz / (2 * dz * h[L])) / (2 * par1);
			C[L] = (par2 + 3 * hx / (2 * dx * h[L])) / (2 * par1);
			D[L] = (par3 + 3 * hz / (2 * dz * h[L])) / (2 * par1);
			E[L] = -3*mi*U*hx/(pow(h[L],3)*par1);

		}
	}

}

void bound() {

	// oil inlet region 
	for (i = 2; i <= nx-1; i++) {
		for (j = 2; j <= nz-1; j++) {
			L = ii[i][j];
			if (inlet[L] == 1) {
				A[L] = 0;
				B[L] = 0;
				C[L] = 0;
				D[L] = 0;
				E[L] = -pps;
			}
		}
	}
	// i=1 Boundary
		i = 1;
		for (j = 2; j <= nz-1; j++) {
			L = ii[i][j];
			in = ii[i + 1][j];
			is = ii[nx-1][j];
			ie = ii[i][j + 1];
			iw = ii[i][j - 1];

			hx = (h[in] - h[is]) / (2 * dx);		//derivative dh/dx
			hz = (h[ie] - h[iw]) / (2 * dz);		//derivative dh/dz

			A[L] = (par2 - 3 * hx / (2 * dx * h[L])) / (2 * par1);
			B[L] = (par3 - 3 * hz / (2 * dz * h[L])) / (2 * par1);
			C[L] = (par2 + 3 * hx / (2 * dx * h[L])) / (2 * par1);
			D[L] = (par3 + 3 * hz / (2 * dz * h[L])) / (2 * par1);
			E[L] = -3 * mi * U * hx / (pow(h[L], 3) * par1);
			
		}
	//i=nx Boundary
		i = nx;
		for (j = 1; j <= nz; j++) {
			L = ii[i][j];
			in = ii[2][j];
			is = ii[i - 1][j];
			ie = ii[i][j + 1];
			iw = ii[i][j - 1];

			hx = (h[in] - h[is]) / (2 * dx);		//derivative dh/dx
			hz = (h[ie] - h[iw]) / (2 * dz);		//derivative dh/dz

			A[L] = (par2 - 3 * hx / (2 * dx * h[L])) / (2 * par1);
			B[L] = (par3 - 3 * hz / (2 * dz * h[L])) / (2 * par1);
			C[L] = (par2 + 3 * hx / (2 * dx * h[L])) / (2 * par1);
			D[L] = (par3 + 3 * hz / (2 * dz * h[L])) / (2 * par1);
			E[L] = -3 * mi * U * hx / (pow(h[L], 3) * par1);
		}
	
	//j=1 Boundary
		j = 1;
		for (i = 1; i <= nx; i++) {
			L = ii[i][j];
			A[L] = 0;
			B[L] = 0;
			C[L] = 0;
			D[L] = 0;
			E[L] = -pps;
		}
	
	//j=nz Boundary
		j = nz;
		for (i = 1; i <= nx; i++) {
			L = ii[i][j];
			A[L] =0;
			B[L] = 0;
			C[L] = 0;
			D[L] = 0;
			E[L] = -pps;
		}
}

void solver() {
	
	//Gauss-Seidel with over-relaxation coefficient

	iter = 1;
	while ( iter < maxiter && error>tol  ) {

		error = 0;

		for (i = 1; i <= nx; i++) {
			for (j = 1; j <= nz; j++) {

				L = ii[i][j];
				in = ii[i + 1][j];
				is = ii[i - 1][j];
				ie = ii[i][j + 1];
				iw = ii[i][j - 1];
				
				sol[L] = (1 - relax) * dsol[L] + relax * (A[L] * dsol[is] + B[L] * dsol[iw] + C[L] * dsol[in] + D[L] * dsol[ie] + E[L]);
				
				//Check if in cavitation region
				if (sol[L]/pps < 1 ) {
					sol[L] = -pps;
				}
				if (sol[L] == 0) {
					res = 0;
				}
				else {
					res = (sol[L] - dsol[L]) / sol[L];
				}
				error = error + pow(res, 2);
			}
		}
		iter++;
		for (i = 1; i <= nx ; i++) {
			for (j = 1; j <= nz; j++) {
				L = ii[i][j];
				dsol[L] = sol[L];
			}
		}
		error = sqrt(error);
	}

}

void force() {

	Ftotalx = 0;
	Ftotaly = 0;
	for (i = 1; i <= nx; i++) {
		for (j = 1; j <= nz; j++) {
			L = ii[i][j];
			if (sol[L] >= maxp) {
				maxp = sol[L];
				posmaxp = L;
			}
			if (i == 1 || i == nx || j == 1 || j == nz) {
				if ((i == 1 && j == 1) || (i == 1 && j == nz) || (i == nx && j == nz) || (i == nx && j == 1)) {
					Ftotalx = Ftotalx - sol[L] * cos(fi[L] - angle + pi) * dx * dz/4;
					Ftotaly = Ftotaly - sol[L] * sin(fi[L] - angle + pi) * dx * dz/4;
				}
				else {
					Ftotalx = Ftotalx - sol[L] * cos(fi[L] - angle + pi) * dx * dz / 2;
					Ftotaly = Ftotaly - sol[L] * sin(fi[L] - angle + pi) * dx * dz / 2;
				}
			}
			else {
				Ftotalx = Ftotalx - sol[L] * cos(fi[L] - angle + pi) * dx * dz;
				Ftotaly = Ftotaly - sol[L] * sin(fi[L] - angle + pi) * dx * dz;
			}
		}	
	}
}

void print() {

	

	std::ofstream xfile("x.dat");
	for (i = 1; i <= nx; i++) {
		for (j = 1; j <= nz; j++) {
			L = ii[i][j];
			xfile << fi[L] << " ";
		}
		xfile << std::endl;
	}
	xfile.close();


	std::ofstream zfile("y.dat");
	for (i = 1; i <= nx; i++) {
		for (j = 1; j <= nz; j++) {
			L = ii[i][j];
			zfile << z[L] << " ";
		}
		zfile << std::endl;
	}
	zfile.close();

	std::ofstream rfile("r.dat");
	for (i = 1; i <= nx; i++) {
		j = nz / 2;
		L = ii[i][j];
		rfile << cos(fi[L] - angle + pi) * ((sol[L]) + 1) << " " << sin(fi[L] - angle + pi) * ((sol[L]) + 1) << std::endl;
		
	}
	rfile.close();

	std::ofstream sfile("s.dat");
	for (i = 1; i <= nx; i++) {
		j = nz / 2;
		L = ii[i][j];
		sfile << cos(fi[L] - angle + pi) << " " << sin(fi[L] - angle + pi) << std::endl;
		
	}
	sfile.close();

	std::ofstream tfile("t.dat");
	for (i = 1; i <= nx; i++) {
		j = nz / 2;
		L = ii[i][j];
		tfile << fi[L] << " " << sol[L] << std::endl;

		tfile << std::endl;
	}
	tfile.close();

	std::ofstream pressfile("z.dat");
	for (i = 1; i <= nx; i++) {
		for (j = 1; j <= nz; j++) {
			L = ii[i][j];
			pressfile << sol[L] << " " ;
		}
		pressfile << std::endl;
	}
	pressfile.close();

}

int main() {

	time_before = time(NULL);

	init();
	discretize_space();
	coefficients();
	bound();
	solver();
	force();
	print();

	time_after = time(NULL);

	std::cout <<"Speed of revolution (rad/s) " << U/jr ;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Re :  " << Re;
	std::cout << std::endl;
	std::cout << "hmin :  " << hmin;
	std::cout << std::endl;
	std::cout << "p_max  :  " << maxp;
	std::cout << std::endl;
	std::cout << "p_max position x,z,fi :  " << x[posmaxp]<<" , "<<z[posmaxp]<<" , "<<fi[posmaxp]*180/pi;
	std::cout << std::endl;
	std::cout << "Fmax  :  " << maxp*area;
	std::cout << std::endl;
	std::cout << "error  :  " << error;
	std::cout << std::endl;
	std::cout << "Fx  :  " << Ftotalx;
	std::cout << std::endl;
	std::cout << "Fy  :  " << Ftotaly;
	std::cout << std::endl;
	std::cout << "simulation time  :  " << time_after-time_before;
	std::cout << std::endl;
	
	return 0;
}