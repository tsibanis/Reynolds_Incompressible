clear; clc;
%rawdata = readmatrix('pressure.csv');
x = readmatrix('x.dat');       % Reshape the column matrix into 51 columns
%x=x';
y = readmatrix('y.dat');  
y=y*1000;
%y=y';
z = readmatrix('z.dat'); 

%z=z';
%N = 50000;

%[X, Y] = ndgrid(x, y);
%[nrow,ncol] = size(z)


surf(x, y, z);
shading interp
colorbar
colormap(jet(256));
%xlim([0 7]);
%ylim([0 300]);

xlabel('Theta (rad)')
ylabel('Bearing Axis (mm)')
legend('Pressure (Pa)')


%h.Limits = [-100 1.4E7]
%clim([-100, 10000]);