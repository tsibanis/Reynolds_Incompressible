x = readmatrix('r.dat');
plot_vectors(x,0);

p = nsidedpoly(1000, 'Center', [-0.1 0.1], 'Radius', 1.2);
plot(p);
axis equal

plot(0,0,".",'Color',"black")
plot(-0.1,0.1,".",'Color',"black")

e=-1.5:0.1:1.5;

plot(e,-e)

legend('Pressure Distribution Profile (dimensionless) ')


