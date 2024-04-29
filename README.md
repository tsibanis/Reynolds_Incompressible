This code is developed by George Kotsidis.

The program consists of several classes that work together
to perform the desired task. The main steps involved in the
simulation are:
Initialization of variables, such as fluid properties, cavity dimensions, and simulation parameters, in the init() function.
Discretization of the simulation domain into a grid of cells
and determination of the coordinates and geometry of each cell
in the discretize space() function.
Calculation of coefficients in the discretized domain for
use in the numerical solution of the flow equations in the
coefficients() function.
Application of boundary conditions in the bound() function.
Solution of the flow equations using an iterative solver in
the solver() function.
Calculation of the force on the fluid in the cavity due to the
rotating motion in the force() function.
Printing of the results in the print() function.
The variables used in the code include:
• jr: the radius of the shaft.
• bl: the length of the shaft.
• cb: the radial clearance between the rotor and stator.
• mi: the dynamic viscosity of the fluid.
• rho: the density of the fluid.
• omega: the angular velocity of the rotor.
• beta: the angle between the eccentricity axis and the
x-axis.
• epsilon: the eccentricity of the rotor.
• pps: the pressure at the supply pocket.
• pla: the angle of the left side of the supply pocket
relative to the real x-axis.
• pa: pocket angle.
• pdh: pocket depth.
• nx: the number of cells in the x-direction.
• nz: the number of cells in the z-direction.
• L: a counter variable used to track the total number of
cells.
• i, j: indices used to loop through the cells.
• ismax: a counter variable used to assign an index to each
cell.
• is, ie, in, iw: indices used to track the location of the
stencil cells.
• tol: the tolerance used for the iterative solver.
• relax: the relaxation factor used for the iterative solver.
• U: the linear velocity of the rotor.
• par1, par2, par3: variables used to calculate the
coefficients in the discretized domain.
• hmin: the minimum value of film thickness.
• thetal, thetar: the angle of the rotor position relative
to the x-axis, with an offset for the inlet.
• maxp: the maximum value of the pressure.
• angle: the angle between the eccentricity axis and the
x-axis.
• hminn: the minimum value of the film thickness.
• Re: the Reynolds number.
• Remin: the minimum Reynolds number at the inlet.
• x: an array of the x-coordinates of the cells.
• z: an array of the z-coordinates of the cells.
• f i: an array of the angles of the cells relative to the
eccentricity axis.
• h: an array of the film thickness for each cell.
• A, B, C, D, E: arrays used to store the coefficients in
the discretized domain.
• dsol: an array used to store the change in the solution
during each iteration of the solver.
• sol: an array used to store the solution.
• ii: a 2D array used to store the index of each cell.
• inlet: an array used to track which cells are part of the
inlet.
