// Gmsh project created on Fri Mar 29 12:33:33 2013

// characteristic size of elements and number of elements in height
dx = 100;
nH = 20;
H = nH*dx;
nL = 2*nH;
L = nL*dx;

// basic line (including future PML on each side)
Point(1) = {-dx, -dx, 0, dx};
Point(2) = {0, -dx, 0, dx};
Point(3) = {L, -dx, 0, dx};
Point(4) = {L+dx, -dx, 0, dx};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

// characteristic length of the elements
Characteristic Length {2, 1, 3, 4} = dx;

// Extrude line to surface air-soil
out1[] = Extrude{0,dx,0}{Line{1,2,3}; Layers{1};Recombine;};
out2[] = Extrude{0,L,0}{Line{out1[4],out1[8],out1[0]}; Layers{nL};Recombine;};
out3[] = Extrude{0,dx,0}{Line{out2[4],out2[8],out2[0]}; Layers{1};Recombine;};

// Extrude surface air-soil to Volume
out4[] = Extrude{0,0,-H}{Surface{out1[1],out1[5],out1[9]}; Layers{nH};Recombine;};
out5[] = Extrude{0,0,-H}{Surface{out2[1],out2[5],out2[9]}; Layers{nH};Recombine;};
out6[] = Extrude{0,0,-H}{Surface{out3[1],out3[5],out3[9]}; Layers{nH};Recombine;};
out7[] = Extrude{0,0,-dx}{Surface{out4[0],out4[6],out4[12]}; Layers{1};Recombine;};
out8[] = Extrude{0,0,-dx}{Surface{out5[0],out5[6],out5[12]}; Layers{1};Recombine;};
out9[] = Extrude{0,0,-dx}{Surface{out6[0],out6[6],out6[12]}; Layers{1};Recombine;};

// output message: stats about the created mesh
Printf("Size of the domain: %g -by-  %g -by-  %g", L, L, H);
Printf("Number of Vertices: %g",  (nL+2)*(nL+2)*(nH+1));
Printf("Number of Elements: %g",  (nL+1)*(nL+1)*nH);
Printf("Expected number of DOFs (8GLL per direction per element): %g",
                                  (7*(nL+1)+1)*(7*(nL+1)+1)*(7*nH+1));

