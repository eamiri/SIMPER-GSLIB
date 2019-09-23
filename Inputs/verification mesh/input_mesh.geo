//Inputs
width=0.01;
depth=10;
gridsize=0.5;


Point(1)={0,0,0,gridsize};
Point(2)={width,0,0,gridsize};
Point(3)={width,0.9*depth,0,0.35*gridsize};
Point(4)={0,0.9*depth,0,0.35*gridsize};
Point(5)={0,depth,0,gridsize};
Point(6)={width,depth,0,0.35*gridsize};


Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line(5)={3,6};
Line(6)={6,5};
Line(7)={5,4};

Transfinite Line {1,3,6} = Ceil(3) Using Progression 1;
Transfinite Line {7,5} = Ceil(201) Using Progression 1;
Transfinite Line {4,2} = Ceil(15) Using Progression 1;

Line Loop(1)={1,2,3,4};
Plane Surface(1) = {1};
Line Loop(2)={-3,5,6,7};
Plane Surface(2) = {2};
Transfinite Surface{1,2};
Recombine Surface{1,2};

Physical Surface("Fine") = {1};
Physical Surface("Coarse") = {2};