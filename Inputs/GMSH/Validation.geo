//Inputs
width=0.01;
depth=1;
gridsize=0.001;


Point(1)={0,0,0,gridsize};
Point(2)={width,0,0,gridsize};
Point(3)={width,depth,0,gridsize};
Point(4)={0,depth,0,gridsize};


Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Transfinite Line {1,3} = 3 Using Progression 1;
Transfinite Line {4,2} = 101 Using Progression 1;

Line Loop(1)={1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1};

Physical Surface("Fine") = {1};