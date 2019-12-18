//Inputs
width=50.0;
depth=5.0;
gridsize=2.50;


Point(1)={0,0,0,gridsize};
Point(2)={width,0,0,gridsize};
Point(3)={width,0.4*depth,0,gridsize};
Point(4)={0,0.4*depth,0,gridsize};
Point(5)={0,depth,0,gridsize};
Point(6)={width,depth,0,gridsize};


Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line(5)={3,6};
Line(6)={6,5};
Line(7)={5,4};

Transfinite Line {1,3,6} = Ceil(30) Using Progression 1;
Transfinite Line {7,5} = Ceil(4) Using Progression 1;
Transfinite Line {4,2} = Ceil(3) Using Progression 1;

Line Loop(1)={1,2,3,4};
Plane Surface(1) = {1};
Line Loop(2)={-3,5,6,7};
Plane Surface(2) = {2};
Transfinite Surface{1,2};
Recombine Surface{1,2};

Physical Surface("Peat") = {1,2};