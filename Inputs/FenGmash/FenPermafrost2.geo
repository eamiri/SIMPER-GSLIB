//Inputs
width=20;
depth=5;
gridsize=depth/8;
pthik=0.6*depth;
pcent=0.6*depth;
phtop=pcent+pthik/2;
ptbot=pcent-pthik/2;

Point(1)={0,0,0,gridsize};
Point(2)={width,0,0,gridsize};
Point(3)={width,depth,0,0.35*gridsize};
Point(4)={0,depth,0,0.35*gridsize};
Point(5)={0,3*depth/5,0,gridsize};
Point(6)={0.2*width,depth,0,0.35*gridsize};
Point(7)={width,ptbot,0,gridsize};
Point(8)={width,pcent,0,gridsize};
//Point(9)={0.4*width,pcent,0,0.35*gridsize};
Point(10)={width,phtop,0,0.35*gridsize};
//Point(11)={0.42*width,0.7*depth,0,0.35*gridsize};
//Point(12)={0.47*width,0.8*depth,0,0.35*gridsize};
//Point(13)={0.405*width,0.65*depth,0,0.35*gridsize};
//Point(14)={0.5*width,0.9*depth,0,0.35*gridsize};

Ellipse(1)={5,4,6,6};
//Ellipse(2)={7,8,9,9};
//Ellipse(3)={10,8,9,9};
//BSpline(3)={9,13,11,12,14,10};
Line(4)={1,2};
Line(5)={2,7};
Line(6)={8,10};
Line(7)={10,3};
Line(8)={3,6};
Line(9)={6,4};
Line(10)={4,5};
Line(11)={5,1};
Line(12)={7,8};

Line Loop(12)={4,5,12,6,7,8,-1,11};
Plane Surface(1)=12;
Line Loop(14)={1,9,10};
Plane Surface(2)=14;
//Line Loop(12)={4,5,2,-3,7,8,-1,11};
//Line Loop(13)={3,-2,12,6};
//Plane Surface(2)=13;


Recombine Surface{1,2,3};//+
Physical Surface("Fen") = {2};
//+
//Physical Surface("Permafrost") = {2};
//+
Physical Surface("Soil") = {1};
