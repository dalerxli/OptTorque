//////////////////////////////////////////////////
//////////////////////////////////////////////////
////  GEO file for flat disc monitor,
////  displaced in z 
//////////////////////////////////////////////////
//////////////////////////////////////////////////
L=2.0; 	// diameter in MICRONS.
Z=0.25;         // monitor displacement from z=0
R=L/2;       // radius 
l=L/30;      // meshing fineness 
//////////////////////////////////////////////////
//////////////////////////////////////////////////
CX = 0.0; //center location 
CY = 0.0; 
CZ = Z; 
Point(1)= { CX, CY, CZ, l }; // center 

Point(10) = {     CX,     CY,   CZ,   l};
Point(11) = {     CX, CY + R,   CZ,   l};
Point(12) = { CX - R,     CY,   CZ,   l};
Point(13) = {     CX, CY - R,   CZ,   l};
Point(14) = { CX + R,     CY,   CZ,   l};

Circle(15) = { 11, 10, 12 };
Circle(16) = { 12, 10, 13 };
Circle(17) = { 13, 10, 14 };
Circle(18) = { 14, 10, 11 };

Line Loop (6) = { 15, 16, 17, 18 };
Plane Surface(1) = { 6 }; 
Physical Surface(1) = {1};   
