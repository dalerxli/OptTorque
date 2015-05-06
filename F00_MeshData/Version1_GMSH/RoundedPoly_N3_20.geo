//////////////////////////////////////////////////
//////////////////////////////////////////////////
//  GEO file for equilateral flat polygon with rounded corners
//  rightmost edge is aligned with y axis 
//  Last edit on: 2014/11/21
//////////////////////////////////////////////////
//////////////////////////////////////////////////

N  = 3; 	// Symmetry of polygon (starting from 3). 
A  = 2*Pi/N; 	// Angle for each edge
L  = 0.4; 	// Edge length in MICRONS. (0.4um)
T  = 0.04;  	// Thickness (0.04um)
CR = L/(2*Sin(A/2)); // Centoid radius
//////////////////////////////////////////////////
// meshing fineness
l=L/20; 

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Point and Line definition to form a single side 
Point(100)  = {0,0,0, l }; // origin

Point(101)  = { CR*Cos(A/2),     CR*Sin(-A/2),  0.5*T, l }; // vertex 
Point(102)  = { CR*Cos(A/2),     CR*Sin(-A/2),      0, l }; //center of rounding arc 
Point(103)  = { CR*Cos(A/2)+T/2, CR*Sin(-A/2),      0, l }; // rounding 
Point(105)  = { CR*Cos(A/2),     CR*Sin(-A/2), -0.5*T, l }; // vertex 

Circle(101) = { 101, 102, 103 }; 
Circle(102) = { 103, 102, 105 };

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Extrude one side
Extrude{0,L,0} { Line{101, 102}; } 
Extrude{ {0,0,1}, {CR*Cos(A/2), CR*Sin(-A/2), 0}, -A} { Line{101, 102}; } 

//////////////////////////////////////////////////
//////////////////////////////////////////////////
top_poly[] = {104}; // Top polygon group 
bot_poly[] = {109}; // Bottom polygon group 
For nn In {1:(N-1)}
  top_poly[] = {top_poly[], Rotate{{0,0,1},{0,0,0},A*nn}{Duplicata{ Line{104}; } } };
  bot_poly[] = {bot_poly[], Rotate{{0,0,1},{0,0,0},A*nn}{Duplicata{ Line{109}; } } };
  Rotate{ {0,0,1}, {0,0,0}, A*nn}{ Duplicata { Line{109, 104}; } } 
  Rotate{ {0,0,1}, {0,0,0}, A*nn}{ Duplicata { Surface{106, 110, 113, 116}; } }	
EndFor
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Line Loop(1000) = {top_poly[]};
Ruled Surface(1001) = {1000};

Line Loop (1002) = {bot_poly[]};
Ruled Surface(1003) = {1002}; 

Show "*";
