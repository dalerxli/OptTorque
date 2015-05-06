/// https://geuz.org/trac/gmsh/wiki/STLRemeshing

Mesh.CharacteristicLengthFactor=0.08;
//Mesh.CharacteristicLengthMin=0.05;
//Mesh.CharacteristicLengthMax=0.09;

Merge "N0_400nm_MICRON_fine.stl";//also possible for .msh, .step, ... files
CreateTopology;

//Compound Surface(100)={2,3}; 
//Compound Surface(200)={1}; 
// We can now define a compound line (resp. surface) for each discrete
// line (resp. surface) in the model
ll[] = Line "*";
For j In {0 : #ll[]-1}
  Compound Line(newl) = ll[j];
EndFor
ss[] = Surface "*";
s = news;
For i In {0 : #ss[]-1}
  Compound Surface(s+i) = ss[i];
EndFor

Physical Surface(1) = {s : s + #ss[]-1};

//uniform = -1;
//If(uniform)
//  // uniform mesh size...
//  Mesh.CharacteristicLengthMin = 2.5;
//  Mesh.CharacteristicLengthMax = 2.5;
//EndIf

Mesh.CharacteristicLengthMin=0.50;
Mesh.CharacteristicLengthMax=0.90;

Mesh.RemeshAlgorithm = 1; // (0) no split (1) automatic (2) automatic only with metis
Mesh.RemeshParametrization = 7; // (0) harmonic (1) conformal spectral (7) conformal finite element
//0=harmonic_circle, 1=conformal_spectral, 2=rbf, 3=harmonic_plane, 4=convex_circle, 5=convex_plane, 6=harmonic square,7=conformal_fe 

Geometry.HideCompounds = 0; // don't hide the compound entities
//Mesh.Algorithm = 6; // Frontal

