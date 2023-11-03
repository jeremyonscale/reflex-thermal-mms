Mesh.Algorithm = 5;
If (Mesh.RecombineAll != 0)
  Mesh.SubdivisionAlgorithm = 2; // Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3: barycentric)
EndIf

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.MeshSizeMin = 1/n;
Mesh.MeshSizeMax = 1/n;
