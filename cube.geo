SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};

Physical Volume("bulk", 1) = {1};
Physical Surface("left", 2) = {1};
Physical Surface("right", 3) = {2};
Physical Surface("bottom", 4) = {3};
Physical Surface("top", 5) = {4};
Physical Surface("front", 6) = {6};
Physical Surface("back", 7) = {5};

n = 1;

