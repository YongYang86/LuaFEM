nx = 1024;
ny = nx;

Point(1) = {0,0,0,0.1};

Extrude {1,0,0} {
  Point{1}; Layers{nx};
}
Extrude {0,1,0} {
  Line{1}; Layers{ny};
}

Physical Line(0) = {3};
Physical Line(1) = {4};
Physical Line(2) = {1};
Physical Line(3) = {2};
Physical Surface(10) = {5};
