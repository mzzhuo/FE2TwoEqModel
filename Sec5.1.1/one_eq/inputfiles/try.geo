
L = 1; 
D = 0.00357143; 

n_x = 280;          // no. of rve in x 
n_y = 1;          // no. of rve in y
n   = n_x * n_y;  // total no. of RVEs

lc1 = L/n_x/2/5;
lc2 = L/n_x/2/5;

x = L/n_x/2;
y = D/n_y/2;
r = L/n_x * 0.3250;


Point(1) = {0.0, 0.0, 0.0, lc1};
Point(2) = {x  , 0.0, 0.0, lc1};
Point(3) = {x  , y  , 0.0, lc1};
Point(4) = {x-0.7071*r, y-0.7071*r, 0.0, lc1};
Point(5) = {x, y-r, 0.0, lc1};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {4, 1};
//
Line(4) = {3, 4};
Line(5) = {3, 5};
Circle(6) = {4,3,5}; 
//
//Transfinite Line { 4,5 } = 4; // Using Progression 1.1;
//Transfinite Line { 3 } = 5; // Using Progression 1.1;
//Transfinite Line { 1 } = 6; // Using Progression 1.1;
Transfinite Line { 2 } = 5 Using Progression 1.1;
Transfinite Line { 3 } = 5 Using Progression 1.1;


linelooprve = newreg;
Line Loop(linelooprve) = {1,2,-6,3};

lineloopcir = newll;
Line Loop(lineloopcir) = {4,6,-5};

psrve = newreg;
Plane Surface(psrve) = {linelooprve};
Printf("Base Matrix Plane Surface: %g.", psrve);

pscircle = newreg;
Plane Surface(pscircle) = {lineloopcir};
Printf("Base Inclusion Plane Surface: %g.", pscircle);

//For i In {1:1}
//  my_new_surfs[] = Translate {i*L/n_x, 0, 0} { Duplicata { Surface{psrve,pscircle}; } };
//  Printf("New surfaces ’%g’ and ’%g’", my_new_surfs[0], my_new_surfs[1]);
//  new_surfs_matrix[2*i-1] = my_new_surfs[0];
//  new_surfs_inclus[2*i  ] = my_new_surfs[1];
//EndFor


Physical Surface("matrix") = {psrve};//,new_surfs_matrix[]};
Physical Surface("inclusion") = {pscircle};//,new_surfs_inclus[]};



//Physical Line("left")  = {3};
//Physical Line("right") = {4};// {2805}; // manually pick

//Mesh 2;
//RefineMesh;
