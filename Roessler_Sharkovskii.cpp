/* ========================================================================
 * From the Sharkovskii theorem to periodic orbits for the Roessler system
 * Anna Gierzkiewicz and Piotr Zgliczyński
 * ========================================================================
 * Program: Roessler_Sharkovskii.cpp
 * Date: 22 VI 2021
 * ========================================================================
 * The program contains the proofs of existence of contracting grids
 * in the Roessler system with four values of parameter a = 4.7, 
 * declared in Section 6.
 * ========================================================================*/ 

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "utils.h"

// Main function

int main()
{
  	cout.precision(16);
	cout << boolalpha;  
	try
	{			
		system3d roessler525(interval(5.25));						// The Roessler system with a=5.25 defined 
		system3d roessler47(interval(4.7));						
		system3d roessler4381(interval(4.381));						
		system3d roessler542(interval(5.42));						
		///===================== variables used in Procedure 1:  =====================
		HSet2D grid3 (IVector ({-6.38401, 0.0327544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({3.63687,0.0004}));			// Attractor's container		
		vector<HSet2D> c3(3);
		c3[0] = HSet2D(IVector ({-3.46642, 0.0346316}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.072,0.00048}));			// cube 1
		c3[1] = HSet2D(IVector ({-6.26401, 0.0326544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.162,0.00066}));			// cube 2
		c3[2] = HSet2D(IVector ({-9.74889, 0.0307529}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.036,0.00072}));			// cube 3

		///===================== variables used in Procedure 2:  =====================	
		HSet2D grid5 (IVector ({-6.1885, 0.0356707}) , IMatrix ({{-1., 0.000778356}, {-0.000778356, -1.}}) , DVector({2.68797,0.0004}));			// Attractor's container		
		vector<HSet2D> c5(5);
		c5[0] = HSet2D(IVector ({-3.86108, 0.0375827}) , IMatrix ({{0.0693366, 1.}, {-0.997593, 0.000984231}}) , DVector({0.0006, 0.00138}));			// cube 1		
		c5[1] = HSet2D(IVector ({-6.82009, 0.0350822}) , IMatrix ({{0.7879108, 1.}, {0.615789, 0.0007307}}) , 
			DVector({0.0012, 0.0024}));		 		// cube 2  
		c5[2] = HSet2D(IVector ({-7.83056, 0.0343732}) , IMatrix ({{0.8138516, 1.}, {0.581073, 0.000671}}) , 
			DVector({0.0012, 0.0042}));				// cube 3  
		c5[3] = HSet2D(IVector ({-5.75153, 0.0359038}) , IMatrix ({{0.9997319, 1.}, {-0.023153, 0.0008062}}) , 
			DVector({0.0228, 0.01116}));			// cube 4		
		c5[4] = HSet2D(IVector ({-8.73615, 0.0337875}) , IMatrix ({{0.8997843, 1.}, {0.436335, 0.00062508}}) , 
			DVector({0.00144, 0.000744}));			// cube 5
	
		
		///===================== variables used in Procedure 3:  =====================		
		IVector pCloseTo6periodic({-7.43681, 0.0363713});							// starting point for INM 
		IVector containsPeriodicPoint;												// will contain periodic point, also in Proc. 5
		
		///===================== variables used in Procedure 4:  =====================		
		HSet2D grid6 (IVector ({-5.99932, 0.0376868}) , IMatrix ({{1., -0.000899679}, {0.000899679, 1.}}) , DVector({2.24683, 0.00022}));			// Attractor's container
		vector<HSet2D> c6(6);	
		c6[0] = HSet2D(IVector ({-7.44827, 0.0363852}) , IMatrix ({{1., 0.8498}, {0.0007825, 0.527106}}) , DVector({0.00225, 0.0005}));			// cube 1			
		c6[1] = HSet2D(IVector ({-5.43268, 0.038121}) , IMatrix ({{1.23042, 0.696746}, {-0.00567154, -0.0240555}}) , 
			DVector({0.00509, 0.015}));			// cube 2  			
		c6[2] = HSet2D(IVector ({-8.14614, 0.0358553}) , IMatrix ({{1., 0.907289}, {0.000736978, 0.420507}}) , 
			DVector({0.000265, 0.00085}));			// cube 3  			
		c6[3] = HSet2D(IVector ({-4.05249, 0.0395383}) , IMatrix ({{0.999999, 0.155044}, {0.00111181, -0.987908}}) , 
			DVector({0.000485, 0.00035}));			// cube 4				
		c6[4] = HSet2D(IVector ({-6.98865, 0.0367524}) , IMatrix ({{1., 0.827066}, {0.000815525, 0.562105}}) , 
			DVector({0.000712, 0.0005}));			// cube 5			
		c6[5] = HSet2D(IVector ({-6.38538, 0.0372585}) , IMatrix ({{1., 0.834783}, {0.000863145, 0.550579}}) , 
			DVector({0.00149, 0.0006}));			// cube 6
		
		///===================== variables used in Procedure 5:  =====================		
		IVector pCloseTo6periodic2({-3.330388727949297,0.03381008102270861});		// starting point for INM 
		
		///===================== variables used in Procedure 6:  =====================		
		HSet2D grid6b (IVector ({-6.60556, 0.0317909}) , IMatrix ({{1., -0.000573253}, {0.000573253, 1.}}) , DVector({3.57445, 0.00035}));			// Attractor's container	
		vector<HSet2D> c6b(6);	
		c6b[0] = HSet2D(IVector ({-3.33039, 0.0338101}) , IMatrix ({{1., 0.0114844}, {0.000763188, -0.999934}}) , DVector({0.0015225, 0.000525}));			// cube 1			
		c6b[1] = HSet2D(IVector ({-6.04388, 0.0319883}) , IMatrix ({{1., 0.566012}, {0.000593828, -0.824397}}) , 
			DVector({0.0029925, 0.0005775}));			// cube 2  						
		c6b[2] = HSet2D(IVector ({-9.93, 0.0299851}) , IMatrix ({{1., 0.866643}, {0.000450065, 0.498928}}) , 
			DVector({0.0021, 0.00105}));			// cube 3  				
		c6b[3] = HSet2D(IVector ({-3.56111, 0.0336361}) , IMatrix ({{1., 0.0148011}, {0.000745026, -0.99989}}) , 
			DVector({0.0043575, 0.000525}));			// cube 4						
		c6b[4] = HSet2D(IVector ({-6.45014, 0.031751}) , IMatrix ({{1., 0.999296}, {0.000574732, -0.0375247}}) , 
			DVector({0.00945, 0.013125}));			// cube 5				
		c6b[5] = HSet2D(IVector ({-10.0618, 0.029926}) , IMatrix ({{1., 0.887687}, {0.000446372, 0.460448}}) , 
			DVector({0.00084, 0.0011025}));			// cube 6
		
		///===================== choice of Procedure:  =====================
		int whichCase;
		
		do{
		cout << "=====================================================================" << endl;
		cout << "|| Choose the number of procedure You wish to perform:     	      " << endl;
		cout << "|| 1. Proof of the contracting grids's existence for a = 5.25	      " << endl;
		cout << "|| 2. Proof of the contracting grids's existence for a = 4.7   	  " << endl;
		cout << "|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 " << endl;
		cout << "|| 4. Proof of the contracting grids's existence for a = 4.381   	  " << endl;
		cout << "|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  " << endl;
		cout << "|| 6. Proof of the contracting grids's existence for a = 5.42   	  " << endl;
		cout << "|| other: exit	" << endl;
		cout << "===========================================================" << endl;
		cout << "Procedure " ;
		cin >> whichCase;
		cout << endl;
		
		switch (whichCase)
			{
				///===================== Procedure 1:  =====================
				case 1:
				{
					cout << "===========================================================" << endl;
					cout << "|| Roessler system, a = 5.25" << endl;
					cout << "===========================================================" << endl;
		
					cout << "Is the grid G3 forward-invariant? ... " << roessler525.inside(grid3,grid3,500,3) << endl;		// check if P(grid) < grid divided into 500x3 pieces
					cout << "----------------------------------------" << endl;				//~ {		
					cout << "P(C1)<C2? ... " << roessler525.inside(c3[0],c3[1],3,1) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C2)<C3? ... " << roessler525.inside(c3[1],c3[2],20,1) << endl;		// true
					cout << "----------------------------------------" << endl;
					cout << "P(C3)<C1? ... " << roessler525.inside(c3[2],c3[0],5,1) << endl;		// true
					cout << "----------------------------------------" <<  endl;

					break;
				}
				
				///===================== Procedure 2:  =====================
				case 2:
				{
					cout << "===========================================================" << endl;
					cout << "|| Roessler system, a = 4.7" << endl;
					cout << "===========================================================" << endl;	
					
					cout << "Is the grid G5 forward-invariant? ... " << roessler47.inside(grid5,grid5,450,1) << endl;		// check if P(grid) < grid divided horizontally into 450  pieces
					cout << "----------------------------------------" << endl;				//~ {		
					cout << "P(C1)<C2? ... " << roessler47.inside(c5[0],c5[1],1,1) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C2)<C3? ... " << roessler47.inside(c5[1],c5[2],20,2) << endl;		// true
					cout << "----------------------------------------" << endl;
					cout << "P(C3)<C4? ... " << roessler47.inside(c5[2],c5[3],5,2) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C4)<C5? ... " << roessler47.inside(c5[3],c5[4],1000,50) << endl;		// true
					cout << "----------------------------------------" << endl;
					cout << "P(C5)<C1? ... " << roessler47.inside(c5[4],c5[0],20,3) << endl;		// true
					cout << "----------------------------------------" <<  endl;
	
					break;
				}			
			
				///===================== Procedure 3:  =====================
				case 3:
				{
					cout << "===========================================================" << endl;
					cout << "|| Roessler system, a = 4.381" << endl;
					cout << "===========================================================" << endl;
					
					containsPeriodicPoint = roessler4381.stationaryPoint(pCloseTo6periodic,6);		// now it contains a 6-periodic point
					cout << "---------------------------------------------------" << endl;
					cout << "6-periodic orbit:" << endl;
					cout << "---------------------------------------------------" << endl;
					cout << "p1 = " << containsPeriodicPoint << endl;
					for(int i=0; i<6; ++i)															// prints all orbit
					{
						containsPeriodicPoint = roessler4381.P2d(containsPeriodicPoint);
						cout << "P^" << i+1 << "(p1) = " << containsPeriodicPoint << endl;
					};	
					cout << "----------------------------------------" << endl << endl;
					break;
				}	
			
				///===================== Procedure 4:  =====================
				case 4:
				{			
					cout << "===========================================================" << endl;
					cout << "|| Roessler system, a = 4.381" << endl;
					cout << "===========================================================" << endl;

					cout << "Is the grid G6 forward-invariant? ... " << roessler4381.inside(grid6,grid6,1300,1) << endl;		// check if P(grid) < grid divided horizontally into 1300 x 1  pieces
					cout << "----------------------------------------" << endl;					
					cout << "P(C1)<C2? ... " << roessler4381.inside(c6[0],c6[1],1,1) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C2)<C3? ... " << roessler4381.inside(c6[1],c6[2],300,300) << endl;	// true 300 x 300
					cout << "----------------------------------------" << endl;
					cout << "P(C3)<C4? ... " << roessler4381.inside(c6[2],c6[3],2,90) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C4)<C5? ... " << roessler4381.inside(c6[3],c6[4],3,7) << endl;		// true 
					cout << "----------------------------------------" << endl;
					cout << "P(C5)<C6? ... " << roessler4381.inside(c6[4],c6[5],1,8) << endl;		// true
					cout << "----------------------------------------" << endl;
					cout << "P(C6)<C1? ... " << roessler4381.inside(c6[5],c6[0],1,9) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					break;
				}
				
				///===================== Procedure 5:  =====================
			
				case 5:
				{		
					cout << "===========================================================" << endl;
					cout << "|| Roessler system, a = 5.42" << endl;
					cout << "===========================================================" << endl;
					
					containsPeriodicPoint = roessler542.stationaryPoint(pCloseTo6periodic2,6);		// now it contains a 6-periodic point
					cout << "---------------------------------------------------" << endl;
					cout << "6-periodic orbit:" << endl;
					cout << "---------------------------------------------------" << endl;
					cout << "p1 = " << containsPeriodicPoint << endl;
					for(int i=0; i<6; ++i)															// prints all orbit
					{
						containsPeriodicPoint = roessler542.P2d(containsPeriodicPoint);
						cout << "P^" << i+1 << "(p1) = " << containsPeriodicPoint << endl;
					};	
					cout << "----------------------------------------" << endl << endl;
					break;
				}
			
				///===================== Procedure 6:  =====================
				case 6:
				{		
					cout << "===========================================================" << endl;
					cout << "|| Roessler system, a = 5.42" << endl;
					cout << "===========================================================" << endl;

					cout << "Is the grid G6 forward-invariant? ... " << roessler542.inside(grid6b,grid6b,800,1) << endl;		// check if P(grid) < grid divided horizontally into 800 x 1  pieces
					cout << "----------------------------------------" << endl;						
					cout << "P(C1)<C2? ... " << roessler542.inside(c6b[0],c6b[1],1,1) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C2)<C3? ... " << roessler542.inside(c6b[1],c6b[2],4,1) << endl;		// true 
					cout << "----------------------------------------" << endl;
					cout << "P(C3)<C4? ... " << roessler542.inside(c6b[2],c6b[3],2,2) << endl;		// true
					cout << "----------------------------------------" <<  endl;
					cout << "P(C4)<C5? ... " << roessler542.inside(c6b[3],c6b[4],1,1) << endl;		// true 
					cout << "----------------------------------------" << endl;
					cout << "P(C5)<C6? ... " << roessler542.inside(c6b[4],c6b[5],60,50) << endl;	// true
					cout << "----------------------------------------" << endl;
					cout << "P(C6)<C1? ... " << roessler542.inside(c6b[5],c6b[0],3,7) << endl;		// true
					cout << "----------------------------------------" <<  endl;

					break;
				}
				
				default:
					cout << "Goodbye!" << endl;
					break;
				}
			}
			while (whichCase>=1 and whichCase <=6);
		
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} 

/* OUTPUT:
time ./Roessler_Sharkovskii 
=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 1

===========================================================
|| Roessler system, a = 5.25
===========================================================
Is the grid G3 forward-invariant? ... true
----------------------------------------
P(C1)<C2? ... true
----------------------------------------
P(C2)<C3? ... true
----------------------------------------
P(C3)<C1? ... true
----------------------------------------
=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 2

===========================================================
|| Roessler system, a = 4.7
===========================================================
Is the grid G5 forward-invariant? ... true
----------------------------------------
P(C1)<C2? ... true
----------------------------------------
P(C2)<C3? ... true
----------------------------------------
P(C3)<C4? ... true
----------------------------------------
P(C4)<C5? ... true
----------------------------------------
P(C5)<C1? ... true
----------------------------------------
=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 3

===========================================================
|| Roessler system, a = 4.381
===========================================================
A stationary point for P^6 in = {[-7.44826514033532, -7.448265140244187],[0.03638524011881493, 0.03638524011973746]} proven.
---------------------------------------------------
6-periodic orbit:
---------------------------------------------------
p1 = {[-7.44826514033532, -7.448265140244187],[0.03638524011881493, 0.03638524011973746]}
P^1(p1) = {[-5.432682771276081, -5.432682771080253],[0.03812100247833106, 0.03812100248150609]}
P^2(p1) = {[-8.146150765219835, -8.146150765118602],[0.03585533157361669, 0.03585533157606319]}
P^3(p1) = {[-4.052482471003891, -4.052482470816507],[0.03953831884313481, 0.03953831884723778]}
P^4(p1) = {[-6.988651597169091, -6.988651596889441],[0.03675237717289087, 0.0367523771731595]}
P^5(p1) = {[-6.385380925198882, -6.385380924637889],[0.03725846245305077, 0.0372584624617641]}
P^6(p1) = {[-7.448265140706449, -7.448265139873075],[0.03638524011216775, 0.03638524012638462]}
----------------------------------------

=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 4

===========================================================
|| Roessler system, a = 4.381
===========================================================
Is the grid G6 forward-invariant? ... true
----------------------------------------
P(C1)<C2? ... true
----------------------------------------
P(C2)<C3? ... true
----------------------------------------
P(C3)<C4? ... true
----------------------------------------
P(C4)<C5? ... true
----------------------------------------
P(C5)<C6? ... true
----------------------------------------
P(C6)<C1? ... true
----------------------------------------
=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 5

===========================================================
|| Roessler system, a = 5.42
===========================================================
A stationary point for P^6 in = {[-3.330388727960296, -3.33038872794934],[0.03381008102270888, 0.03381008102286536]} proven.
---------------------------------------------------
6-periodic orbit:
---------------------------------------------------
p1 = {[-3.330388727960296, -3.33038872794934],[0.03381008102270888, 0.03381008102286536]}
P^1(p1) = {[-6.043878148233535, -6.043878148213811],[0.03198830541026752, 0.03198830541028062]}
P^2(p1) = {[-9.93000468871182, -9.930004688693574],[0.02998512226572283, 0.02998512226583182]}
P^3(p1) = {[-3.561109751505439, -3.561109751469876],[0.03363611142511445, 0.03363611142566204]}
P^4(p1) = {[-6.450138010324274, -6.450138010261234],[0.03175097764903493, 0.03175097764907362]}
P^5(p1) = {[-10.06181179891221, -10.06181179888145],[0.02992604451798922, 0.02992604451853264]}
P^6(p1) = {[-3.330388727981797, -3.330388727927837],[0.03381008102231699, 0.03381008102325724]}
----------------------------------------

=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 6

===========================================================
|| Roessler system, a = 5.42
===========================================================
Is the grid G6 forward-invariant? ... true
----------------------------------------
P(C1)<C2? ... true
----------------------------------------
P(C2)<C3? ... true
----------------------------------------
P(C3)<C4? ... true
----------------------------------------
P(C4)<C5? ... true
----------------------------------------
P(C5)<C6? ... true
----------------------------------------
P(C6)<C1? ... true
----------------------------------------
=====================================================================
|| Choose the number of procedure You wish to perform:     	      
|| 1. Proof of the contracting grids's existence for a = 5.25	      
|| 2. Proof of the contracting grids's existence for a = 4.7   	  
|| 3. Proof of the 6-periodic orbit's existence by INM for a = 4.381 
|| 4. Proof of the contracting grids's existence for a = 4.381   	  
|| 5. Proof of the 6-periodic orbit's existence by INM for a = 5.42  
|| 6. Proof of the contracting grids's existence for a = 5.42   	  
|| other: exit	
===========================================================
Procedure 7

Goodbye!

real	5m36,248s
user	5m21,820s
sys	0m0,017s

*/
