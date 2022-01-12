//////////////////////////////////////////////////////////////////////////////
///
///  @file utils.cpp
///  
///  @author ag  @date   Mar 5, 2020
//////////////////////////////////////////////////////////////////////////////

// Standard libraries and the CAPD library
#include <iostream>
#include "utils.h"
#include "capd/covrel/HSet2D.h"
#include "capd/dynsys/DynSysMap.h"
#include "capd/covrel/HSetMD.h"
#include "capd/intervals/lib.h"
#include "capd/capdlib.h"
//#include <stdexcept>
//#include <sstream>
//#include "capd/vectalg/iobject.hpp"
//#include "capd/dynset/C0DoubletonSet.h"
//#include "capd/geomset/CenteredDoubletonSet.hpp"
//#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;
using namespace std;

typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap;
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;

//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

// cut
// cuts the first coordinate of an interval 3-vector
IVector cut(IVector x)						
{
	if(x.dimension()==2) return x;
	IVector y(2);
	y[0]=x[1];
	y[1]=x[2];
	return y;
}

DVector cut(DVector x)						
{
	if(x.dimension()==2) return x;
	DVector y(2);
	y[0]=x[1];
	y[1]=x[2];
	return y;
}

// cuts an interval 3x3-matrix to 2x2
IMatrix cut(IMatrix M)						
{
	if(M.numberOfColumns()==2 and M.numberOfRows()==2) return M;
	IMatrix y(2,2);
	y[0][0]=M[1][1];
	y[0][1]=M[1][2];
	y[1][0]=M[2][1];
	y[1][1]=M[2][2];
	return y;
}

//expand
//prepends an interval 2-vector to 3-vector with zero 
IVector expand(IVector x)					
{	
	if(x.dimension()==3) return x;										
	IVector y({0.,0.,0.});
	y[1]=x[0];
	y[2]=x[1];
	return y;
}
DVector expand(DVector x)					
{	
	if(x.dimension()==3) return x;										
	DVector y({0.,0.,0.});
	y[1]=x[0];
	y[2]=x[1];
	return y;
}

//expands an interval 2x2 matrix to 3x3 with identity
IMatrix expand(IMatrix M)					
{		
	if(M.numberOfColumns()==3 and M.numberOfRows()==3) return M;									
	IMatrix y(3,3);
	y[1][1]=M[0][0];
	y[1][2]=M[0][1];
	y[2][1]=M[1][0];
	y[2][2]=M[1][1];
	y[0][1]=y[1][0]=y[2][0]=y[0][2]=interval(0.);
	y[0][0]=interval(1.);
	return y;
}

/////////////////////////////////////////////////////////////////////////
/// SecMap class
/////////////////////////////////////////////////////////////////////////

// image of the iteration iter, stores the derivative in DP
IVector SecMap::image(const IVector &x,IMatrix &DP, const int iter) const  
{
	IVector X({0.,0.,0.});
	X[1]=x[0];
	X[2]=x[1];
	
	C1Rect2Set Q(X);
	IMatrix DP3(3,3);
	
	IVector Y=(*P)(Q,DP3,iter);
	DP3=P->computeDP(Y,DP3,iter);
	
	IVector y(2);
	y[0]=Y[1];
	y[1]=Y[2];
	
	for(int i=0;i<2;++i)
	{
		for(int j=0;j<2;++j) DP[i][j]=DP3[i+1][j+1];
	}
	return y;
}

// image of the iteration iter only, no derivative calculation (faster)
IVector SecMap::image(const IVector &x, const int iter) const		
{
	IVector X({0.,0.,0.});
	X[1]=x[0];
	X[2]=x[1];
	
	C0Rect2Set Q(X);
	
	IVector Y=(*P)(Q,iter);
	
	IVector y(2);
	y[0]=Y[1];
	y[1]=Y[2];
	
	return y;
}

// derivative only
IMatrix SecMap::derivative(const IVector &x, const int iter) const		
{
	IMatrix DP(2,2);
	image(x,DP,iter);
	return DP;
}


//////////////////////////////////////////////////////////////////////////////
///	system3d structure methods
//////////////////////////////////////////////////////////////////////////////

// looks for a stationary point close to the starting point x0
IVector system3d::stationaryPoint(const IVector &x0, int iter)
{	  	
	IVector x = x0;
	IVector Px, xPrev, N(2);			
	IMatrix Id = IMatrix({{1.,0.},{0.,1.}});	// Identity 2x2 matrix	
	IMatrix DP(2,2);							
	int a=1;	
        try{
			do
			{
				x=midVector(x);					// Preliminary iteration to get a 'stable' point
				xPrev = x;
				Px=P2d(x,DP,iter);
				x=x-gauss(DP-Id,Px-x);
				++a;
			} while(not(subset(xPrev,x)) && a<50);

			/// if x is a good candidate, then we try a proper interval Newton method:
			if(a<50)
			{		
				const double delta = pow(10.0,-9);			// the size of the box X around x
				IVector X(2);
				X[0]=x[0]+delta*interval(-1.0,1.0);
				X[1]=x[1]+delta*interval(-1.0,1.0);	
				P2d(X,DP,iter);								// calculate DP = [dP^iter(X)]
				N = x-gauss(DP-Id,P2d(x,iter)-x);				// interval Newton
						
				if(subsetInterior(N,X))						// is N contained in the primary neighbourhood?;
				{
					cout << "A stationary point for P^" << iter << " in = " << N << " proven." << endl;
				}
			}
			else cout << "For the starting point " << x0 << " no success." << endl;
			} catch(exception& e){
              cout << "For the starting point " << x0 << " exception caught: " << e.what() << endl;
            }

return N;
}

// checks if the image of hset1 lies inside hset2
bool system3d::inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iteration) 
{
	GridSet gridSet(2);																		// stores grid
	hset1.gridSet(gridSet,howManyPiecesH,howManyPiecesV);		
	HSet2D hset2Base(IVector({0.,0.}),IMatrix::Identity(2),hset2.get_r());
  
	C0Rect2Set Set3d(expand(hset1.center()));
	IVector Pset2d(2);  
	interval rt(0.);
	IMatrix expandedGridSetCoordinateSystem = expand(gridSet.coordinateSystem());			// expanded to 3d, where P is defined
	IVector expandedGridSetBox = expand(gridSet.box());
	IMatrix expandedHSet2InvCoordinateSystem = expand(hset2.invCoordinateSystem());
	IVector expandedHSet2Center = expand(hset2.center());
	
	bool liesInside = 1;
	for(auto i = gridSet.begin();i!=gridSet.end();++i)
	{
		Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
		for(int j=0;j<iteration-1;++j)
			P(Set3d);
		Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));		// the image of the small box
		liesInside=liesInside && hset2Base.inside(Pset2d);									// do all lie inside?
		
		if(!liesInside) 
			{
				cout << "Does not lie inside, for this tiny box: " << Pset2d << " does not lie inside the box " << hset2Base << endl;
				return 0;
			}
	}
    return liesInside;
}
