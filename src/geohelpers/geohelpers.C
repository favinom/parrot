#include "geohelpers.h"

#define myeps 1e-15

bool doesEdgeIntersectElement_2D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & hmin, RealVectorValue const & hmax)
{

	Real const & left  =hmin(0);
	Real const & right =hmax(0);
	Real const & bottom=hmin(1);
	Real const & top   =hmax(1);

	Real max_p_x=std::max(p1(0),p2(0));
	Real min_p_x=std::min(p1(0),p2(0));

	Real max_p_y=std::max(p1(1),p2(1));
	Real min_p_y=std::min(p1(1),p2(1));

	if (right+myeps<min_p_x)
		return false;
	if (max_p_x+myeps<left)
		return false;
	if (top+myeps<min_p_y)
		return false;
	if (max_p_y+myeps<bottom)
		return false;

	Real xb=CalcX_2D(bottom,p1,p2);
	Real xt=CalcX_2D(top   ,p1,p2);
	Real yl=CalcY_2D(left  ,p1,p2);
	Real yr=CalcY_2D(right ,p1,p2);

	//int intersections=0;

	if( left-myeps<xb && xb<right+myeps )
	{
		//intersections=intersections+1;
		return true;
	}
	if( left-myeps<xt && xt<right+myeps )
	{
		//intersections=intersections+1;
		return true;
	}
	if( bottom-myeps<yl && yl<top+myeps )
	{
		//intersections=intersections+1;
		return true;
	}
	if( bottom-myeps<yr && yr<top+myeps )
	{
		//intersections=intersections+1;
		return true;
	}

	if ( hmin(0)-myeps<p1(0) && p1(0)<hmax(0)+myeps && hmin(1)-myeps<p1(1) && p1(1)<hmax(1)+myeps )
		return true;
	if ( hmin(0)-myeps<p2(0) && p2(0)<hmax(0)+myeps && hmin(1)-myeps<p2(1) && p2(1)<hmax(1)+myeps )
		return true;

	return false;
}

Real CalcY_2D(Real const & xval, RealVectorValue const & p1, RealVectorValue const & p2)
{
    Real x0=p1(0);
    Real x1=p2(0);
    Real y0=p1(1);
    Real y1=p2(1);
    if ( std::fabs(x0-x1)<myeps )
    	return NAN;

    return y0 + (xval - x0)*(y1 - y0)/(x1 - x0);
}

Real CalcX_2D(Real const & yval, RealVectorValue const & p1, RealVectorValue const & p2)
{
	Real x0=p1(0);
    Real x1=p2(0);
    Real y0=p1(1);
    Real y1=p2(1);

    if ( std::fabs(y0-y1)<myeps )
    	return NAN;

    return x0 + (yval - y0)*(x1 - x0)/(y1 - y0);
}


