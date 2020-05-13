#include "libmesh/vector_value.h"

using namespace libMesh;

bool doesEdgeIntersectElement_2D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & hmin, RealVectorValue const & hmax);

Real CalcY_2D(Real const & xval, RealVectorValue const & p1, RealVectorValue const & p2);

Real CalcX_2D(Real const & yval, RealVectorValue const & p1, RealVectorValue const & p2);