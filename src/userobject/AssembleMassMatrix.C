//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleMassMatrix.h"

#include "FEProblem.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "FractureUserObject.h"
// #include "MooseVariableFEBase.h"

#include "libmesh/quadrature_gauss.h"
// #include "libmesh/exodusII_io.h"
// #include "libmesh/nonlinear_implicit_system.h"
// #include "libmesh/petsc_matrix.h"
// #include "libmesh/petsc_vector.h"
// #include "libmesh/sparse_matrix.h"
// #include "libmesh/equation_systems.h"
// #include "libmesh/linear_implicit_system.h"
// #include "libmesh/transient_system.h"
// #include "libmesh/dirichlet_boundaries.h"
// #include "libmesh/zero_function.h"
// #include "libmesh/const_function.h"
// #include "libmesh/parsed_function.h"
// #include "libmesh/petsc_matrix.h"
// #include "libmesh/petsc_vector.h"

registerMooseObject("parrotApp", AssembleMassMatrix);


template <>
InputParameters
validParams<AssembleMassMatrix>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<std::vector<int>>("block_id","block_id");
  params.addRequiredParam<std::vector<Real>>("value_p","value_p");
  params.addRequiredParam<bool>("constrain_matrix","constrain_matrix");
  params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  return params;
}

AssembleMassMatrix::AssembleMassMatrix(const InputParameters & parameters) :
GeneralUserObject(parameters),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
userObjectName(getParam<UserObjectName>("operator_userobject")),
_constrainMatrices(getParam<bool>("constrain_matrix")),
_code_dof_map(true),
_hasMeshModifier( isParamValid("fractureMeshModifier") )
{
	if (_hasMeshModifier)
		_meshModifierName=getParam<std::string>("fractureMeshModifier");
}

void AssembleMassMatrix::execute()
{
	assemble_mass_matrix();
};

void AssembleMassMatrix::assemble_mass_matrix(){

   _console << "Assemble_Mass_matrix() begin "  << std::endl;

  MeshModifier       const * _myMeshModifier_ptr;
  FractureUserObject const * _fractureUserObject_ptr;

  if (_hasMeshModifier)
  {
  	MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
  	_myMeshModifier_ptr=&_myMeshModifier;
  	FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
  	_fractureUserObject_ptr=&_fractureUserObject;
  }


   //StoreOperators const & storeOperatorsUO(getUserObject<StoreOperators>(userObjectName));
   //_mass_matrix           = const_cast<StoreOperators&>(storeOperatorsUO).MassMatrix();
   //_poro_mass_matrix      = const_cast<StoreOperators&>(storeOperatorsUO).PoroMassMatrix();
   //_lump_mass_matrix      = const_cast<StoreOperators&>(storeOperatorsUO).LumpMassMatrix();
   //_poro_lump_mass_matrix = const_cast<StoreOperators&>(storeOperatorsUO).PoroLumpMassMatrix();

   DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();

   StoreOperators & storeOperatorsUO=(_fe_problem.getUserObjectTempl<StoreOperators>(userObjectName));
   _interpolator          = storeOperatorsUO.Interpolator();
   _mass_matrix           = storeOperatorsUO.MassMatrix();
   _poro_mass_matrix      = storeOperatorsUO.PoroMassMatrix();
   _lump_mass_matrix      = storeOperatorsUO.LumpMassMatrix();
   _poro_lump_mass_matrix = storeOperatorsUO.PoroLumpMassMatrix();

   if (_code_dof_map)
   {
   	_interpolator->attach_dof_map(dof_map);
   	_mass_matrix->attach_dof_map(dof_map);
   	_poro_mass_matrix->attach_dof_map(dof_map);
   	_lump_mass_matrix->attach_dof_map(dof_map);
   	_poro_lump_mass_matrix->attach_dof_map(dof_map);
   	_interpolator->init();
   	_mass_matrix->init();
   	_poro_mass_matrix->init();
   	_lump_mass_matrix->init();
   	_poro_lump_mass_matrix->init();
   }
   else
   {
   	int m=dof_map.n_dofs();
   	int n=dof_map.n_dofs();
   	int m_l=dof_map.n_local_dofs();
   	int n_l=dof_map.n_local_dofs();

   	_interpolator->init(m,n,m_l,n_l);
   	_mass_matrix->init(m,n,m_l,n_l);
   	_poro_mass_matrix->init(m,n,m_l,n_l);
   	_lump_mass_matrix->init(m,n,m_l,n_l);
   	_poro_lump_mass_matrix->init(m,n,m_l,n_l);

   }


   // Get a constant reference to the mesh object.
   MeshBase     const & mesh = _fe_problem.es().get_mesh();
   unsigned int const   dim  = mesh.mesh_dimension();

    // Get a reference to our system.
    TransientNonlinearImplicitSystem const & _system = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType const & fe_type = _system.get_dof_map().variable_type(0);
	//FEType fe_type = system.variable_type(0);
	UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    //QGauss qrule (dim, fe_type.default_quadrature_order());
	QGauss qrule (dim, TENTH);
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    const std::vector<Real>& JxW      = fe->get_JxW();    
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    //const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
    const std::vector<Point>& q_points = fe->get_xyz();
    //const DofMap& dof_map = system.get_dof_map();

    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_l;
    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_l_p;
    std::vector<dof_id_type> dof_indices_i;

    DenseMatrix<Number> Me;
    DenseMatrix<Number> Me_p;
    DenseMatrix<Number> Me_l;
    DenseMatrix<Number> Me_l_p;
    DenseMatrix<Number> Me_i;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
    	const Elem * elem = *el;
    	fe->reinit (elem);
    	dof_map.dof_indices(elem, dof_indices);
    	dof_map.dof_indices(elem, dof_indices_l);
    	dof_map.dof_indices(elem, dof_indices_p);
    	dof_map.dof_indices(elem, dof_indices_l_p);
    	dof_map.dof_indices(elem, dof_indices_i);

    	int const loc_n=dof_indices.size();

    	Me.resize(loc_n,loc_n);
    	Me_p.resize(loc_n,loc_n);
    	Me_l.resize(loc_n,loc_n);
    	Me_l_p.resize(loc_n,loc_n);
    	Me_i.resize(loc_n,loc_n);
    	Me.zero();
    	Me_p.zero();
    	Me_l.zero();
    	Me_l_p.zero();
    	Me_i.zero();

    	for (unsigned int i=0; i<phi.size(); i++)
    	{
    		for (unsigned int j=0; j<phi.size(); j++)
    		{
    			for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    			{
    				Real poro=ComputeMaterialProprties(elem);
    				if(_hasMeshModifier)
    				{
    					if ( _fractureUserObject_ptr[0].isInside(q_points[qp]) )
    					{
    						poro=_vector_value.at(_vector_value.size()-1);
    					}
    				}


    				Me  (i,j)   +=  JxW[qp] * phi[i][qp] * phi[j][qp];
    				Me_p(i,j)   += poro * JxW[qp] * phi[i][qp] * phi[j][qp];
    				Me_l(i,i)   += JxW[qp] * phi[i][qp] * phi[j][qp];
    				Me_l_p(i,i) += poro * JxW[qp] * phi[i][qp] * phi[j][qp];

    			}
    		}
    	}

    	if (_constrainMatrices)
    	{
    		dof_map.constrain_element_matrix(Me,dof_indices,true);
    		dof_map.constrain_element_matrix(Me_p,dof_indices_p,true);
    		dof_map.constrain_element_matrix(Me_l,dof_indices_l,true);
    		dof_map.constrain_element_matrix(Me_l_p,dof_indices_l_p,true);
    	}
    	{
    		dof_map.constrain_element_matrix(Me_i,dof_indices_i,true);
    		for (int i=0; i<Me_i.m(); ++i)
    		{
    			if (Me_i(i,i)<0.5)
    			{
    				Me_i(i,i)=1;
    			}
    			else
    			{
    				for (int j=0; j<Me_i.n(); ++j)
    				{
    					Me_i(i,j)=-1.0*Me_i(i,j);
    				}
    				Me_i(i,i)=0.0;
    			}
    		}
    	}

    	(*_mass_matrix).add_matrix(Me, dof_indices);
    	(*_poro_mass_matrix).add_matrix (Me_p, dof_indices_p);
    	(*_lump_mass_matrix).add_matrix (Me_l, dof_indices_l);
    	(*_poro_lump_mass_matrix).add_matrix (Me_l_p, dof_indices_l_p);
    	(*_interpolator).add_matrix (Me_i, dof_indices_i);
   }

   if (!_code_dof_map)
   {
   	MatSetOption(_mass_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   	MatSetOption(_poro_mass_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   	MatSetOption(_poro_lump_mass_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   	MatSetOption(_lump_mass_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   	MatSetOption(_interpolator->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   }

   (*_mass_matrix).close();
   (*_poro_mass_matrix).close();
   (*_poro_lump_mass_matrix).close();
   (*_lump_mass_matrix).close();
   (*_interpolator).close();
   
    _console << "Assemble_Mass_matrix() end "  << std::endl;

}

Real
AssembleMassMatrix::ComputeMaterialProprties(const Elem *elem)
{
	Real permeability=0.0;
	for(int ll=0; ll<_vector_p.size(); ll++)
	{
		if (elem->subdomain_id()==_vector_p[ll])
		{
			permeability = _vector_value[ll];
		}
	}
	return permeability;
}



