//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeStatistics.h"

#include "FEProblem.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "FractureUserObject.h"
#include "Assembly.h"
// #include "MooseVariableFEBase.h"

#include "libmesh/quadrature_gauss.h"
#include "libmesh/petsc_vector.h"


registerMooseObject("parrotApp", ComputeStatistics);

// static void getRow(PetscMatrix<Number> & matrix, int const & row, std::vector<Real> & values, std::vector<int> & columns)
// {
//     Mat const & mat=matrix.mat();
//     PetscInt ncol;
//     PetscInt const *col;
//     PetscScalar const *val;
//     MatGetRow(mat,row,&ncol,&col,&val);
//     values.resize(ncol);
//     columns.resize(ncol);
//     for (int i=0; i<ncol; ++i)
//     {
//         values[i] =val[i];
//         columns[i]=col[i];
//     }
//     MatRestoreRow(mat,row,&ncol,&col,&val);
// }


template <>
InputParameters
validParams<ComputeStatistics>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  return params;
}

ComputeStatistics::ComputeStatistics(const InputParameters & parameters) :
GeneralUserObject(parameters),
_hasMeshModifier( isParamValid("fractureMeshModifier") )
{
	if (_hasMeshModifier)
		_meshModifierName=getParam<std::string>("fractureMeshModifier");
}

void ComputeStatistics::execute()
{
  _console << "ComputeStatistics::execute() begin"  << std::endl;

  MeshModifier       const * _myMeshModifier_ptr;
  FractureUserObject const * _fractureUserObject_ptr;

  if (_hasMeshModifier)
  {
    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    _myMeshModifier_ptr=&_myMeshModifier;
    FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
    _fractureUserObject_ptr=&_fractureUserObject;
  }

  DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();

  auto &comm = _fe_problem.es().get_mesh().comm();

  int const n  =dof_map.n_dofs();
  int const n_l=dof_map.n_local_dofs();

  PetscVector<Number> sol(comm,n,n_l);

  PetscMatrix<Number> mat(comm);
  mat.attach_dof_map(dof_map);
  mat.init();
  
  // Get a constant reference to the mesh object.
  MeshBase     const & mesh = _fe_problem.es().get_mesh();
  unsigned int const   dim  = mesh.mesh_dimension();

  // Get a reference to our system.
  TransientNonlinearImplicitSystem const & _system = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

  std::vector<dof_id_type> dof_indices;

  DenseVector<Number> Re;
  DenseMatrix<Number> Ke;
    
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();

  std::vector<int> elementCount(4,0.0);
  std::vector<int> howMany;

  for ( ; el != end_el; ++el)
  {
    const Elem * elem = *el;
    Point point=elem->centroid();
    howMany=_fractureUserObject_ptr->whichIsInside(point);
    int temp=howMany.size();
    elementCount.at(temp)=elementCount.at(temp)+1;

    dof_map.dof_indices(elem, dof_indices);

    int const loc_n=dof_indices.size();

    Re.resize(loc_n);
    Ke.resize(loc_n,loc_n);

    for (int i=0; i<loc_n; ++i)
    {
      Re(i)=1.0;
    }

    for (int i=0; i<loc_n; ++i)
    {
      for (int j=0; j<loc_n; ++j)
        Ke(i,j)=1.0;
    }

    mat.add_matrix(Ke,dof_indices);

    dof_map.constrain_element_vector(Re,dof_indices,true);

    for (int i=0; i<Re.size(); ++i)
    {
      if (Re(i)>0.5)
        Re(i)=0.0;
      else
        Re(i)=1.0;
    }

    sol.add_vector(Re,dof_indices);

  }
  sol.close();
  mat.close();

  std::vector<int> columns;
  std::vector<Number> values;
  int nnz=0;
  for (int i=0; i<mat.m(); ++i)
  {
    getRow(mat, i, values, columns);
    nnz+=values.size();
  }

  int hanging_nodes=0;
  
  std::vector<Number> v_local;
  sol.localize (v_local);

  for (int i=0; i<v_local.size(); ++i)
    if (v_local.at(i)>0.5)
      ++hanging_nodes;

//    mat.print_matlab();

  int regular_nodes=n-hanging_nodes;

  std::cout<<"Tot dofs="<<n<<std::endl;
  std::cout<<"Reg dofs="<<regular_nodes<<std::endl;
  std::cout<<"Han dofs="<<hanging_nodes<<std::endl;
  std::cout<<"nnz     ="<<nnz<<std::endl;
  std::cout<<"Background     el="<<elementCount.at(0)<<std::endl;
  std::cout<<"One   fracture el="<<elementCount.at(1)<<std::endl;
  std::cout<<"Two   fracture el="<<elementCount.at(2)<<std::endl;
  std::cout<<"Three fracture el="<<elementCount.at(3)<<std::endl;

  std::cout<<elementCount.at(3)<<","<<elementCount.at(2)<<","<<elementCount.at(1)<<","<<elementCount.at(0)<<","<<regular_nodes<<","<<nnz<<std::endl;

   _console << "ComputeStatistics::execute() end"  << std::endl;

};
