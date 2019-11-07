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
#include "libmesh/quadrature_gauss.h"
#include "FEProblem.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/parsed_function.h"
#include <iostream>
#include "MooseVariableFEBase.h"
#include "NonlinearSystemBase.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include <string>
#include "StoreOperators.h"
using namespace std;

registerMooseObject("parrotApp", AssembleMassMatrix);

// void
// assembly_diffusion(EquationSystems & es, const std::string & system_name)
// {
//     std::cout<<"CIAO\n";
//     SolveDiffusion * pippo = es.parameters.get< SolveDiffusion*>("pippo");
//     pippo->AssembleDiffusionOP(es, system_name);
// }


template <>
InputParameters
validParams<AssembleMassMatrix>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<std::vector<int>>("block_id",
                                                     "The name of the nodeset to create");
  params.addRequiredParam<std::vector<Real>>("value_p",
                                                     "The name of the nodeset to create");
  params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
  return params;
}

AssembleMassMatrix::AssembleMassMatrix(const InputParameters & parameters) :
GeneralUserObject(parameters),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_operator_storage(getUserObject<StoreOperators>("operator_userobject"))

{
    
    // auto &comm = _fe_problem.es().get_mesh().comm();

    // auto mat_1 = std::make_shared<PetscMatrix<Number>>(comm);

    // auto mat_2 = std::make_shared<PetscMatrix<Number>>(comm);

    // auto mat_3 = std::make_shared<PetscMatrix<Number>>(comm);

}

 void AssembleMassMatrix::execute(){

  assemble_mass_matrix();
};


 void AssembleMassMatrix::initialize()
{
}

 void AssembleMassMatrix::finalize()
{
}

 void AssembleMassMatrix::threadJoin(const UserObject & uo)
{
    
};


void AssembleMassMatrix::assemble_mass_matrix(){

   _console << "Assemble_Mass_matrix() begin "  << std::endl;

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _fe_problem.es().get_mesh();

    const DofMap & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();

    int m=dof_map.n_dofs();
    
    int n=dof_map.n_dofs();
    
    int m_l=dof_map.n_local_dofs();
    
    int n_l=dof_map.n_local_dofs();
    
    int nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());


    _mass_matrix      = const_cast<StoreOperators&>(_operator_storage).MassMatrix();

    // std::cout << "[MASSMATRIX-AssembleMassMatrix]" << _mass_matrix.get() << std::endl;

    _poro_mass_matrix = const_cast<StoreOperators&>(_operator_storage).PoroMassMatrix();


    _lump_mass_matrix = const_cast<StoreOperators&>(_operator_storage).LumpMassMatrix();

    _mass_matrix->init(m,n,m_l,n_l,nnz_x_row);

    _poro_mass_matrix->init(m,n,m_l,n_l,nnz_x_row);

    _lump_mass_matrix->init(m,n,m_l,n_l,nnz_x_row);

    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to our system.
    auto & _system = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);

    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

    DenseMatrix<Number> Me;

    DenseMatrix<Number> Me_p;

    DenseMatrix<Number> Me_l;

    QGauss qrule (dim, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // first we need to manually zero the matrix
    _mass_matrix->zero();

    _poro_mass_matrix->zero();

    _lump_mass_matrix->zero();

    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;

        Elem * ele = *el;

        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);

        Me.resize (dof_indices.size(), dof_indices.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                    Me(i,j)   += JxW[qp] * phi[i][qp] * phi[j][qp];

                 }
             }
         }


         dof_map.constrain_element_matrix(Me,dof_indices,true);

        (* _mass_matrix).add_matrix (Me, dof_indices);


        dof_map.dof_indices(elem, dof_indices);

        Me_p.resize (dof_indices.size(), dof_indices.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            for (unsigned int i=0; i<phi.size(); i++)
            {

                for (unsigned int j=0; j<phi.size(); j++){

                    Me_p(i,j) += ComputeMaterialProprties(elem) * JxW[qp] * phi[i][qp] * phi[j][qp];

                }
            }
        }

        dof_map.constrain_element_matrix(Me_p,dof_indices,true);

        (*_poro_mass_matrix).add_matrix (Me_p, dof_indices);



        dof_map.dof_indices(elem, dof_indices);

        Me_l.resize (dof_indices.size(), dof_indices.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            for (unsigned int i=0; i<phi.size(); i++)
            {

                for (unsigned int j=0; j<phi.size(); j++){

                    Me_l(i,j) += ComputeMaterialProprties(elem) * JxW[qp] * phi[i][qp] * phi[j][qp];

                }
            }
        }



        for (int i=0; i<Me_l.m(); ++i)
        {
            for (int j=0; j<Me_l.n(); ++j)
            {
                if (i!=j)
                {
                    Me_l(i,i)+=Me_l(i,j);
                    Me_l(i,j)=0.0;
                }
            }
        }
                       


        dof_map.constrain_element_matrix(Me_l,dof_indices,true);

        (*_lump_mass_matrix).add_matrix (Me_l, dof_indices);
    }


    (*_mass_matrix).close();
   
    // _mass_matrix->print_matlab("final.m");

    (*_poro_mass_matrix).close();

    (*_lump_mass_matrix).close();

   
    _console << "Assemble_Mass_matrix() end "  << std::endl;


}

Real
AssembleMassMatrix::ComputeMaterialProprties(const Elem *elem){
    
   // _console << "_vector_p.size()"  << _vector_p.size() <<std::endl;

Real permeability=0.0;

    for(int ll=0; ll<_vector_p.size(); ll++){
        if (elem->subdomain_id()==_vector_p[ll]) {

            permeability = _vector_value[ll];
        }
    }
  
    return permeability;
}



