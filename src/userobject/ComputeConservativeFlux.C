//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeConservativeFlux.h"
#include <iostream>
#include <string>
#include "FEProblem.h"
#include "MooseVariableFEBase.h"
#include "NonlinearSystemBase.h"
#include "Assembly.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh_base.h"

#include "libmesh/linear_implicit_system.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"


#include "libmesh/quadrature_gauss.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/quadrature.h"

using namespace std;

registerMooseObject("parrotApp", ComputeConservativeFlux);


template <>
InputParameters
validParams<ComputeConservativeFlux>()
{
    InputParameters params = validParams<GeneralUserObject>();
    
    params.addRequiredParam<std::vector<int>>("block_id","block_id");
      params.addRequiredParam<std::vector<Real>>("value_p","value_p");
    params.addRequiredParam<std::vector<AuxVariableName>>("aux_variable", "The auxiliary variable to store the transferred values in.");
    params.addRequiredParam<std::vector<boundary_id_type>>("boundary_N_bc","boundary_N_bc");
    params.addRequiredParam<std::vector<Real>>("value_N_bc", "The value of Neumann");
    
//    params.addParam<std::string>("output_file", "the file name of the output");
    params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
    params.addRequiredParam<bool>("conservative","use a conservative scheme?");
    
    return params;
}

ComputeConservativeFlux::ComputeConservativeFlux(const InputParameters & parameters) :
GeneralUserObject(parameters),
_aux_var_names(getParam<std::vector<AuxVariableName>>("aux_variable")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_boundary_N_ids(getParam<std::vector<boundary_id_type>>("boundary_N_bc")),
_value_N_bc(getParam<std::vector<Real>>("value_N_bc")),
_pp_comm(_fe_problem.es().get_mesh().comm()),
_stiff_matrix(_pp_comm),
//_has_output_file( isParamValid("output_file") ),
_hasMeshModifier( isParamValid("fractureMeshModifier") ),
_conservativeScheme( getParam<bool>("conservative") ),
_qrule(_assembly.qRule())
{
    
//    if (_has_output_file)
//        _output_filename=getParam<std::string>("output_file");
//
    if (_hasMeshModifier)
    {
        _meshModifierName=getParam<std::string>("fractureMeshModifier");
    }
    
    if (_aux_var_names.size() == 1)
        _aux_var_name = _aux_var_names.at(0);
    else
        paramError("variable", "You need to specify one and only one variable");
    
}



void ComputeConservativeFlux::initialize()
{
    _console<<"BEGIN initialize\n";
    
    EquationSystems & mooseEquationSystems=_fe_problem.es();
    
    NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    const DofMap & dof_map = _nl.dofMap();
    
    _stiff_matrix.attach_dof_map(dof_map);
    _stiff_matrix.init();
    
    PetscVector<Number> _rhs(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    
    _rhs.zero();
    
    AssembleDiffusionOP(_rhs);
    
    
    MooseVariableFEBase  & aux_var = _fe_problem.getVariable(0, _aux_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
    
   
    System & aux_sys = aux_var.sys().system();
    
    NumericVector<Number> * aux_solution = aux_sys.solution.get();
    
    aux_solution->print();
    
    PetscVector<Number> _ris(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _ris.zero();
    
    _stiff_matrix.vector_mult(_ris,*aux_solution);
    
    _ris.print_matlab("sol.m");
    
    _rhs.print_matlab("rhs.m");
    
    _console<<"END initialize\n";
}



void ComputeConservativeFlux::AssembleDiffusionOP(PetscVector<Number> &_rhs)
{
    _console << "BEGIN Assemble_Diffusion"  << std::endl;
    
    MeshModifier       const * _myMeshModifier_ptr;
    FractureUserObject const * _fractureUserObject_ptr;
    
    if (_hasMeshModifier)
    {
        MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
        _myMeshModifier_ptr=&_myMeshModifier;
        FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
        _fractureUserObject_ptr=&_fractureUserObject;
        
    }
    
    EquationSystems & _es =_fe_problem.es();
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _es.get_mesh();
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    // Get a reference to our system.
    NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _nl.system().get_dof_map().variable_type(0);
    
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    
    //QGauss qrule (dim, NINTH );
    //QGrid qrule (dim, NINTH  ); //TENTH
    std::unique_ptr<QBase> qrule( QBase::build (_qrule->type(),dim,_qrule->get_order()));
    
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (qrule.get());

    std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
    
    // Boundary integration requires another quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    // In 1D, the Clough and Gauss quadrature rules are identical.
    std::unique_ptr<QBase> qface(fe_type.default_quadrature_rule(dim-1));
    
    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule (qface.get());
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
    const std::vector<Point> & q_points = fe->get_xyz();
    
    const DofMap & dof_map = _nl.system().get_dof_map();
    
    
    
    std::vector<dof_id_type> dof_indices;
    
    DenseMatrix<Number> ke;
    DenseVector<Number> re;

    
    MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
    MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();
    
    std::vector<Number> permeability;
    
    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;
        fe->reinit (elem);
        
        dof_map.dof_indices(elem, dof_indices);
        const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());
        
        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to (n_dofs, phi.size());
        
        ke.resize (n_dofs , n_dofs);
        ke.zero();
        
        Real localPermeability=ComputeMaterialProprties(elem);
        permeability.resize( qrule->n_points());
        
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
            permeability.at(qp)=localPermeability;
            if(_hasMeshModifier)
            {
                if ( _fractureUserObject_ptr[0].isInside(q_points[qp]) )
                {
                    permeability.at(qp)=_vector_value.at(_vector_value.size()-1);
                }
            }
        }
        
        if(_conservativeScheme)
        {
            if (1)
            {
                Real ook=0.0;
                Real vol=0.0;
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                {
                    ook+=(JxW[qp]/permeability.at(qp));
                    vol+=JxW[qp];
                }
                Real conservativePermeability=vol/ook;
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                {
                    permeability.at(qp)=conservativePermeability;
                }
            }
            else
            {
                Real myk=0.0;
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                {
                    myk+=permeability.at(qp);
                }
                Real conservativePermeability=myk/qrule->n_points();
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                {
                    permeability.at(qp)=conservativePermeability;
                }
            }
        }
        
        
        for (unsigned int i=0; i<phi.size(); i++)
            for (unsigned int j=0; j<phi.size(); j++)
                for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                {
                    //std::cout<<permeability.at(qp)<<std::endl;
                    ke(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
                }
        
        re.resize(n_dofs);
        re.zero();
        
        
        for (auto side : elem->side_index_range())
        {
            if (elem->neighbor_ptr(side) == nullptr)
            {
                const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
                const std::vector<Real> & JxW_face = fe_face->get_JxW();
                
                fe_face->reinit(elem, side);
                
                for(int k =0; k < _boundary_N_ids.size(); k++){
                    
                    if (mesh.get_boundary_info().has_boundary_id (elem, side, _boundary_N_ids[k])) // Apply a traction on the right side
                    {
                        for (unsigned int qp=0; qp<qface->n_points(); qp++)
                            for (unsigned int i=0; i != n_dofs; i++)
                                re(i) += JxW_face[qp] * -1.0 * _value_N_bc.at(k) * phi_face[i][qp];
                    }
                }
            }
        }
        
        
        dof_map.constrain_element_matrix_and_vector(ke, re, dof_indices);
        
        _stiff_matrix.add_matrix (ke, dof_indices);
        _rhs.add_vector(re, dof_indices);
        
    }
    
    
    
    
    _stiff_matrix.close();
    _rhs.close();
    
    _console << "END Assemble_Diffusion"  << std::endl;
}


Real
ComputeConservativeFlux::ComputeMaterialProprties(const Elem *elem)
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




