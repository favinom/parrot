//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParrotProblem3.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "petscmat.h"
#include "petscmat.h"
#include "chrono"



registerMooseObject("parrotApp", ParrotProblem3);

template <>
InputParameters
validParams<ParrotProblem3>()
{
    InputParameters params = validParams<FEProblem>();
    params.addRequiredParam<bool>("use_AFC","use_AlgFluxCorr");
    params.addParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    params.addParam<UserObjectName>("antidiffusive_fluxes","The userobject that computes antidiffusive_fluxes");
    params.addRequiredParam<int>("solver_type","solver_type");
    return params;
}


ParrotProblem3::ParrotProblem3(const InputParameters & parameters) :
FEProblem(parameters),
_pp_comm(_mesh.comm()),
_equationSystems(this->es()),
_nl_libMesh( _equationSystems.get_system<TransientNonlinearImplicitSystem>("nl0") ),
_mat_SM(*_nl_libMesh.matrix),
_rhs_NV(*_nl_libMesh.rhs),
_sol_NV(_nl->solution()),
_stab_matrix(_pp_comm),
_res_m(_pp_comm),
_use_afc(getParam<bool>("use_AFC")),
_solverType(getParam<int>("solver_type")),
_regularNodes(_pp_comm),
_ones(_pp_comm),
_scale_interpolator(_pp_comm)
{
 
    if(parameters.isParamValid("operator_userobject"))
    {
        _hasStoreOperatorsUO=true;
        userObjectName=new UserObjectName(getParam<UserObjectName>("operator_userobject"));
    }

    else
    {
        _hasStoreOperatorsUO=false;
        _storeOperatorsUO=NULL;
    }

    if(parameters.isParamValid("antidiffusive_fluxes"))
    {
        _ComputeAntidiffusiveFluxes=true;
        userObjectNameFluxes=new UserObjectName(getParam<UserObjectName>("antidiffusive_fluxes"));
    }
    else
    {
          _ComputeAntidiffusiveFluxes=false;
    }

    _const_jacobian=true;
    _is_stab_matrix_assembled=false;
    _is_jac_matrix_assembled=false;

    if ( 0<_solverType && _solverType <4)
        _solve=false;

    parrotSolver=new ParrotSolver(_solverType,_pp_comm);

}

void ParrotProblem3::initialSetup()
{
    FEProblem::initialSetup();
    
    if (_hasStoreOperatorsUO)
    {
        _storeOperatorsUO=&getUserObject<StoreOperators>(userObjectName[0]);
        PetscMatrix<Number> * & localmatrix=_storeOperatorsUO[0].StabMatrix();
        localmatrix=&_stab_matrix;

        DofMap const & dof_map = _nl->dofMap();
        _res_m.init(dof_map.n_dofs(), dof_map.n_local_dofs());


                      _ones.init(dof_map.n_dofs(), dof_map.n_local_dofs());
              _regularNodes.init(dof_map.n_dofs(), dof_map.n_local_dofs());
        _scale_interpolator.init(dof_map.n_dofs(), dof_map.n_local_dofs());


        _poro_lumped = _storeOperatorsUO->PoroLumpMassMatrix();

        _interpolator          = _storeOperatorsUO->Interpolator();

        _ones=1;
        _interpolator->vector_mult_add(_scale_interpolator,_ones);
        _scale_interpolator.reciprocal();

        //_interpolator->print_matlab("int.m");
        //_scale_interpolator.print_matlab("scale.m");

        _dirichlet_bc = _storeOperatorsUO->BcVec();
        _value_dirichlet_bc = _storeOperatorsUO->ValueBcVec();

        _value_dirichlet_bc->close();

        _ones=1;
        _ones.close();

        _poro_lumped->vector_mult_add(_regularNodes,_ones);

        for (int i=_regularNodes.first_local_index(); i<_regularNodes.last_local_index(); ++i)
        {
            Real v=_regularNodes(i);
            if (std::fabs(v)<1e-10)
            {
                _regularNodes.set(i,0.0);
            }
            else
            {
                _regularNodes.set(i,1.0);
            }
        }
        _regularNodes.close();
        //_regularNodes.print_matlab("regular.m");
    }


    const DofMap & dof_map = _nl->dofMap(); 
    // _sol_vec                = _storeOperatorsUO->SolVec();
    // _sol_vec->init(dof_map.n_dofs(), dof_map.n_local_dofs());


    if (_ComputeAntidiffusiveFluxes)
    {

        _ComputeAF=&getUserObject<AntidiffusiveFluxes>(userObjectNameFluxes[0]);

        _JMatrix  = _storeOperatorsUO[0].JacMatrix();

        _JMatrix->attach_dof_map(dof_map);

        _JMatrix->init();

        _JMatrix->zero();

        _JMatrix->close();
    }    
};


void
ParrotProblem3::solve()
{
    if ( 0<_solverType && _solverType <4)
    {
            SparseMatrix <Number> * mat_SM = _nl_libMesh.matrix;
            NumericVector<Number> * rhs_NV = _nl_libMesh.rhs;

            std::unique_ptr<NumericVector<Number>> & solutionPointer(_nl_libMesh.solution);
            NumericVector<Number> * sol_NV=solutionPointer.get();

        if (!_is_jac_matrix_assembled)
        {
            std::cout<<"\n\nAssembling Jacobian and stabilization matrix...";
            auto t_start = std::chrono::high_resolution_clock::now();
            computeJacobianSys(_nl_libMesh, *_nl->currentSolution(), _mat_SM);
            auto t_stop = std::chrono::high_resolution_clock::now();
            auto diff=std::chrono::duration<double, std::milli>(t_stop-t_start).count();
            std::cout<<" Done!\n It took ";
            std::cout<<diff<<" ms\n";

            _is_jac_matrix_assembled=true;


            parrotSolver->setMatrixAndVectors(mat_SM,rhs_NV,sol_NV);
            parrotSolver->setConstantMatrix(_const_jacobian);

        }
        NumericVector<Number> & solOld=_nl->solutionOld() ;
        //computeResidualSys(_nl_libMesh, *_nl->currentSolution(),_rhs_NV);
        computeResidualSys(_nl_libMesh, solOld ,_rhs_NV);
        parrotSolver->solve();

        //_interpolator->vector_mult_add(_ones,_sol_NV);
        //_sol_NV.pointwise_mult(_ones,_scale_interpolator);

    }
    else
    {
        FEProblemBase::solve();
    }

}


void
ParrotProblem3::computeResidualSys(NonlinearImplicitSystem & /*sys*/,
                                   const NumericVector<Number> & soln,
                                   NumericVector<Number> & residual)
{

    _res_m.zero();
    _poro_lumped->vector_mult_add(_res_m,soln);
    Real inv_dt = 1.0/this->dt();
    _res_m.scale(inv_dt);

    _ones.pointwise_mult(_res_m,*_dirichlet_bc);
    residual.pointwise_mult(_res_m,_regularNodes);
    residual.add(*_value_dirichlet_bc);

    // _res_m.zero();
    // _poro_lumped->vector_mult_add(_res_m,soln);
    // Real inv_dt = 1.0/this->dt();
    // _res_m.scale(inv_dt);
    // residual.pointwise_mult(_res_m,*_dirichlet_bc);
    // 

    // residual.add(*_value_dirichlet_bc);
}

void
ParrotProblem3::computeJacobianSys(NonlinearImplicitSystem & /*sys*/,
                                  const NumericVector<Number> & soln,
                                  SparseMatrix<Number> & jacobian)
{
    FEProblemBase::computeJacobian(soln, jacobian);
    
    if (_use_afc==true && _is_stab_matrix_assembled==false)
    {
        computeStabilizationMatrix(jacobian);
        jacobian.add(1.0,_stab_matrix);

        if (_ComputeAntidiffusiveFluxes){
             _poro_lumped = _storeOperatorsUO->PoroLumpMassMatrix();
            _JMatrix->add(1.0, jacobian);
            _JMatrix->add(-1.0, *_poro_lumped);
        }

    }
}

void
ParrotProblem3::computeStabilizationMatrix(SparseMatrix<Number> & jacobian)
{
    PetscMatrix<Number> & jac_PM=dynamic_cast<PetscMatrix<Number> &> (jacobian);
   // Declare a PetscMatrix that will contain the transpose
    PetscMatrix<Number> jac_tr_PM(_pp_comm);
    // Transpose of the Jacobian
    jacobian.get_transpose (jac_tr_PM);
        
    int rb=jac_PM.row_start();
    int re=jac_PM.row_stop();

    DofMap const & dof_map = this->getNonlinearSystemBase().dofMap();

    _stab_matrix.attach_dof_map(dof_map);

    _stab_matrix.init();//m,n,m_l,n_l,30);

    _stab_matrix.zero();
    //MatSetOption(_stab_matrix_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    
    std::vector<Real> values;
    std::vector<Real> values_tr;
    std::vector<int> columns;
    std::vector<int> columns_tr;

    for (int row=rb; row<re; ++row)
    {
        getRow(jac_PM,row,values,columns);
        getRow(jac_tr_PM,row,values_tr,columns_tr);
        
        if (columns.size()!=columns_tr.size())
        {
            std::cout<<"ncols!=ncols_tr\n";
            exit(1);
        }
        
        for (int p=0; p<columns.size(); ++p)
        {
            if (columns[p]!=columns_tr[p])
            {
                std::cout<<"cols[p]!=cols_tr[p]"<<std::endl;
                exit(1);
            }
            
            int col=columns[p];
            if (row!=col)
            {
                // we are in a extradiagonal
                Real Aij=values[p];
                Real Aji=values_tr[p];
                
                Real maxEntry=std::max(Aij,Aji);
                
                if (maxEntry>-1e-15)
                {
                    Real value=-1.0*maxEntry;

                    _stab_matrix.set(row,col,value);
                    _stab_matrix.add(row,row,maxEntry);
                }
            }
        }
    }
     _stab_matrix.close();

    _is_stab_matrix_assembled=true;
}

// void
// ParrotProblem3::update_sol()
// {
//     std::cout<<"ParrotProblem3::update_sol BEGIN"<<std::endl;

//     MooseVariableFEBase  & main_var = this->getVariable(0, "CM", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

//     System & main_sys = main_var.sys().system();

//     NumericVector<Number> * main_solution = main_sys.solution.get();


//     _sol_vec                = _storeOperatorsUO->SolVec();

//     //_sol_vec->print_matlab("sol_moose.m");

//     { 
        
//         for (const auto & node : this->es().get_mesh().local_node_ptr_range())

//         {
//             for (unsigned int comp = 0; comp < node->n_comp(main_sys.number(), main_var.number()); comp++)

//             {

//                 //const dof_id_type proj_index = node->dof_number(_nl.number(), sol_var.number(), comp);

//                 const dof_id_type to_index = node->dof_number(main_sys.number(), main_var.number(), comp);

//                 //const dof_id_type to_index_c = node->dof_number(aux_sys.number(), aux_var_c.number(), comp);


//                 main_solution->set(to_index, (*_sol_vec)(to_index));

//                 //aux_solution->set(to_index_c, correction(proj_index));
//             }

//         }
//     }



//     main_solution->close();
//     //aux_solution->close();
//     main_sys.update();

//     std::cout<<"ParrotProblem3::update_sol END"<<std::endl;

// }
