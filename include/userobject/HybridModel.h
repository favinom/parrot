//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
//#include "UserObject.h"

#include "GeneralUserObject.h"
#include "StoreUtopiaOperators.h"
#include "utopia.hpp"
#include "utopia_fe.hpp"
#include "libmesh/quadrature.h"
#include "ParrotProblem.h"

// Forward declarations
class HybridModel;

template <>
InputParameters validParams<HybridModel>();

/**
 * Base class for user-specific data

 */

class HybridModel : public GeneralUserObject
{

    
    
public:
    HybridModel(const InputParameters & params);

    virtual void execute() override;

    /**
    * Called before execute() is ever called so that data can be cleared.
    */
    virtual void initialize() override;

    /**
    * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where
    * you want to do MPI communication!
    */
    virtual void finalize() override {};

    /**
    * Must override.
    *
    * @param uo The UserObject to be joined into _this_ object.  Take the data from the uo object and
    * "add" it into the data for this object.
    */
    virtual void threadJoin(const UserObject &) override {};


    typedef utopia::USparseMatrix SparseMatT;

    typedef utopia::UVector VecT;


    bool solve();

    int  solve_pressure_monolithic(bool & ok);

    int solve_transport_monolithic(bool & ok);

    


  protected:

    std::string _dc_var_m;

    std::vector<int> _dc_boundary_id_m;

    std::vector<std::vector<int> > _dc_variables_id_m;

    std::vector<int> zero_rows;

    VariableName _f_var_name, _m_var_name;

    // const StoreOperators & _operator_storage;

    std::vector<int> _vector_p;

    std::vector<Real> _vector_value;

    bool _pressure, _transport;

    std::string _multiapp_name;

    std::shared_ptr<SparseMatT> T_ = NULL; 

    std::shared_ptr<SparseMatT> _A_t_store = NULL;

    utopia::UVector sol_f, sol_m;

    void CopyMatrixSolution(utopia::UVector to_sol);
    
    void CopyFractureSolution(utopia::UVector to_sol);

    void constraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat);

    void constraint_vec(utopia::UVector &boundary, utopia::UVector &vec);

    void boundary_conditions(libMesh::DofMap &dof_map, utopia::USparseMatrix &mat, utopia::UVector &vec, utopia::UVector &vec2);

    void set_zero_at_constraint_rows(DofMap &dof_map, utopia::USparseMatrix &mat);

    void find_boundary(std::vector<int> &zero_rows, 
                       std::vector<int> &_dc_boundary_id);

    void determine_dc_bnd_var_id(const std::vector<std::string> & BC_var);

    std::vector<std::string> split_string(const std::string & s, char delim);

    std::string _userobject_name_A_tot = "Matrix_A_tot";

    std::string _userobject_name_M_f   = "Matrix_M_f";

    std::string _userobject_name_M_m   = "Matrix_M_m";

    void assemble_poro_mass_matrix(FEProblemBase & problem, std::string & _userobject_name);

    void stabilize_matrix(utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix);

    Real ComputeMaterialProprties(const Elem *elem);

    std::unique_ptr<ExodusII_IO> _ex_writer;

    QBase const * const & _qrule;


    KSP ksp;
    SNESType ttype;
    KSPType type;

    KSP_PARROT * _ksp_ptr;

    utopia::UVector  rhs_m_c, rhs_f_c;



    //KSP_PARROT * _ksp_ptr;





};

