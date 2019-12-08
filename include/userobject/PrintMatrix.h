#ifndef PrintMatrix_H
#define PrintMatrix_H

#include "MultiAppTransfer.h"
#include <memory>
#include "UserObjectInterface.h"
#include "GeneralUserObject.h"

#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "StoreOperators.h"
class PrintMatrix
;

/*
 *                                  _.._
 *                                .' .-'`
 *                               /  /
 *                               |  |
 *                               \  '.___.;
 *                                '._  _.'
 *
 *                 ) \     /_(
 *                  )_`-)-_(_
 *                   `,' `__,)
 *                  _/   ((
 *         ______,-'    )
 *        (            ,
 *         \  )    (   |
 *        /| /`-----` /|
 *        \ \        / |
 *        |\|\      /| |\
 */

template <>
InputParameters validParams<PrintMatrix>();

/**
 * Project values from one domain to another
 */
class PrintMatrix : public GeneralUserObject

{
public:
    PrintMatrix(const InputParameters & parameters);
    
    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;
    void constraint_matrix(EquationSystems & _es, const std::string & system_name) ;

    Parallel::Communicator const & _pp_comm;
    
    std::string _dc_var;

    std::shared_ptr<PetscMatrix<Number>> jmat, smat;

    std::vector<int> _dc_boundary_id;

    std::vector<std::vector<int> > _dc_variables_id;

    std::vector<int> zero_rows;

    void stabilize_coeffiecient();


    UserObjectName const _userObjectName;
    StoreOperators * _storeOperatorsUO;



    void find_boundary(std::vector<int> &zero_rows, 
                       std::vector<int> &_dc_boundary_id);

    void determine_dc_bnd_var_id(const std::vector<std::string> & BC_var);

    void set_solution(PetscVector<Number> &correction, PetscVector<Number> &sol);

    std::vector<std::string> split_string(const std::string & s, char delim);

    
};

#endif /* FractureAppConforming */
