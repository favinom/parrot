#ifndef AntidiffusiveFluxes3_H
#define AntidiffusiveFluxes3_H

#include "GeneralUserObject.h"
#include "StoreOperators.h"

#include "libmesh/petsc_vector.h"

class AntidiffusiveFluxes3;


template <>
InputParameters validParams<AntidiffusiveFluxes3>();

/**
 * Project values from one domain to another
 */
class AntidiffusiveFluxes3 : public GeneralUserObject

{
public:
    AntidiffusiveFluxes3(const InputParameters & parameters);
    
    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override;

    Parallel::Communicator const & _pp_comm;
    
    std::string _dc_var;

    std::shared_ptr<PetscMatrix<Number>> jmat, smat;

    std::vector<int> _dc_boundary_id;

    std::vector<std::vector<int> > _dc_variables_id;

    std::vector<int> zero_rows_fluxes;

    void stabilize_coeffiecient();


    UserObjectName const _userObjectName;
    StoreOperators * _storeOperatorsUO;

    bool _write_correction;

    
    void find_boundary(std::vector<int> &zero_rows_fluxes, 
                       std::vector<int> &_dc_boundary_id);

    void determine_dc_bnd_var_id(const std::vector<std::string> & BC_var);

    void set_solution(PetscVector<Number> &correction);

    std::vector<std::string> split_string(const std::string & s, char delim);

    
};

#endif /* FractureAppConforming */
