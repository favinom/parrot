/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*    Immersed_Boundary- ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/


#include "AntidiffusiveFluxes.h"

#include <algorithm>    // std::max

#include "NonlinearSystemBase.h"
#include "Transient.h"




registerMooseObject("parrotApp",AntidiffusiveFluxes);

template <>
InputParameters
validParams<AntidiffusiveFluxes>()
{
    
    InputParameters params = validParams<GeneralUserObject>();
    params.addParam<std::string>("dc_boundaries", "-1", "Dirichlet Boundary ID");
    params.addParam<std::string>("dc_variables" , "-1", "Variable to which given BC_id applies");
    params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    params.addRequiredParam<bool>("WriteCorrection","boolean to write correction");

    return params;
    
}

static void getRow(PetscMatrix<Number> & matrix, int const & row, std::vector<Real> & values, std::vector<int> & columns)
{
    Mat const & mat=matrix.mat();
    PetscInt ncol;
    PetscInt const *col;
    PetscScalar const *val;
    MatGetRow(mat,row,&ncol,&col,&val);
    values.resize(ncol);
    columns.resize(ncol);
    for (int i=0; i<ncol; ++i)
    {
        values[i] =val[i];
        columns[i]=col[i];
    }
    MatRestoreRow(mat,row,&ncol,&col,&val);
}


AntidiffusiveFluxes::AntidiffusiveFluxes(const InputParameters & parameters):
GeneralUserObject(parameters),
_pp_comm(_fe_problem.es().get_mesh().comm()),
_dc_var(getParam<std::string>("dc_variables")),
_userObjectName(getParam<UserObjectName>("operator_userobject")),
_write_correction(getParam<bool>("WriteCorrection"))
{
 
    std::vector<std::string> tmp = split_string(parameters.get<std::string>("dc_boundaries"), ' ');
    for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
    {
        _dc_boundary_id.push_back(atoi(str_tmp->c_str()));
    }

   

}



void
AntidiffusiveFluxes::initialize()
{
  
}
  

void
AntidiffusiveFluxes::execute()
{
  stabilize_coeffiecient();
}


void
AntidiffusiveFluxes::finalize()
{
}

void
AntidiffusiveFluxes::stabilize_coeffiecient()
{
    auto &_sys = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    const DofMap & dof_map = _nl.dofMap();    
    PetscMatrix<Number> *petsc_mat = dynamic_cast<libMesh::PetscMatrix<Number>* >(_sys.matrix);

    //petsc_mat->print_matlab("petsc_mat.txt");
    //_stab_matrix.print_matlab("original.txt");    

    //    int m=dof_map.n_dofs();
    //
    //    int n=dof_map.n_dofs();
    //
    //    int m_l=dof_map.n_local_dofs();
    //
    //    int n_l=dof_map.n_local_dofs();


    // int nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());

    StoreOperators & storeOperatorsUO=_fe_problem.getUserObjectTempl<StoreOperators>(_userObjectName);

    auto _M  = storeOperatorsUO.MassMatrix();
    auto _L  = storeOperatorsUO.LumpMassMatrix();
    auto _J  = storeOperatorsUO.JacMatrix();
    auto _PL = storeOperatorsUO.PoroLumpMassMatrix();
    auto _PM = storeOperatorsUO.PoroMassMatrix();
    auto _D  = storeOperatorsUO.StabMatrix();

    //_M->print_matlab("porous_matrix.m");

    if(_fe_problem.timeStep()==1)
    {
        
        AntidiffusiveFluxes::determine_dc_bnd_var_id(AntidiffusiveFluxes::split_string(_dc_var, ' '));
        find_boundary(zero_rows, _dc_boundary_id);
    }




    PetscVector<Number> _inv(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _inv.zero();

    PetscVector<Number> _inv_p(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _inv_p.zero();

    PetscVector<Number> _u_dot(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _u_dot.zero();

    PetscVector<Number> _tmp(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _tmp.zero();


    _sol_vec_fluxes  = storeOperatorsUO.SolVec();

    //_sol_vec_fluxes->print_matlab("_sol_vec_fluxes.m");

    
    //NumericVector<Number> &ghosted_solution = *_sys.current_local_solution.get();
    
  
    std::vector<double>_vec_localize;
    std::vector<double>_vec_localize_dot;

    _sol_vec_fluxes->localize(_vec_localize);
    _J->vector_mult(_tmp, *_sol_vec_fluxes);

    _L->get_diagonal(_inv); 
    _inv.reciprocal();

    _PL->get_diagonal(_inv_p); 
    _inv_p.reciprocal();


    _u_dot.pointwise_mult(_inv_p,_tmp);
    _u_dot.localize(_vec_localize_dot);

    
    // NumericVector<Number> * _u_dot_moose = _nl.solutionUDot();
    // PetscVector<Number> &_u_dot = dynamic_cast<libMesh::PetscVector<libMesh::Number>& >(*_u_dot_moose);
    //_u_dot.print_matlab("_u_dot.m");

    PetscMatrix<Number> _I_tr(_pp_comm);
    //_J->get_transpose(_J_tr);
    //Mat L_petsc     = _L->mat();
    
    int r_start = _M->row_start();
    int r_stop  = _M->row_stop();


    Mat M_petsc    = _M->mat();
    Mat D_petsc    = _D->mat();


    PetscVector<Number> _R_p(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _R_p.zero();

    PetscVector<Number> _R_m(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _R_m.zero();

    PetscVector<Number> _a_bar(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _a_bar.zero();

    PetscVector<Number> _m_i(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _m_i.zero();

    PetscVector<Number> _ones(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _ones.zero();


    _ones.add(1.0);


    
    PetscMatrix<Number> _f_mat(_pp_comm);
    _f_mat.attach_dof_map(dof_map);
    _f_mat.init();

    PetscMatrix<Number> _alpha_mat(_pp_comm);
    _alpha_mat.attach_dof_map(dof_map);
    _alpha_mat.init();

    {
        std::vector<Real> v1;
        std::vector<int> c1;
        int rb=petsc_mat->row_start();
        int re=petsc_mat->row_stop();
        int rb1=_f_mat.row_start();
        //int re1=_f_mat.row_stop();
        //int rb2=_alpha_mat.row_start();
        int re2=_alpha_mat.row_stop();

        if (rb!=rb1)
        {
            std::cout<<"rb\n";
            exit(1);
        }
        if (re!=re2)
        {
            std::cout<<"re\n";
            exit(1);
        }
        if (rb!=rb1)
        {
            std::cout<<"rb\n";
            exit(1);
        }
        if (re!=re2)
        {
            std::cout<<"re\n";
            exit(1);
        }
        for (int row=rb; row<re; ++row)
        {
            getRow(*petsc_mat,row,v1,c1);
            for (int ii=0; ii<v1.size(); ++ii)
            {
                int col=c1.at(ii);
                _f_mat.set(row,col,0.0);
                _alpha_mat.set(row,col,0.0);
            }
        }
    }


    _f_mat.close();
    _alpha_mat.close();

    std::unique_ptr<NumericVector<Number>> f = _sol_vec_fluxes->zero_clone();

    
    _PL->vector_mult(_m_i,_ones);
    
    //Real volume = _m_i.dot(_ones);
    //std::cout<<"domian_sum"<<volume<<std::endl;


    Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();



    for (int row=r_start; row<r_stop; ++row)
    {
        PetscInt ncols_m, ncols_d; //, ncols_l;

        PetscInt const *cols_m; 
        PetscInt const *cols_d; 
        //PetscInt const *cols_l;


        PetscScalar const *val_m; 
        PetscScalar const *val_d; 
        //PetscScalar const *val_l;
        
        Real _P_p=0.0; 
        Real _Q_p=0.0; 

        Real _P_m=0.0; 
        Real _Q_m=0.0; 

        Real f_ij = 0.0;
        //Real f_ij_sum = 0.0;

        MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
        MatGetRow(D_petsc,row,&ncols_d,&cols_d,&val_d);


        
        for (int p=0; p<ncols_m; p++)
        {
            int col=cols_m[p];

            if (row!=col)
            {
                
                Real Dij  = val_d[p];

                Real PMij = val_m[p];
                
                auto v_i = (*_sol_vec_fluxes)(row);

                auto v_d_i = _u_dot(row);

                f_ij = PMij * (v_d_i - _vec_localize_dot.at(col)) - 1.0 * Dij * (v_i - _vec_localize.at(col));


                double check = (_vec_localize.at(col) - v_i) * f_ij;

                //if(check>1.0e-16)std::cout<<"check"<<check<<std::endl;

                if(check > 1.0e-16) _f_mat.set(row, col,0.0);
               
                else _f_mat.set(row, col, f_ij);

                _P_p+=std::max(0.0,f_ij);

                _P_m+=std::min(0.0,f_ij);             

            }

        }



        MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
        
        MatRestoreRow(D_petsc,row,&ncols_d,&cols_d,&val_d);

        
        //Real u_max = ghosted_solution(row);
        //Real u_min = ghosted_solution(row);

        std::vector<Real> values;
        values.clear();


        MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m);

        for (int p=0; p<ncols_m; ++p){
                int col=cols_m[p];
                values.push_back(_vec_localize.at(col));
            }

        MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);   

        double v_i = (*_sol_vec_fluxes)(row);

        double m_i = _m_i(row);

        _Q_p = m_i/dt * (*std::max_element(std::begin(values), end(values)) - v_i);

        _Q_m = m_i/dt * (*std::min_element(std::begin(values), end(values)) - v_i);   
            
        
        Real value_p = 0.0;

        if(std::abs(_P_p)>1e-16) {
            
            value_p = _Q_p/_P_p;
            //std::cout<<"value_p==>"<<_P_p<<std::endl;
        }

        Real value_m = 0.0;

        if(std::abs(_P_m)>1e-16) {

            value_m = _Q_m/_P_m;
        }

        Real entry_p = std::min(1.0,value_p);

        Real entry_m = std::min(1.0,value_m);

        auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

        //std::cout<<"zero_rows.size()"<<zero_rows.size()<<std::endl;

        if(it == zero_rows.end()){

            _R_p.set(row,entry_p);

            _R_m.set(row,entry_m);
         

        }
        else

        {
            _R_p.set(row,1.0);
            _R_m.set(row,1.0);
            
        }
    }

     
    _f_mat.close();
   
    _R_m.close();

    _R_p.close();

    //_f_mat.print_matlab("f_d.m");

    std::vector<double>_vec_localize_r_p;
    _R_p.localize(_vec_localize_r_p);



    std::vector<double>_vec_localize_r_m;
    _R_m.localize(_vec_localize_r_m);


    // _R_m.print_matlab("R_m.m");

    // _R_p.print_matlab("R_p.m");



    for (int row=r_start; row<r_stop; ++row)
    {
        PetscInt ncols_m, ncols_f; //, ncols_l;

        PetscInt const *cols_m, *cols_f; //, *cols_l;
     
        PetscScalar const *val_m, *val_f; //, *val_l;

        std::vector<double> values;
        
        values.clear();
        

        MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m); 

        for (int p=0; p<ncols_m; ++p){
                int col=cols_m[p];
                values.push_back(_vec_localize.at(col));
            }
        MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);   
   
        MatGetRow(_f_mat.mat(),row,&ncols_f,&cols_f,&val_f); 

        MatGetRow(M_petsc,row,&ncols_m,&cols_m,&val_m); 
 

        for (int p=0; p<ncols_m; ++p)
        {
            int col=cols_m[p];



            if(row!=col) {

                Real ris = 0.0;

                if(val_f[p]>=0.0) {
                                   
                    double entry_p = std::min(_R_p(row),_vec_localize_r_m.at(col));
                    ris = val_f[p] * entry_p;
                    _alpha_mat.set(row, col, ris);


                 
                    
                }

                else 
                {

                    double entry_m = std::min(_R_m(row),_vec_localize_r_p.at(col));
                    ris = val_f[p] * entry_m;
                    _alpha_mat.set(row, col, ris);

                }

            }
                
        }

        MatRestoreRow(M_petsc,row,&ncols_m,&cols_m,&val_m);  
 
        MatRestoreRow(_f_mat.mat(),row,&ncols_f,&cols_f,&val_f); 


    }


    _alpha_mat.close();


    //_alpha_mat.print_matlab("alpha.m");

    _alpha_mat.vector_mult(_a_bar,_ones);


    for (int row=r_start; row<r_stop; ++row)
    {
          
        auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

        if(it == zero_rows.end()){

            //std::cout<<"_a_bar(row)"<<_a_bar(row)<<std::endl;
            //std::cout<<"_inv(row)"<<_inv(row)<<std::endl;
            

            auto f_bar = _a_bar(row) * _inv(row) * dt;

            (*f).set(row,f_bar);
        }
        else{

           (*f).set(row,0.0);
        }


    }


   

             
        
    (*f).close();

    PetscVector<Number> _f_I(_pp_comm, dof_map.n_dofs(), dof_map.n_local_dofs());
    _f_I.zero();

    auto _I = storeOperatorsUO.H_Interpolator();

    _I->get_transpose(_I_tr);

    //_I->print_matlab("Interp.m");

    _I_tr.vector_mult(_f_I,*f);

    //(*f).print_matlab("f.m");

    //if(_fe_problem.timeStep()>1) 

    _sol_vec_fluxes->add(_f_I);

    // ghosted_solution.print_matlab("sol.m");


    PetscVector<Number> &f_c= dynamic_cast<libMesh::PetscVector<libMesh::Number>& >(*f);


    //PetscVector<Number> &sol = dynamic_cast<libMesh::PetscVector<libMesh::Number>& >(ghosted_solution);
    

    if(_write_correction) set_solution(f_c);

}

void
AntidiffusiveFluxes::find_boundary(std::vector<int> &zero_rows, std::vector<int> &_dc_boundary_id){

  ConstBndNodeRange & bnd_nodes = *_fe_problem.mesh().getBoundaryNodeRange();
  NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();
    unsigned int i = 0;
    // std::cout<<"_dc_variables_id"<< _dc_variables_id[0].size()<<std::endl;
    // std::cout<<"_dc_boundary_id"<< _dc_boundary_id.size()<<std::endl;
    for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
      {
        // iterate just over boundary nodes
            for (const auto & bnode : bnd_nodes)
            {
                  libMesh::Node * current_node = bnode->_node;

                  // check if node is in active boundary list
                  if (_fe_problem.mesh().isBoundaryNode(current_node->id(), *boundary))
                  {
                    // loop over all variables at this node

                    for (auto v = 0; v < _fe_problem.getNonlinearSystemBase().nVariables(); v++)
                    {
                      const Variable & var = _nl.system().variable(v);
                      unsigned int var_num = var.number();
                        //std::cout<<"nnnnnnn"<< var_num <<std::endl;

                      // see if this variable has any dofs at this node
                      if (current_node->n_dofs(_fe_problem.getNonlinearSystemBase().number(), var_num) > 0)
                      {
                        // check if given variable has BC on node

                        if(std::find(_dc_variables_id[i].begin(), _dc_variables_id[i].end(), var_num) != _dc_variables_id[i].end())
                        {

                          // different components are not supported by moose at the moment...
                          //std::cout<<"kkkkkkkk"<< std::endl;
                          zero_rows.push_back(
                              current_node->dof_number(_fe_problem.getNonlinearSystemBase().number(), var_num, 0));
                        }
                    }
                }
            } 
        }
    }

    //std::cout<<"zero_rows"<< zero_rows.size()<<std::endl;
}


void 
AntidiffusiveFluxes::determine_dc_bnd_var_id(const std::vector<std::string> & BC_var){
    // automatic fill-in
     NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();

    std::vector<int> vec(_nl.nVariables());

    std::iota(vec.begin(), vec.end(), 0);

    unsigned int i;

    auto str_tmp = BC_var.begin();

    PetscFunctionBegin;
    // going over all BC_ids
    for(i = 0; str_tmp != BC_var.end(); i++, str_tmp++)
    {
        std::vector<std::string> tmp = AntidiffusiveFluxes::split_string(*str_tmp, '-');

        // check if variable assigned in the input file exists for given simulation
        bool var_flg = 1;
        for(auto t = tmp.begin(); t != tmp.end(); ++t)
        {
            if(atoi(t->c_str()) >= _nl.nVariables())
                var_flg = 0;
        }

        // in case u havent put anything into input file, or u put too much
        if(*str_tmp == "-1" || var_flg == 0)
        {
            //std::cout<<"no_si"<<_nl.nVariables()<<std::endl;
            _dc_variables_id.push_back(vec);
        }
        else
        {
            unsigned int j;
            std::vector<int > one_BC_id;
            auto str_in = tmp.begin();
            for(j = 0; str_in != tmp.end(); j ++, str_in++)
            {
                one_BC_id.push_back(atoi(str_in->c_str()));
            }
            _dc_variables_id.push_back(one_BC_id);
        }
    }

    // check if u have same number of BC_ids in both parameters
    if(_dc_variables_id.size() != _dc_boundary_id.size())
    {
        _dc_variables_id.clear();
        for(auto i = 0; i != _dc_boundary_id.size(); i++)
        {
            _dc_variables_id.push_back(vec);
        }
    }

    // print out what is considered for zero-ing
    std::cout<<" ------ BC CONDITIONS  ------ \n";
    unsigned int t = 0;
    //std::cout<<"_dc_variables_id.begin()"<<_dc_variables_id.size()<<std::endl;
    for(auto i = _dc_variables_id.begin(); i != _dc_variables_id.end();  t++, i++)
    {
        std::cout<<"\n BC_id:  "<< _dc_boundary_id[t] << "   var_ids:  ";
        std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
    }

}


    

     
std::vector<std::string>
AntidiffusiveFluxes::split_string(const std::string & s, char delim)
{

      std::vector<std::string> v;

      if (s.length() == 0)
        std::cerr << "Got an empty string. Split_string(...) is confused. \n";

      auto i = 0;
      auto pos = s.find(delim);
      while (pos != std::string::npos)
      {
        v.push_back(s.substr(i, pos - i));
        i = ++pos;
        pos = s.find(delim, pos);

        if (pos == std::string::npos)
          v.push_back(s.substr(i, s.length()));
      }

      if (v.size() == 0) // if only one word is in the string
        v.push_back(s);

      return v;
}

void AntidiffusiveFluxes::set_solution(PetscVector<Number> &correction)
{
    // copy projected solution into target es
    
    MooseVariableFEBase  & main_var = _fe_problem.getVariable(0, "CM", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

    MooseVariableFEBase  & aux_var_c = _fe_problem.getVariable(0, "correction", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

    MooseVariableFEBase  & sol_var = _fe_problem.getVariable(0, "CM", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
    
    // solution of the original system
    System & main_sys = main_var.sys().system();

    //NumericVector<Number> * main_solution = main_sys.solution.get();

    System & aux_sys = aux_var_c.sys().system();

    NumericVector<Number> * aux_solution = aux_sys.solution.get();

    NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();


  
    { 
        
        for (const auto & node : _fe_problem.es().get_mesh().local_node_ptr_range())

        {
            for (unsigned int comp = 0; comp < node->n_comp(main_sys.number(), main_var.number()); comp++)

            {

                const dof_id_type proj_index = node->dof_number(_nl.number(), sol_var.number(), comp);

                //const dof_id_type to_index = node->dof_number(main_sys.number(), main_var.number(), comp);

                const dof_id_type to_index_c = node->dof_number(aux_sys.number(), aux_var_c.number(), comp);

                //main_solution->set(to_index, sol(proj_index));

                aux_solution->set(to_index_c, correction(proj_index));
            }

        }
    }



  //main_solution->close();
  aux_solution->close();
  //main_sys.update();
  aux_sys.update();

  //ExodusII_IO (_fe_problem.es().get_mesh()).write_equation_systems("matrix_c.e", _fe_problem.es());

    
}



