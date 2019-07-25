    //* This file is part of the MOOSE framework
    //* This file is part of the MOOSE framework
    //* https://www.mooseframework.org
    //*
    //* All rights reserved, see COPYRIGHT for full restrictions
    //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
    //*
    //* Licensed under LGPL 2.1, please see LICENSE for details
    //* https://www.gnu.org/licenses/lgpl-2.1.html

    #include "SolveDiffusion.h"
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
    #include <string> 
    using namespace std;

    registerMooseObject("parrotApp", SolveDiffusion);

    void
    assembly_diffusion(EquationSystems & es, const std::string & system_name)
    {
        std::cout<<"CIAO\n";
        SolveDiffusion * pippo = es.parameters.get< SolveDiffusion*>("pippo");
        pippo->AssembleDiffusionOP(es, system_name);
    }


    template <>
    InputParameters
    validParams<SolveDiffusion>()
    {
      InputParameters params = validParams<GeneralUserObject>();
      
      params.addRequiredParam<std::vector<int>>("block_id",
                                                         "The name of the nodeset to create");
      params.addRequiredParam<std::vector<Real>>("value_p",
                                                         "The name of the nodeset to create");

      params.addRequiredParam<std::vector<boundary_id_type>>("boundary_D_bc",
                                                         "The name of the boundary to create");

      params.addRequiredParam<std::vector<boundary_id_type>>("boundary_N_bc",
                                                         "The name of the boundary to create");

      params.addRequiredParam<std::vector<std::string>>("function_D_bc",
                                                         "The value of the dirichlet_bc");

      params.addRequiredParam<std::vector<Real>>("value_N_bc", "The value of Neumann");

      params.addRequiredParam<std::vector<Real>>("value_D_bc", "The value of Dirichlet");

      params.addRequiredParam<std::vector<AuxVariableName>>(
          "aux_variable", "The auxiliary variable to store the transferred values in.");




      return params;
    }

    SolveDiffusion::SolveDiffusion(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _aux_var_names(getParam<std::vector<AuxVariableName>>("aux_variable")),
    _vector_p(getParam<std::vector<int>>("block_id")),
    _vector_value(getParam<std::vector<Real>>("value_p")),
    _boundary_D_ids(getParam<std::vector<boundary_id_type>>("boundary_D_bc")),
    _boundary_N_ids(getParam<std::vector<boundary_id_type>>("boundary_N_bc")),
    _string_expr(getParam<std::vector<std::string>>("function_D_bc")),
    _value_N_bc(getParam<std::vector<Real>>("value_N_bc")),
    _value_D_bc(getParam<std::vector<Real>>("value_D_bc"))

    {
      if (_aux_var_names.size() == 0)
        paramError("variable", "You need to specify at least one variable");

      /* Right now, most of transfers support one variable only */
      if (_aux_var_names.size() == 1)
        _aux_var_name = _aux_var_names[0];
    }

     void SolveDiffusion::execute(){};


     void SolveDiffusion::initialize()
    {
        std::cout<<"initialize\n";

       // _fe_problem.es().get_mesh().elem_ref(0).set_refinement_flag (Elem::REFINE);

       // Create an equation systems object.
       EquationSystems equation_systems (_fe_problem.es().get_mesh());

       equation_systems.parameters.set<bool>("adaptivity") = true;

       LinearImplicitSystem & _system = equation_systems.add_system<LinearImplicitSystem> ("Diffusion");

       equation_systems.parameters.set<SolveDiffusion *>("pippo") = this;

       _p_var = _system.add_variable ("pressure", FIRST);

       _system.attach_assemble_function(assembly_diffusion);

       //add_bc(_system);

       equation_systems.init();

       equation_systems.print_info();

       solve(equation_systems, _system);

       set_solution(equation_systems);


       //    _p_var = _fe_problem.es().add_system<LinearImplicitSystem>("Diffusion").add_variable("pressure",FIRST),
        
       //    _fe_problem.es().reinit();
     
        
    }

     void SolveDiffusion::finalize()
    {
        std::cout<<"finalize\n";
    }

     void SolveDiffusion::threadJoin(const UserObject & uo)
    {
    	
    };



    int SolveDiffusion::solve(EquationSystems & _es, LinearImplicitSystem & _system){


    	   std::cout<<"solve\n";

    	     //MeshBase & _mesh_m = _fe_problem.es().get_mesh();

           //LinearImplicitSystem & _system = equation_systems.get_system<LinearImplicitSystem> ("Diffusion");
           
           //libMesh::Parallel::Communicator const & comm_in=_mesh_m.comm();
           
           AssembleDiffusionOP(_es,"Diffusion");

           NumericVector<Number> & sol_NV=*_system.solution;
           NumericVector<Number> & rhs_NV=*_system.rhs;

           PetscVector<Number> & sol_PV= dynamic_cast<PetscVector<Number> &>(sol_NV);
           PetscVector<Number> & rhs_PV= dynamic_cast<PetscVector<Number> &>(rhs_NV);

           libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_system.matrix);

    
           //_system.solve()

           PetscErrorCode ierr;
           PC _diff_problem;

           ierr = PCCreate(PETSC_COMM_WORLD, &_diff_problem);
           CHKERRQ(ierr);

           std::cout<<"solve::Begin1\n";
           
           ierr = PCSetType(_diff_problem,PCLU);
           CHKERRQ(ierr);

           std::cout<<"solve::Begin2\n";

           ierr = PCSetOperators(_diff_problem, petsc_mat_m->mat(),petsc_mat_m->mat());
           CHKERRQ(ierr);

           std::cout<<"solve::Begin3\n";
          
           ierr = PCFactorSetMatSolverPackage(_diff_problem,MATSOLVERMUMPS);
           CHKERRQ(ierr);
           std::cout<<"solve::Begin4\n";

           


           ierr = PCApply(_diff_problem,rhs_PV.vec(),sol_PV.vec()); CHKERRQ(ierr);

           std::cout<<"solve::END\n";

        

           //solution->print_matlab();

           ExodusII_IO (_es.get_mesh()).write_equation_systems("matrix_c.e", _es);


            return 0 ;
           
           //NumericVector<Number> & sol_ref = solution[0];

    }

    void SolveDiffusion::AssembleDiffusionOP(EquationSystems & _es, const std::string & system_name){

    	_console << "Assemble_Diffusion::Begin "  << std::endl;

        // Get a constant reference to the mesh object.
        const MeshBase & mesh = _es.get_mesh();

        // The dimension that we are running.
        const unsigned int dim = mesh.mesh_dimension();

        // Get a reference to our system.
        LinearImplicitSystem & _system = _es.get_system<LinearImplicitSystem>(system_name);

        // Get a constant reference to the Finite Element type
        // for the first (and only) variable in the system.
        FEType fe_type = _system.get_dof_map().variable_type(0);

        UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

        SparseMatrix<Number> & matrix_K = *_system.matrix;
     
        DenseMatrix<Number> ke;

        DenseVector<Number> re;

        QGauss qrule (dim, fe_type.default_quadrature_order());

        // Tell the finite element object to use our quadrature rule.
        fe->attach_quadrature_rule (&qrule);

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

        const DofMap & dof_map = _system.get_dof_map();

        std::vector<dof_id_type> dof_indices;

        MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

        const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

        // first we need to manually zero the matrix
        matrix_K.zero();

        for ( ; el != end_el; ++el)
        {
            const Elem * elem = *el;

            Elem * ele = *el;

            fe->reinit (elem);

            dof_map.dof_indices(elem, dof_indices);

            const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());

            // With one variable, we should have the same number of degrees
            // of freedom as shape functions.
            libmesh_assert_equal_to (n_dofs, phi.size());

            //std::cout<<"n_dofs "<<n_dofs<<std::endl;

            ke.resize (n_dofs , n_dofs);
            ke.zero();

            for (unsigned int qp=0; qp<qrule.n_points(); qp++){

                for (unsigned int i=0; i<phi.size(); i++){

                    for (unsigned int j=0; j<phi.size(); j++){

                        ke(i,j) += ComputeMaterialProprties(elem) * JxW[qp] * dphi[i][qp] * dphi[j][qp];
                    }
                }

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

                    std::cout<<_boundary_N_ids[k] <<" "<<_value_N_bc.at(k) << std::endl;

                    if (mesh.get_boundary_info().has_boundary_id (elem, side, _boundary_N_ids[k])) // Apply a traction on the right side
                      {
                         for (unsigned int qp=0; qp<qface->n_points(); qp++)
                         for (unsigned int i=0; i != n_dofs; i++)
                         		re(i) += JxW_face[qp] * -1.0 * _value_N_bc.at(k) * phi_face[i][qp];
                      }
                   }
                 }
              }
         
           const Real penalty = 1.e10;

           for (auto s : elem->side_index_range())
           {

            if (elem->neighbor_ptr(s) == nullptr)
            {
                  fe_face->reinit(elem, s);

                  const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
                  const std::vector<Real> & JxW_face = fe_face->get_JxW();

                  std::vector<Real> whichOnBoundary(phi_face.size(),0.0);

                  for (int i=0; i<phi_face.size(); ++i)
                  {
                    for (int qp=0; qp<phi_face.at(i).size(); ++qp)
                      whichOnBoundary.at(i)+=phi_face.at(i).at(qp);
                  }

                  Real boundaryValue=0;

                  for(int k =0; k < _boundary_D_ids.size(); k++)
                  {

                    if (mesh.get_boundary_info().has_boundary_id (elem, s, _boundary_D_ids[k]))
                    {
                       boundaryValue=_value_D_bc.at(k);


                       _console << "Assemble_Diffusion:: Add penalty BC on Id"  << _boundary_D_ids[k] <<std::endl;
                    

                    for (int k=0; k<whichOnBoundary.size(); ++k)
                      {
                        if (whichOnBoundary.at(k)>0.5){
                       re(k) += penalty * boundaryValue;
                       ke(k,k) += penalty;
                     }
                   }
                 }
               }
             }
           }
          


            std::cout<<ke<<std::endl;
            dof_map.constrain_element_matrix_and_vector (ke, re, dof_indices,true);
            //dof_map.constrain_element_matrix_and_vector (ke, re, dof_indices);
            std::cout<<ke<<std::endl;

            _system.matrix->add_matrix (ke, dof_indices);

            _system.rhs->add_vector(re, dof_indices);

            //exit(1);



        }



        _system.matrix->close();

        _system.rhs->close();

        _system.matrix->print_matlab("Diff.m");
        _system.rhs->print_matlab("rhs.m");

        _console << "Assemble_Diffusion::end "  << std::endl;

        
    }


    Real
    SolveDiffusion::ComputeMaterialProprties(const Elem *elem){
        
       // _console << "_vector_p.size()"  << _vector_p.size() <<std::endl;

    Real permeability=0.0;

        for(int ll=0; ll<_vector_p.size(); ll++){
            if (elem->subdomain_id()==_vector_p[ll]) {

                permeability = _vector_value[ll]; 
            }
        }
      
        return permeability;   
    }

    void SolveDiffusion::add_bc(LinearImplicitSystem & _system){

    	_console << "SolveDiffusion::add_bc"  << std::endl;

      std::vector<ParsedFunction<Number> > pf;

      pf.clear();

      DofMap & dof_map = _system.get_dof_map();

      // ParsedFunction<Number> pf_d2(_string_expr);

      // Most DirichletBoundary users will want to supply a "locally
      // indexed" functor

      std::vector<unsigned int> variables(1);

      variables[0] = _p_var; 

      for(int k =0; k < _boundary_D_ids.size(); k++){

      	    std::set<boundary_id_type> boundary_dr;

      	    boundary_dr.clear();

            boundary_dr.insert(_boundary_D_ids[k]);

            ParsedFunction<Number> temp(_string_expr.at(k));

            pf.push_back(temp);

            std::cout<<_boundary_D_ids[k] <<" "<<_string_expr.at(k) << std::endl;

            dof_map.add_dirichlet_boundary(DirichletBoundary(boundary_dr,
                                                             variables,
                                                             pf.at(k), LOCAL_VARIABLE_ORDER));

      }



      //  DirichletBoundary dirichlet_bc(boundary_dr, variables, pf,
                                     //LOCAL_VARIABLE_ORDER);


      //LinearImplicitSystem & _system = _fe_problem.es().get_system<LinearImplicitSystem>("Diffusion");

      // We must add the Dirichlet boundary condition _before_
      // we call equation_systems.init()
      //_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

    }



    void SolveDiffusion::set_solution(EquationSystems & _es)
    {
        // copy projected solution into target es

        return;
        MeshBase & _mesh = _fe_problem.es().get_mesh();


        
      //  MooseVariableFEBase  & _var = _fe_problem.getVariable(0, _aux_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
        
        // solution of the original system
        // System & _sys = _var.sys().system();

        // NumericVector<Number> * _solution = _sys.solution.get();
        
      
        LinearImplicitSystem & ls = _es.get_system<LinearImplicitSystem>("Diffusion");

        NumericVector<Number> & from_solution = *ls.solution;
        // { // loop through our local elements and set the solution from the projection
            
        //     for (const auto & node : _mesh.local_node_ptr_range())

        //     {
        //     	for (unsigned int comp = 0; comp < node->n_comp(_sys.number(), _var.number()); comp++)

        //     	{

        //         std::cout<<"uno"<<std::endl;

        //     		const dof_id_type proj_index = node->dof_number(ls.number(), _p_var, comp);
                
        //         std::cout<<"due"<<std::endl;

        //     		const dof_id_type to_index = node->dof_number(_sys.number(), _var.number(), comp);

        //         std::cout<<"tre"<<std::endl;

        //          //

        //     		_solution->set(to_index, (*ls.solution)(proj_index));
        //     	}

        //     }
        // }

      //auto &aux_sys = _fe_problem.getAuxiliarySystem();



      MooseVariableFEBase & aux_var = _fe_problem.getVariable(
          0, _aux_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

      //_nl.setSolution((*ls.solution));

      System & aux_sys = aux_var.sys().system();

      NumericVector<Number> * aux_solution = aux_sys.solution.get();

      // std::set<dof_id_type> from_var_indices;
      // ls.local_dof_indices(_p_var, from_var_indices);

      // unsigned int aux_vn = aux_var.number();
      // std::set<dof_id_type> aux_var_indices;
      // _fe_problem.es().get_system("aux0").local_dof_indices(aux_vn, aux_var_indices);

      //  // copy the values from from solution vector to to solution vector
      //  std::set<dof_id_type>::iterator from_it = from_var_indices.begin();
      //  std::set<dof_id_type>::iterator from_it_end = from_var_indices.end();
      //  std::set<dof_id_type>::iterator aux_it = aux_var_indices.begin();

      //  for (; from_it != from_it_end; ++from_it, ++aux_it){

      //    aux_sys.solution->set(*aux_it, from_solution(*from_it));
      //  }

      //  std::cout<<"CIAO sono fuori"<<std::endl;
      //  aux_sys.solution->close();
      //  aux_sys.update();


      // NumericVector<Number> & from_solution = *from_sys->solution;

      // for (; from_it != from_it_end; ++from_it, ++to_it)
      //   to_sys->solution->set(*to_it, from_solution(*from_it));


      // ;

       //MeshBase::const_node_iterator it = _mesh.local_nodes_begin();
         //   const MeshBase::const_node_iterator end_it = _mesh.local_nodes_end();
           // for (; it != end_it; ++it)
            for (const auto & node : _mesh.local_node_ptr_range())
            {
                // set solution for variable 1
                //const Node * node = *it;

                 //unsigned int comp = 0;
                  for (unsigned int comp = 0;
                       comp < node->n_comp(aux_sys.number(), aux_var.number()); comp++)
                 {
                    std::cout<<"comp"<<comp<<" "<<aux_sys.number()<<" "<<aux_var.number()<<std::endl;
                    const dof_id_type pre_index =
                    node->dof_number(ls.number(), _p_var, comp);
               
                    const dof_id_type index =
                    node->dof_number(aux_sys.number(), aux_var.number(), comp);
                    //std::cout<<"index"<<index<<std::endl;
                    aux_solution->set(index, (*ls.solution)(pre_index));
                  
                }
            }
        

      // for (const auto & node : _mesh.local_node_ptr_range())
      // {
      //   for (unsigned int comp = 0; comp < node->n_comp(ls.number(), _p_var); comp++)
      //   {
      //     const dof_id_type pre_index = node->dof_number(ls.number(), _p_var, comp);
      //     std::cout<<comp<<std::endl;
      //     const dof_id_type aux_index = node->dof_number(aux_sys.number(), aux_var.number(), comp);
      //     aux_solution->set(aux_index, (*ls.solution)(pre_index));
      //   }
      // }


      aux_solution->close();
      aux_sys.update();

        
    }
