//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleVolumeVectors.h"
#include "FractureUserObject.h"

registerMooseObject("parrotApp", AssembleVolumeVectors);


template <>
InputParameters
validParams<AssembleVolumeVectors>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addParam<int>("NRegions",0,"NRegions");
  params.addParam<std::vector<int>>("block_id","block_id");
  // params.addRequiredParam<bool>("constrain_matrix","constrain_matrix");
  //params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
  params.addParam<bool>("FractureRegions",false,"FractureRegions");
  params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  return params;
}

AssembleVolumeVectors::AssembleVolumeVectors(const InputParameters & parameters) :
GeneralUserObject(parameters),
_hasMeshModifier( isParamValid("fractureMeshModifier") ),
_r_fracture(getParam<bool>("FractureRegions")),
_r_fractureUO(isParamValid("fractureMeshModifier")),
_n_regions(getParam<int>("NRegions"))
// _vector_value(getParam<std::vector<Real>>("value_p")),
// userObjectName(getParam<UserObjectName>("operator_userobject")),
// _constrainMatrices(getParam<bool>("constrain_matrix")),
// _code_dof_map(true),

// _qrule(_assembly.qRule())
{
  _vectorAllocated=false;
  if (_hasMeshModifier)
    _meshModifierName=getParam<std::string>("fractureMeshModifier");

  if(_r_fracture)
    _block_id=getParam<std::vector<int>>("block_id");
}

void AssembleVolumeVectors::execute()
{
  _console << "AssembleVolumeVectors::execute() begin"  << std::endl;

  //MeshModifier     const * _myMeshModifier_ptr;
  RegionUserObject const * _myRegionUserObject_ptr;
  FractureUserObject const * _fractureUserObject_ptr;

  if (_hasMeshModifier)
  {
    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    //_myMeshModifier_ptr=&_myMeshModifier;
    RegionUserObject const & _regionUserObject( dynamic_cast<RegionUserObject const &>(_myMeshModifier) );
    
    if(_r_fractureUO) 
    {
        FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
        
        _fractureUserObject_ptr=&_fractureUserObject;
    }
    else{

      _myRegionUserObject_ptr=&_regionUserObject;
    
    }
  
  }

  DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();
  auto &comm = _fe_problem.es().get_mesh().comm();

  if (_r_fracture==false) _fn = _myRegionUserObject_ptr[0].get_fn();
  else _fn = _n_regions;

  volumes.resize(_fn);

  for (int i=0; i<_fn; ++i)
  {
    volumes[i]=new PetscVector<Number>(comm);
    volumes[i][0].init( dof_map.n_dofs() , dof_map.n_local_dofs() );
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

  QBase const * const & qbase(_assembly.qRule());

  std::unique_ptr<QBase> qrule( QBase::build (qbase->type(),dim,qbase->get_order()));

  fe->attach_quadrature_rule (qrule.get());
  
  const std::vector<Real>& JxW      = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  //const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  const std::vector<Point>& q_points = fe->get_xyz();

  std::vector<dof_id_type> dof_indices;
  // this has to be resized
  std::vector< DenseVector<Number> > Re;
  std::vector< std::vector<Number> > coeff;

  Re.resize(_fn);
  coeff.resize(_fn);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem * elem = *el;
    fe->reinit (elem);
    dof_map.dof_indices(elem, dof_indices);

    int const loc_n=dof_indices.size();
    for (int i=0;i<_fn;++i)
    {
      Re.at(i).resize(loc_n);
      Re.at(i).zero();
      coeff.at(i).resize( qrule->n_points() );
    }




    Real ze=0.0;
    for (int qp=0; qp<qrule->n_points(); ++qp)
    {
      for (int i=0; i<_fn; ++i)
      {
        if (_r_fracture == false && _r_fractureUO==false && _myRegionUserObject_ptr[0].isInsideRegion(q_points[qp],i,ze) ){
          
          coeff.at(i).at(qp)=1.0;
        }
        else if (_r_fractureUO && _fractureUserObject_ptr[0].isInsideRegion(q_points[qp],i,ze)){
          
          coeff.at(i).at(qp)=1.0;
        }
        else if(_r_fracture  && elem->subdomain_id() == _block_id[i]){
          
          coeff.at(i).at(qp)=1.0;
        }
        else{
          coeff.at(i).at(qp)=0.0;
        }
      }
    }
    
    for (unsigned int i=0; i<phi.size(); i++)
    {
      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      {
        for (unsigned int j=0; j<_fn; ++j )
          Re.at(j)(i) += coeff.at(j).at(qp) * JxW[qp] * phi[i][qp];
      }
    }

    for (int i=0; i<_fn; ++i)
    {
      volumes[i][0].add_vector(Re.at(i), dof_indices);
    }
  }

  for (int i=0; i<_fn; ++i)
  {
    volumes[i][0].close();
    //std::cout<<volumes.at(i)[0].sum()<<std::endl<<std::endl;
  }

   _console << "AssembleVolumeVectors::execute() end"  << std::endl;

}

AssembleVolumeVectors::~AssembleVolumeVectors()
{
 
  RegionUserObject const * _myRegionUserObject_ptr;

  if (_hasMeshModifier)
  {
    MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
    
    RegionUserObject const & _regionUserObject( dynamic_cast<RegionUserObject const &>(_myMeshModifier) );

    if(!_r_fracture) _myRegionUserObject_ptr=&_regionUserObject;
  }

  int fn = 0;

  if(_r_fracture) {
    
    fn = _n_regions;
  }
  else{
    
    fn = _myRegionUserObject_ptr[0].get_fn();
  }

  for (int i=0; i<fn; ++i) delete volumes[i];

}
