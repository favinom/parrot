#include "FlowAndTransport.h"
#include "FractureUserObject.h"

  // MOOSE includes
  //#include "Assembly.h"
  //#include "MooseMesh.h"

registerMooseObject("parrotApp", FlowAndTransport);

  template <>
InputParameters
validParams<FlowAndTransport>()
{
	InputParameters params = validParams<Material>();
	params.addRequiredParam<Real>("phi","phi");
	params.addRequiredParam<Real>("k","k");
	params.addRequiredParam<bool>("conservative","use a conservative scheme?");
	params.addCoupledVar("pressure","The gradient of this variable will be used as the velocity vector.");

	params.addParam<std::string>("fractureMeshModifier","fractureMeshModifier");
	params.addParam<Real>("phiFrac","phiFrac");
	params.addParam<Real>("kFrac","kFrac");

	return params;
}

FlowAndTransport::FlowAndTransport(const InputParameters &parameters) :
Material(parameters),
_condInput(getParam<Real>("k")),
_poroInput(getParam<Real>("phi")),
_dim(_mesh.dimension()),
_isPressureValid(parameters.isParamValid("pressure")),
_conservativeScheme(getParam<bool>("conservative")),
_gradP(_isPressureValid ? coupledGradient("pressure"): _grad_zero),
_poro(declareProperty<Real>("Porosity")),
_K(declareProperty<RealTensorValue>("PermeabilityTensor")),
_Kscalar(declareProperty<Real>("Permeability")),
_U(declareProperty<RealVectorValue>("VelocityVector")),
_hasMeshModifier(parameters.isParamValid("fractureMeshModifier"))
{
	_id=RealTensorValue(0.0,0.0,0.0,
		0.0,0.0,0.0,
		0.0,0.0,0.0);

	for (int i=0; i<_dim; ++i)
		_id(i,i)=1.0;

	if ( _hasMeshModifier )
	{
		if ( !isParamValid("phiFrac") || !isParamValid("kFrac") )
		{
			mooseError("you provided a mesh meodifier but not phiFrac or poroFrac");
		}
	}

	if ( isParamValid("phiFrac") )
	{
		if ( !_hasMeshModifier || !isParamValid("kFrac") )
		{
			mooseError("you provided phiFrac but not a mesh meodifier or poroFrac");
		}
	}

	if ( isParamValid("kFrac") )
	{
		if ( !_hasMeshModifier || !isParamValid("phiFrac") )
		{
			mooseError("you provided poroFrac but not a mesh meodifier or phiFrac");
		}
	}

	if (_hasMeshModifier)
	{
		_meshModifierName=getParam<std::string>("fractureMeshModifier");
  		_condFracture=getParam<Real>("kFrac");
  		_poroFracture=getParam<Real>("phiFrac");
	}

}

void
FlowAndTransport::computeQpProperties()
{
	_poro[_qp]   =_poroInput;
	_Kscalar[_qp]=_condInput;
	_K[_qp]      =_condInput * _id;

	if (_hasMeshModifier)
	{
		MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
		FractureUserObject const & _fractureUserObject( dynamic_cast<FractureUserObject const &>(_myMeshModifier) );
		if (_fractureUserObject.isInside(_q_point[_qp]))
		{
			_poro[_qp]   =_poroFracture;
			_Kscalar[_qp]=_condFracture;
			_K[_qp]      =_condFracture * _id;
		}
	}

	if (_isPressureValid)
	{
		_U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];
	}
	else
	{
		for (int i=0; i<3; ++i)
			_U[_qp](i)=0.0/0.0;
	}

	if (_conservativeScheme && _qp==_qrule->n_points()-1)
	{
		RealVectorValue b;
		RealTensorValue A;
		b.zero();
		A.zero();
		for (int i=0; i<_dim; ++i)
		{
			for (int qp=0; qp<_qrule->n_points(); ++qp)
			{
				b(i)+=_JxW[qp]*_gradP[qp](i);
				A(i,i)+=-1.0*_JxW[qp]/_Kscalar[qp];
			}
			_u_elem(i)=b(i)/A(i,i);
		}
		for (int qp=0; qp<_qrule->n_points(); ++qp)
		{
			_U[qp]=_u_elem;
		}
	}
}
