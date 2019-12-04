  #include "FlowAndTransport.h"

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
	params.addCoupledVar("pressure",
		"The gradient of this variable will be used as "
		"the velocity vector.");
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
_U(declareProperty<RealVectorValue>("VelocityVector"))
{
	_id=RealTensorValue(0.0,0.0,0.0,
		0.0,0.0,0.0,
		0.0,0.0,0.0);

	for (int i=0; i<_dim; ++i)
		_id(i,i)=1.0;

}

void
FlowAndTransport::computeQpProperties()
{
	_poro[_qp]=_poroInput;
	_K[_qp]=  _condInput * _id;
	//std::cout<<_qrule->n_points()<<std::endl;
	if (_isPressureValid)
	{
		if (_conservativeScheme)
		{
			if (_qp==0)
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
						A(i,i)+=-1.0*_JxW[qp]/_condInput;
					}
					_u_elem(i)=b(i)/A(i,i);
				}
			}
			//RealVectorValue temp=-1.0 * _K[_qp] * _gradP[_qp];
			//std::cout<<temp<<std::endl;
			//std::cout<<_u_elem<<std::endl;
			//std::cout<<std::endl;
			_U[_qp]=_u_elem;
		}
		else
		{
			_U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];
		}
	}
	else
	{
		for (int i=0; i<3; ++i)
			_U[_qp](i)=0.0/0.0;
	}
}
