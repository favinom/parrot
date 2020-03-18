//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParrotProblem.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "petscmat.h"
#include "libmesh/petsc_matrix.h"

registerMooseObject("parrotApp", ParrotProblem);

template <>
InputParameters
validParams<ParrotProblem>()
{
    InputParameters params = validParams<FEProblem>();
    params.addRequiredParam<bool>("use_AFC","use_AlgFluxCorr");
    params.addParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
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


ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters),
_pp_comm(_mesh.comm()),
_stab_matrix(_pp_comm),
_use_afc(getParam<bool>("use_AFC"))
{
    
    if(parameters.isParamValid("operator_userobject"))
    {
        _hasStoreOperatorsUO=true;
        userObjectName=new UserObjectName(getParam<UserObjectName>("operator_userobject"));
        
        // std::vector<UserObject *> userobjs;
        // theWarehouse().query().condition<AttribSystem>("UserObject").queryInto(userobjs);
        // std::cout<<userobjs.size()<<std::endl;
        // exit(1);
    }
    else
    {
        _hasStoreOperatorsUO=false;
        _storeOperatorsUO=NULL;
    }
    
    PCCreate(PETSC_COMM_WORLD, &_problem_PC);
    PCSetType(_problem_PC,PCLU);
    
    _const_jacobian=true;
    _is_stab_matrix_assembled=false;
}

void ParrotProblem::initialSetup()
{
    std::cout<<"BEGIN ParrotProblem::initialSetup"<<std::endl;
    
    FEProblem::initialSetup();
    
    Moose::PetscSupport::petscSetOptions(*this);
    
    PetscErrorCode ierr;
    
    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());
    
    SNES snes = petsc_solver->snes();
    KSP ksp;
    SNESType ttype;
    KSPType type;
    ierr = SNESGetKSP(snes,&ksp);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    std::cout<< "SNES type from PARROTPROBLEM exec:"<<ttype<<std::endl;
    ierr = KSPGetType(ksp, &type);
    std::cout<< "KSP  type from PARROTPROBLEM exec:"<<type<<std::endl;
    _factorized = 0 ;
    
    if (_hasStoreOperatorsUO)
    {
        _storeOperatorsUO=&getUserObject<StoreOperators>(userObjectName[0]);
        PetscMatrix<Number> * & localmatrix=_storeOperatorsUO[0].StabMatrix();
        localmatrix=&_stab_matrix;
    }
    
    std::cout<<"END ParrotProblem::initialSetup"<<std::endl;
    
};

void ParrotProblem::timestepSetup()
{
    std::cout<<"BEGIN ParrotProblem::timestepSetup"<<std::endl;
    
    FEProblem::timestepSetup();
    
    PetscErrorCode ierr;
    Moose::PetscSupport::petscSetOptions(*this);
    
    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());
    SNES snes = petsc_solver->snes();
    KSP ksp;
    PC pc;
    SNESType ttype;
    KSPType type;
    PCType ptype;
    
    ierr = SNESGetKSP(snes,&ksp);
    ierr = KSPGetPC(ksp,&pc);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    ierr = KSPGetType(ksp, &type);
    ierr = PCGetType(pc, &ptype);
    
    std::cout<< "SNES type da passo steady exec:"<<ttype<<std::endl;
    std::cout<< "KSP type da passo steady exec:"<<type<<std::endl;
    std::cout<< "PC type da passo steady exec:"<<ptype<<std::endl;
    
    _ksp_ptr = (KSP_PARROT *)ksp->data;
    (_ksp_ptr[0].local_pc)=&_problem_PC;
    (*_ksp_ptr).factorized=&_factorized;
    //_ksp_ptr->_fe_problem = this;
    
    std::cout<<"END ParrotProblem::timestepSetup"<<std::endl;
};

void
ParrotProblem::computeJacobianSys(NonlinearImplicitSystem & /*sys*/,
                                  const NumericVector<Number> & soln,
                                  SparseMatrix<Number> & jacobian)
{
    
    
    
    
    
    NonlinearSystemBase & nl_sys = this->getNonlinearSystemBase();
    
    DofMap const & dof_map = nl_sys.dofMap();
    
    
    FEProblemBase::computeJacobian(soln, jacobian);
    
    
    //libMesh::PetscMatrix<libMesh::Number> & petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>&>(jacobian);
    
    //FEProblemBase::computeJacobian(soln, jacobian);
    
    if (_use_afc==true && _is_stab_matrix_assembled==false)
    {
        computeStabilizationMatrix(jacobian);
        jacobian.add(1.0,_stab_matrix);
        _jac_matrix->add(1.0,jacobian);
    }
    
    
}

void
ParrotProblem::computeResidualSys(NonlinearImplicitSystem & /*sys*/,
                                  const NumericVector<Number> & soln,
                                  NumericVector<Number> & residual)
{
    
    NonlinearImplicitSystem * a;
    //FEProblemBase::computeResidualSys(a[0],soln,residual);
    
    NonlinearSystemBase & nl_sys = this->getNonlinearSystemBase();
    
    DofMap const & dof_map = nl_sys.dofMap();
    
    NumericVector<Number> & old_solution = nl_sys.solutionOld();
    
    PetscVector<Number> residual_mio(this->es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
    
    
    if(this->timeStep()==1){
        //    _poro_lumped = _storeOperatorsUO->PoroLumpMassMatrix();
        FEProblemBase::computeResidualSys(a[0],soln,residual);
        //    _jac_matrix = _storeOperatorsUO->JacMatrix();
    }
    else
    {
        _poro_lumped = _storeOperatorsUO->PoroLumpMassMatrix();
        
        PetscVector<Number> res_j(this->es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
        
        PetscVector<Number> res_m(this->es().get_mesh().comm(), dof_map.n_dofs(), dof_map.n_local_dofs());
        
        _dirichlet_bc = _storeOperatorsUO->BcVec();
        
        _dirichlet_bc->print_matlab("bc.m");
        
        FEProblemBase::computeResidualSys(a[0],soln,residual);
        
        residual.print_matlab("residual_moose_prima.m");
        
        soln.print_matlab("soln.m");
        
        //PetscMatrix()
        
        // FEProblemBase::computeJacobian(soln, jacobian);
        
        
        //_jac_matrix->vector_mult_add(res_j,soln);
        
        
        //old_solution.print_matlab("old_solution.m");
        
        _poro_lumped->vector_mult(res_m, soln);
        
        // res_m.add(res_j);
        
        residual_mio.pointwise_mult(res_m,*_dirichlet_bc);
        
        Real inv_dt = -1.0/this->dt();
        
        residual_mio.scale(inv_dt);
        
        residual_mio.print_matlab("residual_mio.m");
        
        //_poro_lumped->print_matlab("poro.m");
        
        
        
        //soln.print_matlab("old_solution_1.m");
    }
    
    
    
    if (_use_afc && _is_stab_matrix_assembled)
    {
        _stab_matrix.vector_mult_add(residual,soln);
        
        
        
        //_stab_matrix.vector_mult_add(residual,soln);
        
        
        //if(this->timeStep()>1) residual.scale(-1.0);
        
        
        
        
        residual.print_matlab("residual_moose_dopo.m");
        
    }
    
}


void
ParrotProblem::computeStabilizationMatrix(SparseMatrix<Number> & jacobian)
{
    std::cout<<"START ParrotProblem::computeStabilizationMatrix\n";
    PetscMatrix<Number> & jac_PM=dynamic_cast<PetscMatrix<Number> &> (jacobian);
    // Declare a PetscMatrix that will contain the transpose
    PetscMatrix<Number> jac_tr_PM(_pp_comm);
    // Transpose of the Jacobian
    jacobian.get_transpose (jac_tr_PM);
    
    int rb=jac_PM.row_start();
    int re=jac_PM.row_stop();
    
    _jac_matrix      = _storeOperatorsUO->JacMatrix();
    
    std::cout<<"get dof map\n";
    DofMap const & dof_map = this->getNonlinearSystemBase().dofMap();
    std::cout<<"done\n";
    
    std::cout<<"attach dofmap\n";
    _stab_matrix.attach_dof_map(dof_map);
    _jac_matrix->attach_dof_map(dof_map);
    std::cout<<"done\n";
    
    std::cout<<"init matrix\n";
    _stab_matrix.init();//m,n,m_l,n_l,30);
    _jac_matrix->init();
    std::cout<<"done\n";
    
    std::cout<<"zeroing\n";
    _stab_matrix.zero();
    _jac_matrix->zero();
    std::cout<<"done\n";
    //MatSetOption(_stab_matrix_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    
    std::vector<Real> values;
    std::vector<Real> values_tr;
    std::vector<int> columns;
    std::vector<int> columns_tr;
    
    std::cout<<"looping\n";
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
                    
                    _jac_matrix->set(row,col,value);
                    _jac_matrix->add(row,row,maxEntry);
                }
            }
        }
    }
    std::cout<<"done\n";
    
    std::cout<<"closing\n";
    _stab_matrix.close();
    _jac_matrix->close();
    std::cout<<"done\n";
    
    _is_stab_matrix_assembled=true;
    
    std::cout<<"STOP ParrotProblem::computeStabilizationMatrix\n";
}







%Vec Object: Vec_0x84000000_0 1 MPI processes
%  type: mpi
Vec_0x84000000_0 = [
-3.7685736594262652e-18
-1.0553853506137496e-17
-1.1859008958523291e-17
-8.1982625312439527e-19
-2.3564138128078803e-17
-1.5047036050916415e-16
-2.4778410770934545e-16
-1.6794429527040446e-15
-1.5406814086771128e-15
-1.5681600030869331e-14
-1.0332280057292315e-14
-1.2492308201088681e-13
-5.1242154109999265e-14
-8.2209460335563023e-13
-1.8906297007442471e-13
-4.0929601942861514e-12
-5.0151620887274102e-13
-1.3943371271904943e-11
-9.7554213701502833e-13
-3.0981412977626988e-11
-1.4642635186152567e-12
-4.6072971050625920e-11
-1.8174269068840661e-12
-4.8722136192675690e-11
-2.0018917503089401e-12
-3.8852213826906937e-11
-2.0797015367239044e-12
-2.4408243072517163e-11
-2.1265916265334017e-12
-1.2746337647028207e-11
-2.1834077901443710e-12
-5.5660747554845296e-12
-2.2612625043587659e-12
-2.0522799332422689e-12
-2.3569944490796523e-12
-6.4140669383451500e-13
-2.4642147518494034e-12
-1.6980350633803213e-13
-2.5780289169359433e-12
-3.7909511598629565e-14
-2.6957793390221202e-12
-7.0834991540940628e-15
-2.8165087688168909e-12
-1.0966312227659422e-15
-2.9399842711796694e-12
-1.3891905199539111e-16
-3.0676193980511732e-12
-1.4178436691451719e-17
-3.1991651599645832e-12
-1.1499461077082629e-18
-3.3547550538963795e-12
-7.4900319006593152e-20
-3.4783169199301590e-12
-3.7589521103280870e-21
-1.7799735228882363e-12
-1.4097345686149967e-22
-6.2038545941477076e-25
-2.0679515313825692e-25
-2.1765189867801541e-23
-2.7409792819680975e-21
-1.9767798393633393e-18
-3.8306427378319978e-16
-2.2826240907454679e-14
-4.4650983001895369e-13
-3.5230590878064074e-12
-1.3906832375919469e-11
-3.2406022921727898e-11
-5.0095093386878322e-11
-5.5624602876490127e-11
-4.6817501568190206e-11
-3.0966335347333632e-11
-1.6479861747846830e-11
-7.1594986614136901e-12
-2.5592064958104247e-12
-7.5487313255436889e-13
-1.8346070492490042e-13
-3.6562397796353596e-14
-5.9284396035400644e-15
-7.7319946811898357e-16
-7.9827982255022525e-17
-6.3792981684033034e-18
-3.8172355249778930e-19
-1.6286690660861793e-20
-4.5492856456589561e-22
-2.0679515313825692e-25
2.0679515313825692e-25
2.5849394142282115e-26
-4.4784075351503764e-24
-1.0836628248767257e-20
-5.6615916422300106e-18
-7.4136226855041550e-16
-3.1648462693480385e-14
-5.1327709239318908e-13
-3.7373236091976572e-12
-1.4512934614097821e-11
-3.4432749116903893e-11
-5.4785652964971250e-11
-6.2352294411302017e-11
-5.3085415123947946e-11
-3.4862626366471372e-11
-1.8024104995323031e-11
-7.4302077813147905e-12
-2.4590638899425493e-12
-6.5467110866358157e-13
-1.3992961718092685e-13
-2.3860864690078011e-14
-3.2088892080387152e-15
-3.3410906239022299e-16
-2.6126389368059249e-17
-1.4640610845946793e-18
-5.4337070823724301e-20
-1.1208938468432610e-21
1.1632227364026952e-25
-1.5509636485369269e-25
-3.4896682092080855e-25
0.0000000000000000e+00
-2.2175498512712335e-23
-5.2664279311649806e-20
-2.4369304039487070e-17
-2.3785764148898339e-15
-7.1065001653401930e-14
-8.6186361601093114e-13
-5.2375168243801863e-12
-1.8382514685613574e-11
-4.0971262263305069e-11
-6.2107330001590180e-11
-6.7297311415291786e-11
-5.4022029771536806e-11
-3.2950173696268866e-11
-1.5539016309518900e-11
-5.7297373864142979e-12
-1.6617689200444188e-12
-3.7945085074935627e-13
-6.7921518900436326e-14
-9.4278490071097759e-15
-9.9470997324259042e-16
-7.7065146112911813e-17
-4.1538676954105394e-18
-1.3893470017196441e-19
-2.1502122981343026e-21
0.0000000000000000e+00
0.0000000000000000e+00
-2.5849394142282115e-25
5.1698788284564230e-26
-3.3862706326389570e-24
-4.0851090032755540e-21
-1.8570338651677128e-18
-2.3311129490131520e-16
-9.9276749279573446e-15
-1.7557691848847116e-13
-1.5116024233366302e-12
-7.1728751806078411e-12
-2.0633792301967106e-11
-3.8664147508369729e-11
-4.9830198388278084e-11
-4.6019520939978729e-11
-3.1402336167455652e-11
-1.6190836542296646e-11
-6.4063865862530476e-12
-1.9640072471104419e-12
-4.6831474128953851e-13
-8.6668668508115556e-14
-1.2332315233953541e-14
-1.3234417244648633e-15
-1.0340363372190791e-16
-5.5109085753322053e-18
-1.7566842146595996e-19
-2.5706077059416370e-21
-2.7832771856383164e-13
-6.3248074254864817e-13
-1.2633194140714804e-12
-5.5413273003057093e-13
-7.4468969988927893e-13
-1.5004741948985093e-12
-8.5064491943627224e-13
-1.7388212309091448e-12
-9.5732652619562771e-13
-1.9892051774874567e-12
-1.0653348441174173e-12
-2.2655508133436376e-12
-1.1735352283026046e-12
-2.5855703760588692e-12
-1.2819908202974612e-12
-3.0201771092130748e-12
-1.3904896587607970e-12
-3.8807023014997127e-12
-1.4982986192202764e-12
-5.8239194751194636e-12
-1.6030195356394429e-12
-9.2751499154817685e-12
-1.6982728085119965e-12
-1.3820439591498434e-11
-1.7697289657221356e-12
-1.8560877298399600e-11
-1.7905396260983762e-12
-2.2744925855410720e-11
-1.7207141378242611e-12
-2.6119897987092446e-11
-1.5197663755352287e-12
-2.8689271926608843e-11
-1.1810586991885139e-12
-3.0551367302391035e-11
-7.6626719729979382e-13
-3.1861162850701227e-11
-3.9227530229115109e-13
-3.2792919346346315e-11
-1.5081310172625204e-13
-3.3512555340343344e-11
-4.2235641602970062e-14
-3.4145754917079013e-11
-8.5537026473868190e-15
-3.4758059968368519e-11
-1.2017355119693712e-15
-3.5356713249802868e-11
-1.0981435149635802e-16
-3.5883964042127466e-11
-6.9797088222365258e-17
-3.6181963261568818e-11
-5.7288827778314470e-17
-3.6247301644164624e-11
-5.3663232499813807e-17
-3.6010540399638375e-11
-2.6059686413101854e-17
-1.8732291844858077e-11
-1.2753103535798405e-12
-5.5832626725075874e-13
-1.5229581377897430e-12
-1.7865456574482337e-12
-2.0808685012510603e-12
-2.4447422980331011e-12
-3.0324103579407127e-12
-4.4668518013916057e-12
-7.9971448441652242e-12
-1.4272775413270138e-11
-2.1983170140371754e-11
-2.8784108825029719e-11
-3.3162198799660696e-11
-3.4998953022708943e-11
-3.4998047664190249e-11
-3.3910067585155982e-11
-3.2348405409850339e-11
-3.0763017366871025e-11
-2.9332361495750625e-11
-2.8072448886730338e-11
-2.7029093410374085e-11
-2.6192640747765520e-11
-2.5517433767611256e-11
-2.4938204130858894e-11
-2.4364569213138318e-11
-2.3673302307068437e-11
-2.3219292543027419e-11
-1.1833127873428219e-11
-1.2910655705456469e-12
-5.6335039763241950e-13
-1.5561068787951254e-12
-1.8668669366781576e-12
-2.3072908904002682e-12
-3.1717021291106853e-12
-5.2506207882838616e-12
-9.7816902444605611e-12
-1.7115042351204424e-11
-2.5016219485296981e-11
-3.0319874650615993e-11
-3.1623586735704594e-11
-2.9681642032285272e-11
-2.6104504703950581e-11
-2.2206599262603699e-11
-1.8687493532575639e-11
-1.5820391906761982e-11
-1.3600669907174179e-11
-1.1864575662293771e-11
-1.0518553551491878e-11
-9.5330815192665099e-12
-8.8048016693076273e-12
-8.2616476451440268e-12
-7.8525681510655396e-12
-7.4935548151345360e-12
-7.1789461945081739e-12
-7.0136688339576108e-12
-3.5403459661109855e-12
-1.3244570463030400e-12
-5.7207267097389823e-13
-1.6727193679972337e-12
-2.2657827972041579e-12
-3.5405941067007970e-12
-6.3511635097270977e-12
-1.1203083954812966e-11
-1.7524365848047033e-11
-2.3033661669702230e-11
-2.5191472866531425e-11
-2.3635305996163410e-11
-1.9823416536101177e-11
-1.5453721356348602e-11
-1.1586048818223853e-11
-8.5775565660584911e-12
-6.3909650222554233e-12
-4.8516018204029691e-12
-3.7687926164626605e-12
-2.9960558576666658e-12
-2.4435325640356912e-12
-2.0469388302962289e-12
-1.7551788144884908e-12
-1.5360785693276425e-12
-1.3668267563329920e-12
-1.2332168423334807e-12
-1.1670542299705343e-12
-1.1455996035145239e-12
-5.7332953772303471e-13
-1.8228416262139308e-12
-8.6632565877381057e-13
-2.5693643318688185e-12
-3.8845805078708027e-12
-6.3358759738182989e-12
-9.9196194776230000e-12
-1.3706622444287737e-11
-1.6277372889166473e-11
-1.6383713494055438e-11
-1.4132111460977005e-11
-1.0772218360962939e-11
-7.5277227145867139e-12
-4.9968586000399096e-12
-3.2475596773830461e-12
-2.1146486649401470e-12
-1.4012657789275116e-12
-9.5257873940174825e-13
-6.6562127549031219e-13
-4.7843089831183527e-13
-3.5403070318422996e-13
-2.6914494524640939e-13
-2.0939039056645504e-13
-1.6606115892511582e-13
-1.3379431549319886e-13
-1.1040976607635412e-13
-9.4045618341107918e-14
-8.1971115554679983e-14
-3.9429874513991953e-14
-5.2889739230267155e-22
-8.3189883314799716e-20
-1.6255909677803483e-19
-1.5424407156411768e-21
-3.1122892922750172e-18
-5.9778134987967324e-18
-6.2551973368381268e-17
-1.1931678368236370e-16
-8.1752044982440108e-16
-1.5544400655770116e-15
-7.5462435169808158e-15
-1.4340374369229721e-14
-5.1322661720519181e-14
-9.7675404503477423e-14
-2.6295395374079938e-13
-5.0208075949217899e-13
-1.0253294410869650e-12
-1.9674305324266123e-12
-3.0483308625808112e-12
-5.8879488171381874e-12
-6.8773695293534717e-12
-1.3395136462355938e-11
-1.1650705451744441e-11
-2.2925496305049023e-11
-1.4579686068647369e-11
-2.9044860618525125e-11
-1.3176088802897450e-11
-2.6638718270741709e-11
-8.3458975406556145e-12
-1.7173629240826855e-11
-3.5636777659614778e-12
-7.4904918777255684e-12
-9.7521531655904441e-13
-2.1034146197387053e-12
-1.6009987788638821e-13
-3.5641534642301208e-13
-1.4454234606174691e-14
-3.3447074891699875e-14
-6.3978363890483905e-16
-1.5504881215571712e-15
-1.1963637883907085e-17
-3.0551294847715369e-17
-7.9325948659587655e-20
-2.1352866828072852e-19
-1.6742652585956126e-22
-4.6741520722928077e-22
-1.8094575899597480e-25
-5.1698788284564230e-25
-1.8094575899597480e-25
-1.9387045606711586e-25
5.1698788284564230e-26
1.0339757656912846e-25
-1.0339757656912846e-25
0.0000000000000000e+00
0.0000000000000000e+00
0.0000000000000000e+00
-1.6343851109030828e-19
-1.9861074167182318e-21
-5.2974105772324465e-18
-1.0142276658353003e-16
-1.3001457879684126e-15
-1.1948435842219785e-14
-8.1772863021760984e-14
-4.2529442357204749e-13
-1.6967284444519903e-12
-5.2008180949085864e-12
-1.2191065329122092e-11
-2.1627984138333808e-11
-2.8574473719467136e-11
-2.7488450079448130e-11
-1.8686293670099461e-11
-8.6324972355058043e-12
-2.5764934806888103e-12
-4.6527865644330989e-13
-4.6658153439611531e-14
-2.3191227389989671e-15
-4.9154869506544953e-17
-3.6899798173075752e-19
-8.4811862180827619e-22
-7.2378303598389922e-25
-1.0339757656912846e-25
1.0339757656912846e-25
0.0000000000000000e+00
0.0000000000000000e+00
-4.0010445223033656e-12
-1.8679895038577459e-12
-4.9764142093021272e-12
-6.7087916785057658e-12
-8.7189324901918231e-12
-1.0385256004687698e-11
-1.0921980039324451e-11
-9.9711647465737347e-12
-7.8746603145518951e-12
-5.4578392831494760e-12
-3.4144887263061388e-12
-1.9933811340950116e-12
-1.1217421519035634e-12
-6.2542048396307723e-13
-3.5267114960528420e-13
-2.0380273141863641e-13
-1.2144422109995123e-13
-7.4763404809609267e-14
-4.7595481798396984e-14
-3.1327041579305995e-14
-2.1263102138455447e-14
-1.4830734689953167e-14
-1.0596694985853375e-14
-7.7587558069705010e-15
-5.8613066730588393e-15
-4.6382641886433239e-15
-4.0181463044711665e-15
-1.9108959243375319e-15
-6.6234252179462119e-12
-3.1145584902134395e-12
-7.3058905132733872e-12
-8.1292627581605790e-12
-8.5448356569179264e-12
-8.1549400762013312e-12
-6.9045936363879981e-12
-5.1284998317572369e-12
-3.3471431811423739e-12
-1.9507653327642342e-12
-1.0430389137763689e-12
-5.2795130914760827e-13
-2.6077112590240281e-13
-1.2890844206225011e-13
-6.4953401361091238e-14
-3.3724190632350514e-14
-1.8129133934359842e-14
-1.0107304190361316e-14
-5.8463516607834254e-15
-3.5034979627897236e-15
-2.1682644290832583e-15
-1.3810112519675214e-15
-9.0329059401326657e-16
-6.0769723596712039e-16
-4.2358494540539086e-16
-3.1169974760057667e-16
-2.4957369079988771e-16
-1.2052415838744988e-16
-3.8514525244501664e-12
-1.8408932730246956e-12
-3.9390370981375040e-12
-3.8357849705450997e-12
-3.4210273265050039e-12
-2.7218569270842137e-12
-1.9065691361829281e-12
-1.1722361803036330e-12
-6.3870770316155488e-13
-3.1502489163739940e-13
-1.4479250432225227e-13
-6.3986634465890245e-14
-2.7975391767958118e-14
-1.2374504865015496e-14
-5.6200325731474393e-15
-2.6405577888921199e-15
-1.2871138353720266e-15
-6.5114382368914132e-16
-3.4134327682142199e-16
-1.8475191896518190e-16
-1.0275377522556377e-16
-5.8453210727817758e-17
-3.3910773744637054e-17
-2.0063168645788494e-17
-1.2168011483772661e-17
-7.6303897801073886e-18
-5.0411251538880844e-18
-1.6613478025389064e-18
];
