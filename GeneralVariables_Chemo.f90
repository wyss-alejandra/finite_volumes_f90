MODULE GeneralVariables
IMPLICIT NONE 
PUBLIC 
!CHARACTER*3 WithManufSol
INTEGER :: icase,Ncells,neq,ngauss,nt,nx,ny,typeBC,NNX,NNY,WithManufSol
INTEGER :: yderivatve_upperB=0 ! 1 means non-zero derivative. 0 means Neumann homogeneous.
INTEGER :: rr=3 ! Stencil: i-rr,i-(rr-1),...,i,i+1,...,i+rr
REAL    :: ax, bx, ay, by,CFL, cfl_difx,cfl_dify,cfl_advx,cfl_advy, cfd, dt, dx, dy,t,tmax
REAL,ALLOCATABLE :: xn(:),yn(:),x(:),y(:),xx(:),yy(:),w(:,:),cg(:), pg(:)

! WENO Reconstruction
REAL                :: M1(3,3), ML(3,3), MR(3,3), Sigma(4,4), epsilon, lambda(0:3), r
REAL                :: pi=acos(-1.)
REAL                :: alpha1, alpha2, lambdaa1, lambdaa2, tau1, tau2, beta1
REAL                :: d_u, d_w, alphav, mus, kappas, fm, xi, kv, lambdav
REAL, ALLOCATABLE   :: xig(:), wg(:), coef(:,:,:,:), WWX0(:,:,:), WWX1(:,:,:), WWX2(:,:,:), WWX3(:,:,:), WWX(:,:,:,:), LO(:,:,:), OI(:,:), omegah(:,:), omegas(:,:), wweno(:,:)
REAL, ALLOCATABLE   :: FFO(:,:,:,:), fR(:,:,:,:,:), fL(:,:,:,:,:), fSO(:,:,:,:,:), fNO(:,:,:,:,:), FDifus(:,:,:,:), GDifus(:,:,:,:)
REAL, ALLOCATABLE   :: fDR(:,:,:,:,:), fDL(:,:,:,:,:), fDSO(:,:,:,:,:), fDNO(:,:,:,:,:)
REAL, ALLOCATABLE   :: source(:,:,:,:,:)!, MLL(:,:), MRR(:,:), MC1(:,:), iMLL(:,:), iMRR(:,:), iMC1(:,:)
REAL, ALLOCATABLE   :: GFO(:,:,:,:), FFLUX(:,:,:), GFLUX(:,:,:), FDFLUX(:,:,:), GDFLUX(:,:,:), SSource(:,:,:)
REAL, ALLOCATABLE   :: FDFLUX1(:,:,:), GDFLUX1(:,:,:), PolL(:,:,:,:), PolR(:,:,:,:), PolD(:,:,:,:), PolU(:,:,:,:)
REAl, ALLOCATABLE   :: DPolR_x(:,:,:,:), DPolL_x(:,:,:,:), DPolD_y(:,:,:,:), DPolU_y(:,:,:,:)
REAL, ALLOCATABLE   :: iMLL(:,:), iMRR(:,:), iMC0(:,:), iMC1(:,:)
REAL :: mm=1.
!INTEGER :: neq=1 ! Number of PDEs + 1ODE
end MODULE GeneralVariables