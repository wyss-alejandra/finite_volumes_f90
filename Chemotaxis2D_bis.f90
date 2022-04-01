USE GeneralVariables
IMPLICIT NONE
character*80 text
INTEGER nx1,i,j,k,l,NP1,kdib,numfile,imetODE
REAL    			:: tolE, difE(3)
REAL    			:: dt1, dt2, FV, k1, k2, k3, k4, taux, u00, uold, vold, ux, uy
REAL				:: x1, x2, Vx, Vy, y1,y2, vvel, Difusx, Difusy, vvelx, vvely, kappa, p, xr, yr, QV, wg1, wg2
REAL,ALLOCATABLE    :: u(:,:,:),u0(:,:,:),uRK(:,:,:),coeff(:,:),Pol(:),DPolx(:),DPoly(:),uu(:),f(:),g(:),df(:),dg(:),Q(:),sou(:)
REAL,ALLOCATABLE    :: vprey(:,:,:,:),vprey0(:,:,:,:),VAver(:,:)
!REAL,ALLOCATABLE    ::  vprey(:,:),vprey0(:,:),u(:,:),unorm(:),u0(:,:),uRK(:,:),coeff(:,:),Pol(:),DPolx(:),uu(:),f(:),df(:),Q(:),sou(:),ddd(:),SourceGauss(:,:,:),dvall(:,:),vall(:,:),dvalr(:,:),valr(:,:)
real :: XG,YG,TestRec,xtest,ytest,DxTestRec,DyTestRec,xvg(2,2),normM,normA,tol,NormL2,d1,auxnorm(3),media
!
! DATA INPUT
! ----------

! Data for upper boundary
! -----------------------

OPEN(0,FILE='Data.dat',STATUS='UNKNOWN')
READ(0,*) text   ! Coordinates of domain
READ(0,*) ax
READ(0,*) bx
READ(0,*) ay
READ(0,*) by
READ(0,*) text   ! Number of cells
READ(0,*) nx
READ(0,*) ny
READ(0,*) text   ! Number of Gaussian quadrature points
READ(0,*) NP1
READ(0,*) text   ! Output time
READ(0,*) tmax
print *,tmax
pause
READ(0,*) text
READ(0,*)  cfl_difx
READ(0,*) text
READ(0,*)  cfl_dify
READ(0,*) text
READ(0,*)  cfl_advx
READ(0,*) text
READ(0,*)  cfl_advy
CLOSE(0)
print *, tmax
pause 6
ngauss=3 !Number of Gaussian points
!Type of problem: icase=1 (Advection-diffusion), icase=2 (Euler)

WithManufSol=1 !0: NO, 1: YES
PRINT *, NP1
!number of variables (equations)
PRINT *, 'CHEMOTAXIS PROBLEM'
PAUSE
! allocating of variables
ALLOCATE(LO(neq,0:nx,0:ny),x(-2:nx+3),y(-2:ny+3),u(neq,-2:nx+2,-2:ny+2),u0(neq,-2:nx+2,-2:ny+2),uRK(neq,-2:nx+2,-2:ny+2))
ALLOCATE(coeff(neq,ngauss*ngauss),xig(NP1),wg(NP1),OI(neq,0:2),omegah(neq,0:2),omegas(neq,0:2),wweno(neq,0:2))
ALLOCATE(coef(neq,ngauss*ngauss,-2:nx+2,-2:ny+2),WWX0(neq,-2:nx+2,-2:ny+2),WWX1(neq,-2:nx+2,-2:ny+2),WWX2(neq,-2:nx+2,-2:ny+2),w(neq,0:8))
ALLOCATE(FFO(neq,0:nx+1,0:ny+1,NP1),FDifus(neq,0:nx+2,0:ny+2,NP1),fR(neq,0:nx+1,0:ny+1,NP1),fL(neq,0:nx+1,0:ny+1,ngauss),FFLUX(neq,0:nx+1,0:ny+1),FDFLUX(neq,0:nx+1,0:ny+1))
!ALLOCATE(FDFLUX1(neq,0:nx+2,0:ny+2),GDFLUX1(neq,0:nx+1,0:ny+1))
ALLOCATE(GFO(neq,0:nx+1,0:ny+1,NP1),GDifus(neq,0:nx+2,0:ny+2,NP1),GFLUX(neq,0:nx+1,0:ny+1),GDFLUX(neq,0:nx+1,0:ny+1),gD(neq,0:nx+1,0:ny+1,NP1),gU(neq,0:nx+1,0:ny+1,NP1))
ALLOCATE(fDR(neq,0:nx+2,0:ny+2,NP1),fDL(neq,0:nx+2,0:ny+2,NP1),gDD(neq,0:nx+2,0:ny+2,NP1),gDU(neq,0:nx+2,0:ny+2,NP1))
ALLOCATE(Pol(neq),DPolx(neq),DPoly(neq),Q(neq),f(neq),g(neq),df(neq),dg(neq),source(neq,0:nx+1,1:ny+1,NP1,NP1),sou(neq),SSource(neq,0:nx+1,0:ny+1))
ALLOCATE(PolL(neq,0:nx+1,0:ny+1,NP1),PolR(neq,0:nx+1,0:ny+1,NP1),PolD(neq,0:nx+1,0:ny+1,NP1),PolU(neq,0:nx+1,0:ny+1,NP1))
ALLOCATE(DPolR_x(neq,0:nx+1,0:ny+1,NP1),DPolL_x(neq,0:nx+1,0:ny+1,NP1),DPolD_y(neq,0:nx+1,0:ny+1,NP1),DPolU_y(neq,0:nx+1,0:ny+1,NP1))
!
! Mesh generation
! ---------------
!
x=0.
y=0.
CALL MALLADO
!DO i = 1, nx
!    print *, I,x(i)
!END DO
!PAUSE
!DO i = 1, nY
!    print *, I,Y(i)
!END DO
!PAUSE
!stop
! For WENO reconstruction.
CALL MATRICES1D
epsilon=1.d-15
r=12.
DO i = 0, 2
	IF (i==0) THEN
		lambda(i) = 1000.
	ELSE
		lambda(i) = 1.
	END IF
END DO
!
! Gaussian parameters
print *, y
pause
CALL Gauss_quadrature (NP1) ! xig, wg: ngauss components.
!
! FILES
!

!!!! TEST RECONSTRUCTION
!CALL TEST_RECONSTRUCTION
!STOP

!
! Initial condition
! -----------------
x(0)=x(1)-dx
y(0)=y(1)-dy
x(-1)=x(0)-dx
x(-2)=x(-1)-dx
y(-1)=y(0)-dy
y(-2)=y(-1)-dy
x(nx)=x(nx-1)+dx
y(ny)=y(ny-1)+dy
x(nx+1)=x(nx)+dx
y(ny+1)=y(ny)+dy
x(nx+2)=x(nx+1)+dx
y(ny+2)=y(ny+1)+dy
u0=0.

    OPEN(0,FILE='INITIALU.DAT',STATUS='UNKNOWN')
    OPEN(1,FILE='INITIALW.DAT',STATUS='UNKNOWN')
    OPEN(2,FILE='INITIALV.DAT',STATUS='UNKNOWN')
DO i = 1, nx-1
	x1=x(i)
	x2=x(i+1)
	DO j = 1, ny-1
		y1=y(j)
		y2=y(j+1)
		IF (WithManufSol==0) THEN 
		    CALL solinic(NP1,x1,x2,y1,y2,Q)
		ELSE
		    CALL SolManuf(0.,x1,x2,y1,y2,Q)
		END IF
		u0(:,i,j) = Q(:)
		WRITE(0,*) (x1+x2)/2.,(y1+y2)/2.,u0(1,i,j)
        WRITE(1,*) (x1+x2)/2.,(y1+y2)/2.,u0(2,i,j)
        WRITE(2,*) (x1+x2)/2.,(y1+y2)/2.,u0(3,i,j)
        k=k+1
	END DO
END DO
DO i = 0, neq-1
    CLOSE(i)
END DO

t =0.
uRK=0.
kdib=1.
numfile=44    
CALL Gauss_quadrature (Np1) ! Gaussian quadrature parameters in (0, 1)
tol=1.e-8
tolE=1.e-8
normM=2.*tol; normA=2.*tol 
difE(:)=2.*tolE
!DO WHILE (difE(1)>tolE.or.difE(2)>tolE)
!DO WHILE (difE(1)>tolE)
DO WHILE (t<tmax)
    CALL calctime(u0)
     
	t = t+dt
	IF (WithManufSol.eq.1) THEN
	    IF (t>tmax) THEN
		    taux = t-tmax
		    t = t - dt
		    dt = dt-taux
	        t = t + dt
		    write(*,*) 'time=',t, 'step size=',dt
	    END IF
	    write(*,*) 'time=',t, 'step size=',dt
	END IF
    nx1=nx-1
!
! Operator L reconstruction from cell averages.
!
	CALL BoundaryConditions(u0)
    coeff=0.
    Pol=0.
    DPolx=0.
    DPoly=0.
    CALL OpL(coeff,NP1,u0,vprey,Pol,DPolx,DPoly)
	DO i = 1, nx-1
		DO j = 1, ny-1
			uRK(:,i,j) = u0(:,i,j) + dt*LO(:,i,j)
		END DO
	END DO
!pause 

    CALL BoundaryConditions(uRK)
	CALL OpL(coeff,NP1,uRK,vprey,Pol,DPolx,DPoly)

	DO i = 1, nx-1
		DO j = 1, ny-1
			uRK(:,i,j) = (3.*u0(:,i,j) + (uRK(:,i,j) + dt*LO(:,i,j)))/4.
		END DO
	END DO


    CALL BoundaryConditions(uRK)

	CALL OpL(coeff,NP1,uRK,vprey,Pol,DPolx,DPoly)
	DO i = 1, nx-1
		DO j = 1, ny-1
			u(:,i,j) = (u0(:,i,j) + 2.*(uRK(:,i,j) + dt*LO(:,i,j)))/3.
			!print *, 'sol',u(2,i,j)
		END DO
	END DO
	IF (WithManufSol.EQ.0) THEN
    	!difE=0.
    	difE(:)=abs(u0(:,1,1)-u(:,1,1))
	    DO i = 1, nx-1
		    DO j = 1, ny-1
			!print *, 'sol',u(2,i,j)
			    !difE(:)=difE(:)+abs(u0(:,i,j)-u(:,i,j))**2
			    auxnorm(:)=abs(u0(:,i,j)-u(:,i,j))
			    IF (difE(1)<auxnorm(1)) THEN
			        difE(1)=auxnorm(1)
			    END IF
			        
		    END DO
	    END DO
        !difE(:)=sqrt(difE(:)/(nx*ny))
        print *, 'Time=',t,dt
        print *, difE(1), difE(2)
	END IF
	DO i = 1, nx-1
		DO j = 1, ny-1
			u0(:,i,j) = u(:,i,j)
!			print *, u(2,i,j)
		END DO
	END DO

END DO
CALL OUTPUT(u,VAver)
media=0.
DO i = 1, nx-1
		DO j = 1, ny-1
			media = media+u(1,i,j)
!			print *, u(2,i,j)
		END DO
	END DO
	media=media/((nx-1)*(ny-1))
	print *, 'media=',media
IF (WithManufSol==1) THEN
    NormL2=0.
    DO i = 1, nx-1
	    x1=x(i)
	    x2=x(i+1)
	    DO j = 1, ny-1
		    y1=y(j)
		    y2=y(j+1)
    	        CALL SolManuf(t,x1,x2,y1,y2,Q)
    	        print *, Q(1),u(1,i,j)
    		    NormL2=NormL2+(Q(1)-u(1,i,j))**2
        END DO
    END DO
    NormL2=sqrt(NormL2/(nx*ny))
    print *, 'Error',NormL2

END IF
STOP
END


FUNCTION Vx(x,y)
IMPLICIT NONE
REAL    :: Vx,x,y
Vx=0.1
END FUNCTION Vx

FUNCTION Vy(x,y)
IMPLICIT NONE
REAL    :: Vy,x,y
Vy=0.1
END FUNCTION Vy

FUNCTION Difusx(x,y)
IMPLICIT NONE
REAL    :: Difusx,x,y
Difusx=0.01
END FUNCTION Difusx

FUNCTION Difusy(x,y)
IMPLICIT NONE
REAL    :: Difusy,x,y
Difusy=0.01
END FUNCTION Difusy


function derf(x)
implicit REAL*8 (a - h, o - z)
dimension a(0 : 64), b(0 : 64)
data (a(i), i = 0, 12) / 0.00000000005958930743d0, -0.00000000113739022964d0, 0.00000001466005199839d0, -0.00000016350354461960d0, 0.00000164610044809620d0, -0.00001492559551950604d0, 0.00012055331122299265d0, -0.00085483269811296660d0, 0.00522397762482322257d0, -0.02686617064507733420d0, 0.11283791670954881569d0, -0.37612638903183748117d0, 1.12837916709551257377d0 / 
data (a(i), i = 13, 25) / 0.00000000002372510631d0, -0.00000000045493253732d0, 0.00000000590362766598d0, -0.00000006642090827576d0, 0.00000067595634268133d0, -0.00000621188515924000d0, 0.00005103883009709690d0, -0.00037015410692956173d0, 0.00233307631218880978d0, -0.01254988477182192210d0, 0.05657061146827041994d0, -0.21379664776456006580d0, 0.84270079294971486929d0 / 
data (a(i), i = 26, 38) / 0.00000000000949905026d0, -0.00000000018310229805d0, 0.00000000239463074000d0, -0.00000002721444369609d0, 0.00000028045522331686d0, -0.00000261830022482897d0, 0.00002195455056768781d0, -0.00016358986921372656d0, 0.00107052153564110318d0, -0.00608284718113590151d0, 0.02986978465246258244d0, -0.13055593046562267625d0, 0.67493323603965504676d0 / 
data (a(i), i = 39, 51) / 0.00000000000382722073d0, -0.00000000007421598602d0, 0.00000000097930574080d0, -0.00000001126008898854d0, 0.00000011775134830784d0, -0.00000111992758382650d0, 0.00000962023443095201d0, -0.00007404402135070773d0, 0.00050689993654144881d0, -0.00307553051439272889d0, 0.01668977892553165586d0, -0.08548534594781312114d0, 0.56909076642393639985d0 / 
data (a(i), i = 52, 64) / 0.00000000000155296588d0, -0.00000000003032205868d0, 0.00000000040424830707d0, -0.00000000471135111493d0, 0.00000005011915876293d0, -0.00000048722516178974d0, 0.00000430683284629395d0, -0.00003445026145385764d0, 0.00024879276133931664d0, -0.00162940941748079288d0, 0.00988786373932350462d0, -0.05962426839442303805d0, 0.49766113250947636708d0 / 
data (b(i), i = 0, 12) / -0.00000000029734388465d0, 0.00000000269776334046d0, -0.00000000640788827665d0, -0.00000001667820132100d0, -0.00000021854388148686d0, 0.00000266246030457984d0, 0.00001612722157047886d0, -0.00025616361025506629d0, 0.00015380842432375365d0, 0.00815533022524927908d0, -0.01402283663896319337d0, -0.19746892495383021487d0, 0.71511720328842845913d0 / 
data (b(i), i = 13, 25) / -0.00000000001951073787d0, -0.00000000032302692214d0, 0.00000000522461866919d0, 0.00000000342940918551d0, -0.00000035772874310272d0, 0.00000019999935792654d0, 0.00002687044575042908d0, -0.00011843240273775776d0, -0.00080991728956032271d0, 0.00661062970502241174d0, 0.00909530922354827295d0, -0.20160072778491013140d0, 0.51169696718727644908d0 / 
data (b(i), i = 26, 38) / 0.00000000003147682272d0, -0.00000000048465972408d0, 0.00000000063675740242d0, 0.00000003377623323271d0, -0.00000015451139637086d0, -0.00000203340624738438d0, 0.00001947204525295057d0, 0.00002854147231653228d0, -0.00101565063152200272d0, 0.00271187003520095655d0, 0.02328095035422810727d0, -0.16725021123116877197d0, 0.32490054966649436974d0 / 
data (b(i), i = 39, 51) / 0.00000000002319363370d0, -0.00000000006303206648d0, -0.00000000264888267434d0, 0.00000002050708040581d0, 0.00000011371857327578d0, -0.00000211211337219663d0, 0.00000368797328322935d0, 0.00009823686253424796d0, -0.00065860243990455368d0, -0.00075285814895230877d0, 0.02585434424202960464d0, -0.11637092784486193258d0, 0.18267336775296612024d0 / 
data (b(i), i = 52, 64) / -0.00000000000367789363d0, 0.00000000020876046746d0, -0.00000000193319027226d0, -0.00000000435953392472d0, 0.00000018006992266137d0, -0.00000078441223763969d0, -0.00000675407647949153d0, 0.00008428418334440096d0, -0.00017604388937031815d0, -0.00239729611435071610d0, 0.02064129023876022970d0, -0.06905562880005864105d0, 0.09084526782065478489d0 / 
w = abs(x)
IF (w .lt. 2.2d0) THEN
   t = w * w
   k = int(t)
   t = t - k
   k = k * 13
   y = ((((((((((((a(k) * t + a(k + 1)) * t + a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + a(k + 11)) * t + a(k + 12)) * w
ELSE IF (w .lt. 6.9d0) THEN
   k = int(w)
   t = w - k
   k = 13 * (k - 2)
   y = (((((((((((b(k) * t + b(k + 1)) * t + b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + b(k + 11)) * t + b(k + 12)
   y = y * y
   y = y * y
   y = y * y
   y = 1 - y * y
ELSE
   y = 1
END IF
IF (x .lt. 0) y = -y
      derf = y
END

SUBROUTINE calctime(ua)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,j,k
REAL    :: velx,vely,p,aa,dta,dt1,dt2,kappa,lambdaa,lambda1(6),ua(neq,-2:nx+2,-2:ny+2),vvel,vvelx,vvely
REAL    :: Vx,Vy,Difusx,Difusy,dtt(4),dtaux
dt=1.e20
  DO i = 1, nx
        DO j = 1, ny
            dta=(ua(1,i+1,j)-ua(1,i,j))/dx+1.e-20
            aa=(ua(1,i,j+1)-ua(1,i,j))/dy+1.e-20
            !print *, aa
            dtt(1)=min(cfl_difx*dx*dx/d_u,cfl_difx*dx*dx/d_w)
            dtt(2)=min(cfl_dify*dy*dy/d_u,cfl_dify*dy*dy/d_w)
            dtt(3)=min(cfl_advy*dy/Xi/ABS(aa),cfl_advy*dy/Xi/ABS(aa))
!            print *, dtt(3)
            dtt(4)=min(cfl_advx*dx/Xi/ABS(dta),cfl_advx*dx/Xi/ABS(dta))
! print *, dtt
            dtaux=minval(dtt)
            IF (dtaux<dt) THEN
               dt=dtaux
            END IF
!            print *, dt
!            pause
        END DO
  END DO

!dt=1.e-5

END SUBROUTINE calctime



SUBROUTINE flux_hyperb(xg,yg,Deriv,Q,f)
USE GeneralVariables
IMPLICIT NONE
INTEGER i
REAL    :: f(neq),g(neq),grav,H,u,Vx,Vy,xg,yg,p,kappa,xxx,yyy,xL,xR,yD,yU,Q(neq),ux,uy,Deriv,DPoly(neq)
f(1) = Xi*Deriv*Q(1)/dx
!f(1) = Xi*Q(1)!/dx
!print *, 'der',f(1),Deriv
f(2) = 0.
f(3) = 0.
END SUBROUTINE flux_hyperb

SUBROUTINE flux_diffusive(idifu,xg,yg,Pol,DPolx,df)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,n,idifu
REAL    :: Pol(neq),DPolx(neq),u,Difusx,Difusy,xg,yg,df(neq)
df(1) = d_u*DPolx(1) ! - Pol(1)*Xi *DPolx(2)
df(2) = d_w*DPolx(2)
df(3) = 0.
END SUBROUTINE flux_diffusive


SUBROUTINE SourceTerm(j,Pol,sou)
USE GeneralVariables
IMPLICIT NONE
INTEGER j
REAL    :: sou(neq),FU,mu,alpha,Pol(neq),rel,kappa,gamma,dvt
FU = Fm*Pol(1)/(kappas+Pol(1))
sou(1)=0.
sou(2)=-mus*Pol(2)+alphav*Pol(3)*FU
!print *, Pol(3)
sou(3)=lambdav*Pol(3)*(1.-Pol(3)/kv)-Pol(3)*FU
END SUBROUTINE SourceTerm

! Functions for Manufactured solutions

SUBROUTINE DUman(xr,yr,tman,DUmanX,DUmanY)
USE GeneralVariables
IMPLICIT NONE
REAL    ::  Uman,expo,xr,yr,tman,DUmanX,DUmanY,bby
!DUmanX=80.*(1.+tman)*(2.*(xr-ax))
!DUmanY=80.*(1.+tman)*(2.*(yr-ay))
bby=by
expo=exp(-xr**2-yr**2)
DUmanX=200.*(1.+tman)*(xr-ax)*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*expo+&
200.*(1.+tman)*(xr-ax)**2*(xr-bx)*(yr-ay)**2*(yr-bby)**2*expo-200.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*&
xr*expo
DUmanY=200.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)*(yr-bby)**2*expo+200.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)*&
expo-200.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*yr*expo
END SUBROUTINE DUman


FUNCTION Uman(xr,yr,tman)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,j
REAL    ::  Uman,expo,xr,yr,tman,u1
Uman=100.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-by)**2*exp(-xr**2-yr**2)
!Uman=80.*(1.+tman)*((xr-ax)**2+(yr-ay)**2)
END FUNCTION Uman

FUNCTION Wman(xr,yr,tman)
USE GeneralVariables
IMPLICIT NONE
REAL    ::  Wman,expo,xr,yr,tman 
Wman = 80.*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-by)**2*exp(-3.*xr**2-3.*yr**2)
!Wman=100.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-by)**2*exp(-xr**2-yr**2)
END FUNCTION Wman

FUNCTION Vman(xr,yr,tman)
USE GeneralVariables
IMPLICIT NONE
REAL    ::  Vman,expo,xr,yr,tman    
Vman = 10.*tman*(xr-ax)**2*(xr-bx)**2.*(yr-ay)**2*(yr-by)**2*exp(-2.*xr**2-3.*yr**2)
END FUNCTION Vman

SUBROUTINE SourceManufact(xL,yD,ForcingTerm)
USE GeneralVariables
IMPLICIT NONE
INTEGER ::  i,j,k
REAL    :: dux,duy,d2ux,d2uy,dwx,dwy,d2wx,d2wy,dut,dvt,dwt,expo,expo1,ForcingTerm(neq)
REAL    ::  bby,t1,t2,t3,t5,t6,t7,t9,t11,t12,t15,t19,t30,t33,tu,tv,tw 
REAL ::   t4,t10,t8,t13,t14,t16,t17,t18,t20,t21,t23,t24,t26,t27,t28,t34,t37
REAL    :: Uman,Vman,Wman,xL,xxig(9),wwg(9),wg1,wg2,xr,yr,yD,Vx,Vy,gama(9),gi(9),time,tman
	   gi(1)=-0.9681602395076261
	   gama(1)=0.08127438836157451
	   gi(2)=-0.8360311073266358
	   gama(2)=0.1806481606948574
	   gi(3)=-0.6133714327005904
	   gama(3)=0.2606106964029359
	   gi(4)=-0.3242534234038089
	   gama(4)=0.3123470770400022
	   gi(5)=0.
	   gama(5)=0.3302393550012596
	   gi(6)=-gi(4)
	   gama(6)=gama(4)
	   gi(7)=-gi(3)
	   gama(7)=gama(3)
	   gi(8)=-gi(2)
	   gama(8)=gama(2)
	   gi(9)=-gi(1)
	   gama(9)=gama(1)

!Gauss quadrature (0, 1)
DO k=1, 9
	xxig(k)=1./2.+1./2.*gi(k)
	wwg(k)=1./2.*gama(k)
END DO
                  time=t-dt 
                  tman=t-dt


ForcingTerm(:)=0.
DO i = 1, 9
   xr   = xL+dx*xxig(i)
   wg1  = wwg(i)*dx
   DO j = 1, 9
              
              wg2 = wwg(j)*dy
              yr  = yD+dy*xxig(j) 
              expo=exp(-xr**2-yr**2)  
              expo1=exp(-3.*xr**2-3.*yr**2)
                  dut = 100.*xr**2*(xr-1.)**2*yr**2*(yr-1.)**2*expo
                  dux = 200.*(1.+tman)*xr*(xr-1.)**2*yr**2*(yr-1.)**2*expo+200.*(1.+tman)*xr**2*(xr-1.)*yr**2*(yr-1.)**2*expo-200.*(1.+tman)*xr**3*(xr-1.)**2*yr**2*(yr-1)**2*expo
                  duy = 200.*(1.+tman)*xr**2*(xr-1.)**2*yr*(yr-1.)**2*expo+200.*(1.+tman)*xr**2*(xr-1.)**2*yr**2*(yr-1)*expo-200.*(1.+tman)*xr**2*(xr-1.)**2*yr**3*(yr-1.)**2*expo
                  t1 = 1.+tman
                  t2 = xr-1.
                  t3 = t2**2
                  t5 = yr**2
                  t7 = (yr-1.)**2
                  t9 = xr**2
                  t11 = exp(-t9-t5)
                  t12 = t5*t7*t11
                  t19 = 100.*t1*t9
                  t30 = t9**2
                  d2ux = -800.*t1*t9*xr*t2*t12+800.*t1*xr*t2*t12+400.*t1*t30*t3*t12+200.*t1*t3*t12-10.*t19*t3*t12+2.*t19*t12
                  t2 = xr**2
                  t3 = 100.*(1.+tman)*t2
                  t5 = (xr-1.)**2
                  t6 = yr-1.
                  t7 = t6**2
                  t9 = yr**2
                  t11 = exp(-t2-t9)
                  t15 = t3*t5
                  t33 = t9**2
                  d2uy = -8.*t15*t9*yr*t6*t11+4.*t15*t33*t7*t11+8.*t15*yr*t6*t11-10.*t15*t9*t7*t11+2.*t3*t5*t7*t11+2.*t3*t5*t9*t11
          
                  bby=by
                  dwt = 80.D0*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)
                  !dvt = 10.*xr**2*(xr-1.)**2*yr**2*(yr-1.)**2*exp(-2.*xr**2-3.*yr**2)
                  dvt = 10.*(xr-ax)**2*(xr-bx)**2.*(yr-ay)**2*(yr-by)**2*exp(-2.*xr**2-3.*yr**2)
                  dwx = 160.*tman*xr*(xr-1.)**2*yr**2*(yr-1.)**2*expo1+160.*tman*xr**2*(xr-1.)*yr**2*(yr-1.)**2*expo1-480.*tman*xr**3*(xr-1.)**2*yr**2*(yr-1.)**2*expo1
                  dwy = 160.*tman*xr**2*(xr-1.)**2*yr*(yr-1.)**2*expo1+160.*tman*xr**2*(xr-1.)**2*yr**2*(yr-1.)*expo1-480.*tman*xr**2*(xr-1.)**2*yr**3*(yr-1.)**2*expo1
                  d2wx = 160.D0*tman*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)+640.D0*tman*(xr-ax)*(xr-bx)*(yr-ay)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)-1920.D0*tman*(xr-ax)*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*xr*exp(-3.D0*xr**2-3.D0*yr**2)+160.D0*tman*(xr-ax)**2*(yr-ay)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)-1920.D0*tman*(xr-ax)**2*(xr-bx)*(yr-ay)**2*(yr-bby)**2*xr*exp(-3.D0*xr**2-3.D0*yr**2)-480.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)+2880.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*xr**2*exp(-3.D0*xr**2-3.D0*yr**2)

                  d2wy = 160.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)+640.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)*(yr-bby)*exp(-3.D0*xr**2-3.D0*yr**2)-1920.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)*(yr-bby)**2*yr*exp(-3.D0*xr**2-3.D0*yr**2)+160.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*exp(-3.D0*xr**2-3.D0*yr**2)-1920.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)*yr*exp(-3.D0*xr**2-3.D0*yr**2)-480.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*exp(-3.D0*xr**2-3.D0*yr**2)+2880.D0*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-bby)**2*yr**2*exp(-3.D0*xr**2-3.D0*yr**2)
                  tw=dwt-d_w*(d2wx+d2wy)+mus*Wman(xr,yr,tman)-alphav*Vman(xr,yr,tman)*(Fm*Uman(xr,yr,tman))/(kappas+Uman(xr,yr,tman))           
!                  tw=0.
                  tu=dut-d_u*(d2ux+d2uy)+Xi*(dux*dwx+Uman(xr,yr,tman)*d2wx+duy*dwy+Uman(xr,yr,tman)*d2wy)
!                  tu=dut-d_u*(d2ux+d2uy)+Xi*(dux+duy)
!                 tu=dut-d_u*(d2ux+d2uy)+Xi*dux+Xi*duy
                  tv = dvt-lambdav*Vman(xr,yr,tman)*(1.-Vman(xr,yr,tman)/kv)+&
                    Vman(xr,yr,tman)*(Fm*Uman(xr,yr,tman))/(kappas+Uman(xr,yr,tman))

              tu = wg1*wg2*tu
              ForcingTerm(1)=ForcingTerm(1)+tu
              tw = wg1*wg2*tw
              ForcingTerm(2)=ForcingTerm(2)+tw
              tv = wg1*wg2*tv
              ForcingTerm(3)=ForcingTerm(3)+tv     
   END DO
END DO
DO i = 1, neq
   ForcingTerm(i)=ForcingTerm(i)/(dx*dy)
END DO
END SUBROUTINE SourceManufact
              

subroutine Gauss_quadrature (NP1)
USE GeneralVariables
implicit none
parameter m=9
integer k,NP1
REAL    :: cg(m),gi(m),gama(m)


! Gauss quadrature (-1, 1)
if (NP1==1) then
		gi(1)=0.
		gama(1)=2.
end if
if (NP1==2) then
		gi(1)=-0.5773502692
		gama(1)=1.
		gi(2)=-gi(1)
		gama(2)=gama(1)
end if
if (NP1==3) then
	   gi(1)=-0.7745966692
	   gama(1)=0.5555555556
	   gi(2)=0.
	   gama(2)=0.8888888889
	   gi(3)=-gi(1)
	   gama(3)=gama(1)
end if
if (NP1==4) then
	   gi(1)=-0.8611361159
	   gama(1)=0.3478548451
	   gi(2)=-0.3399810436
	   gama(2)=0.6521451549
	   gi(3)=-gi(2)
	   gama(3)=gama(2)
	   gi(4)=-gi(1)
	   gama(4)=gama(1)
end if
if (NP1==8) then
	   gi(1)=-0.9602898564975362
	   gama(1)=0.1012285362903763
	   gi(2)=-0.7966664774136267
	   gama(2)=0.2223810344533743
	   gi(3)=-0.5255324099163290
	   gama(3)=0.3137066458778873
	   gi(4)=-0.1834346424956498
	   gama(4)=0.3626837833783625
	   gi(5)=-gi(4)
	   gama(5)=gama(4)
	   gi(6)=-gi(3)
	   gama(6)=gama(3)
	   gi(7)=-gi(2)
	   gama(7)=gAma(2)
	   gi(8)=-gi(1)
	   gama(8)=gama(1)
end if
if (NP1==9) then
	   gi(1)=-0.9681602395076261
	   gama(1)=0.08127438836157451
	   gi(2)=-0.8360311073266358
	   gama(2)=0.1806481606948574
	   gi(3)=-0.6133714327005904
	   gama(3)=0.2606106964029359
	   gi(4)=-0.3242534234038089
	   gama(4)=0.3123470770400022
	   gi(5)=0.
	   gama(5)=0.3302393550012596
	   gi(6)=-gi(4)
	   gama(6)=gama(4)
	   gi(7)=-gi(3)
	   gama(7)=gama(3)
	   gi(8)=-gi(2)
	   gama(8)=gama(2)
	   gi(9)=-gi(1)
	   gama(9)=gama(1)
end if
!Gauss quadrature (0, 1)
DO k=1, NP1
	xig(k)=1./2.+1./2.*gi(k)
	wg(k)=1./2.*gama(k)
END DO
END SUBROUTINE Gauss_quadrature

SUBROUTINE MALLADO
USE GeneralVariables
IMPLICIT NONE
INTEGER i


Ncells = (nx-1)*(ny-1)
dx = (bx-ax)/(nx-1)
dy = (by-ay)/(ny-1)
DO i = 1, nx
	x(i)=ax+(i-1)*dx
END DO
DO i = 1, ny
	y(i)=ay+(i-1)*dy
END DO

RETURN
END



SUBROUTINE OpL(coeff,NP1,ua,vprey,Pol,DPolx,DPoly)
USE GeneralVariables
IMPLICIT NONE
INTEGER NP1,i,j,k,l,m,kappa,iesq,kk
REAL    :: alpha,xL,xR,yD,yU,ua(neq,-2:nx+2,-2:ny+2),coeff(neq,9),Vx,Vy,Difusx,Difusy,slope1(2),slope2(2),slopex,slopey
REAL    :: sumfk,sumgk,sumsk,xg,yg,wg1,wg2,Pol(neq),DPolx(neq),DPoly(neq), vprey(-2:nx+2,-2:ny+2,NP1,NP1), preyy
REAL    :: ForcingTerm(neq)
REAL    :: FLF(neq),GLF(neq),ULW(neq),FLW(neq),GLW(neq),FLU(neq),u00(neq),DUmanX,DUmanY,Uman,Wman
REAL    :: vall(neq,-1:nx+2),valr(neq,-1:nx+2),dvall(neq,-1:nx+2),dvalr(neq,-1:nx+2),SourceGauss(neq,0:nx+1,2),dddf,dddg,derw
REAL    :: um1(neq),um2(neq),ui(neq),u1(neq),u2(neq),d1,slope(2),DDD(2),Pol1(neq),sou(neq)
REAL    :: Pol1L(neq,-2:nx+2,-2:ny+2),Pol1R(neq,-2:nx+2,-2:ny+2),Pol1D(neq,-2:nx+2,-2:ny+2),Pol1U(neq,-2:nx+2,-2:ny+2),f(neq),di(2)
REAL,ALLOCATABLE    :: fk(:),gk(:),Dfk(:),Dgk(:),sk(:),QQQ(:)
ALLOCATE(fk(neq),gk(neq),Dfk(neq),Dgk(neq),sk(neq),QQQ(neq))
alpha=1.
kappa=1.4
iesq=1
coef=0.
DDD(1)=d_u; DDD(2)=d_w
IF (iesq==0) THEN
    fL=0.; fR=0.
    DO i = 0, nx
	    DO j = 0, ny
	        slope1=0; slope2=0; slope=0
            slope1(:)=(ua(:,i+1,j)-ua(:,i,j))/dx   
            slope2(:)=(ua(:,i,j)-ua(:,i-1,j))/dx
            DO k = 1, 2
                IF (abs(slope1(k)) > abs(slope2(k))) THEN
                    slope(k)=slope2(k)
                ELSE
                    slope(k)=slope1(k)
                END IF
            END DO
            fDL(:,i,j,1)=DDD(:)*slope(:)
            fDR(:,i,j,1)=DDD(:)*slope(:)
           
            Pol1L(:,i,j)=ua(:,i,j)-slope(:)*dx/2.
            Pol1R(:,i,j)=ua(:,i,j)+slope(:)*dx/2.
            CALL flux_hyperb(x(i),y(j),slope(2)*dx,Pol1L(1,i,j),f)
            fL(1,i,j,1)=f(1)
            CALL flux_hyperb(x(i),y(j),slope(2)*dx,Pol1R(1,i,j),f)
            fR(1,i,j,1)=f(1)
            
!            f(1) = Xi*Deriv*Q(1)/dx
!!f(1) = Xi*Q(1)
!!print *, 'der',f(1),Deriv
!f(2) = 0.
!f(3) = 0.
            slope1=0; slope2=0; slope=0
            slope1(:)=(ua(:,i,j+1)-ua(:,i,j))/dy   
            slope2(:)=(ua(:,i,j)-ua(:,i,j-1))/dy
            
            DO k = 1, 2
                IF (abs(slope1(k)) > abs(slope2(k))) THEN
                    slope(k)=slope2(k)
                ELSE
                    slope(k)=slope1(k)
                END IF
            END DO

!            GDFLUX(1,i,j)=d_u*slope(1)
!            GDFLUX(2,i,j)=d_w*slope(2)
            gDD(:,i,j,1)=DDD(:)*slope(:)
            gDU(:,i,j,1)=DDD(:)*slope(:)
            Pol1D(:,i,j)=ua(:,i,j)-slope(:)*dy/2.
            Pol1U(:,i,j)=ua(:,i,j)+slope(:)*dy/2.
            CALL flux_hyperb(x(i),y(j),slope(2)*dy,Pol1D(1,i,j),f)
            gD(1,i,j,1)=f(1)
            CALL flux_hyperb(x(i),y(j),slope(2)*dy,Pol1U(1,i,j),f)
            gU(1,i,j,1)=f(1)
 
            Pol1(:)=ua(:,i,j)
            CALL SourceTerm(0,Pol1,sou)
            SSource(:,i,j)=sou(:)
        END DO
    END DO
    DO i = 1, nx
	    DO j = 1, ny
			FDFLUX(:,i,j)=(fDL(:,i,j,1)+fDR(:,i-1,j,1))/2. 
			GDFLUX(:,i,j)=(gDD(:,i,j,1)+gDU(:,i,j-1,1))/2.
            dddf = FDFLUX(2,i,j)/d_w ! (DPolL_x(2,i,j,k)+DPolR_x(2,i-1,j,k))/2.
!                dddf = 1.
            IF (dddf*Xi>0.) THEN
		        FFlux(:,i,j)=fR(:,i-1,j,1)
		    ELSE
		        FFlux(:,i,j)=fL(:,i,j,1)
		    END IF
		    
             ! (DPolL_x(2,i,j,k)+DPolR_x(2,i-1,j,k))/2.
            !FFlux(:,i,j)=0.5*(fR(:,i,j-1,1)+fL(:,i,j,1) - abs(Xi*dddf)* ( Pol1L(:,i,j) - Pol1R(:,i,j-1) ) )
!                dddf = 1.
            dddf = GDFLUX(2,i,j)/d_w

            IF (dddf*Xi>0.) THEN
		        GFlux(:,i,j)=gU(:,i,j-1,1)
		    ELSE
		        GFlux(:,i,j)=gD(:,i,j,1)
		    END IF
            !GFlux(:,i,j)=0.5*(gU(:,i,j-1,1)+gD(:,i,j,1) - abs(Xi*dddf)* ( Pol1D(:,i,j) - Pol1U(:,i,j-1) ) )      

        END DO
    END DO
ELSE
    CALL ReconsWenoDimDim(ua)
    fL=0.
    fR=0.
    gU=0.
    gD=0.
    fDL=0.
    fDR=0.
    gDD=0.
    gDU=0.
    fK=0.
    gk=0.
    sk=0.
    Dfk=0.
    Dgk=0.
    source=0.
    coeff(:,:)=0.
    DO i = 0, nx
	    DO j = 0, ny
    !	           print *, 'Volume',i,j

            DO k = 1,9
                coeff(:,k) = coef(:,k,i,j)   
            END DO
            xg = 1.
		    DO k = 1, NP1
		       yg = xig(k)
               CALL PolRecons(0,2,xg,yg,coeff,Pol,DPolx)
               d1=Dpolx(2)
               CALL flux_hyperb(xg,yg,d1,Pol,fk)
 		       fR(:,i,j,k) = fk(:)
               Dfk=0.
               Dgk=0.           
               CALL flux_diffusive(0,xg,yg,Pol,DPolx,Dfk)
!               DO kk = 1, nx
!                  CALL DUman(x(kk),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
!                  print*, 'Vol',i,j
!                  print*, kk,Dpolx(1)/dx,DUmanX
!               END DO
!               pause 333
!               !----------------------------------------------------------
!               print *, 'SOLUTION EAST SIDE OF CONTROL VOLUME'
!               print *, k,Pol(1),Uman(x(i+1),y(j)+ xig(k)*dy,t-dt)
!               print *, 'DERIVATIVES EAST SIDE OF CONTROL VOLUME'
!               CALL DUman(x(i+1),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
!               print *, 'x-derivative'
!               print *, k,Dpolx(1)/dx,DUmanX
!               !------------------------------------------------------------  
!               
!    print *, k,Pol(1),Uman(x(i+1),y(j)+ xig(k)*dy,t-dt)    
               fDR(:,i,j,k) = Dfk(:)/dx

		    END DO
    !		pause 222
    !		pause
		    xg=0.
		    DO k = 1, NP1
		       yg=xig(k)
               CALL PolRecons(0,4,xg,yg,coeff,Pol,DPolx)
               DPolR_x(:,i,j,k)=DPolx(:)
                    d1=Dpolx(2)
             
                CALL flux_hyperb(xg,yg,d1,Pol,fk)
		       fL(:,i,j,k) = fk(:)
               Dfk=0.
               Dgk=0.
               CALL flux_diffusive(0,xg,yg,Pol,DPolx,Dfk)
               !----------------------------------------------------------
!               print *, 'SOLUTION WEST SIDE OF CONTROL VOLUME'
!               print *, k,Pol(1),Uman(x(i),y(j)+ xig(k)*dy,t-dt)
!               print *, 'DERIVATIVES WEST SIDE OF CONTROL VOLUME'
!               CALL DUman(x(i),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
!               print *, 'x-derivative'
!               print *, k,Dpolx(1)/dx,DUmanX
               !------------------------------------------------------------ 
               fDL(:,i,j,k) = Dfk(:)/dx
    !                      print *, i,j,fDL(1,i,j,k)
               PolR(:,i,j,k)=Pol(:)
               PolL(:,i,j,k)=Pol(:)
               DPolL_x(:,i,j,k)=DPolx(:)
		    END DO
!		    pause
		    yg=0.
		    DO k = 1, NP1
		       xg=xig(k)
               CALL PolRecons(1,1,xg,yg,coeff,Pol,DPolx)
                    d1=Dpolx(2)
            
                CALL flux_hyperb(xg,yg,d1,Pol,fk)
		       gD(:,i,j,k) = fk(:)
               Dfk=0.
               Dgk=0.
               CALL flux_diffusive(1,xg,yg,Pol,DPolx,Dfk)
               !----------------------------------------------------------
!               print *, 'SOLUTION SOUTH SIDE OF CONTROL VOLUME'
!               print *, k,Pol(1),Uman(x(i)+ xig(k)*dx,y(j),t-dt)
!               print *, 'DERIVATIVES SOUTH SIDE OF CONTROL VOLUME'
!               CALL DUman(x(i)+ xig(k)*dx,y(j),t-dt,DUmanX,DUmanY)
!               print *, 'x-derivative'
!               print *, k,Dpolx(1)/dy,DUmanY
!               !------------------------------------------------------------
!                CALL DUman(x(i)+ xig(k)*dx,y(j),t-dt,DUmanX,DUmanY)
!               print *, k,Dpolx(1)/dx,DUmanY
               gDD(:,i,j,k) = Dfk(:)/dy
               PolD(:,i,j,k)=Pol(:)
               DPolD_y(:,i,j,k)=DPolx(:)
		    END DO
!    		pause
		    yg=1.
		    DO k = 1, NP1
		       xg=xig(k)
               CALL PolRecons(1,3,xg,yg,coeff,Pol,DPolx)
               d1=Dpolx(2)       
               CALL flux_hyperb(xg,yg,d1,Pol,fk)
		       gU(:,i,j,k) = fk(:)
               Dfk=0.
               Dgk=0.
               CALL flux_diffusive(1,xg,yg,Pol,DPolx,Dfk)
               !----------------------------------------------------------
!               print *, 'SOLUTION NORTH SIDE OF CONTROL VOLUME'
!               print *, k,Pol(1),Uman(x(i)+ xig(k)*dx,y(j+1),t-dt)
!               print *, 'DERIVATIVES NORTH SIDE OF CONTROL VOLUME'
!               CALL DUman(x(i)+ xig(k)*dx,y(j+1),t-dt,DUmanX,DUmanY)
!               print *, 'y-derivative'
!               print *, k,Dpolx(1)/dy,DUmanY
               !------------------------------------------------------------           
               gDU(:,i,j,k) = Dfk(:)/dy
               PolU(:,i,j,k)=Pol(:)
               DPolU_y(:,i,j,k)=DPolx(:)
		    END DO
    !			           print *, 'Volume',i,j

!    pause 6666
          END DO
	    END DO
    !	print *, FDR(1,1,29,1)
    !	CALL DUman(x(2),y(29)+ xig(1)*dy,t-dt,DUmanX,DUmanY)
    !	print *, DUmanx*d_u
    	
    FFO=0.
    GFO=0.
    ! Intercell Diffusive Flux
    FDifus=0.
    GDifus=0.
    CALL Gauss_quadrature (NP1)
    !print *, 'valizq malo',FDL(1,19,29,1)
    !print *, 'valder malo',FDR(1,20,29,1)
    CALL ReconsFluxDifus(NP1,ua)
    CALL Forceold(ua,NP1)
    DO i = 0, nx
	    DO j = 1, ny
         DO l = 1,neq
	        sumfk = 0.
		    DO k = 1, NP1
                wg1=wg(k)
			    sumfk=sumfk + wg1*FFO(l,i,j,k)
		    END DO
            FFLUX(l,i,j)=sumfk
         END DO
         !FFLUX(1,i,j)=0.
         FFLUX(2,i,j)=0.
         FFLUX(3,i,j)=0.
	    END DO
    END DO
    DO i = 0, nx
	    DO j = 0, ny
          DO l = 1,neq
		    sumgk = 0.
		    DO k = 1, NP1
                wg2=wg(k)
			    sumgk=sumgk + wg2*GFO(l,i,j,k)
		    END DO
		    GFLUX(l,i,j)=sumgk
          END DO
          !GFLUX(1,i,j)=0.
          GFLUX(2,i,j)=0.
          GFLUX(3,i,j)=0.
	    END DO
    END DO



    CALL Gauss_quadrature (NP1)

    ! Source terms
    DO i = 1, nx
       DO j = 1, ny
          DO k = 1,ngauss*ngauss
             coeff(:,k) = coef(:,k,i,j)
          END DO
          DO k = 1, NP1
             xg = xig(k)  
             DO l = 1, NP1
                yg = xig(l)
                CALL PolRecons(0,5,xg,yg,coeff,Pol,DPolx)  
    !            print *, Pol
                CALL SourceTerm(j,Pol,sk)
 		        source(:,i,j,k,l) = sk(:)
		     END DO
          END DO
       END DO
    END DO
    SSource = 0.
    sumsk = 0.
    DO i = 0, nx
       DO j = 1, ny
          DO l = 1,neq
             sumsk = 0.
             DO k = 1, NP1
                wg1 = wg(k)
                DO m = 1, NP1
                    wg2 = wg(m)
                    sumsk = sumsk + wg1*wg2*source(l,i,j,k,m)
                END DO
            END DO
            SSource(l,i,j)=sumsk
          END DO
       END DO
    END DO
END IF
LO=0.
!di(1)=d_u; di(2)=d_w
!DO i = 1, nx-1
!	DO j = 1, ny-1
!!!FDFLUX(:,i,j)=di(:)*(ua(:,i,j)-ua(:,i-1,j))/dx
!!!FDFLUX(:,i+1,j)=di(:)*(ua(:,i+1,j)-ua(:,i,j))/dx	
!!!GDFLUX(:,i,j)=di(:)*(ua(:,i,j)-ua(:,i,j-1))/dy
!!!GDFLUX(:,i,j+1)=di(:)*(ua(:,i,j+1)-ua(:,i,j))/dy	
!        derw=(ua(2,i,j)-ua(2,i-1,j))/dx
!        !derw=1.
!        IF (Xi*derw>0.) THEN
!            FFLUX(:,i,j)=Xi*derw*ua(:,i-1,j)
!        ELSE
!            FFLUX(:,i,j)=Xi*derw*ua(:,i,j)
!        END IF
!        derw=(ua(2,i+1,j)-ua(2,i,j))/dx
!        !derw=1.
!        IF (Xi*derw>0.) THEN
!            FFLUX(:,i+1,j)=Xi*derw*ua(:,i,j)
!        ELSE
!            FFLUX(:,i+1,j)=Xi*derw*ua(:,i+1,j)
!        END IF
!        derw=(ua(2,i,j)-ua(2,i,j-1))/dy
!        !derw=1.
!        IF (Xi*derw>0.) THEN
!            GFLUX(:,i,j)=Xi*derw*ua(:,i,j-1)
!        ELSE
!            GFLUX(:,i,j)=Xi*derw*ua(:,i,j)
!        END IF
!        derw=(ua(2,i,j+1)-ua(2,i,j))/dy
!       ! derw=1.
!        IF (Xi*derw>0.) THEN
!            GFLUX(:,i,j+1)=Xi*derw*ua(:,i,j)
!        ELSE
!            GFLUX(:,i,j+1)=Xi*derw*ua(:,i,j+1)
!        END IF
!    END DO
!END DO
      !SSource=0.
DO i = 1, nx-1
	DO j = 1, ny-1
        IF (WithManufSol==0) THEN
            ForcingTerm=0.
        ELSE
            CALL SourceManufact(x(i),y(j),ForcingTerm)
        END IF
        SSource(:,i,j)=SSource(:,i,j)+ForcingTerm(:)   
        LO(:,i,j) = (FFLUX(:,i,j)-FFLUX(:,i+1,j))/dx + (GFLUX(:,i,j)-GFLUX(:,i,j+1))/dy
        LO(:,i,j) =LO(:,i,j)-(GDFLUX(:,i,j)-GDFLUX(:,i,j+1))/dy&
         - (FDFLUX(:,i,j)-FDFLUX(:,i+1,j))/dx   +  SSource(:,i,j)
    END DO
END DO
END SUBROUTINE OpL




SUBROUTINE BoundaryConditions(ua)
USE GeneralVariables
IMPLICIT NONE
INTEGER iBC(4),i,j,k,kk
REAL    :: ua(neq,-2:nx+2,-2:ny+2),cdirE(4),cdirW(4),cdirN(4),cdirS(4),eps,gA,time
REAL    :: u1,u2,um,um1,um2,um3,um4
! 1: South; 2: East; 3: North; 4: West

DO i = 1, 4
    iBC(i)=0 !No flow
!    iBC(i)=-1 !Periodic
END DO

IF (iBC(2)==-1) THEN ! Periodic in East side of the domain
    DO j = -2, ny+2
        ua(:,nx,j)     = ua(:,1,j)
        ua(:,nx+1,j)   = ua(:,2,j)
        ua(:,nx+2,j)   = ua(:,3,j)
    END DO
ELSE IF (iBC(2)==0) THEN ! No flow in East side of the domain
    DO j = -2, ny+2
        ua(:,nx,j)     = ua(:,nx-1,j)
        ua(:,nx+1,j)   = ua(:,nx-2,j)
        ua(:,nx+2,j)   = ua(:,nx-3,j)
    END DO
ELSE IF (iBC(2)==1) THEN  !Dirichlet
    DO j = -2, ny+2
        ua(:,nx,j)     = 2.*cdirE(:)-ua(:,nx-1,j)
        ua(:,nx+1,j)   = 2.*cdirE(:)-ua(:,nx-2,j)
        ua(:,nx+2,j)   = 2.*cdirE(:)-ua(:,nx-3,j)
    END DO
ELSE IF (iBC(2)==3) THEN  !Solid walls
    DO j = -2, ny+2
        ua(1,nx,j)     = ua(1,nx-1,j)
        ua(1,nx+1,j)   = ua(1,nx-2,j)
        ua(1,nx+2,j)   = ua(1,nx-3,j)
    END DO
    DO j = -2, ny+2
        ua(2,nx,j)     = -ua(2,nx-1,j)
        ua(2,nx+1,j)   = -ua(2,nx-2,j)
        ua(2,nx+2,j)   = -ua(2,nx-3,j)
    END DO
    DO j = -2, ny+2
        ua(3,nx,j)     = ua(3,nx-1,j)
        ua(3,nx+1,j)   = ua(3,nx-2,j)
        ua(3,nx+2,j)   = ua(3,nx-3,j)
    END DO
    DO j = -2, ny+2
        ua(4,nx,j)     = ua(4,nx-1,j)
        ua(4,nx+1,j)   = ua(4,nx-2,j)
        ua(4,nx+2,j)   = ua(4,nx-3,j)
    END DO
ELSEIF (iBC(2)==4) THEN ! Absorbing
    DO j = -2, ny+2
        ua(:,nx,j)     = ua(:,nx-1,j)
        ua(:,nx+1,j)   = ua(:,nx-1,j)
        ua(:,nx+2,j)   = ua(:,nx-1,j)
    END DO

END IF
IF (iBC(4)==-1) THEN ! Periodic in West side of the domain
    DO j = -2, ny+2
        ua(:,0,j)     = ua(:,nx-1,j)
        ua(:,-1,j)    = ua(:,nx-2,j)
        ua(:,-2,j)    = ua(:,nx-3,j)
    END DO
ELSE IF (iBC(4)==0) THEN ! No flow in West side of the domain
    DO j = -2, ny+2
        ua(:,0,j)     = ua(:,1,j)
        ua(:,-1,j)    = ua(:,2,j)
        ua(:,-2,j)    = ua(:,3,j)
    END DO
ELSE IF (iBC(4)==1) THEN ! Dirichlet
!print *, ua(:,1,ny+2)
!pause
        DO j = -2, ny+2
           ua(:,0,j)     = 2.*cdirW(:)-ua(:,1,j)
           ua(:,-1,j)    = 2.*cdirW(:)-ua(:,2,j)
          ua(:,-2,j)    = 2.*cdirW(:)-ua(:,3,j)
!          print *, i,j,ua(1,-1,j)
        END DO
!       pause
ELSE IF (iBC(4)==3) THEN ! Reflecting
    DO j = -2, ny+2
        ua(1,0,j)     = ua(1,1,j)
        ua(1,-1,j)    = ua(1,2,j)
        ua(1,-2,j)    = ua(1,3,j)
    END DO
    DO j = -2, ny+2
        ua(2,0,j)     = -ua(2,1,j)
        ua(2,-1,j)    = -ua(2,2,j)
        ua(2,-2,j)    = -ua(2,3,j)
    END DO
    DO j = -2, ny+2
        ua(3,0,j)     = ua(3,1,j)
        ua(3,-1,j)    = ua(3,2,j)
        ua(3,-2,j)    = ua(3,3,j)
    END DO
    DO j = -2, ny+2
        ua(4,0,j)     = ua(4,1,j)
        ua(4,-1,j)    = ua(4,2,j)
        ua(4,-2,j)    = ua(4,3,j)
    END DO
ELSE IF (iBC(4)==4) THEN ! Absorbing
    DO j = -2, ny+2
        ua(:,0,j)     = ua(:,1,j)
        ua(:,-1,j)    = ua(:,1,j)
        ua(:,-2,j)    = ua(:,1,j)
    END DO

END IF
IF (iBC(1)==-1) THEN ! Periodic in South side of the domain
    DO i = -2, nx+2
        ua(:,i,0)     = ua(:,i,ny-1)
        ua(:,i,-1)    = ua(:,i,ny-2)
        ua(:,i,-2)    = ua(:,i,ny-3)
    END DO
ELSE IF (iBC(1)==0) THEN ! No flow in South side of the domain
    DO i = -2, nx+2
        ua(:,i,0)     = ua(:,i,1)
        ua(:,i,-1)    = ua(:,i,2)
        ua(:,i,-2)    = ua(:,i,3)
    END DO
ELSE IF (iBC(1)==1) THEN  !Dirichlet
    DO i = -2, nx+2
        ua(:,i,0)     = 2.*cdirS(:)-ua(:,i,1)
        ua(:,i,-1)   = 2.*cdirS(:)-ua(:,i,2)
        ua(:,i,-2)   = 2.*cdirS(:)-ua(:,i,3)
    END DO
ELSE IF (iBC(1)==2) THEN ! Periodic
    DO i = -2, nx+2
        ua(:,i,0)     = ua(:,i,nx-1)
        ua(:,i,-1)    = ua(:,i,nx-2)
        ua(:,i,-2)    = ua(:,i,nx-3)
    END DO

END IF
IF (iBC(3)==-1) THEN ! Periodic in North side of the domain
    DO i = -2, nx+2
        ua(:,i,ny)    = ua(:,i,1)
        ua(:,i,ny+1)  = ua(:,i,2)
        ua(:,i,ny+2)  = ua(:,i,3)
    END DO
ELSE IF (iBC(3)==0) THEN ! No flow in North side of the domain
    DO i = -2, nx+2
        ua(:,i,ny)    = ua(:,i,ny-1)
        ua(:,i,ny+1)  = ua(:,i,ny-2)
        ua(:,i,ny+2)  = ua(:,i,ny-3)
    END DO
ELSE IF (iBC(3)==1) THEN  !Dirichlet
    DO i = -2, nx+2
        ua(:,i,ny)     = 2.*cdirN(:)-ua(:,i,ny-1)
        ua(:,i,ny+1)   = 2.*cdirN(:)-ua(:,i,ny-2)
        ua(:,i,ny+2)   = 2.*cdirN(:)-ua(:,i,ny-3)
    END DO
ELSE IF (iBC(3)==2) THEN ! Periodic
    DO i = -2, nx+2
        ua(:,i,ny)    = ua(:,i,1)
        ua(:,i,ny+1)  = ua(:,i,2)
        ua(:,i,ny+2)  = ua(:,i,3)
    END DO
END IF
END SUBROUTINE BoundaryConditions


SUBROUTINE FORCEOld(ua,NP1)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,j,k,l,NP1
REAL    ::  fk,gk,alpha,ua(neq,-2:nx+2,-2:ny+2),ddx,ddy,dddf,dddg
REAL, ALLOCATABLE   :: FLF(:),GLF(:),ULW(:),FLW(:),GLW(:)
ALLOCATE(FLF(neq),GLF(neq),ULW(neq),FLW(neq),GLW(neq))
! SORT OF FLUX: 1 (FORCE), 2 (UPWIND), 3 (RUSANOV)
!IFLUX=1
alpha=2.
ddx=dx
ddy=dy
DO i = 1, nx
	DO j = 1, ny-1	
        DO k = 1, NP1
		    FLF(:) =(fR(:,i-1,j,k)+fL(:,i,j,k))/2.+ddx*(PolR(:,i-1,j,k)-PolL(:,i,j,k))/(2.*alpha*dt)	
			ULW(:) = (PolL(:,i,j,k)+PolR(:,i-1,j,k))/2.+dt*alpha*(fR(:,i-1,j,k)-fL(:,i,j,k))/(2.*ddx)

!            dddf=dx/d_w*(fDL(2,i-1,j,k)+fDR(2,i,j,k))/2. 
!			dddg=dy/d_w*(gDD(2,i,j-1,k)+gDU(2,i,j,k))/2.

!        dddf=1./d_w*(fDR(2,i-1,j,k)+fDL(2,i,j,k))/2. 
!			dddg=1./d_w*(gDU(2,i,j-1,k)+gDD(2,i,j,k))/2.



                dddf = (DPolL_x(2,i,j,k)+DPolR_x(2,i-1,j,k))/2.
!                dddf = 1.
            CALL flux_hyperb(1.,1.,dddf,ULW,FLW)
!			FFO(:,i,j,k)=(FLW(:)+FLF(:))/(2.)
			FFO(:,i,j,k)=FLw(:)
!            FFO(:,i,j,k)=fR(:,i-1,j,k)
!            IF (dddf*Xi>0.) THEN
!		        FFO(:,i,j,k)=fR(:,i-1,j,k)
!		    ELSE
!		        FFO(:,i,j,k)=fL(:,i,j,k)
!		    END IF
!FFO(:,i,j,k) = 0.5*(fL(:,i,j,k)+fR(:,i-1,j,k) - abs(dddf*Xi) * ( PolL(:,i,j,k) - PolR(:,i-1,j,k) ) )
FFO(:,i,j,k) = 0.5*(fL(:,i,j,k)+fR(:,i-1,j,k) - abs(Xi) * ( PolL(:,i,j,k) - PolR(:,i-1,j,k) ) )
		END DO
      !END DO
	END DO
END DO


DO i = 1, nx-1
	DO j = 1, ny
     ! DO l = 1,neq
		DO k = 1, NP1
			GLF(:)         = (gU(:,i,j-1,k)+gD(:,i,j,k))/2.+ddy*(PolU(:,i,j-1,k)-PolD(:,i,j,k))/(2.*alpha*dt)	
			ULW(:)            = (PolD(:,i,j,k)+PolU(:,i,j-1,k))/2.+dt*alpha*(gU(:,i,j-1,k)-gD(:,i,j,k))/(2.*ddy)	

            
                dddg = (DPolU_y(2,i,j-1,k)+DPolD_y(2,i,j,k))/2.
!  dddg=1.        
			CALL flux_hyperb(1.,1.,dddg,ULW,GLW) 
!			GFO(:,i,j,k)   =(Glw(:)+GLF(:))/(2.)  
			   GFO(:,i,j,k)   =GLw(:)
!            IF (dddg*Xi>0.) THEN
!		        GFO(:,i,j,k)   =gU(:,i,j-1,k)
!		    ELSE
!		        GFO(:,i,j,k)   =gD(:,i,j,k)
!		    END IF
!GFO(:,i,j,k) = 0.5*(gU(:,i,j-1,k)+gD(:,i,j,k) - abs(dddg*Xi) * ( PolD(:,i,j,k) - PolU(:,i,j-1,k) ) )
GFO(:,i,j,k) = 0.5*(gU(:,i,j-1,k)+gD(:,i,j,k) - abs(Xi)* ( PolD(:,i,j,k) - PolU(:,i,j-1,k) ) )
        END DO
    ! END DO
	END DO
END DO
!print *, ffo
!DEALLOCATE(FLF,GLF,ULW,FLW,GLW)
END SUBROUTINE FORCEOld


!SUBROUTINE ReconsFluxDifus(NP1,u0)
!USE GeneralVariables
!IMPLICIT NONE
!INTEGER i,j,k,l,NP1
!REAL    :: u0(neq,-2:nx+2,-2:ny+2),sumfk,sumgk,val,val1,theta
!!FDifus=0.
!DO l = 1,neq-1
!    DO i = 1, nx
!        DO j = 1, ny
!            sumfk=0.; sumgk=0.
!            DO k = 1, NP1
!                val1=(gDD(l,i,j,k)+gDU(l,i,j-1,k))/2.
!                sumgk=sumgk + wg(k)*val1
!            END DO
!            GDFLUX(l,i,j)=sumgk
!         END DO
!    END DO
!    DO j = 1, ny
!        DO i = 1, nx
!            sumfk=0.; sumgk=0.
!            DO k = 1, NP1
!                val=(fDR(l,i-1,j,k)+fDL(l,i,j,k))/2.
!                sumfk=sumfk + wg(k)*val
!            END DO
!            FDFLUX(l,i,j)=sumfk
!         END DO
!    END DO   
!END DO            
!            
!!            
!!            
!!DO i = 1, nx
!!	DO j = 1, ny
!!	  DO l = 1,neq-1	
!!	    sumfk=0.; sumgk=0.
!!        DO k = 1, NP1
!!			val=(fDR(l,i-1,j,k)+fDL(l,i,j,k))/2. ! -0.1*abs(PolL(l,i,j,k)-PolR(l,i-1,j,k))&
!!			val1=(gDD(l,i,j+1,k)+gDU(l,i,j,k))/2. ! -0.1*abs(PolD(l,i,j,k)-PolU(l,i,j-1,k))&
!!			    sumfk=sumfk + wg(k)*val
!!			    sumgk=sumgk + wg(k)*val1
!!		END DO
!!       
!!        FDFLUX(l,i,j)=sumfk 
!!        GDFLUX(l,i,j)=sumgk
!!      END DO
!!      FDFLUX(3,i,j)=0. 
!!      GDFLUX(3,i,j)=0.
!!	END DO
!!END DO
!END SUBROUTINE ReconsFluxDifus



SUBROUTINE ReconsFluxDifus(NP1,u0)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,j,k,l,NP1
REAL    :: u0(neq,-2:nx+2,-2:ny+2),sumfk,sumgk,val,val1,theta,dd(2)
REAL    ::  derf,x1,x2,xr,y1,y2,yr
REAL    ::  tman,t4,t5,t7,t8,t9,t12,t13,t14,t16,t17,t18,t19,t20,t21,t22,t23,t25,t26,t30,t32,t34,t35,t36,t37,t38
REAL    ::  t41,t43,t54,t55,t56,t62,t65,t79,t82,t85,t98,t105,valint,DUmanX,DUmanY,dif1,dif2

!FDifus=0.
tman=t-dt
!print *, tman
!print *, 'Valor aproximado en izqdo de VC 20'
!print *, fDR(1,19,29,1), fDL(1,20,29,1)
!valint=(fDR(1,19,29,1)+fDL(1,20,29,1))/2.
!print *, valint 
!print *, 'Valor exacto en izqdo de VC 20'
!CALL DUman(x(20),y(29)+ xig(1)*dy,t-dt,DUmanX,DUmanY)
!print *, DUmanX*d_u
!dif1=abs(DUmanX*d_u-valint)
!print *, '*******************************************'
!print *, 'VALORES QUE NO DEBERIAN SER'
!print *, fDL(1,19,29,1), fDR(1,20,29,1)
!valint=(fDL(1,19,29,1)+fDR(1,20,29,1))/2.
!print *, valint 
!print *, '*******************************************'
dif2=abs(DUmanX*d_u-valint)
!IF (dif1<dif2) THEN
!  print *, 'Gana el bueno'
!ELSE
!  print *, 'Gana el malo'
!END IF
!print *, FDR(1,1,29,1)
!pause 22
!pause 201
dd(1)=d_u; dd(2)=d_w
DO i = 1, nx
	DO j = 1, ny
	  DO l = 1,neq	
	    sumfk=0.; sumgk=0.
        DO k = 1, NP1
!        IF (l==1) print *,i,j,k,fDR(l,i-1,j,k),fDL(l,i,j,k)
!        pause 414
!         IF (l==1) THEN
!         ! Polinomio por la izquierda de x(i)
!             PRINT*,i,j,k,fDR(l,i-1,j,k)
!             PRINT*,i,j,k,fDL(l,i-1,j,k)
!         ! Polinomio por la derecha de x(i)
!             PRINT*,i,j,k,fDL(l,i,j,k)
!             PRINT*,i,j,k,fDL(l,i-1,j,k)
!             PRINT*,i,j,k,fDR(l,i,j,k)
!         ! Media
!             PRINT*, 'Media',i,j,k,(fDR(l,i-1,j,k)+fDL(l,i,j,k))/2.
!             PRINT*, 'La otra Media',i,j,k,(fDL(l,i-1,j,k)+fDR(l,i,j,k))/2.
!         ! Derivada exacta en x(i)
!             CALL DUman(x(i),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
!             PRINT*, 'val exacto',DUmanX*d_u
!         END IF
!         print*, 'ERRORES' 
!         print*, 'fDR(l,i-1,j,k)',fDR(l,i-1,j,k)-DUmanX*d_u
!         print*, 'fDL(l,i,j,k)',fDL(l,i,j,k)-DUmanX*d_u
!         print*, 'fDL(l,i-1,j,k)',fDL(l,i-1,j,k)-DUmanX*d_u
!         print*, 'fDR(l,i,j,k)',fDR(l,i,j,k)-DUmanX*d_u
!         print*, 'Media 1:',(fDR(l,i-1,j,k)+fDL(l,i,j,k))/2.-DUmanX*d_u
!         print*, 'Media 2:',(fDL(l,i-1,j,k)+fDR(l,i,j,k))/2.-DUmanX*d_u
!         PAUSE 33
		  val=(fDR(l,i-1,j,k)+fDL(l,i,j,k))/2. !  -d_u/2.*(PolL(l,i,j,k)-PolR(l,i-1,j,k)) !+2./(3.)*(PolL(l,i,j,k)-PolR(l,i-1,j,k))/dx
			
		!val=(fDL(l,i-1,j,k)+fDR(l,i,j,k))/2. 
!			val=d_u*(PolL(l,i-1,j,k)-PolR(l,i,j,k))/dx
!			val=(u0(l,i,j)-u0(l,i-1,j))/dx
!			-0.5*abs(fDR(l,i-1,j,k)-fDL(l,i,j,k))
			!+(u0(l,i,J)-u0(l,i-1,j))
			val1=(gDD(l,i,j,k)+gDU(l,i,j-1,k))/2. ! -d_u/2.*(PolD(l,i,j,k)-PolU(l,i,j-1,k)) !+2./(3.)*(PolD(l,i,j,k)-PolU(l,i,j-1,k))/dy
!			val1=(gDD(l,i,j-1,k)+gDU(l,i,j,k))/2.
!		    CALL DUman(x(i),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
!		    print *, i,j
!            print *, val,d_u*DUmanx
!			val1=d_w*(PolD(l,i,j-1,k)-PolU(l,i,j,k))/dy
!			-0.5*abs(gDD(l,i,j,k)-gDU(l,i,j-1,k)) 
!            val1=(u0(l,i,j)-u0(l,i,j-1))/dy    
!            theta=0.5
!print *, val,val1
			    sumfk=sumfk + wg(k)*val !*dx
			    sumgk=sumgk + wg(k)*val1 !*dy
!			    print *, i,j,fDL(1,i,j,k)/d_u
		END DO
!		pause
!        FDFLUX(1,i,j)=d_u*(u0(l,i,j)-u0(l,i-1,j))/dx !sumfk 
!        GDFLUX(1,i,j)=d_u*(u0(l,i,J)-u0(l,i,j-1))/dy !sumgk
!        FDFLUX(2,i,j)=d_w*(u0(l,i,j)-u0(l,i-1,j))/dx !sumfk 
!        GDFLUX(2,i,j)=d_w*(u0(l,i,J)-u0(l,i,j-1))/dy !sumgk!        
        FDFLUX(l,i,j)=sumfk
        GDFLUX(l,i,j)=sumgk
      END DO
     y1=y(j)
      y2=y(j+1)
      x1=x(i)
      x2=x(i+1)
      xr=x(i)
      yr=y(j)
      t4 = y2**2      
      t5 = exp(t4)
      t8 = y1**2
      t9 = t8*y1
      t12 = t5*xr
      t16 = xr**2
      t17 = t16*xr
      t18 = t17*t5
      t20 = exp(t8)
      t21 = t20*xr
      t23 = t4*y2
      t30 = t17*t20
      t32 = t16*t5
      t34 = t16*t20
      t37 = exp(t4+t8)
      t38 = derf(y1)
      t41 = derf(y2)
      t56 = -2000000000.D0*t5*y1-800000000.D0*t5*t9-3200000000.D0*t12+1600000000.D0*t5*t8+1600000000.D0*t18+3200000000.D0*t21+800000000.D0*t20*t23+2000000000.D0*t20*y2-1600000000.D0*t20*t4-1600000000.D0*t30-1600000000.D0*t32+1600000000.D0*t34+1772453851.D0*t37*t38-1772453851.D0*t37*t41+1600000000.D0*t12*t9+800000000.D0*t32*t9+2000000000.D0*t32*y1-1600000000.D0*t32*t8+2000000000.D0*t30*y2+4000000000.D0*t12*y1
      t79 = t37*xr
      t82 = t17*t37
      t85 = t16*t37
      t98 = -3200000000.D0*t12*t8-800000000.D0*t18*t9+1600000000.D0*t18*t8-1600000000.D0*t21*t23-4000000000.D0*t21*y2+3200000000.D0*t21*t4+800000000.D0*t30*t23-1600000000.D0*t30*t4-800000000.D0*t34*t23-2000000000.D0*t34*y2+1600000000.D0*t34*t4-3544907702.D0*t79*t38+1772453851.D0*t82*t38-1772453851.D0*t85*t38+3544907702.D0*t79*t41-1772453851.D0*t82*t41+1772453851.D0*t85*t41-2000000000.D0*t18*y1+1600000000.D0*t5-1600000000.D0*t20
      t105 = exp(-1.D0*t16-1.D0*t8-1.D0*t4)
      valint = 0.125D-6*(1.D0+tman)*xr*(xr-1.D0)*(t56+t98)*t105

      
!      print *, i,j,FDFLUX(1,i,j),d_u*valint/dy
       
      
      t4 = x2**2      
      t5 = x1**2
      t7 = exp(t4+t5)
      t8 = t7*yr
      t9 = derf(x2)
      t12 = yr**2
      t13 = t12*t7
      t14 = derf(x1)
      t17 = exp(t4)
      t18 = t17*yr
      t19 = t5*x1
      t22 = t12*t17
      t25 = t12*yr
      t26 = t25*t7
      t35 = exp(t5)
      t36 = t12*t35
      t43 = t25*t17
      t54 = t35*yr
      t55 = t4*x2
      t62 = t25*t35
      t65 = 0.1772453851D2*t8*t9-0.8862269255D1*t13*t14+8.D0*t18*t19+10.D0*t22*x1-0.8862269255D1*t26*t9+0.8862269255D1*t13*t9-0.1772453851D2*t8*t14+0.8862269255D1*t26*t14-10.D0*t36*x2+20.D0*t18*x1-16.D0*t18*t5-4.D0*t43*t19-10.D0*t43*x1+8.D0*t43*t5+4.D0*t22*t19-8.D0*t22*t5-8.D0*t54*t55-20.D0*t54*x2+16.D0*t54*t4+10.D0*t62*x2
      t98 = -8.D0*t62*t4-4.D0*t36*t55+8.D0*t36*t4+4.D0*t62*t55+0.8862269255D1*t7*t14-0.8862269255D1*t7*t9-16.D0*t18+8.D0*t36-8.D0*t22-4.D0*t17*t19-10.D0*t17*x1+8.D0*t17*t5+8.D0*t43+16.D0*t54+4.D0*t35*t55+10.D0*t35*x2-8.D0*t35*t4-8.D0*t62+8.D0*t17-8.D0*t35
      t105 = exp(-1.D0*t12-1.D0*t5-1.D0*t4)
      valint = 25.D0*(1.D0+tman)*yr*(yr-1.D0)*(t65+t98)*t105

      
      
!      print *, i,j,GDFLUX(1,i,j),d_u*valint/dx
!      pause
      FDFLUX(3,i,j)=0. 
      GDFLUX(3,i,j)=0.
	END DO
END DO
!print *, FDFLUX
END SUBROUTINE ReconsFluxDifus

SUBROUTINE SOLINIC(NP1,xL,xR,yD,yU,Q)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,NP1,j,k
REAL    :: derf,xg,yg,xL,xR,yD,yU,Q(neq),wwg,XXX,YYY,H,eepsilon
REAL    :: kappa,p,rho,ux,uy,xreal,yreal,wg1,wg2,sum1,sum2,sum3,gama(9),gi(9),sss,aa,bb,Rad,Uman

gi(1)=-0.9681602395076261
gama(1)=0.08127438836157451
gi(2)=-0.8360311073266358
gama(2)=0.1806481606948574
gi(3)=-0.6133714327005904
gama(3)=0.2606106964029359
gi(4)=-0.3242534234038089
gama(4)=0.3123470770400022
gi(5)=0.
gama(5)=0.3302393550012596
gi(6)=-gi(4)
gama(6)=gama(4)
gi(7)=-gi(3)
gama(7)=gama(3)
gi(8)=-gi(2)
gama(8)=gama(2)
gi(9)=-gi(1)
gama(9)=gama(1)
pi=acos(-1.)
sum1 = 0.
sum3 = 0.
eepsilon=5.e-2
DO i = 1, 9
   xreal   = (xL+xR)/2+dx*gi(i)/2.
   wg1  = dx*gama(i)/2.
   DO j = 1, 9
      wg2 = dy*gama(j)/2.
      yreal  = (yD+yU)/2+dy*gi(j)/2.       
      sum2 = 0.
!      sum1 = sum1 + wg1*wg2*(0.15+eepsilon*cos(-2.*pi*xreal)*cos(-2.*pi*yreal))
!      sum3 = sum3 + wg1*wg2*(0.1+eepsilon*cos(-2.*pi*xreal)*cos(-2.*pi*yreal))
      
      sum1 = sum1 + wg1*wg2*(0.15+eepsilon*cos(2.*pi*xreal)*cos(2.*pi*yreal))
      sum3 = sum3 + wg1*wg2*(0.1+eepsilon*cos(pi*xreal)*cos(pi*yreal))
!      sum1 = sum1 + wg1*wg2*(10.+cos(pi*xreal)*cos(pi*yreal))
!      sum3 = sum3 + wg1*wg2*(21.+20.*cos(pi*xreal)*cos(pi*yreal))

   END DO
END DO 
Q(1) = sum1/(dx*dy)
Q(2) = sum2/(dx*dy)
Q(3) = sum3/(dx*dy)
END SUBROUTINE SOLINIC

SUBROUTINE SolManuf(tman,xL,xR,yD,yU,Q)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,NP1,j,k
REAL    :: derf,xg,yg,xL,xR,yD,yU,Q(neq),wwg,XXX,YYY,tman
REAL    :: kappa,p,rho,ux,uy,xreal,yreal,wg1,wg2,sum1,sum2,sum3,gama(9),gi(9),Uman,Wman,Vman,SumaU,SumaV,SumaW
gi(1)=-0.9681602395076261
gama(1)=0.08127438836157451
gi(2)=-0.8360311073266358
gama(2)=0.1806481606948574
gi(3)=-0.6133714327005904
gama(3)=0.2606106964029359
gi(4)=-0.3242534234038089
gama(4)=0.3123470770400022
gi(5)=0.
gama(5)=0.3302393550012596
gi(6)=-gi(4)
gama(6)=gama(4)
gi(7)=-gi(3)
gama(7)=gama(3)
gi(8)=-gi(2)
gama(8)=gama(2)
gi(9)=-gi(1)
gama(9)=gama(1)
pi=acos(-1.)
SumaU=0.
SumaW=0.
SumaV=0.
DO i = 1, 9
   xreal   = (xL+xR)/2.+dx*gi(i)/2.
   wg1  = dx*gama(i)/2.
   DO j = 1, 9
      wg2 = dy*gama(j)/2.
      yreal  = (yD+yU)/2.+dy*gi(j)/2.       
      SumaU = SumaU + wg1*wg2*Uman(xreal,yreal,tman)
      SumaW = SumaW + wg1*wg2*Wman(xreal,yreal,tman)
      SumaV = SumaV + wg1*wg2*Vman(xreal,yreal,tman)
   END DO
END DO
Q(1) = SumaU  /(dx*dy)
Q(2) = SumaW /(dx*dy)
Q(3) = SumaV/(dx*dy)   
END SUBROUTINE SolManuf

SUBROUTINE TEST_RECONSTRUCTION
USE GeneralVariables
INTEGER i,j,k
REAL    :: poltest(neq,-2:nx+2,-2:ny+2),ipoltest(neq,-2:nx+2,-2:ny+2),mu,nu,xL,xR,yD,yU,coeff(neq,9),Pol(neq),DPolx(neq),DPoly(neq)
REAL    :: err(4),dpoltestx(neq,-2:nx+2,-2:ny+2),dpoltesty(neq,-2:nx+2,-2:ny+2)


ipoltest=0.
poltest=0.
Pol=0.
DO i = 1, nx-1
	xL=x(i)
	xR=x(i+1)
	DO j = 1, ny-1
		yD=y(j)
		yU=y(j+1)
       
        ipoltest(1,i,j) = ((xR-xL)*(yU**3-yD**3)/3.+(-xR**2+xL**2)*(yU**2-yD**2)/2+xR**3*(yU-yD)/3.-&
                          xL**3*(yU-yD)/3.)/(dx*dy)
!        ipoltest(2,i,j) = ((-3.*xR+3*xL)*(yU**3.-yD**3)/3.+(xL-xR)*(yU**2-yD**2)/2.+xR**3*(yU-yD)/3-xL**3*(yU-yD)/3.+&
!                          xR**2*(yU-yD)/2-xL**2*(yU-yD)/2)/(dx*dy)
!        ipoltest(3,i,j) = ((-4.*xR+4.*xL)*(yU**3-yD**3)/3+(5./2.*xR**2-5./2.*xL**2)*(yU**2-yD**2)/2.+&
!                          xR**3*(yU-yD)-xL**3*(yU-yD))/(dx*dy)
!        ipoltest(4,i,j) = ((8.*xR-8.*xL)*(yU**3-yD**3)/3+(-3./2.*xR**2+3./2.*xL**2-xR+xL)*(yU**2-yD**2)/2.+&
!                          2.*xR**3*(yU-yD)-2.*xL**3*(yU-yD)+xR**2*(yU-yD)/2.-xL**2*(yU-yD)/2.)/(dx*dy)
    END DO
END DO

CALL BoundaryConditions(ipoltest)
CALL ReconsWenoDimDim(ipoltest)
!print *, coef
!pause
DO i = 1, nx-1
	!mu=x(i)+dx/2.
!	mu=x(i)+dx
	DO j = 1, ny-1
		!nu=y(j)+dy/2.
!        nu=y(j)+dy/2.
        DO k = 1,ngauss*ngauss
            coeff(:,k) = coef(:,k,i,j)
        END DO
        mu=x(i)+dx; nu=y(j)+dy/2.
        CALL PolRecons(0,5,1.,0.5,coeff,Pol,DPolx)
        poltest(1,i,j) = mu**2+nu**2-2*mu*nu
        dpoltestx(1,i,j) = 2.*mu-2.*nu
        DO k = 1, neq
           err(k)=poltest(k,i,j)-Pol(k)
           print *, i,j,'err',err(k)
        END DO
         DO k = 1, neq
           err(k)=dpoltestx(k,i,j)-DPolx(k)/dx
           print *, i,j,'derrx',err(k)
        END DO
        pause 
         
        mu=x(i)+dx/2.; nu=y(j)+dy
        CALL PolRecons(1,5,0.5,1.,coeff,Pol,DPoly) 
        poltest(1,i,j) = mu**2+nu**2-2*mu*nu

!        poltest(2,i,j) = mu**2-3*nu**2+mu-nu
!        poltest(3,i,j) = 3*mu**2-4*nu**2+5*mu*nu
!        poltest(4,i,j) = 6*mu**2+8*nu**2-3*mu*nu+mu-nu
!        dpoltestx(1,i,j) = 2.*mu-2.*nu
!        dpoltestx(2,i,j) = 2*mu+1
!        dpoltestx(3,i,j) = 6.*mu+5.*nu
!        dpoltestx(4,i,j) = 12.*mu-3.*nu+1.
        dpoltesty(1,i,j) = 2.*nu-2.*mu
!        dpoltesty(2,i,j) = -6.*nu-1
!        dpoltesty(3,i,j) = -8.*nu+5.*mu
!        dpoltesty(4,i,j) = 16.*nu-3*mu-1
!        
        DO k = 1, neq
           err(k)=poltest(k,i,j)-Pol(k)
           print *, i,j,'err',err(k)
        END DO
        DO k = 1, neq
           err(k)=dpoltesty(k,i,j)-DPoly(k)/dy
           print *, i,j,'derry',err(k)
        END DO
        pause
    END DO
END DO


END SUBROUTINE TEST_RECONSTRUCTION

SUBROUTINE OUTPUT(u,VAver)
USE GeneralVariables
INTEGER  i,j
REAL    :: x1,x2,y1,y2,ux,uy,vvel,kappa,p,u(neq,-2:nx+2,-2:ny+2),Q(neq),VAver(-2:nx+2,-2:ny+2),QV
!
! FILES OPENING
!

    OPEN (1,  FILE='SolU03.DAT',STATUS='UNKNOWN')
    OPEN (2,FILE='SolW03.DAT',STATUS='UNKNOWN')
    OPEN (3,FILE='SolV03.DAT',STATUS='UNKNOWN')
    OPEN (4,FILE='UCuty05Num.DAT',STATUS='UNKNOWN')
    OPEN (41,FILE='WCuty05Num.DAT',STATUS='UNKNOWN')
    OPEN (42,FILE='VCuty05Num.DAT',STATUS='UNKNOWN')
    OPEN (43,FILE='UCutx05Num.DAT',STATUS='UNKNOWN')
    DO j = 1, ny-1
	    y1=y(j)
	    y2=y(j+1)
	    DO i = 1, nx-1
		    x1=x(i)
		    x2=x(i+1)
	        WRITE(1,*) (x1+x2)/2.,(y1+y2)/2.,u(1,i,j)
            WRITE(2,*) (x1+x2)/2.,(y1+y2)/2.,u(2,i,j)
            WRITE(3,*) (x1+x2)/2.,(y1+y2)/2.,u(3,i,j)
        END DO
    END DO
    DO i = 1, nx-1
		    x1=x(i)
		    x2=x(i+1)
	        WRITE(4,*) (x1+x2)/2.,u(1,i,ny/2)
	        WRITE(41,*) (x1+x2)/2.,u(2,i,ny/2)
	        WRITE(42,*) (x1+x2)/2.,u(3,i,ny/2)
!            WRITE(3,*) (x1+x2)/2.,(y1+y2)/2.,u(3,i,j)
    END DO
    DO i = 1, ny-1
		    y1=y(i)
		    y2=y(i+1)
	        WRITE(43,*) (y1+y2)/2.,u(1,nx/2,i)
    END DO
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    IF (WithManufSol==1) THEN
        OPEN(2,FILE='EXACT.dat',STATUS='UNKNOWN')
        OPEN(3,FILE='UEXACTy05.dat',STATUS='UNKNOWN')
        OPEN(31,FILE='WEXACT05.dat',STATUS='UNKNOWN')
        OPEN(32,FILE='VEXACT05.dat',STATUS='UNKNOWN')
        OPEN(33,FILE='UEXACTx05.dat',STATUS='UNKNOWN')
        DO i = 1, nx-1
            x1=x(i)
		    x2=x(i+1)
            DO j = 1, ny-1
                y1=y(j)
	            y2=y(j+1)
                CALL SolManuf(t,x1,x2,y1,y2,Q)
                IF (j==ny/2) THEN
                    WRITE(3,*)  (x1+x2)/2.,Q(1)
                    WRITE(31,*) (x1+x2)/2.,Q(2)
                    WRITE(32,*) (x1+x2)/2.,Q(3)
                END IF
                IF (i==nx/2) THEN
                    WRITE(33,*)  (y1+y2)/2.,Q(1)
                END IF
                WRITE(2,*) (x1+x2)/2.,(y1+y2)/2.,Q(1)
      
            END DO
        END DO
        CLOSE(2)
        CLOSE(3)
        CLOSE(31)
        CLOSE(32)
END IF
CLOSE(1)
CLOSE(2)
END SUBROUTINE OUTPUT




SUBROUTINE ReconsWenoDimDim(ua)
USE GeneralVariables
IMPLICIT NONE
INTEGER ::  i,idir,ileg,j,jleg,kleg,l
REAL    :: Pol,phi(9),ua(neq,-2:nx+2,-2:ny+2),oo
REAL,ALLOCATABLE    :: uu(:,:)
ALLOCATE(uu(neq,5))
! idir=1 Reconstruction in x-direction, idir=2 Reconstruction in y-direction
!coef=0.
!print *, ua
!pause 44
!print *, ua
WWX0=0.
WWX1=0.
WWX2=0.
DO idir = 1, 2
    IF (idir==1) THEN
      !DO l = 1,neq
        DO i = 0, nx
!                    print *, ua(1,i,ny+2)

            DO j = -2, ny+2
!            pause 444
                uu(:,1)                 = ua(:,i-2, j)
                uu(:,2)                 = ua(:,i-1, j)
                uu(:,3)                 = ua(:,i, j)
                uu(:,4)                 = ua(:,i+1, j)
                uu(:,5)                 = ua(:,i+2, j)
!                print *, 'u1',i,j,uu(:,1)
!                print *, 'u2',i,j,uu(:,2)
!                print *, 'u3',i,j,uu(:,3)
!                print *, 'u4',i,j,uu(:,4)
!                print *, 'u5',i,j,uu(:,5)
                CALL WENO1D(uu)
!print *,wweno
                WWX0(:,i, j)            = wweno(:,0)
                WWX1(:,i, j)            = wweno(:,1)
                WWX2(:,i, j)            = wweno(:,2)
!    print *, i,j,WWX1(:,i, j)
            END DO
        END DO 
!        print *, WWX1  
      !END DO
    ELSE
! Reconstruction in y-direction
        DO i = 0, nx
            DO j = 0, ny
    
              !DO l = 1,neq
		        uu(:,1)              =  WWX0(:, i, j-2)
		        uu(:,2)              =  WWX0(:, i, j-1)
		        uu(:,3)              =  WWX0(:, i, j)
		        uu(:,4)              =  WWX0(:, i, j+1)
		        uu(:,5)              =  WWX0(:, i, j+2)
 
!                PRINT *, 'uu',i,j,uu
!                pause
!print *, 'uu1',uu

		        CALL WENO1D(uu)
	
                coef(:,1,i,j)        = wweno(:,0)
		        coef(:,2,i,j)        = wweno(:,1)
		        coef(:,3,i,j)        = wweno(:,2)

                uu(:,1)              = WWX1(:, i, j-2)
				uu(:,2)              = WWX1(:, i, j-1)
				uu(:,3)              = WWX1(:, i, j)
				uu(:,4)              = WWX1(:, i, j+1)
				uu(:,5)              = WWX1(:, i, j+2)
				
!print *, 'uu2',uu
				
				CALL WENO1D(uu)
                coef(:,4,i,j)        = wweno(:,0)
				coef(:,5,i,j)        = wweno(:,1)
				coef(:,6,i,j)        = wweno(:,2)
		        uu(:,1)              = WWX2(:, i, j-2)
	            uu(:,2)              = WWX2(:, i, j-1)
		        uu(:,3)              = WWX2(:, i, j)
		        uu(:,4)              = WWX2(:, i, j+1)
		        uu(:,5)              = WWX2(:, i, j+2)
!print *, 'uu3',uu
		        
		        CALL WENO1D(uu)
		        coef(:, 7,i,j)       = wweno(:,0)
		        coef(:, 8,i,j)       = wweno(:,1)
		        coef(:, 9,i,j)       = wweno(:,2)
   
              END DO
            END DO
        !END DO
     END IF
END DO
!DEALLOCATE(uu)
END SUBROUTINE ReconsWenoDimDim

SUBROUTINE WENO1D(uu)
USE GeneralVariables
IMPLICIT NONE
INTEGER l
REAL    :: uu(neq,5),WX0(3,1),WX1(3,1),WX2(3,1),oo,u00(3,1),c1,c2
!REAL,ALLOCATABLE    :: u00(:,:,:)

	! Centred stencil OI(neq,0:2)
!ALLOCATE(oo(neq),u00(neq,3,1))
 c1=13./12.; c2=0.25

OI=0.
DO l = 1,neq
    u00(1,1)    = uu(l,2)
    u00(2,1)    = uu(l,3)
    u00(3,1)    = uu(l,4)
   
    WX0         = matmul(M1,u00)
    CALL OSCIND(WX0,oo)
!    print *, u00
    !OI(l,0)     = oo
    OI(l,0) = c1*(u00(1,1)-2.*u00(2,1)+u00(3,1))**2+c2*(u00(1,1)-u00(3,1))**2
!    print *, 'OI1',OI
! print *, l,u00
!pause 777
!END DO
	! Left-sided stencil
!DO l = 1,neq
    u00(1,1)    = uu(l,1)
    u00(2,1)    = uu(l,2)
    u00(3,1)    = uu(l,3)
    
    WX1         = matmul(ML,u00)
    CALL OSCIND(WX1,oo)
    !OI(l,1)     = oo
    OI(l,1) = c1*(u00(1,1)-2.*u00(2,1)+u00(3,1))**2+c2*(u00(1,1)-4.*u00(2,1)+3.*u00(3,1))**2
!    print *, 'OI2',OI

!END DO
	! Right-sided stencil
!DO l = 1,neq
    u00(1,1)    = uu(l,3)
    u00(2,1)    = uu(l,4)
    u00(3,1)    = uu(l,5)
    WX2         = matmul(MR,u00)
!    print *, u00
    CALL OSCIND(WX2,oo)
    !OI(l,2)     = oo
    OI(l,2) = c1*(u00(1,1)-2.*u00(2,1)+u00(3,1))**2+c2*(3.*u00(1,1)-4.*u00(2,1)+u00(3,1))**2

!    print *, 'OI3',OI

  
!end do
!DO l = 1,neq
    w(l,0)      = WX0(1,1)
    w(l,1)      = WX1(1,1)
    w(l,2)      = WX2(1,1)

    w(l,3)      = WX0(2,1)
    w(l,4)      = WX1(2,1)
    w(l,5)      = WX2(2,1)

    w(l,6)      = WX0(3,1)
    w(l,7)      = WX1(3,1)
    w(l,8)      = WX2(3,1)

END DO
CALL COEFF_WENO

END SUBROUTINE WENO1D

SUBROUTINE COEFF_WENO
USE GeneralVariables
IMPLICIT NONE
INTEGER i,l
REAL    :: suma
DO l = 1,neq
    suma=0.
    DO i = 0, 2
!    print *, lambda(i),OI(l,i),epsilon,r
	    omegah(l,i)=lambda(i)/((OI(l,i)+epsilon)**r)
	    suma = suma + omegah(l,i)
    END DO

    DO i = 0, 2
        omegas(l,i) = omegah(l,i)/suma
    END DO
END DO
!omegas(:,0)=1.; 
!omegas(:,1)=0.; 
!omegas(:,2)=0.;
!do l =1,neq
!    print *, 'cero',l,omegas(:,0)
!    print *, 'uno',l,omegas(:,1)
!    print *, 'dos',l,omegas(:,1)
!end do
!pause
wweno(:,0) = omegas(:,0)*w(:,0) + omegas(:,1)*w(:,1) + omegas(:,2)*w(:,2)
wweno(:,1) = omegas(:,0)*w(:,3) + omegas(:,1)*w(:,4) + omegas(:,2)*w(:,5)
wweno(:,2) = omegas(:,0)*w(:,6) + omegas(:,1)*w(:,7) + omegas(:,2)*w(:,8)
END SUBROUTINE COEFF_WENO

SUBROUTINE OSCIND(WW,oo)
USE GeneralVariables
IMPLICIT NONE
REAL    ::  A(1,3), OI1(1,1), oo,WW(3,1)
!DO l = 1:neq
    !WWWW(:,:) = WW(l,:,:)
!    print *, 'ww',ww
    A=matmul(transpose(WW),Sigma)
    OI1=matmul(A,WW)
    oo=OI1(1,1) 
!    print *, 'oo',oo

!END DO
END SUBROUTINE OSCIND

SUBROUTINE MATRICES1D
USE GeneralVariables
IMPLICIT NONE
INTEGER i,j

! Universal Oscillation matrix

	DO i = 1, 3
	   DO j = 1, 3
		  Sigma(i,j)=0.
	   END DO
	END DO
	Sigma(2,2) = 4.
	Sigma(3,3) = 156.

! Inverse of WENO matrices
!LEGENDRE
	M1(1,2) = 1.      
	M1(2,1) = -0.25
	M1(2,3) = 0.25
	M1(3,1) = 1./12.
	M1(3,2) = -1./6.
	M1(3,3) = 1./12.
!LAGRANGE
! M1(1,1) = (45.+2.*sqrt(15.))*sqrt(15.)/900.     
!  M1(1,2) = 14./15.
!      M1(1,3) = (-45.+2.*sqrt(15.))*sqrt(15.)/900.
!      M1(2,1) = -1./24.
!      M1(2,2) = 13./12.
!      M1(2,3) = -1./24.
!      M1(3,1) = (-45.+2.*sqrt(15.))*sqrt(15.)/900.
!      M1(3,2) = 14./15.
!      M1(3,3) = (45.+2.*sqrt(15.))*sqrt(15.)/900.



!LEGENDRE
	MR(1,1) = 1.     
	MR(2,1) = -0.75
	MR(2,2) = 1.
	MR(2,3) = -0.25
	MR(3,1) = 1./12.
	MR(3,2) = -1./6.
	MR(3,3) = 1./12.
!LAGRANGE
! MR(1,1) = (135.+62.*sqrt(15.))*sqrt(15.)/900.    
! MR(1,2) = -(45.+sqrt(15.))*sqrt(15.)/225.
!      MR(1,3) = (45.+2.*sqrt(15.))*sqrt(15.)/900.
!      MR(2,1) = 23./24.
!      MR(2,2) = 1./12.
!      MR(2,3) = -1./24.
!      MR(3,1) = (62.*sqrt(15.)-135.)*sqrt(15.)/900.
!      MR(3,2) = -(-45.+sqrt(15.))*sqrt(15.)/225.
!      MR(3,3) = (-45.+2.*sqrt(15.))*sqrt(15.)/900.
	
	
! LEGENDRE
	ML(1,3) = 1.      
	ML(2,1) = 0.25
	ML(2,2) = -1.
	ML(2,3) = 0.75
	ML(3,1) = 1./12.
	ML(3,2) = -1./6.
	ML(3,3) = 1./12.
! LAGRANGE
!  ML(1,1) = (-45+2*sqrt(15.))*sqrt(15.)/900.      
!  ML(1,2) = -(-45+sqrt(15.))*sqrt(15.)/225.
!      ML(1,3) = (62*sqrt(15.)-135)*sqrt(15.)/900.
!      ML(2,1) = -1./24.
!      ML(2,2) = 1./12.
!      ML(2,3) = 23./24.
!      ML(3,1) = (45+2*sqrt(15.))*sqrt(15.)/900.
!      ML(3,2) = -(45+sqrt(15.))*sqrt(15.)/225.
!      ML(3,3) = (135+62*sqrt(15.))*sqrt(15.)/900.



	
	
END SUBROUTINE MATRICES1D

SUBROUTINE LEGENDRE(icamb,xi,eta,Lx,Ly,DLx,DLy)
IMPLICIT NONE
INTEGER :: icamb
REAL    :: DLx(3),DLy(3),Lx(3),Ly(3),xi,eta,p(3)
DLx(1) = 0.
DLx(2) = 2.
DLx(3) = 12.*xi - 6.
DLy(1) = 0.
DLy(2) = 2.
DLy(3) = 12.*eta - 6.
Lx(1)  = 1.
Lx(2)  = 2.*xi - 1.
Lx(3)  = 6.*xi**2 - 6.*xi + 1.
Ly(1)  = 1.
Ly(2)  = 2.*eta - 1.
Ly(3)  = 6.*eta**2 - 6.*eta + 1.

END SUBROUTINE LEGENDRE


SUBROUTINE PolRecons(iside,icamb,xg,yg,coeff,Pol,DPolx)!(0,xg,yg,coeff,Pol,DPolx)
USE GeneralVariables
IMPLICIT NONE
INTEGER i,j,k,iside,icamb
REAL    :: DPolx(neq),DPoly(neq),Pol(neq),Lx(3),Ly(3),DLx(3),DLy(3),xg,yg,phi(9),Dphix(9),Dphiy(9),coeff(neq,9)
Lx=0.; Ly=0.; DLx=0.; DLy=0.
phi=0.; Dphix=0.
IF (iside==0) THEN ! East-West
    CALL LEGENDRE(icamb,xg,yg,Lx,Ly,DLx,DLy)
ELSE               ! North-South
    CALL LEGENDRE(icamb,xg,yg,Lx,Ly,DLx,DLy)
END IF
k = 1
DO i = 1, 3
	DO j = 1, 3
		phi(k)     =  Lx(i)*Ly(j)
        IF (iside==0) THEN
		    Dphix(k)   =  DLx(i)*Ly(j)
        ELSE
            Dphix(k)   =  DLy(j)*Lx(i)
        END IF
        k = k+1
	END DO
END DO
Pol(:)    = 0.
DPolx(:)  = 0.
DPoly(:)  = 0.
        
DO k = 1, 9
    Pol(:)     = Pol(:)   +  phi(k)    *coeff(:,k)    
    DPolx(:)   = DPolx(:) +  Dphix(k)  *coeff(:,k)
END DO
END SUBROUTINE PolRecons

! Chemotaxis routines

SUBROUTINE PREY(NP1,vprey0)
USE GeneralVariables
PARAMETER gg=3
INTEGER :: i, j, k, nng
REAL    :: xiv(gg),xpp,xp,yp,ypp,VINIT,vprey0(-2:nx+2,-2:ny+2,NP1,NP1),xmed0,ymed0,dx2,dy2,xmed,ymed
nng=3  ! Number of Gaussian points in (-1,1)
IF (nng==2) THEN
    xiv(1)=-1./sqrt(3.)
    xiv(2)=-xiv(1)
ELSE
    xiv(1)=-sqrt(3./5.)
    xiv(2)=0.
    xiv(3)=-xiv(1)
END IF
xpp =  ax-dx
xmed0 = (xpp+ax)/2.
dx2 = dx/2.
xpp = ay-dy
ymed0 = (xpp+ay)/2.
dy2 = dy/2.

DO i = 0, nx
    xmed = (x(i)+x(i+1))/2.
    DO j = 0, ny
        ymed = (y(j)+y(j+1))/2.
        DO k = 1, nng
           IF (i==0) THEN        
              xp = xmed0+dx2*xiv(k)
           ELSE
              xp = xmed+dx2*xiv(k)
           END IF
!           xp=1./2.+xiv(k)*dx/2. !Prueba      
           DO l = 1, nng
              IF (j==0) THEN        
                 yp = ymed0+dy2*xiv(l)
              ELSE
                 yp = ymed+dy/2.*xiv(l)
              END IF     
!              yp=1./2.+xiv(l)/2. !Prueba
              vprey0(i,j,k,l) = VINIT(xp,yp)
              
            !WRITE(2,*)  (x(i)+x(i+1))/2.+dx/2.*xiv(1),vprey0(i,j,1)
           END DO
        END DO
    END DO
END DO
END SUBROUTINE PREY

FUNCTION VINIT(xp,yp)
USE GeneralVariables
IMPLICIT NONE
INTEGER :: i
REAL    :: VINIT,xp,yp,xleft,xright,xg1,xg2,xii(3),ww(3),xg(3),uuu,Vman,eepsilon
xii(1) = -sqrt(3./5.)
xii(2) = 0.
xii(3) = -xii(1)
ww(1) = 5./9.
ww(2) = 8./9.
ww(3) = 5./9.
eepsilon=1.
IF (WithManufSol==0) THEN
          VINIT = (120.+eepsilon*cos(60*pi*xp))*(cos(58.*pi*yp))
ELSE
    VINIT = Vman(xp,yp,0.)
END IF
END FUNCTION VINIT


FUNCTION FV(v,u,xr,yr,tman)
USE GeneralVariables
IMPLICIT NONE
REAL    :: dvt,FU,FV,u,v,xs,tman,ForcingVmanufactured,Vman,Uman,xr,yr
pi = acos(-1.)
FU = Fm*u/(kappas+u)
FV = lambdav*v*(1.-v/Kv)-v*FU
IF (WithManufSol==1) THEN
     dvt = 10.*xr**2*(xr-1.)**2*yr**2*(yr-1.)**2*exp(-2.*xr**2-3.*yr**2)
     ForcingVmanufactured = dvt-lambdav*Vman(xr,yr,tman)*(1.-Vman(xr,yr,tman)/kv)+Vman(xr,yr,tman)*(Fm*Uman(xr,yr,tman))/(kappas+Uman(xr,yr,tman))
     FV = FV + ForcingVmanufactured
     
END IF
END FUNCTION FV

FUNCTION TestRec(x,y) 
implicit none
real    :: TestRec,x,y
TestRec=2.*x**2+3.*y**2-5.*x*y-2*x+4.*y+8.
END FUNCTION

FUNCTION DxTestRec(x,y) 
implicit none
real    :: DxTestRec,x,y
DxTestRec=4.*x-5.*y-2.
END FUNCTION DxTestRec

FUNCTION DyTestRec(x,y) 
implicit none
real    :: DyTestRec,x,y
DyTestRec=6.*y-5.*x+4.
END FUNCTION DyTestRec

SUBROUTINE TestRecinit(xL,xR,yD,yU,uu)
IMPLICIT NONE
REAL    :: xL, xR, yD, yU, uu
uu=dble((3.*xR - 3. * xL) * (yU ** 3 - yD ** 3)) / 3. + (dble(4 * xR) - dble(4 * xL) - 5./2. * dble(xR ** 2) + 5./2. * dble(xL ** 2)) * dble(yU ** 2 - yD ** 2) / 2. + 2. / 3. * dble(xR ** 3) * dble(yU - yD) - 2. / 3. * dble(xL ** 3) * dble(yU - yD) - dble(xR ** 2 * (yU - yD)) + dble(xL ** 2 * (yU - yD)) - dble(8 * xL * (yU - yD)) + dble(8 * xR * (yU - yD)) 
END SUBROUTINE TestRecinit

! OLD WENO

SUBROUTINE RecWENO1D(ua,vall,valr,dvall,dvalr,SourceGauss)
USE GeneralVariables
IMPLICIT NONE
INTEGER :: i,j
REAL    :: alphaw(0:2,2),omegaw(0:2,2),sumalph(2),betaw(0:2,2),dW(0:2),dS(0:2),eps,c1,c2,uim1(2),ui(2),uim2(2),ui1(2),ui2(2)
REAL    :: ul(2,0:2),dul(2,0:2),ur(2,0:2),dur(2,0:2),vall(neq,-1:nx+2),valr(neq,-1:nx+2),dvall(neq,-1:nx+2),dvalr(neq,-1:nx+2),suma(2),suma1(2)
REAL    :: uq1(2),uq2(2),ua(neq,-2:nx+2),uq(2,0:2),SourceGauss(neq,0:nx+1,2)


dW(0) = 1.
dW(1) = 1000.
dW(2) = 1.

eps = 1.e-20
c1 = 13./12.
c2 = 1./4.
!dW(0) = 0.3
!dW(1) = 0.6
!dW(2) = 0.1
!PRINT *, UA
!PAUSE
DO i = 0, nx

!DO j = 0, 2
!        print *, d(j)
!END DO
!PAUSE 33
   uim2(:)  = ua(:,i-2)
   uim1(:)  = ua(:,i-1)
   ui(:)    = ua(:,i)
   ui1(:)   = ua(:,i+1)
   ui2(:)   = ua(:,i+2)
!   PRINT *, I,UIM2(:)
!if (ireac==0.and.dabs(ui1)>6) then       
 !      write(*,*) ui1,ui2
!pause 77
!end if

  ul(:,2)  = ui(:)/3.+5./6.*uim1(:)-uim2(:)/6.     !!Left stencil
  ur(:,2)  = 11./6.*ui(:)-7./6.*uim1(:)+uim2(:)/3.
  dul(:,2) = (ui(:)-uim1(:))/dx
  dur(:,2) = (2.*ui(:)-3.*uim1(:)+uim2(:))/dx

  ul(:,1)  = 5./6.*ui(:)+uim1(:)/3.-ui1(:)/6.     !!Centred stencil
  ur(:,1)  = 5./6.*ui(:)-uim1(:)/6.+ui1(:)/3.
  dul(:,1) = (ui(:)-uim1(:))/dx
  dur(:,1) = (-ui(:)+ui1(:))/dx

  ul(:,0)  = 11./6.*ui(:)-7./6.*ui1(:)+ui2(:)/3.  !!Right stencil   
  ur(:,0)  = ui(:)/3.+5./6.*ui1(:)-ui2(:)/6.
  dul(:,0) = (-2.*ui(:)+3.*ui1(:)-ui2(:))/dx
  dur(:,0) = (-ui(:)+ui1(:))/dx
!
!
!		ur(2,:) = -7./6. * uim1(:) + 11. / 6. * ui(:) + uim2(:)/ 3. ! {i-2,i-1,i}
!		dur(2,:) = (-3. * uim1(:) + 2. * ui(:) + uim2(:)) / h ! {i-2,i-1,i}
!					
!		ur(1,:) = -uim1(:) / 6. + 5. / 6. * ui + ui1 /3. ! {i-1,i,i+1}
!		dur(1,:) = (ui1(:) - ui(:)) / h ! {i-1,i,i+1}
!
!		ur(0,:) = 5./6. * ui1(:) + ui(:) / 3. - ui2(:) / 6. ! {i,i+1,i+2}
!		dur(0,:) = (ui1(:) - ui(:)) / h ! {i,i+1,i+2}
!


		betaw(0,:) = c1*(ui(:)-2.*ui1(:)+ui2(:))**2+c2*(3.*ui(:)-4.*ui1(:)+ui2(:))**2
		betaw(1,:) = c1*(uim1(:)-2.*ui(:)+ui1(:))**2+c2*(uim1(:)-ui1(:))**2
		betaw(2,:) = c1*(uim2(:)-2.*uim1(:)+ui(:))**2+c2*(uim2(:)-4.*uim1(:)+3.*ui(:))**2
		sumalph(:)=0.
		DO j = 0, 2
!        PRINT *, betaw(j,:)
        !print *, d(j)!/((eps+betaw(j,:))**2)
!            PRINT *, betaw(j,:)
			alphaw(j,:) = dW(j)/((eps+betaw(j,:))**2)
			sumalph(:) = sumalph(:) + alphaw(j,:)
		END DO
		DO j = 0, 2
			omegaw(j,:) = alphaw(j,:)/sumalph(:)
		END DO
!        omegaw(0,:)=0.
!        omegaw(1,:)=1.
!        omegaw(2,:)=0.
		suma(:) = 0.
		suma1(:) = 0.
		DO  j = 0, 2
			suma(:)  = suma(:) + omegaw(j,:)*ur(:,j)
			suma1(:) = suma1(:) + omegaw(j,:)*dur(:,j) 
		END DO
		valR(:,i)  = suma(:)
		dvalR(:,i) = suma1(:)  

		suma(:) = 0.
		suma1(:) = 0.
		DO  j = 0, 2
			suma(:)  = suma(:) + omegaw(j,:)*ul(:,j)
			suma1(:) = suma1(:) + omegaw(j,:)*dul(:,j) 
		END DO
		valL(:,i)  = suma(:)
		dvalL(:,i) = suma1(:)  


!if (ireac==0.and.dabs(dvalr)>6.) then       
!write(*,*) valR,dvalR
!pause 77
!end if

		! DIFFUSIVE TERMS

	!	dvalr=dble((25. * uim1 - 245. * ui + 245. * ui1 - 25. * ui2 - 2. * uim2 + 2. * ui3) / h) / 0.180D3

!dvalr=dur(1)
		! SOURCE
		




!dvalr=dur(0)
		!
		! WENO FOR REACTIVE TERMS 
		! ---------------------------------------------------------
		!
		! ********************
		! First Gaussian point
		! ********************
		!


		dS(0)=(210.-sqrt(3.))/(1080.)
		dS(1)=11./18.
		dS(2)=(210.+sqrt(3.))/(1080.)
		sumalph(:)=0.
		DO j = 0, 2
			alphaw(j,:)=dS(j)/((eps+betaw(j,:))**2)
			sumalph(:) = sumalph(:) + alphaw(j,:)
		END DO
		DO j = 0, 2
			omegaw(j,:) = alphaw(j,:)/sumalph(:)
		END DO
		uq(:,2)= ui(:) - (-4. * uim1(:) + 3. * ui(:) + uim2(:)) * sqrt(3.) / 12. ! {i-2,i-1,i}		
		uq(:,1) = ui(:) - (-uim1 + ui1) * sqrt(3.) / 12. ! {i-1,i,i+1}
		uq(:,0) = ui(:) + (ui2(:) + 3. * ui(:) - 4. * ui1(:)) * sqrt(3.) / 12. ! {i,i+1,i+2}
		suma(:) = 0.
		DO  j = 0, 2
			suma(:) = suma(:) + omegaw(j,:)*uq(:,j)
		END DO
		SourceGauss(:,i,1) = suma(:)
		!
		! ********************
		! Second Gaussian point
		! ********************
		!
		dS(0)=(210.+sqrt(3.))/(1080.)
		dS(1)=11./18.
		dS(2)=(210.-sqrt(3.))/(1080.)
		sumalph(:)=0.
		DO j = 0, 2
			alphaw(j,:)=dS(j)/((eps+betaw(j,:))**2)
			sumalph(:) = sumalph(:) + alphaw(j,:)
		END DO
		DO j = 0, 2
			omegaw(j,:) = alphaw(j,:)/sumalph
		END DO

		uq(:,2)= ui(:) + (-4. * uim1(:) + 3. * ui(:) + uim2(:)) * sqrt(3.) / 12. ! {i-2,i-1,i}		
		uq(:,1) = ui(:) + (-uim1(:) + ui1(:)) * sqrt(3.) / 12. ! {i-1,i,i+1}
		uq(:,0) = ui(:) - (ui2(:) + 3. * ui(:) - 4. * ui1(:)) * sqrt(3.) / 12.
		 ! {i,i+1,i+2}

		suma(:) = 0.
		DO  j = 0, 2
			suma(:) = suma(:) + omegaw(j,:)*uq(:,j)
		END DO
		SourceGauss(:,i,2) = suma(:)
    end do
END SUBROUTINE RecWENO1D

