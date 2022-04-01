! 11:38
program ch2_main

    use ch2_general_vars
    use ch2_subroutines

    use, intrinsic:: iso_fortran_env, only: stdin=>input_unit

    character*80 text

    integer             :: nx1, i, j, k, l, NP1, kdib, numfile, imetODE, i_type
    real    			:: tolE, difE(3)
    real    			:: dt1, dt2, FV, k1, k2, k3, k4, taux, u00, uold, vold, ux, uy
    real				:: x1, x2, Vx, Vy, y1,y2, vvel, Difusx, Difusy, vvelx, vvely, kappa, p, xr, yr, QV, wg1, wg2
    real                :: XG, YG, TestRec, xtest, ytest, DxTestRec, DyTestRec, xvg(2,2), normM, normA, tol, NormL2, d1
    real                :: auxnorm(3), media
    real, allocatable   :: u(:,:,:), u0(:,:,:), uRK(:,:,:), coeff(:,:), Pol(:), DPolx(:), DPoly(:), uu(:), f(:), g(:)
    real, allocatable   :: df(:), dg(:), Q(:), sou(:)
    real, allocatable   :: vprey(:,:,:,:), vprey0(:,:,:,:), VAver(:,:)
    !REAL,ALLOCATABLE   ::  vprey(:,:),vprey0(:,:),u(:,:),unorm(:),u0(:,:),uRK(:,:),coeff(:,:),Pol(:),DPolx(:),uu(:),f(:),df(:),Q(:),sou(:),ddd(:),SourceGauss(:,:,:),dvall(:,:),vall(:,:),dvalr(:,:),valr(:,:)

    character(80)       :: userinp

    !!! ------------------------------------------------------------------------------------------------------------ !!!
    !!!                                                 DATA INPUT                                                   !!!
    !!! ------------------------------------------------------------------------------------------------------------ !!!

    !   -------------- Data for upper boundary ---------------------------------------------------------------------   !


    print *, 'Input variables:'
    open(unit=9, file='Data.dat', status='old')
    read(9,*, IOSTAT=iostat) text   ! Coordinates of domain                    ! 1
    read(9,*, IOSTAT=iostat) ax                                                ! 2
    print *, 'ax=', ax
    read(9,*, IOSTAT=iostat) bx                                                ! 3
    print *, 'bx=', bx
    read(9,*, IOSTAT=iostat) ay                                                ! 4
    print *, 'ay=', ay
    read(9,*, IOSTAT=iostat) by                                                ! 5
    print *, 'by', by
    read(9,*, IOSTAT=iostat) text   ! Number of cells                          ! 6
    read(9,*, IOSTAT=iostat) nx                                                ! 7
    print *, 'nx=', nx
    read(9,*, IOSTAT=iostat) ny                                                ! 8
    print *, 'ny=', ny
    read(9,*, IOSTAT=iostat) text   ! Number of Gaussian quadrature points     ! 9
    read(9,*, IOSTAT=iostat) NP1                                               ! 10
    print *, 'NP1=', NP1
    read(9,*, IOSTAT=iostat) text   ! Output time                              ! 11
    read(9,*, IOSTAT=iostat) tmax                                              ! 12
    print *, 'tmax=', tmax
    read(9,*, IOSTAT=iostat) text                                              ! 13
    read(9,*, IOSTAT=iostat) cfl_difx                                          ! 14
    print *, 'cfl_difx=', clf_difx
    read(9,*, IOSTAT=iostat) text                                              ! 15
    read(9,*, IOSTAT=iostat) cfl_dify                                          ! 16
    print *, 'cfl_dify=', cfl_dify
    read(9,*, IOSTAT=iostat) text                                              ! 17
    read(9,*, IOSTAT=iostat) cfl_advx                                          ! 18
    print *, 'cfl_advx=', cfl_advx
    read(9,*, IOSTAT=iostat) text                                              ! 19
    read(9,*, IOSTAT=iostat) cfl_advy                                          ! 20
    print *, 'cfl_advy=', cfl_advy

    close(9)

    ngauss=3 !Number of Gaussian points
    !Type of problem: icase=1 (Advection-diffusion), icase=2 (Euler)

    WithManufSol=1 !0: NO, 1: YES
    !number of variables (equations)
    print *, 'CHEMOTAXIS PROBLEM'
    !PAUSE

    ! allocating of variables
    allocate(LO(neq,0:nx,0:ny), x(-2:nx+3), y(-2:ny+3))
    allocate(u(neq,-2:nx+2,-2:ny+2), u0(neq,-2:nx+2,-2:ny+2), uRK(neq,-2:nx+2,-2:ny+2))
    allocate(coeff(neq,ngauss*ngauss), xig(NP1), wg(NP1))
    allocate(OI(neq,0:2), omegah(neq,0:2), omegas(neq,0:2), wweno(neq,0:2))
    allocate(coef(neq, ngauss*ngauss, -2:nx+2, -2:ny+2))
    allocate(WWX0(neq,-2:nx+2,-2:ny+2), WWX1(neq,-2:nx+2,-2:ny+2), WWX2(neq,-2:nx+2,-2:ny+2),w(neq,0:8))
    allocate(FFO(neq,0:nx+1,0:ny+1,NP1), FDifus(neq,0:nx+2,0:ny+2,NP1))
    allocate(fR(neq,0:nx+1,0:ny+1,NP1), fL(neq,0:nx+1,0:ny+1,ngauss))
    allocate(FFLUX(neq,0:nx+1,0:ny+1), FDFLUX(neq,0:nx+1,0:ny+1))
    !ALLOCATE(FDFLUX1(neq,0:nx+2,0:ny+2),GDFLUX1(neq,0:nx+1,0:ny+1))
    allocate(GFO(neq,0:nx+1,0:ny+1,NP1), GDifus(neq,0:nx+2,0:ny+2,NP1))
    allocate(GFLUX(neq,0:nx+1,0:ny+1), GDFLUX(neq,0:nx+1,0:ny+1), gD(neq,0:nx+1,0:ny+1,NP1), gU(neq,0:nx+1,0:ny+1,NP1))
    allocate(fDR(neq,0:nx+2,0:ny+2,NP1), fDL(neq,0:nx+2,0:ny+2,NP1))
    allocate(gDD(neq,0:nx+2,0:ny+2,NP1), gDU(neq,0:nx+2,0:ny+2,NP1))
    allocate(Pol(neq), DPolx(neq), DPoly(neq), Q(neq), f(neq), g(neq), df(neq), dg(neq))
    allocate(source(neq,0:nx+1,1:ny+1,NP1,NP1), sou(neq), SSource(neq,0:nx+1,0:ny+1))
    allocate(PolL(neq,0:nx+1,0:ny+1,NP1), PolR(neq,0:nx+1,0:ny+1,NP1))
    allocate(PolD(neq,0:nx+1,0:ny+1,NP1), PolU(neq,0:nx+1,0:ny+1,NP1))
    allocate(DPolR_x(neq,0:nx+1,0:ny+1,NP1), DPolL_x(neq,0:nx+1,0:ny+1,NP1))
    allocate(DPolD_y(neq,0:nx+1,0:ny+1,NP1), DPolU_y(neq,0:nx+1,0:ny+1,NP1))

    !!! ------------------------------------------------------------------------------------------------------------ !!!
    !!!                                             MESH GENERATION                                                  !!!
    !!! ------------------------------------------------------------------------------------------------------------ !!!

    x=0.
    y=0.

    call build_mesh ! aka MALLADO

    !!! ------------------------------------------------------------------------------------------------------------ !!!
    !!!                                BUILD MATRICES FOR WENO RECONSTURCTION                                        !!!
    !!! ------------------------------------------------------------------------------------------------------------ !!!

    i_type = 2  ! 1=legendre, 2=lagrange
    call build_weno_matrices_1d(i_type) ! aka MATIRCES1D

    epsilon=1.d-15
    r=12.
    DO i = 0, 2
        IF (i==0) THEN
            lambda(i) = 1000.
        ELSE
            lambda(i) = 1.
        END IF
    END DO

    !!! ------------------------------------------------------------------------------------------------------------ !!!
    !!!                                               GAUSSIAN PARAMETERS                                            !!!
    !!! ------------------------------------------------------------------------------------------------------------ !!!

    print *, 'y='
    print *, y

    call gauss_quadrature(NP1) ! xig, wg: ngauss components.


    !!! ------------------------------------------------------------------------------------------------------------ !!!
    !!!                                               INITIAL CONDITIONS                                             !!!
    !!! ------------------------------------------------------------------------------------------------------------ !!!

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

    open(10,file='initial_u.dat',status='unknown')
    open(11,file='initial_w.dat',status='unknown')
    open(12,file='initial_v.dat',status='unknown')

    DO i = 1, nx-1
        x1=x(i)
        x2=x(i+1)
        DO j = 1, ny-1
            y1=y(j)
            y2=y(j+1)
            if (WithManufSol==0) then
                call sol_inic(NP1,x1,x2,y1,y2,Q) ! aka SOLINIC
            else
                call sol_manuf(0.,x1,x2,y1,y2,Q) ! aka SolManuf
            end if
            u0(:,i,j) = Q(:)
            write(10,*) (x1+x2)/2.,(y1+y2)/2.,u0(1,i,j)
            write(11,*) (x1+x2)/2.,(y1+y2)/2.,u0(2,i,j)
            write(12,*) (x1+x2)/2.,(y1+y2)/2.,u0(3,i,j)
            k=k+1
        END DO
    END DO
    DO i = 0, neq-1
        CLOSE(i)
    END DO

    t = 0.
    uRK = 0.
    kdib = 1.
    numfile = 44

    call gauss_quadrature (Np1) ! Gaussian quadrature parameters in (0, 1)

    tol = 1.e-8
    tolE = 1.e-8
    normM = 2.*tol; normA=2.*tol
    difE(:)=2.*tolE
    !DO WHILE (difE(1)>tolE.or.difE(2)>tolE)
    !DO WHILE (difE(1)>tolE)

    !print *, 'Waiting for Enter.'
    !read(stdin,*)

    do while (t<tmax)
        call calctime(u0) ! aka calctime

        t = t+dt
        if (WithManufSol.eq.1) then
            if (t>tmax) then
                taux = t-tmax
                t = t - dt
                dt = dt-taux
                t = t + dt
                write(*,*) 'time=',t, 'step size=',dt
            end if
            write(*,*) 'time=',t, 'step size=',dt
        end if
        nx1=nx-1
        !!! -------------------------------------------------------------------------------------------------------- !!!
        !!!                                OPERATOR L RECONSTRUCTION FROM CELL AVERAGES                              !!!
        !!! -------------------------------------------------------------------------------------------------------- !!!

        call boundary_conditions(u0) ! aka BoundaryConditions

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

        call boundary_conditions(uRK)
        CALL OpL(coeff,NP1,uRK,vprey,Pol,DPolx,DPoly)

        DO i = 1, nx-1
            DO j = 1, ny-1
                uRK(:,i,j) = (3.*u0(:,i,j) + (uRK(:,i,j) + dt*LO(:,i,j)))/4.
            END DO
        END DO


        call boundary_conditions(uRK)

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

end program ch2_main