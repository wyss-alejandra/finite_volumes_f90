module ch2_subroutines

    use ch2_general_vars

    implicit none

contains

    subroutine build_mesh  ! aka MALLADO
        implicit none
        integer :: i

        Ncells = (nx-1)*(ny-1)
        dx = (bx-ax)/(nx-1)
        dy = (by-ay)/(ny-1)
        do i = 1, nx
            x(i)=ax+(i-1)*dx
        end do

        do i = 1, ny
            y(i)=ay+(i-1)*dy
        end do
    end subroutine build_mesh


    subroutine build_weno_matrices_1d(i_type) ! aka MATRICES1D
        ! i_type: 1=legendre, 2=lagrange
        implicit none
        integer, intent(in) :: i_type

        integer :: i,j

        ! Universal Oscillation matrix

        DO i = 1, 3
            DO j = 1, 3
                Sigma(i,j)=0.
            END DO
        END DO
        Sigma(2,2) = 4.
        Sigma(3,3) = 156.

        ! Inverse of WENO matrices

        if (i_type == 1) then
            !LEGENDRE
            M1(1,2) = 1.
            M1(2,1) = -0.25
            M1(2,3) = 0.25
            M1(3,1) = 1./12.
            M1(3,2) = -1./6.
            M1(3,3) = 1./12.

            MR(1,1) = 1.
            MR(2,1) = -0.75
            MR(2,2) = 1.
            MR(2,3) = -0.25
            MR(3,1) = 1./12.
            MR(3,2) = -1./6.
            MR(3,3) = 1./12.

            ML(1,3) = 1.
            ML(2,1) = 0.25
            ML(2,2) = -1.
            ML(2,3) = 0.75
            ML(3,1) = 1./12.
            ML(3,2) = -1./6.
            ML(3,3) = 1./12.
        else if (i_type == 2) then
            !LAGRANGE
            M1(1,1) = (45.+2.*sqrt(15.))*sqrt(15.)/900.
            M1(1,2) = 14./15.
            M1(1,3) = (-45.+2.*sqrt(15.))*sqrt(15.)/900.
            M1(2,1) = -1./24.
            M1(2,2) = 13./12.
            M1(2,3) = -1./24.
            M1(3,1) = (-45.+2.*sqrt(15.))*sqrt(15.)/900.
            M1(3,2) = 14./15.
            M1(3,3) = (45.+2.*sqrt(15.))*sqrt(15.)/900.

            MR(1,1) = (135.+62.*sqrt(15.))*sqrt(15.)/900.
            MR(1,2) = -(45.+sqrt(15.))*sqrt(15.)/225.
            MR(1,3) = (45.+2.*sqrt(15.))*sqrt(15.)/900.
            MR(2,1) = 23./24.
            MR(2,2) = 1./12.
            MR(2,3) = -1./24.
            MR(3,1) = (62.*sqrt(15.)-135.)*sqrt(15.)/900.
            MR(3,2) = -(-45.+sqrt(15.))*sqrt(15.)/225.
            MR(3,3) = (-45.+2.*sqrt(15.))*sqrt(15.)/900.

            ML(1,1) = (-45+2*sqrt(15.))*sqrt(15.)/900.
            ML(1,2) = -(-45+sqrt(15.))*sqrt(15.)/225.
            ML(1,3) = (62*sqrt(15.)-135)*sqrt(15.)/900.
            ML(2,1) = -1./24.
            ML(2,2) = 1./12.
            ML(2,3) = 23./24.
            ML(3,1) = (45+2*sqrt(15.))*sqrt(15.)/900.
            ML(3,2) = -(45+sqrt(15.))*sqrt(15.)/225.
            ML(3,3) = (135+62*sqrt(15.))*sqrt(15.)/900.
        end if
    end subroutine build_weno_matrices_1d

    subroutine gauss_quadrature(NP1)
        implicit none
        integer, parameter  :: m=9
        integer             :: k,NP1
        real                :: cg(m),gi(m),gama(m)

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
    end subroutine gauss_quadrature

    subroutine sol_inic(NP1,xL,xR,yD,yU,Q) ! aka SOLINIC
        implicit none
        integer     :: i, NP1, j, k
        real        :: derf, xg, yg, xL, xR, yD, yU, Q(neq), wwg, XXX, YYY, H, eepsilon
        real        :: kappa, p, rho, ux, uy, xreal, yreal, wg1, wg2
        real        :: sum1, sum2, sum3, gama(9), gi(9), sss, aa, bb, Rad, Uman

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
    end subroutine sol_inic

    function U_man(xr,yr,tman)
        implicit none
        integer         :: i,j
        real            :: U_man,expo,xr,yr,tman,u1

        U_man=100.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-by)**2*exp(-xr**2-yr**2)
        !Uman=80.*(1.+tman)*((xr-ax)**2+(yr-ay)**2)
    end function U_man

    function W_man(xr,yr,tman)
        implicit none
        real            ::  W_man,expo,xr,yr,tman
        W_man = 80.*tman*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-by)**2*exp(-3.*xr**2-3.*yr**2)
        !Wman=100.*(1.+tman)*(xr-ax)**2*(xr-bx)**2*(yr-ay)**2*(yr-by)**2*exp(-xr**2-yr**2)
    end function W_man

    function V_man(xr,yr,tman)
        implicit none
        real            ::  V_man,expo,xr,yr,tman
        V_man = 10.*tman*(xr-ax)**2*(xr-bx)**2.*(yr-ay)**2*(yr-by)**2*exp(-2.*xr**2-3.*yr**2)
    end function V_man

    subroutine sol_manuf(tman,xL,xR,yD,yU,Q) ! aka SolManuf
        implicit none
        integer     :: i, NP1, j, k
        real        :: derf, xg, yg, xL, xR, yD, yU, Q(neq), wwg, XXX, YYY, tman
        real        :: kappa, p, rho, ux, uy, xreal, yreal, wg1, wg2, sum1, sum2, sum3
        real        :: gama(9), gi(9), sum_U, sum_V, sum_W

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
        sum_U=0.
        sum_W=0.
        sum_V=0.
        DO i = 1, 9
            xreal   = (xL+xR)/2.+dx*gi(i)/2.
            wg1  = dx*gama(i)/2.
            DO j = 1, 9
                wg2 = dy*gama(j)/2.
                yreal  = (yD+yU)/2.+dy*gi(j)/2.
                sum_U = sum_U + wg1*wg2*U_man(xreal,yreal,tman)
                sum_W = sum_W + wg1*wg2*W_man(xreal,yreal,tman)
                sum_V = sum_V + wg1*wg2*V_man(xreal,yreal,tman)
            END DO
        END DO
        Q(1) = sum_U / (dx*dy)
        Q(2) = sum_W / (dx*dy)
        Q(3) = sum_V / (dx*dy)
    end subroutine sol_manuf

    subroutine calc_time(ua) ! aka calctime
        implicit none
        integer     :: i, j, k
        real        :: velx, vely, p, aa, dta, dt1, dt2, kappa, lambdaa, lambda1(6), ua(neq,-2:nx+2,-2:ny+2)
        real        :: vvel, vvelx, vvely, Vx, Vy, Difusx, Difusy, dtt(4), dtaux

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
                if (dtaux<dt) then
                    dt=dtaux
                end if
            END DO
        END DO
    end subroutine calc_time

    subroutine boundary_conditions(ua) ! aka BoundaryConditions
        implicit none
        integer     :: iBC(4),i,j,k,kk
        real        :: ua(neq,-2:nx+2,-2:ny+2), cdirE(4), cdirW(4), cdirN(4), cdirS(4), eps, gA, time
        real        :: u1, u2, um, um1, um2, um3, um4
        ! 1: South; 2: East; 3: North; 4: West

        DO i = 1, 4
            iBC(i)=0 !No flow
            !    iBC(i)=-1 !Periodic
        END DO

        if (iBC(2) == -1) then ! Periodic in East side of the domain
            DO j = -2, ny+2
                ua(:,nx,j)     = ua(:,1,j)
                ua(:,nx+1,j)   = ua(:,2,j)
                ua(:,nx+2,j)   = ua(:,3,j)
            END DO
        else if (iBC(2) == 0) then ! No flow in East side of the domain
            DO j = -2, ny+2
                ua(:,nx,j)     = ua(:,nx-1,j)
                ua(:,nx+1,j)   = ua(:,nx-2,j)
                ua(:,nx+2,j)   = ua(:,nx-3,j)
            END DO
        else if (iBC(2) == 1) then  !Dirichlet
            DO j = -2, ny+2
                ua(:,nx,j)     = 2.*cdirE(:)-ua(:,nx-1,j)
                ua(:,nx+1,j)   = 2.*cdirE(:)-ua(:,nx-2,j)
                ua(:,nx+2,j)   = 2.*cdirE(:)-ua(:,nx-3,j)
            END DO
        else if (iBC(2) == 3) then  !Solid walls
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
        else if (iBC(2) == 4) then ! Absorbing
            DO j = -2, ny+2
                ua(:,nx,j)     = ua(:,nx-1,j)
                ua(:,nx+1,j)   = ua(:,nx-1,j)
                ua(:,nx+2,j)   = ua(:,nx-1,j)
            END DO

        end if
        if (iBC(4) == -1) then ! Periodic in West side of the domain
            DO j = -2, ny+2
                ua(:,0,j)     = ua(:,nx-1,j)
                ua(:,-1,j)    = ua(:,nx-2,j)
                ua(:,-2,j)    = ua(:,nx-3,j)
            END DO
        else if (iBC(4) == 0) then ! No flow in West side of the domain
            DO j = -2, ny+2
                ua(:,0,j)     = ua(:,1,j)
                ua(:,-1,j)    = ua(:,2,j)
                ua(:,-2,j)    = ua(:,3,j)
            END DO
        else if (iBC(4) == 1) then ! Dirichlet
            !print *, ua(:,1,ny+2)
            DO j = -2, ny+2
                ua(:,0,j)     = 2.*cdirW(:)-ua(:,1,j)
                ua(:,-1,j)    = 2.*cdirW(:)-ua(:,2,j)
                ua(:,-2,j)    = 2.*cdirW(:)-ua(:,3,j)
                !          print *, i,j,ua(1,-1,j)
            END DO
        else if (iBC(4) == 3) then ! Reflecting
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
        else if (iBC(4) == 4) then ! Absorbing
            DO j = -2, ny+2
                ua(:,0,j)     = ua(:,1,j)
                ua(:,-1,j)    = ua(:,1,j)
                ua(:,-2,j)    = ua(:,1,j)
            END DO

        end if
        if (iBC(1) == -1) then ! Periodic in South side of the domain
            DO i = -2, nx+2
                ua(:,i,0)     = ua(:,i,ny-1)
                ua(:,i,-1)    = ua(:,i,ny-2)
                ua(:,i,-2)    = ua(:,i,ny-3)
            END DO
        else if (iBC(1) == 0) then ! No flow in South side of the domain
            DO i = -2, nx+2
                ua(:,i,0)     = ua(:,i,1)
                ua(:,i,-1)    = ua(:,i,2)
                ua(:,i,-2)    = ua(:,i,3)
            END DO
        else if (iBC(1) == 1) then  !Dirichlet
            DO i = -2, nx+2
                ua(:,i,0)     = 2.*cdirS(:)-ua(:,i,1)
                ua(:,i,-1)    = 2.*cdirS(:)-ua(:,i,2)
                ua(:,i,-2)    = 2.*cdirS(:)-ua(:,i,3)
            END DO
        else if (iBC(1) == 2) then ! Periodic
            DO i = -2, nx+2
                ua(:,i,0)     = ua(:,i,nx-1)
                ua(:,i,-1)    = ua(:,i,nx-2)
                ua(:,i,-2)    = ua(:,i,nx-3)
            END DO

        end if

        if (iBC(3) == -1) then ! Periodic in North side of the domain
            DO i = -2, nx+2
                ua(:,i,ny)    = ua(:,i,1)
                ua(:,i,ny+1)  = ua(:,i,2)
                ua(:,i,ny+2)  = ua(:,i,3)
            END DO
        else if (iBC(3) == 0) then ! No flow in North side of the domain
            DO i = -2, nx+2
                ua(:,i,ny)    = ua(:,i,ny-1)
                ua(:,i,ny+1)  = ua(:,i,ny-2)
                ua(:,i,ny+2)  = ua(:,i,ny-3)
            END DO
        else if (iBC(3) == 1) then  !Dirichlet
            DO i = -2, nx+2
                ua(:,i,ny)     = 2.*cdirN(:)-ua(:,i,ny-1)
                ua(:,i,ny+1)   = 2.*cdirN(:)-ua(:,i,ny-2)
                ua(:,i,ny+2)   = 2.*cdirN(:)-ua(:,i,ny-3)
            END DO
        else if (iBC(3) == 2) then ! Periodic
            DO i = -2, nx+2
                ua(:,i,ny)    = ua(:,i,1)
                ua(:,i,ny+1)  = ua(:,i,2)
                ua(:,i,ny+2)  = ua(:,i,3)
            END DO
        end if
    end subroutine boundary_conditions

    subroutine flux_hyperb(xg,yg,Deriv,Q,f)
        implicit none
        integer     :: i
        real        :: f(neq), g(neq), grav, H, u, Vx, Vy, xg, yg, p, kappa, xxx, yyy
        real        :: xL, xR, yD, yU, Q(neq), ux, uy, Deriv, DPoly(neq)
        f(1) = Xi*Deriv*Q(1)/dx
        !f(1) = Xi*Q(1)!/dx
        !print *, 'der',f(1),Deriv
        f(2) = 0.
        f(3) = 0.
    end subroutine flux_hyperb

    subroutine recons_weno_dim_dim(ua) ! aka ReconsWenoDimDim
        implicit none
        integer             :: i, idir, ileg, j, jleg, kleg, l
        real                :: Pol, phi(9), ua(neq,-2:nx+2,-2:ny+2), oo
        real, allocatable   :: uu(:,:)

        allocate(uu(neq,5))
        ! idir=1 Reconstruction in x-direction, idir=2 Reconstruction in y-direction
        !coef=0.
        !print *, ua
        !pause 44
        !print *, ua
        WWX0=0.
        WWX1=0.
        WWX2=0.
        DO idir = 1, 2
            if (idir==1) then
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
                        call WENO1D(uu)
                        !print *,wweno
                        WWX0(:,i, j)            = wweno(:,0)
                        WWX1(:,i, j)            = wweno(:,1)
                        WWX2(:,i, j)            = wweno(:,2)
                        !    print *, i,j,WWX1(:,i, j)
                    END DO
                END DO
                !        print *, WWX1  
                !END DO
            else
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

                        call WENO1D(uu)

                        coef(:,1,i,j)        = wweno(:,0)
                        coef(:,2,i,j)        = wweno(:,1)
                        coef(:,3,i,j)        = wweno(:,2)

                        uu(:,1)              = WWX1(:, i, j-2)
                        uu(:,2)              = WWX1(:, i, j-1)
                        uu(:,3)              = WWX1(:, i, j)
                        uu(:,4)              = WWX1(:, i, j+1)
                        uu(:,5)              = WWX1(:, i, j+2)

                        !print *, 'uu2',uu

                        call WENO1D(uu)
                        coef(:,4,i,j)        = wweno(:,0)
                        coef(:,5,i,j)        = wweno(:,1)
                        coef(:,6,i,j)        = wweno(:,2)
                        uu(:,1)              = WWX2(:, i, j-2)
                        uu(:,2)              = WWX2(:, i, j-1)
                        uu(:,3)              = WWX2(:, i, j)
                        uu(:,4)              = WWX2(:, i, j+1)
                        uu(:,5)              = WWX2(:, i, j+2)
                        !print *, 'uu3',uu

                        call WENO1D(uu)
                        coef(:, 7,i,j)       = wweno(:,0)
                        coef(:, 8,i,j)       = wweno(:,1)
                        coef(:, 9,i,j)       = wweno(:,2)

                    END DO
                END DO
                !END DO
            end if
        END DO
        !DEALLOCATE(uu)
    end subroutine recons_weno_dim_dim

    SUBROUTINE pol_recons(iside,icamb,xg,yg,coeff,Pol,DPolx)!(0,xg,yg,coeff,Pol,DPolx) ! aka PolRecons
        implicit none
        integer     :: i, j, k, iside, icamb
        real        :: DPolx(neq), DPoly(neq), Pol(neq), Lx(3), Ly(3), DLx(3), DLy(3), xg, yg
        real        :: phi(9), Dphix(9), Dphiy(9), coeff(neq,9)

        Lx=0.; Ly=0.; DLx=0.; DLy=0.
        phi=0.; Dphix=0.

        if (iside==0) then ! East-West
            call LEGENDRE(icamb,xg,yg,Lx,Ly,DLx,DLy)
        else               ! North-South
            call LEGENDRE(icamb,xg,yg,Lx,Ly,DLx,DLy)
        end if
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

    end subroutine pol_recons

    subroutine legendre(icamb,xi,eta,Lx,Ly,DLx,DLy)
        implicit none
        integer     :: icamb
        real        :: DLx(3), DLy(3), Lx(3), Ly(3), xi, eta, p(3)

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

    end subroutine legendre

    subroutine flux_diffusive(idifu,xg,yg,Pol,DPolx,df)
        implicit none
        integer     :: i, n, idifu
        real        :: Pol(neq), DPolx(neq), u, Difusx, Difusy, xg, yg, df(neq)

        df(1) = d_u*DPolx(1) ! - Pol(1)*Xi *DPolx(2)
        df(2) = d_w*DPolx(2)
        df(3) = 0.
    end subroutine flux_diffusive

    subroutine recons_flux_difus(NP1,u0) ! aka ReconsFluxDifus
        implicit none
        integer     :: i, j, k, l, NP1
        real        :: u0(neq,-2:nx+2,-2:ny+2), sumfk, sumgk,val, val1, theta, dd(2)
        real        :: derf, x1, x2, xr, y1, y2, yr
        real        :: tman, t4, t5, t7, t8, t9, t12, t13, t14, t16, t17, t18, t19, t20, t21, t22, t23, t25, t26
        real        :: t30, t32, t34, t35, t36, t37, t38
        real        :: t41, t43, t54, t55, t56, t62, t65, t79, t82, t85, t98, t105, valint, DUmanX, DUmanY, dif1,dif2

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
                t56 = -2000000000.D0*t5*y1 - 800000000.D0*t5*t9 - 3200000000.D0*t12 + 1600000000.D0*t5*t8 &
                        + 1600000000.D0*t18 + 3200000000.D0*t21 + 800000000.D0*t20*t23 + 2000000000.D0*t20*y2 &
                        - 1600000000.D0*t20*t4 - 1600000000.D0*t30 - 1600000000.D0*t32 + 1600000000.D0*t34 &
                        + 1772453851.D0*t37*t38 - 1772453851.D0*t37*t41 + 1600000000.D0*t12*t9 &
                        + 800000000.D0*t32*t9 + 2000000000.D0*t32*y1 - 1600000000.D0*t32*t8 + 2000000000.D0*t30*y2 &
                        + 4000000000.D0*t12*y1
                t79 = t37*xr
                t82 = t17*t37
                t85 = t16*t37
                t98 = -3200000000.D0*t12*t8 - 800000000.D0*t18*t9 + 1600000000.D0*t18*t8 - 1600000000.D0*t21*t23 &
                        - 4000000000.D0*t21*y2 + 3200000000.D0*t21*t4 + 800000000.D0*t30*t23 - 1600000000.D0*t30*t4 &
                        - 800000000.D0*t34*t23 - 2000000000.D0*t34*y2 + 1600000000.D0*t34*t4 - 3544907702.D0*t79*t38 &
                        + 1772453851.D0*t82*t38 - 1772453851.D0*t85*t38 + 3544907702.D0*t79*t41 &
                        - 1772453851.D0*t82*t41 + 1772453851.D0*t85*t41 - 2000000000.D0*t18*y1 + 1600000000.D0*t5 &
                        - 1600000000.D0*t20
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
                t65 = 0.1772453851D2*t8*t9 - 0.8862269255D1*t13*t14 + 8.D0*t18*t19 + 10.D0*t22*x1 &
                        - 0.8862269255D1*t26*t9 + 0.8862269255D1*t13*t9 - 0.1772453851D2*t8*t14 &
                        + 0.8862269255D1*t26*t14 - 10.D0*t36*x2 + 20.D0*t18*x1 - 16.D0*t18*t5 - 4.D0*t43*t19 &
                        - 10.D0*t43*x1 + 8.D0*t43*t5 + 4.D0*t22*t19 - 8.D0*t22*t5-8.D0*t54*t55 - 20.D0*t54*x2 &
                        + 16.D0*t54*t4 + 10.D0*t62*x2
                t98 = -8.D0*t62*t4 - 4.D0*t36*t55 + 8.D0*t36*t4 + 4.D0*t62*t55 + 0.8862269255D1*t7*t14 &
                        - 0.8862269255D1*t7*t9 - 16.D0*t18 + 8.D0*t36-8.D0*t22 - 4.D0*t17*t19 - 10.D0*t17*x1 &
                        + 8.D0*t17*t5 + 8.D0*t43+16.D0*t54 + 4.D0*t35*t55 + 10.D0*t35*x2 - 8.D0*t35*t4 - 8.D0*t62 &
                        + 8.D0*t17 - 8.D0*t35
                t105 = exp(-1.D0*t12-1.D0*t5-1.D0*t4)
                valint = 25.D0*(1.D0+tman)*yr*(yr-1.D0)*(t65+t98)*t105



                !      print *, i,j,GDFLUX(1,i,j),d_u*valint/dx
                !      pause
                FDFLUX(3,i,j)=0.
                GDFLUX(3,i,j)=0.
            END DO
        END DO
        !print *, FDFLUX
    end subroutine recons_flux_difus

    subroutine OpL(coeff,NP1,ua,vprey,Pol,DPolx,DPoly)
        implicit none
        integer             :: NP1, i, j, k, l, m, kappa, iesq, kk
        real                :: alpha, xL, xR, yD, yU, ua(neq,-2:nx+2,-2:ny+2), coeff(neq,9), Vx, Vy
        real                :: Difusx, Difusy, slope1(2), slope2(2), slopex, slopey
        real                :: sumfk, sumgk, sumsk, xg, yg, wg1, wg2, Pol(neq), DPolx(neq), DPoly(neq)
        real                :: vprey(-2:nx+2,-2:ny+2,NP1,NP1), preyy
        real                :: ForcingTerm(neq)
        real                :: FLF(neq), GLF(neq), ULW(neq), FLW(neq), GLW(neq), FLU(neq), u00(neq)
        real                :: DUmanX, DUmanY, Uman, Wman
        real                :: vall(neq,-1:nx+2), valr(neq,-1:nx+2), dvall(neq,-1:nx+2), dvalr(neq,-1:nx+2) 
        real                :: SourceGauss(neq,0:nx+1,2), dddf, dddg, derw
        real                :: um1(neq), um2(neq), ui(neq), u1(neq), u2(neq), d1, slope(2), DDD(2), Pol1(neq), sou(neq)
        real                :: Pol1L(neq,-2:nx+2,-2:ny+2), Pol1R(neq,-2:nx+2,-2:ny+2)
        real                :: Pol1D(neq,-2:nx+2,-2:ny+2), Pol1U(neq,-2:nx+2,-2:ny+2), f(neq), di(2)
        real, allocatable   :: fk(:), gk(:), Dfk(:), Dgk(:), sk(:), QQQ(:)

        allocate(fk(neq), gk(neq), Dfk(neq), Dgk(neq), sk(neq), QQQ(neq))
        
        alpha=1.
        kappa=1.4
        iesq=1
        coef=0.
        DDD(1)=d_u; DDD(2)=d_w
        if (iesq==0) then
            fL=0.; fR=0.
            DO i = 0, nx
                DO j = 0, ny
                    slope1=0; slope2=0; slope=0
                    slope1(:)=(ua(:,i+1,j)-ua(:,i,j))/dx
                    slope2(:)=(ua(:,i,j)-ua(:,i-1,j))/dx
                    DO k = 1, 2
                        if (abs(slope1(k)) > abs(slope2(k))) then
                            slope(k)=slope2(k)
                        else
                            slope(k)=slope1(k)
                        end if
                    END DO
                    fDL(:,i,j,1)=DDD(:)*slope(:)
                    fDR(:,i,j,1)=DDD(:)*slope(:)

                    Pol1L(:,i,j)=ua(:,i,j)-slope(:)*dx/2.
                    Pol1R(:,i,j)=ua(:,i,j)+slope(:)*dx/2.
                    call flux_hyperb(x(i),y(j),slope(2)*dx,Pol1L(1,i,j),f)
                    fL(1,i,j,1)=f(1)
                    call flux_hyperb(x(i),y(j),slope(2)*dx,Pol1R(1,i,j),f)
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
                        if (abs(slope1(k)) > abs(slope2(k))) then
                            slope(k)=slope2(k)
                        else
                            slope(k)=slope1(k)
                        end if
                    END DO

                    !            GDFLUX(1,i,j)=d_u*slope(1)
                    !            GDFLUX(2,i,j)=d_w*slope(2)
                    gDD(:,i,j,1)=DDD(:)*slope(:)
                    gDU(:,i,j,1)=DDD(:)*slope(:)
                    Pol1D(:,i,j)=ua(:,i,j)-slope(:)*dy/2.
                    Pol1U(:,i,j)=ua(:,i,j)+slope(:)*dy/2.
                    call flux_hyperb(x(i),y(j),slope(2)*dy,Pol1D(1,i,j),f)
                    gD(1,i,j,1)=f(1)
                    call flux_hyperb(x(i),y(j),slope(2)*dy,Pol1U(1,i,j),f)
                    gU(1,i,j,1)=f(1)

                    Pol1(:)=ua(:,i,j)
                    call SourceTerm(0,Pol1,sou)
                    SSource(:,i,j)=sou(:)
                END DO
            END DO
            DO i = 1, nx
                DO j = 1, ny
                    FDFLUX(:,i,j)=(fDL(:,i,j,1)+fDR(:,i-1,j,1))/2.
                    GDFLUX(:,i,j)=(gDD(:,i,j,1)+gDU(:,i,j-1,1))/2.
                    dddf = FDFLUX(2,i,j)/d_w ! (DPolL_x(2,i,j,k)+DPolR_x(2,i-1,j,k))/2.
                    !                dddf = 1.
                    if (dddf*Xi>0.) then
                        FFlux(:,i,j)=fR(:,i-1,j,1)
                    else
                        FFlux(:,i,j)=fL(:,i,j,1)
                    end if

                    ! (DPolL_x(2,i,j,k)+DPolR_x(2,i-1,j,k))/2.
                    !FFlux(:,i,j)=0.5*(fR(:,i,j-1,1)+fL(:,i,j,1) - abs(Xi*dddf)* ( Pol1L(:,i,j) - Pol1R(:,i,j-1) ) )
                    !                dddf = 1.
                    dddf = GDFLUX(2,i,j)/d_w

                    if (dddf*Xi>0.) then
                        GFlux(:,i,j)=gU(:,i,j-1,1)
                    else
                        GFlux(:,i,j)=gD(:,i,j,1)
                    end if
                    !GFlux(:,i,j)=0.5*(gU(:,i,j-1,1)+gD(:,i,j,1) - abs(Xi*dddf)* ( Pol1D(:,i,j) - Pol1U(:,i,j-1) ) )      

                END DO
            END DO
        else
            call recons_weno_dim_dim(ua)
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
                        call pol_recons(0,2,xg,yg,coeff,Pol,DPolx)
                        d1=Dpolx(2)
                        call flux_hyperb(xg,yg,d1,Pol,fk)
                        fR(:,i,j,k) = fk(:)
                        Dfk=0.
                        Dgk=0.
                        call flux_diffusive(0,xg,yg,Pol,DPolx,Dfk)
                        !               DO kk = 1, nx
                        !                  call DUman(x(kk),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
                        !                  print*, 'Vol',i,j
                        !                  print*, kk,Dpolx(1)/dx,DUmanX
                        !               END DO
                        !               pause 333
                        !               !----------------------------------------------------------
                        !               print *, 'SOLUTION EAST SIDE OF CONTROL VOLUME'
                        !               print *, k,Pol(1),Uman(x(i+1),y(j)+ xig(k)*dy,t-dt)
                        !               print *, 'DERIVATIVES EAST SIDE OF CONTROL VOLUME'
                        !               call DUman(x(i+1),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
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
                        call pol_recons(0,4,xg,yg,coeff,Pol,DPolx)
                        DPolR_x(:,i,j,k)=DPolx(:)
                        d1=Dpolx(2)

                        call flux_hyperb(xg,yg,d1,Pol,fk)
                        fL(:,i,j,k) = fk(:)
                        Dfk=0.
                        Dgk=0.
                        call flux_diffusive(0,xg,yg,Pol,DPolx,Dfk)
                        !----------------------------------------------------------
                        !               print *, 'SOLUTION WEST SIDE OF CONTROL VOLUME'
                        !               print *, k,Pol(1),Uman(x(i),y(j)+ xig(k)*dy,t-dt)
                        !               print *, 'DERIVATIVES WEST SIDE OF CONTROL VOLUME'
                        !               call DUman(x(i),y(j)+ xig(k)*dy,t-dt,DUmanX,DUmanY)
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
                        call pol_recons(1,1,xg,yg,coeff,Pol,DPolx)
                        d1=Dpolx(2)

                        call flux_hyperb(xg,yg,d1,Pol,fk)
                        gD(:,i,j,k) = fk(:)
                        Dfk=0.
                        Dgk=0.
                        call flux_diffusive(1,xg,yg,Pol,DPolx,Dfk)
                        !----------------------------------------------------------
                        !               print *, 'SOLUTION SOUTH SIDE OF CONTROL VOLUME'
                        !               print *, k,Pol(1),Uman(x(i)+ xig(k)*dx,y(j),t-dt)
                        !               print *, 'DERIVATIVES SOUTH SIDE OF CONTROL VOLUME'
                        !               call DUman(x(i)+ xig(k)*dx,y(j),t-dt,DUmanX,DUmanY)
                        !               print *, 'x-derivative'
                        !               print *, k,Dpolx(1)/dy,DUmanY
                        !               !------------------------------------------------------------
                        !                call DUman(x(i)+ xig(k)*dx,y(j),t-dt,DUmanX,DUmanY)
                        !               print *, k,Dpolx(1)/dx,DUmanY
                        gDD(:,i,j,k) = Dfk(:)/dy
                        PolD(:,i,j,k)=Pol(:)
                        DPolD_y(:,i,j,k)=DPolx(:)
                    END DO
                    !    		pause
                    yg=1.
                    DO k = 1, NP1
                        xg=xig(k)
                        call pol_recons(1,3,xg,yg,coeff,Pol,DPolx)
                        d1=Dpolx(2)
                        call flux_hyperb(xg,yg,d1,Pol,fk)
                        gU(:,i,j,k) = fk(:)
                        Dfk=0.
                        Dgk=0.
                        call flux_diffusive(1,xg,yg,Pol,DPolx,Dfk)
                        !----------------------------------------------------------
                        !               print *, 'SOLUTION NORTH SIDE OF CONTROL VOLUME'
                        !               print *, k,Pol(1),Uman(x(i)+ xig(k)*dx,y(j+1),t-dt)
                        !               print *, 'DERIVATIVES NORTH SIDE OF CONTROL VOLUME'
                        !               call DUman(x(i)+ xig(k)*dx,y(j+1),t-dt,DUmanX,DUmanY)
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
            !	call DUman(x(2),y(29)+ xig(1)*dy,t-dt,DUmanX,DUmanY)
            !	print *, DUmanx*d_u

            FFO=0.
            GFO=0.
            ! Intercell Diffusive Flux
            FDifus=0.
            GDifus=0.
            call gauss_quadrature (NP1)
            !print *, 'valizq malo',FDL(1,19,29,1)
            !print *, 'valder malo',FDR(1,20,29,1)
            call ReconsFluxDifus(NP1,ua)
            call Forceold(ua,NP1)
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

            call gauss_quadrature (NP1)

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
                            call PolRecons(0,5,xg,yg,coeff,Pol,DPolx)
                            !            print *, Pol
                            call SourceTerm(j,Pol,sk)
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
        end if
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
        !        if (Xi*derw>0.) then
        !            FFLUX(:,i,j)=Xi*derw*ua(:,i-1,j)
        !        else
        !            FFLUX(:,i,j)=Xi*derw*ua(:,i,j)
        !        end if
        !        derw=(ua(2,i+1,j)-ua(2,i,j))/dx
        !        !derw=1.
        !        if (Xi*derw>0.) then
        !            FFLUX(:,i+1,j)=Xi*derw*ua(:,i,j)
        !        else
        !            FFLUX(:,i+1,j)=Xi*derw*ua(:,i+1,j)
        !        end if
        !        derw=(ua(2,i,j)-ua(2,i,j-1))/dy
        !        !derw=1.
        !        if (Xi*derw>0.) then
        !            GFLUX(:,i,j)=Xi*derw*ua(:,i,j-1)
        !        else
        !            GFLUX(:,i,j)=Xi*derw*ua(:,i,j)
        !        end if
        !        derw=(ua(2,i,j+1)-ua(2,i,j))/dy
        !       ! derw=1.
        !        if (Xi*derw>0.) then
        !            GFLUX(:,i,j+1)=Xi*derw*ua(:,i,j)
        !        else
        !            GFLUX(:,i,j+1)=Xi*derw*ua(:,i,j+1)
        !        end if
        !    END DO
        !END DO
        !SSource=0.
        DO i = 1, nx-1
            DO j = 1, ny-1
                if (WithManufSol==0) then
                    ForcingTerm=0.
                else
                    call SourceManufact(x(i),y(j),ForcingTerm)
                end if
                SSource(:,i,j)=SSource(:,i,j)+ForcingTerm(:)
                LO(:,i,j) = (FFLUX(:,i,j)-FFLUX(:,i+1,j))/dx + (GFLUX(:,i,j)-GFLUX(:,i,j+1))/dy
                LO(:,i,j) =LO(:,i,j)-(GDFLUX(:,i,j)-GDFLUX(:,i,j+1))/dy&
                        - (FDFLUX(:,i,j)-FDFLUX(:,i+1,j))/dx   +  SSource(:,i,j)
            END DO
        END DO
    end subroutine OpL
    
end module ch2_subroutines