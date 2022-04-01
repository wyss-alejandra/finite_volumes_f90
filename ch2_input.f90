module ch2_read_input
    implicit none

    !!! ------------------------------------------------------------------------------------------------------------ !!!
    !!!                                                 DATA INPUT                                                   !!!
    !!! ------------------------------------------------------------------------------------------------------------ !!!


    ! ---------------------------------------------------------------------------------------------------------------- !
    ! Data for upper boundary
    ! ---------------------------------------------------------------------------------------------------------------- !

    open(0,FILE='Data.dat', STATUS='UNKNOWN')
    read(0,*) text   ! Coordinates of domain
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


    end module ch2_read_input