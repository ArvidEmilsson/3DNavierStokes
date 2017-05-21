module bc_skew_flex

  use abk
  use parameter
  use stoffgroessen
  use geomdata
  use rhsdata
  use ableitung

  use initialize_time
  use initialize_boundaries

  public  :: set_bc_flux,cornersCorrection
  private :: get_phi, get_targets_nonref


contains
  subroutine set_bc_flux(u,u_sub,dt,kv,kStep,fluxes,first)
    real(kind=rk),dimension(:,:,:,:), intent(in)    :: u
    real(kind=rk),dimension(:,:,:,:), intent(in)    :: u_sub
    real(kind=rk), intent(in)                       :: dt
    real(kind=rk),dimension(:,:,:,:), intent(inout) :: kv
    integer,intent(in)                              :: kStep
    real(kind=rk),dimension(:,:,:,:),intent(inout)  :: fluxes
    integer                                         :: first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer               :: n1,n2,n3,n4
    real(kind=rk)         :: gamma,gasc,molmass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Allgemein
    real(kind=rk),dimension(:,:,:,:),allocatable    :: delta_fluxes    
    real(kind=rk),dimension(:,:,:,:),allocatable    :: delta_u    
    real(kind=rk),dimension(:,:,:,:),allocatable    :: u_test
    real(kind=rk)                                   :: dt_sub
    integer                                         :: ieq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Euler-Slip Wand
    real(kind=rk),dimension(:,:,:),allocatable      :: Ex,Ey,Ez
    real(kind=rk),allocatable,dimension(:,:,:)      :: sqrtRHOU_tangential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Isothermal No-Slip Wall
    real(kind=rk)                                   :: T0
    real(kind=rk)                                   :: TWall
    real(kind=rk),allocatable,dimension(:,:,:)      :: pTarget
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  adiabatic No-Slip Wall
    real(kind=rk),allocatable,dimension(:,:,:)      :: p_Flux
    real(kind=rk),allocatable,dimension(:,:,:,:)    :: phi_loc
    real(kind=rk)                                   :: dxi,deta,dzeta
    real(kind=rk),allocatable,dimension(:,:,:,:)    :: Ws
    real(kind=rk),allocatable,dimension(:)          :: ds
    real(kind=rk)                                   :: bndir, spVnr,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Nonreflecting BCs
    real(kind=rk)                                   :: rho0,p0,u0,relaxation_param
    real(kind=rk),allocatable,dimension(:,:,:,:)    :: uRef
    real(kind=rk),allocatable,dimension(:,:,:)      :: kx,ky,kz
    real(kind=rk),allocatable,dimension(:,:,:,:)    :: TargetVals, DeltaY
    real(kind=rk),allocatable,dimension(:,:,:)      :: u_eff, projectY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 'Fixed' type BC
    real(kind=rk),allocatable,dimension(:,:,:,:)    :: uTarget
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Variables for Flex structure
    integer                                         :: nob                                 ! Number of boundaries
    integer                                         :: is,ie,js,je,ks,ke                   ! local boundary surface
    integer                                         :: bi                                  ! current boundary
    character(len=80)                               :: BC_type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Collect general quantities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call get_param('n1',n1) 
    call get_param('n2',n2) 
    call get_param('n3',n3) 
    n4 = size(u,4)

    call get_param('gamma',gamma,default=1.4_rk)
    call get_param('T0',T0)
    call get_param('p_oo',p0)
    call get_param('gasc',gasc)
    call get_molmass(molmass)
    call get_param('rho_oo',rho0,default=287.0_rk*T0/p0)

!!!! TEST TEST TEST TEST !!!

    call give_boundaries(u,kStep,first)


!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!! 


    nob = boundaries%nob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Set up dt_sub and initial u_test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(delta_fluxes(n1,n2,n3,n4))
    allocate(delta_u(n1,n2,n3,n4))
    allocate(u_test(n1,n2,n3,n4)) ! Nötig?


    dt_sub       = dt

    u_test       = u + dt_sub*kv(:,:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Loop over boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do bi=1,nob

       BC_type = boundaries%bc_type(bi)


       delta_u      = nll
       delta_fluxes = nll





       if (trim(BC_type).eq.'wall_adiabatic') then

          allocate(phi_loc(n1,n2,n3,3))
          call get_phi(u_test,phi_loc)

       end if


       if(boundaries%belong(bi).eqv..true.) then

          is      = boundaries%is(bi)
          ie      = boundaries%ie(bi)
          js      = boundaries%js(bi)
          je      = boundaries%je(bi)
          ks      = boundaries%ks(bi)
          ke      = boundaries%ke(bi)

          relaxation_param = boundaries%relaxation_param(bi)

          select case(trim(BC_type))

          case("leer")

          case("wall_euler")


             allocate(sqrtRHOU_tangential(n1,n2,n3))
             allocate(Ex(n1,n2,n3))
             allocate(Ey(n1,n2,n3))
             allocate(Ez(n1,n2,n3))

             sqrtRHOU_tangential = nll

             Ex = boundaries%Ex(:,:,:,bi)
             Ey = boundaries%Ey(:,:,:,bi)
             Ez = boundaries%Ez(:,:,:,bi)

             sqrtRHOU_tangential(is:ie,js:je,ks:ke) = u_test(is:ie,js:je,ks:ke,2)*Ex(is:ie,js:je,ks:ke)+u_test(is:ie,js:je,ks:ke,3)*Ey(is:ie,js:je,ks:ke)+u_test(is:ie,js:je,ks:ke,4)*Ez(is:ie,js:je,ks:ke)

             delta_u(is:ie,js:je,ks:ke,2)    = delta_u(is:ie,js:je,ks:ke,2) + Ex(is:ie,js:je,ks:ke)*sqrtRHOU_tangential(is:ie,js:je,ks:ke)  - u_test(is:ie,js:je,ks:ke,2)
             delta_u(is:ie,js:je,ks:ke,3)    = delta_u(is:ie,js:je,ks:ke,3) + Ey(is:ie,js:je,ks:ke)*sqrtRHOU_tangential(is:ie,js:je,ks:ke)  - u_test(is:ie,js:je,ks:ke,3)
             delta_u(is:ie,js:je,ks:ke,4)    = delta_u(is:ie,js:je,ks:ke,4) + Ez(is:ie,js:je,ks:ke)*sqrtRHOU_tangential(is:ie,js:je,ks:ke)  - u_test(is:ie,js:je,ks:ke,4) 

             u_test(is:ie,js:je,ks:ke,:) = u_test(is:ie,js:je,ks:ke,:) + delta_u(is:ie,js:je,ks:ke,:)

             deallocate(sqrtRHOU_tangential,Ex,Ey,Ez)


          case("wall_isothermal")


             TWall = boundaries%TWall(bi)

             allocate(pTarget(n1,n2,n3))

             pTarget = TWall*u_test(is:ie,js:je,ks:ke,1)**2*gasc/molmass



             delta_u(is:ie,js:je,ks:ke,2)    = delta_u(is:ie,js:je,ks:ke,2) + nll  - u_test(is:ie,js:je,ks:ke,2)
             delta_u(is:ie,js:je,ks:ke,3)    = delta_u(is:ie,js:je,ks:ke,3) + nll  - u_test(is:ie,js:je,ks:ke,3)
             delta_u(is:ie,js:je,ks:ke,4)    = delta_u(is:ie,js:je,ks:ke,4) + nll  - u_test(is:ie,js:je,ks:ke,4) 
             delta_u(is:ie,js:je,ks:ke,5)    = delta_u(is:ie,js:je,ks:ke,5) + pTarget - u_test(is:ie,js:je,ks:ke,5) 

             u_test(is:ie,js:je,ks:ke,:) = u_test(is:ie,js:je,ks:ke,:) + delta_u(is:ie,js:je,ks:ke,:)

             deallocate(pTarget)

          case("wall_adiabatic")

             delta_u(is:ie,js:je,ks:ke,2)    = delta_u(is:ie,js:je,ks:ke,2) + nll  - u_test(is:ie,js:je,ks:ke,2)
             delta_u(is:ie,js:je,ks:ke,3)    = delta_u(is:ie,js:je,ks:ke,3) + nll  - u_test(is:ie,js:je,ks:ke,3)
             delta_u(is:ie,js:je,ks:ke,4)    = delta_u(is:ie,js:je,ks:ke,4) + nll  - u_test(is:ie,js:je,ks:ke,4) 

             u_test(is:ie,js:je,ks:ke,:) = u_test(is:ie,js:je,ks:ke,:) + delta_u(is:ie,js:je,ks:ke,:)

          case("fixed")


             allocate(uTarget(n1,n2,n3,n4))

             uTarget = boundaries%uRef(:,:,:,:,bi)

!!$             uTarget(:,:,:,1) = sqrt(uTarget(:,:,:,1))
!!$             uTarget(:,:,:,2) = uTarget(:,:,:,2)*uTarget(:,:,:,1)
!!$             uTarget(:,:,:,3) = uTarget(:,:,:,3)*uTarget(:,:,:,1)
!!$             uTarget(:,:,:,4) = uTarget(:,:,:,4)*uTarget(:,:,:,1)


             delta_u(is:ie,js:je,ks:ke,1)    = delta_u(is:ie,js:je,ks:ke,1) +  uTarget(is:ie,js:je,ks:ke,1)  - u_test(is:ie,js:je,ks:ke,1)
             delta_u(is:ie,js:je,ks:ke,2)    = delta_u(is:ie,js:je,ks:ke,2) +  uTarget(is:ie,js:je,ks:ke,2)  - u_test(is:ie,js:je,ks:ke,2)
             delta_u(is:ie,js:je,ks:ke,3)    = delta_u(is:ie,js:je,ks:ke,3) +  uTarget(is:ie,js:je,ks:ke,3)  - u_test(is:ie,js:je,ks:ke,3)
             delta_u(is:ie,js:je,ks:ke,4)    = delta_u(is:ie,js:je,ks:ke,4) +  uTarget(is:ie,js:je,ks:ke,4)  - u_test(is:ie,js:je,ks:ke,4) 
             delta_u(is:ie,js:je,ks:ke,5)    = delta_u(is:ie,js:je,ks:ke,5) +  uTarget(is:ie,js:je,ks:ke,5)  - u_test(is:ie,js:je,ks:ke,5) 
             
	     do ieq = 6,n4
	        delta_u(is:ie,js:je,ks:ke,ieq)    = delta_u(is:ie,js:je,ks:ke,ieq) +  uTarget(is:ie,js:je,ks:ke,ieq)  - u_test(is:ie,js:je,ks:ke,ieq) 
             enddo

             u_test(is:ie,js:je,ks:ke,:) = u_test(is:ie,js:je,ks:ke,:) + delta_u(is:ie,js:je,ks:ke,:)

             deallocate(uTarget)



          case("nonreflecting")


             allocate(TargetVals(n1,n2,n3,n4))
             allocate(uRef(n1,n2,n3,n4))
             allocate(kx(n1,n2,n3))
             allocate(ky(n1,n2,n3))
             allocate(kz(n1,n2,n3))
             allocate(u_eff(n1,n2,n3))
             allocate(projectY(n1,n2,n3))
             allocate(DeltaY(n1,n2,n3,n4))

             uRef   =  boundaries%uRef(:,:,:,:,bi)
             kx     = -boundaries%kx(:,:,:,bi)
             ky     = -boundaries%ky(:,:,:,bi)
             kz     = -boundaries%kz(:,:,:,bi)

             TargetVals = nll
             DeltaY     = nll
             projectY   = nll
             u_eff      = nll

             call get_targets_nonref(TargetVals,u_test,kx,ky,kz,uRef,bi,relaxation_param)

             do ieq = 6,n4

                u_eff(is:ie,js:je,ks:ke) = -(kx(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,2) +  & 
                                             ky(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,3) +  &
                                             kz(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,4) )/u_sub(is:ie,js:je,ks:ke,1) !! Velocity at the boundary
                
                where(u_eff.lt.nll) projectY  = eins                           !! Mask on fluxes for Passive Scalar
 
                TargetVals(:,:,:,ieq)          = boundaries%uRef(:,:,:,ieq,bi)*TargetVals(:,:,:,1)**2  !! Passive Scalar Target: rho*Y
                !DeltaY(is:ie,js:je,ks:ke,ieq) = projectY*(TargetVals(is:ie,js:je,ks:ke,ieq)  - u_test(is:ie,js:je,ks:ke,ieq))
                DeltaY(is:ie,js:je,ks:ke,ieq)  =          (TargetVals(is:ie,js:je,ks:ke,ieq)  - u_test(is:ie,js:je,ks:ke,ieq))

             enddo

	     !! target minus current val
             delta_u(is:ie,js:je,ks:ke,1)    = delta_u(is:ie,js:je,ks:ke,1) + TargetVals(is:ie,js:je,ks:ke,1)  - u_test(is:ie,js:je,ks:ke,1)
             delta_u(is:ie,js:je,ks:ke,2)    = delta_u(is:ie,js:je,ks:ke,2) + TargetVals(is:ie,js:je,ks:ke,2)  - u_test(is:ie,js:je,ks:ke,2)
             delta_u(is:ie,js:je,ks:ke,3)    = delta_u(is:ie,js:je,ks:ke,3) + TargetVals(is:ie,js:je,ks:ke,3)  - u_test(is:ie,js:je,ks:ke,3)
             delta_u(is:ie,js:je,ks:ke,4)    = delta_u(is:ie,js:je,ks:ke,4) + TargetVals(is:ie,js:je,ks:ke,4)  - u_test(is:ie,js:je,ks:ke,4) 
             delta_u(is:ie,js:je,ks:ke,5)    = delta_u(is:ie,js:je,ks:ke,5) + TargetVals(is:ie,js:je,ks:ke,5)  - u_test(is:ie,js:je,ks:ke,5)
                  
             do ieq = 6,n4
                 delta_u(is:ie,js:je,ks:ke,ieq)    = delta_u(is:ie,js:je,ks:ke,ieq) + DeltaY(is:ie,js:je,ks:ke,ieq)
             enddo

             u_test(is:ie,js:je,ks:ke,:) = u_test(is:ie,js:je,ks:ke,:) + delta_u(is:ie,js:je,ks:ke,:)

             deallocate(TargetVals,kx,ky,kz,uRef)
             deallocate(u_eff,projectY,DeltaY) 

          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! Updating delta_fluxes !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          delta_fluxes(is:ie,js:je,ks:ke,1)   = J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1)*delta_u(is:ie,js:je,ks:ke,1)*zwei/dt_sub 
          delta_fluxes(is:ie,js:je,ks:ke,2)   = J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1)*delta_u(is:ie,js:je,ks:ke,2)/dt_sub
          delta_fluxes(is:ie,js:je,ks:ke,3)   = J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1)*delta_u(is:ie,js:je,ks:ke,3)/dt_sub
          delta_fluxes(is:ie,js:je,ks:ke,4)   = J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1)*delta_u(is:ie,js:je,ks:ke,4)/dt_sub
          delta_fluxes(is:ie,js:je,ks:ke,5)   = J_tilde(is:ie,js:je,ks:ke)*delta_u(is:ie,js:je,ks:ke,5)/((gamma-eins)*dt_sub)
          
          do ieq = 6,n4
             delta_fluxes(is:ie,js:je,ks:ke,ieq)   = J_tilde(is:ie,js:je,ks:ke)*delta_u(is:ie,js:je,ks:ke,ieq)/dt_sub
          enddo 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Setting BCs using Fluxes !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Adiabatic wall !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (trim(BC_type).eq.'wall_adiabatic') then
             allocate(p_Flux(n1,n2,n3))
             allocate(Ws(n1,n2,n3,3))
             allocate(ds(3))

             p_FLux = nll

             ds(1) = eins/real(i1*n1-1,rk)
             ds(2) = eins/real(i2*n2-1,rk)
             ds(3) = eins/real(i3*n3-1,rk)

             Ws(:,:,:,1) = W_eta*W_zeta
             Ws(:,:,:,2) = W_xi*W_zeta
             Ws(:,:,:,3) = W_xi*W_eta

             do k=1,3

                bndir = boundaries%boundary_directions(k,bi)
                spVnr = boundaries%space_variables(k,bi)


                if (spVnr.ne.0) then 
                   p_Flux(is:ie,js:je,ks:ke) = p_flux(is:ie,js:je,ks:ke) + ( bndir*Ws(is:ie,js:je,ks:ke,spVnr)*phi_loc(is:ie,js:je,ks:ke,spVnr))/ds(spVnr)
                end if

             end do

             delta_fluxes(is:ie,js:je,ks:ke,5) = - ( p_Flux(is:ie,js:je,ks:ke) + fluxes(is:ie,js:je,ks:ke,5) )

             deallocate(p_Flux,Ws,ds)

          end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  Setting new Fluxes   !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          fluxes(is:ie,js:je,ks:ke,:) = fluxes(is:ie,js:je,ks:ke,:) + delta_fluxes(is:ie,js:je,ks:ke,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!  Setting new kv   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          kv(is:ie,js:je,ks:ke,1)     = kv(is:ie,js:je,ks:ke,1) +  delta_fluxes(is:ie,js:je,ks:ke,1)/(J_tilde(is:ie,js:je,ks:ke)*zwei*u_sub(is:ie,js:je,ks:ke,1))
          kv(is:ie,js:je,ks:ke,2)     = kv(is:ie,js:je,ks:ke,2) +  delta_fluxes(is:ie,js:je,ks:ke,2)/(J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1))
          kv(is:ie,js:je,ks:ke,3)     = kv(is:ie,js:je,ks:ke,3) +  delta_fluxes(is:ie,js:je,ks:ke,3)/(J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1))
          kv(is:ie,js:je,ks:ke,4)     = kv(is:ie,js:je,ks:ke,4) +  delta_fluxes(is:ie,js:je,ks:ke,4)/(J_tilde(is:ie,js:je,ks:ke)*u_sub(is:ie,js:je,ks:ke,1))
          kv(is:ie,js:je,ks:ke,5)     = kv(is:ie,js:je,ks:ke,5) + (delta_fluxes(is:ie,js:je,ks:ke,5)/J_tilde(is:ie,js:je,ks:ke)) * (gamma-1)

          do ieq = 6,n4
             kv(is:ie,js:je,ks:ke,ieq) = kv(is:ie,js:je,ks:ke,ieq) +  delta_fluxes(is:ie,js:je,ks:ke,ieq)/J_tilde(is:ie,js:je,ks:ke)
          enddo


       end if

       if (trim(BC_type).eq.'wall_adiabatic') then

          deallocate(phi_loc)

       end if


    end do


  end subroutine set_bc_flux


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine get_targets_nonref(TargetVals,u,kx,ky,kz,uRef,bi,relaxation_param)
    real(kind=rk),dimension(:,:,:,:),intent(out)  :: TargetVals
    real(kind=rk),dimension(:,:,:,:),intent(in)   :: u,uRef
    real(kind=rk),dimension(:,:,:)                :: kx,ky,kz
    integer                                       :: bi
    real(kind=rk)                                 :: relaxation_param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),allocatable  :: primeU
    real(kind=rk),dimension(:,:,:,:),allocatable  :: du,duN
    real(kind=rk),dimension(:,:,:),allocatable    :: R1,R2,R3,R4,R5
    real(kind=rk),dimension(:,:,:),allocatable    :: rhoR,uR,vR,wR,pR
    real(kind=rk),dimension(:,:,:),allocatable    :: lamS,lamP,lamM
    real(kind=rk),dimension(:,:,:),allocatable    :: c,M
    real(kind=rk)                                 :: c0,gamma,po
    integer                                       :: n1,n2,n3,n4



    n1 = size(u,1)
    n2 = size(u,2)
    n3 = size(u,3)
    n4 = size(u,4)

    call get_param('gamma',gamma,default=1.4_rk)



    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(rhoR(n1,n2,n3),uR(n1,n2,n3),vR(n1,n2,n3),wR(n1,n2,n3),pR(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))
    allocate(du(n1,n2,n3,n4),duN(n1,n2,n3,n4),primeU(n1,n2,n3,n4))
    allocate(c(n1,n2,n3),M(n1,n2,n3))


    du  = nll
    duN = nll

    rhoR = uRef(:,:,:,1)
    uR   = uRef(:,:,:,2)
    vR   = uRef(:,:,:,3)
    wR   = uRef(:,:,:,4)
    pR   = uRef(:,:,:,5)


    primeU(:,:,:,1) = u(:,:,:,1)**2
    primeU(:,:,:,2) = u(:,:,:,2)/u(:,:,:,1)
    primeU(:,:,:,3) = u(:,:,:,3)/u(:,:,:,1)
    primeU(:,:,:,4) = u(:,:,:,4)/u(:,:,:,1)
    primeU(:,:,:,5) = u(:,:,:,5)

    c  = sqrt(gamma*(primeU(:,:,:,5)/primeU(:,:,:,1)))


    du(:,:,:,1) = PrimeU(:,:,:,1) - rhoR
    du(:,:,:,2) = PrimeU(:,:,:,2) - uR
    du(:,:,:,3) = PrimeU(:,:,:,3) - vR
    du(:,:,:,4) = PrimeU(:,:,:,4) - wR
    du(:,:,:,5) = PrimeU(:,:,:,5) - pR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    R1 = kx*du(:,:,:,1)                     +  kz*du(:,:,:,3)  -  ky*du(:,:,:,4)  -            kx*du(:,:,:,5)/(c**2)
    R2 = ky*du(:,:,:,1)  -  kz*du(:,:,:,2)                     +  kx*du(:,:,:,4)  -            ky*du(:,:,:,5)/(c**2)
    R3 = kz*du(:,:,:,1)  +  ky*du(:,:,:,2)  -  kx*du(:,:,:,3)                     -            kz*du(:,:,:,5)/(c**2)
    R4 =                    kx*du(:,:,:,2)  +  ky*du(:,:,:,3)  +  kz*du(:,:,:,4)  +  du(:,:,:,5)/(PrimeU(:,:,:,1)*c)
    R5 =                 -  kx*du(:,:,:,2)  -  ky*du(:,:,:,3)  -  kz*du(:,:,:,4)  +  du(:,:,:,5)/(PrimeU(:,:,:,1)*c)


    lamS = kx*primeU(:,:,:,2) + ky*primeU(:,:,:,3) + kz*primeU(:,:,:,4) 
    lamP = kx*primeU(:,:,:,2) + ky*primeU(:,:,:,3) + kz*primeU(:,:,:,4)  + c
    lamM = kx*primeU(:,:,:,2) + ky*primeU(:,:,:,3) + kz*primeU(:,:,:,4)  - c


    where(lamS.gt.nll) R1 = nll
    where(lamS.gt.nll) R2 = nll
    where(lamS.gt.nll) R3 = nll
    where(lamP.gt.nll) R4 = nll
    where(lamM.gt.nll) R5 = nll



    if(relaxation_param.ne.nll) then                                          ! relaxation_param = dx/L * sigma

!!$       R5 = zehn/(64.0_rk)/(primeU(:,:,:,2)-c) * 0.58_rk*(1-0.57_rk**2)/(zehn) * (primeU(:,:,:,5)-100000.0_rk)/(primeU(:,:,:,1))

       M  = sqrt( primeU(:,:,:,2)**2 + primeU(:,:,:,3)**2 + primeU(:,:,:,4)**2 )/c
       call get_param('po',po)

       R5 = relaxation_param/(primeU(:,:,:,2)-c)*(eins-maxval(M)**2) * (primeU(:,:,:,5)-po)/(primeU(:,:,:,1))

    end if




    duN(:,:,:,1) =   kx*R1  +  ky*R2  +  kz*R3  +  PrimeU(:,:,:,1)/(zwei*c)*R4  +  PrimeU(:,:,:,1)/(zwei*c)*R5
    duN(:,:,:,2) =          -  kz*R2  +  ky*R3  +                   halb*kx*R4  -                   halb*kx*R5
    duN(:,:,:,3) =   kz*R1            -  kx*R3  +                   halb*ky*R4  -                   halb*ky*R5
    duN(:,:,:,4) =  -ky*R1  +  kx *R2           +                   halb*kz*R4  -                   halb*kz*R5
    duN(:,:,:,5) =                              +    halb*PrimeU(:,:,:,1)*c*R4  +    halb*PrimeU(:,:,:,1)*c*R5


    targetVals(:,:,:,1) = rhoR + duN(:,:,:,1)
    targetVals(:,:,:,2) = uR   + duN(:,:,:,2)
    targetVals(:,:,:,3) = vR   + duN(:,:,:,3)
    targetVals(:,:,:,4) = wR   + duN(:,:,:,4)
    targetVals(:,:,:,5) = pR   + duN(:,:,:,5)



    targetVals(:,:,:,1) = sqrt(targetVals(:,:,:,1))
    targetVals(:,:,:,2) = targetVals(:,:,:,2)*targetVals(:,:,:,1)
    targetVals(:,:,:,3) = targetVals(:,:,:,3)*targetVals(:,:,:,1)
    targetVals(:,:,:,4) = targetVals(:,:,:,4)*targetVals(:,:,:,1)



    deallocate(R1,R2,R3,R4,R5,lamS,lamP,lamM)
    deallocate(du,duN,PrimeU,c,M)

  end subroutine get_targets_nonref


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine get_phi(u,phi)
    real(kind=rk),dimension(:,:,:,:),intent(in)   :: u
    real(kind=rk),dimension(:,:,:,:),intent(out)  :: phi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),allocatable  :: phi_loc
    real(kind=rk),dimension(:,:,:),allocatable    :: T,T_x,T_y,T_z
    integer                                       :: n1,n2,n3
    real(kind=rk),dimension(:,:,:),allocatable    :: lambda,mu
    real(kind=rk)                                 :: Pr,Cp

    n1 = size(u,1)
    n2 = size(u,2)
    n3 = size(u,3)

    allocate(T(n1,n2,n3))
    allocate(T_x(n1,n2,n3))
    allocate(T_y(n1,n2,n3))
    allocate(T_z(n1,n2,n3))
    allocate(phi_loc(n1,n2,n3,3))
    allocate(lambda(n1,n2,n3))
    allocate(mu(n1,n2,n3))



    phi_loc = nll

    call get_T(u,T)
    call visk(T,mu)

    call get_param('Pr',Pr,default=0.71_rk)
    call get_Cp(Cp)

    lambda  = Cp*mu/Pr



    call grad_zentral(T,T_x,T_y,T_z)

    T_x = T_x/W_xi
    T_y = T_y/W_eta
    T_z = T_z/W_zeta


    phi_loc(:,:,:,1) = lambda*(m11*T_x + m21*T_y + m31*T_z)
    phi_loc(:,:,:,2) = lambda*(m12*T_x + m22*T_y + m32*T_z)
    phi_loc(:,:,:,3) = lambda*(m13*T_x + m23*T_y + m33*T_z)

    phi_loc(:,:,:,1) = phi_loc(:,:,:,1)/Jacobian
    phi_loc(:,:,:,2) = phi_loc(:,:,:,2)/Jacobian
    phi_loc(:,:,:,3) = phi_loc(:,:,:,3)/Jacobian

    phi(:,:,:,1) =  (m11*phi_loc(:,:,:,1) + m12*phi_loc(:,:,:,2) + m13*phi_loc(:,:,:,3))   
    phi(:,:,:,2) =  (m21*phi_loc(:,:,:,1) + m22*phi_loc(:,:,:,2) + m23*phi_loc(:,:,:,3))   
    phi(:,:,:,3) =  (m31*phi_loc(:,:,:,1) + m32*phi_loc(:,:,:,2) + m33*phi_loc(:,:,:,3))   
!!$    phi = phi_loc



    deallocate(T,T_x,T_y,T_z,phi_loc,lambda)


  end subroutine get_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cornersCorrection(u)  !!! correct wall adiabat + non-reflect corner in species
    real(kind=rk),dimension(:,:,:,:),intent(inout) :: u
    real(kind=rk),dimension(:,:,:,:),allocatable   :: temp
    real(kind=rk),dimension(2)                     :: rhoYs
    integer                                        :: nob     ! Number of boundaries
    integer                                        :: bi      ! current boundary
    integer                                        :: n1, n2, n3, n4, jj, j
    integer                                        :: i1, i2, j1, j2, k1, k2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1    = size(u,1)
    n2    = size(u,2)
    n3    = size(u,3)
    n4    = size(u,4)

    allocate(temp(n1,n2,n3,n4))
    temp  = nll
    nob   = boundaries%nob
    rhoYs = [n4-1,n4]    !! Species Fields to correct -- Products Outflow
    jj    = size(rhoYs) 

    do bi = 1,nob

       if (boundaries%corner_Corr(bi).eq..true.) then

          i1 = boundaries%is(bi)
          i2 = boundaries%ie(bi)
          j1 = boundaries%js(bi)
          j2 = boundaries%je(bi)
          k1 = boundaries%ks(bi)
          k2 = boundaries%ke(bi)

          do j = 1,jj
	    !! Right ---
	    temp(i2,j1,k1:k2,rhoYs(j)) = u(i2-1,j1,k1:k2,rhoYs(j))/u(i2-1,j1,k1:k2,1)**2	    
            u(i2,j1,k1:k2,rhoYs(j))    = temp(i2,j1,k1:k2,rhoYs(j))*u(i2,j1,k1:k2,1)**2

	    !! Left ---
	    temp(i1,j1,k1:k2,rhoYs(j)) = u(i1+1,j1,k1:k2,rhoYs(j))/u(i1+1,j1,k1:k2,1)**2
	    u(i1,j1,k1:k2,rhoYs(j))    = temp(i1,j1,k1:k2,rhoYs(j))*u(i1,j1,k1:k2,1)**2
	  enddo
       endif
    enddo

    deallocate(temp)

  end subroutine cornersCorrection

end module bc_skew_flex
