module iterationstep_skew
!!!
!!! Navier Stokes equation in skew formulation and variables ( sqrt(rho), sqrt(rho)u, p )
!!! on distorted grids  ||||| 
!!!
#include "parameter.h"
!!!
  use abk  
  use ableitung
  use parameter
  use stoffgroessen
  use blas
  use geomdata
  use rhsdata
!!! --- OPTIMIZATION ---
  use OPT_CASE
  use statistik
  use sponge


  use initialize_time
  use bc_skew_flex

!!! For Debugging
!!!  use ioprimitives
!!!  use ieee
!!! TMP


  public  :: iteration_step
  private :: iteration_step_real, iteration_step_cmplx
  private :: get_u_sub

  interface iteration_step
     module procedure iteration_step_real, iteration_step_cmplx
  end interface iteration_step

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine iteration_step_real(u,kv,fluxes,dt,kStep,scheme,scheme_type,w_T)
    real(kind=rk), dimension(:,:,:,:),intent(in)             :: u
    real(kind=rk), dimension(:,:,:,:,:),intent(inout)        :: kv
    real(kind=rk), dimension(:,:,:,:),intent(inout)          :: fluxes
    real(kind=rk)                                            :: dt
    integer                                                  :: kStep
    character(len=80)                                        :: scheme,scheme_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Constants
    integer                                        :: i,j,k,n1,n2,n3,b1,b2,meqn
    real(kind=rk)                                  :: gamma,gamma1, Cv, Cp, gasc, Pr
    integer                                        :: nt,STEP,nunit
    character(len=10)                              :: ntstring
    character(len=80)                              :: fname
    character(len=80)                              :: name
    integer                                        :: avs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Field
    real(kind=rk),dimension(:,:,:,:),allocatable   :: u_sub
    real(kind=rk),dimension(:,:,:,:),allocatable   :: velocities
    real(kind=rk),dimension(:,:,:),allocatable     :: velx,vely,velz
    real(kind=rk),dimension(:,:,:),allocatable     :: velxt,velyt,velzt,divu
    real(kind=rk),dimension(:,:,:),allocatable     :: tau11,tau22,tau33,tau12,tau13,tau23
    real(kind=rk),dimension(:,:,:),allocatable     :: tau11t,tau12t,tau13t
    real(kind=rk),dimension(:,:,:),allocatable     :: tau21t,tau22t,tau23t
    real(kind=rk),dimension(:,:,:),allocatable     :: tau31t,tau32t,tau33t
    real(kind=rk),dimension(:,:,:),allocatable     :: phi1,phi2,phi3
    real(kind=rk),dimension(:,:,:),allocatable     :: mu,mu_d,lambda
    real(kind=rk),dimension(:,:,:),allocatable     :: fric_u,fric_v,fric_w,fric_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Derivative fields
    real(kind=rk),dimension(:,:,:),allocatable     :: velx_x,velx_y,velx_z
    real(kind=rk),dimension(:,:,:),allocatable     :: vely_x,vely_y,vely_z
    real(kind=rk),dimension(:,:,:),allocatable     :: velz_x,velz_y,velz_z
    real(kind=rk),dimension(:,:,:),allocatable     :: p_x,p_y,p_z,T_x,T_y,T_z
    real(kind=rk),dimension(:,:,:),allocatable     :: dummy
    character(len=keyword_length)                  :: Gittertypus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! boundary
    real(kind=rk),dimension(:,:,:),allocatable      :: c    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Chemie
    real(kind=rk),dimension(:,:,:), allocatable,intent(inout),optional    :: w_T

    call get_param('Gittertypus', Gittertypus,default='kartesisch')

    call get_param('nt',nt)

    call allocate_rhs()

    n1   = size(u,1)
    n2   = size(u,2)
    n3   = size(u,3)
    meqn = size(u,4)

    call get_param('gamma',gamma)
    gamma1 = gamma-eins
    call get_param('gasc',gasc)
    call get_Cv(Cv)
    call get_Cp(Cp)
    call get_param('Pr',Pr,default=0.71_rk)
    call get_param('STEP',STEP)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocating variables
    allocate(u_sub(n1,n2,n3,meqn))
    allocate(velocities(n1,n2,n3,3))
    allocate(velx(n1,n2,n3),vely(n1,n2,n3),velz(n1,n2,n3))
    allocate(velxt(n1,n2,n3),velyt(n1,n2,n3),velzt(n1,n2,n3))
    allocate(divu(n1,n2,n3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocating derivatives
    allocate(velx_x(n1,n2,n3),velx_y(n1,n2,n3),velx_z(n1,n2,n3))
    allocate(vely_x(n1,n2,n3),vely_y(n1,n2,n3),vely_z(n1,n2,n3))
    allocate(velz_x(n1,n2,n3),velz_y(n1,n2,n3),velz_z(n1,n2,n3))
    allocate(p_x(n1,n2,n3),p_y(n1,n2,n3),p_z(n1,n2,n3))
    allocate(T_x(n1,n2,n3),T_y(n1,n2,n3),T_z(n1,n2,n3))

    allocate(dummy(n1,n2,n3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocate friction terms
    allocate(mu(n1,n2,n3))
    allocate(mu_d(n1,n2,n3))
    allocate(lambda(n1,n2,n3))

    allocate(tau11(n1,n2,n3),tau22(n1,n2,n3),tau33(n1,n2,n3),tau12(n1,n2,n3),tau13(n1,n2,n3),tau23(n1,n2,n3))
    allocate(tau11t(n1,n2,n3),tau12t(n1,n2,n3),tau13t(n1,n2,n3))
    allocate(tau21t(n1,n2,n3),tau22t(n1,n2,n3),tau23t(n1,n2,n3))
    allocate(tau31t(n1,n2,n3),tau32t(n1,n2,n3),tau33t(n1,n2,n3))
    allocate(fric_u(n1,n2,n3),fric_v(n1,n2,n3),fric_w(n1,n2,n3),fric_p(n1,n2,n3))
    allocate(phi1(n1,n2,n3),phi2(n1,n2,n3),phi3(n1,n2,n3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocate boundary terms

    !    allocate(flux_rho(n1,n2,n3),flux_u(n1,n2,n3),flux_v(n1,n2,n3),flux_w(n1,n2,n3),flux_p(n1,n2,n3))
    !    allocate(c(n1,n2,n3),bc_tau1(n1,n2,n3),bc_tau2(n1,n2,n3),bc_tau3(n1,n2,n3),bc_phi(n1,n2,n3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Setting up the time integration scheme  !!!!!!!!!!!!!!!
!!!!!!!!!        and half values                   !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call get_u_sub(u,kv,dt,scheme,kStep,u_sub)



    call get_rho(u_sub,rho)
    call get_uvw(u_sub,velocities)
    velx = velocities(:,:,:,1)
    vely = velocities(:,:,:,2)
    velz = velocities(:,:,:,3)

    call grad_zentral(velx,velx_x,velx_y,velx_z)
    call grad_zentral(vely,vely_x,vely_y,vely_z)
    call grad_zentral(velz,velz_x,velz_y,velz_z)

    velx_x = velx_x/W_xi
    velx_y = velx_y/W_eta
    velx_z = velx_z/W_zeta

    vely_x = vely_x/W_xi
    vely_y = vely_y/W_eta
    vely_z = vely_z/W_zeta

    velz_x = velz_x/W_xi
    velz_y = velz_y/W_eta
    velz_z = velz_z/W_zeta

!!!!!!!!!!!!!!!!!!!!!!! Computing projected velocities components

    velxt = m11*velx + m12*vely + m13*velz
    velyt = m21*velx + m22*vely + m23*velz
    velzt = m31*velx + m32*vely + m33*velz

!!!!!!! Computing derivatives of p with respect to physical coordinates

    if (trim(Gittertypus).eq.'cylindrical') then

       call diff1x_zentral(u_sub(:,:,:,5),dummy)
       p_x = m11*dummy/W_xi
       call diff1y_zentral(u_sub(:,:,:,5),dummy)
       p_x = p_x + m21*dummy/W_eta
       call diff1z_zentral(u_sub(:,:,:,5),dummy)
       p_x = p_x + m31*dummy/W_zeta
   
       call diff1x_zentral(u_sub(:,:,:,5),dummy)
       p_y = m12*dummy/W_xi
       call diff1y_zentral(u_sub(:,:,:,5),dummy)
       p_y = p_y + m22*dummy/W_eta
       call diff1z_zentral(u_sub(:,:,:,5),dummy)
       p_y = p_y + m32*dummy/W_zeta
   
       call diff1x_zentral(u_sub(:,:,:,5),dummy)
       p_z = m13*dummy/W_xi
       call diff1y_zentral(u_sub(:,:,:,5),dummy)
       p_z = p_z + m23*dummy/W_eta
       call diff1z_zentral(u_sub(:,:,:,5),dummy)
       p_z = p_z + m33*dummy/W_zeta

    else

       call diff1x_zentral(m11*u_sub(:,:,:,5),dummy)
       p_x = dummy/W_xi
       call diff1y_zentral(m21*u_sub(:,:,:,5),dummy)
       p_x = p_x + dummy/W_eta
       call diff1z_zentral(m31*u_sub(:,:,:,5),dummy)
       p_x = p_x + dummy/W_zeta
   
       call diff1x_zentral(m12*u_sub(:,:,:,5),dummy)
       p_y = dummy/W_xi
       call diff1y_zentral(m22*u_sub(:,:,:,5),dummy)
       p_y = p_y + dummy/W_eta
       call diff1z_zentral(m32*u_sub(:,:,:,5),dummy)
       p_y = p_y + dummy/W_zeta
   
       call diff1x_zentral(m13*u_sub(:,:,:,5),dummy)
       p_z = dummy/W_xi
       call diff1y_zentral(m23*u_sub(:,:,:,5),dummy)
       p_z = p_z + dummy/W_eta
       call diff1z_zentral(m33*u_sub(:,:,:,5),dummy)
       p_z = p_z + dummy/W_zeta

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!  Statistics and Stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DONE BEFORE RHS Evaluation to have all statistics available for calculation of BC fluxes (esp. Recycling!!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_param('av_start',avs,default=1)
    if( STEP.eq.1.and.(kStep==1).and.(nt.ge.avs).and.is_true('statistik')) then 
       if (is_true('statistik').or.is_true('statistik_taylor'))  then
          call mittelwerte(u)
       endif
       if (is_true('statistik_kor')) call mittelwerte_kor(u)
       if (is_true('statistik_strich')) call fluctuations(u)
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!! Evaluating the equation of mass   !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)

    call diff1x_zentral(rho*velxt ,dummy)
    kv(:,:,:,1,kStep) = -dummy/W_xi
    call diff1y_zentral(rho*velyt ,dummy)
    kv(:,:,:,1,kStep) = kv(:,:,:,1,kStep) - dummy/W_eta
    call diff1z_zentral(rho*velzt ,dummy)
    kv(:,:,:,1,kStep) = kv(:,:,:,1,kStep) - dummy/W_zeta

    kv(:,:,:,1,kStep) = kv(:,:,:,1,kStep) + fluxes(:,:,:,1)/(W_xi*W_eta*W_zeta)

    kv(:,:,:,1,kStep) = kv(:,:,:,1,kStep)*halb/Jacobian/u_sub(:,:,:,1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Set boundary fluxes for after eq. of mass !!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(scheme_type.eq.'implicit') then
       call set_bc_flux(u,u_sub,dt,kv(:,:,:,:,kStep),kStep,fluxes,1)
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! Updating u_sub !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(scheme_type.eq.'implicit') then
       call get_u_sub(u,kv,dt,scheme,kStep,u_sub)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! FRICTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.not.is_true('Eulerrechnung')) then

       ! Compute mu

       call get_T(u_sub,T)
       call visk(T,mu)


       mu_d = nll

       call get_param('Pr',Pr,default=0.71_rk)


       lambda  = Cp*mu/Pr

       ! Compute tau_ij



!!!!!! tau11

       tau11 = mu*zwei*(m11*velx_x +m21*velx_y + m31*velx_z)

       call diff1x_zentral( velxt , dummy )
       divu = dummy/W_xi
       call diff1y_zentral( velyt , dummy )
       divu = divu + dummy/W_eta
       call diff1z_zentral( velzt , dummy )
       divu = divu + dummy/W_zeta

       tau11 = tau11 + (mu_d-zweidrittel*mu)*divu

       tau11 = tau11/Jacobian

!!!!! tau22

       tau22 = mu*zwei*( m12*vely_x + m22*vely_y + m32*vely_z )

       tau22 = tau22 + (mu_d-zweidrittel*mu)*divu

       tau22 = tau22/Jacobian

!!!!! tau33

       tau33 = mu*zwei*( m13*velz_x + m23*velz_y + m33*velz_z )

       tau33 = tau33 + (mu_d-zweidrittel*mu)*divu

       tau33 = tau33/Jacobian

!!!!!! tau12

       tau12 = mu*( m12*velx_x + m11*vely_x  +  m22*velx_y + m21*vely_y  +  m32*velx_z + m31*vely_z )

       tau12 = tau12/Jacobian

!!!!!! tau13

       tau13 = mu*( m13*velx_x + m11*velz_x  +  m23*velx_y + m21*velz_y  +  m33*velx_z + m31*velz_z )

       tau13 = tau13/Jacobian


!!!!!! tau 23

       tau23 = mu*( m13*vely_x + m12*velz_x  +  m23*vely_y + m22*velz_y  + m33*vely_z + m32*velz_z )

       tau23 = tau23/Jacobian



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute projected Tau_i* components

       tau11t = m11*tau11 + m12*tau12 + m13*tau13  
       tau12t = m21*tau11 + m22*tau12 + m23*tau13  
       tau13t = m31*tau11 + m32*tau12 + m33*tau13  

       tau21t = m11*tau12 + m12*tau22 + m13*tau23  
       tau22t = m21*tau12 + m22*tau22 + m23*tau23  
       tau23t = m31*tau12 + m32*tau22 + m33*tau23  

       tau31t = m11*tau13 + m12*tau23 + m13*tau33  
       tau32t = m21*tau13 + m22*tau23 + m23*tau33  
       tau33t = m31*tau13 + m32*tau23 + m33*tau33  




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!! Friction terms for Momentum equation = div(tau_i*)/(J*srho)

!!!!!! Friction for momentum in u
       call diff1x_zentral( tau11t , dummy )
       fric_u = dummy/W_xi
       call diff1y_zentral( tau12t , dummy )
       fric_u = fric_u + dummy/W_eta
       call diff1z_zentral( tau13t , dummy )
       fric_u = fric_u + dummy/W_zeta


       fric_u = fric_u/(Jacobian*u_sub(:,:,:,1))


!!!!!! Friction for momentum in v

       call diff1x_zentral( tau21t , dummy )
       fric_v = dummy/W_xi
       call diff1y_zentral( tau22t , dummy )
       fric_v = fric_v + dummy/W_eta
       call diff1z_zentral( tau23t , dummy )
       fric_v = fric_v + dummy/W_zeta


       fric_v = fric_v/(Jacobian*u_sub(:,:,:,1))





!!!!!! Friction for momentum in w

       call diff1x_zentral( tau31t , dummy )
       fric_w = dummy/W_xi
       call diff1y_zentral( tau32t , dummy )
       fric_w = fric_w + dummy/W_eta
       call diff1z_zentral( tau33t , dummy )
       fric_w = fric_w + dummy/W_zeta

       fric_w = fric_w/(Jacobian*u_sub(:,:,:,1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Friction terms for the energy equation

!!!!!! Heat Flux

       call grad_zentral(T,T_x,T_y,T_z)

       T_x = T_x/W_xi
       T_y = T_y/W_eta
       T_z = T_z/W_zeta


       phi1 = lambda*(m11*T_x + m21*T_y + m31*T_z)
       phi1 = phi1/Jacobian

       phi2 = lambda*(m12*T_x + m22*T_y + m32*T_z)
       phi2 = phi2/Jacobian

       phi3 = lambda*(m13*T_x + m23*T_y + m33*T_z)
       phi3 = phi3/Jacobian


!!!!!! Friction Terms for energy equations

       ! All simple divergence terms for u_i*tau_ik and phi_k
       call diff1x_zentral( ( velx*tau11t + vely*tau21t + velz*tau31t + (m11*phi1 + m12*phi2 + m13*phi3) )  ,dummy)
       fric_p = dummy/W_xi
       call diff1y_zentral( ( velx*tau12t + vely*tau22t + velz*tau32t + (m21*phi1 + m22*phi2 + m23*phi3) )  ,dummy)
       fric_p = fric_p + dummy/W_eta
       call diff1z_zentral( ( velx*tau13t + vely*tau23t + velz*tau33t + (m31*phi1 + m32*phi2 + m33*phi3) )  ,dummy)
       fric_p = fric_p + dummy/W_zeta


       ! u_i*dx_k (tau_ik) terms
       call diff1x_zentral( tau11t  , dummy)         !( m11*tau11 + m12*tau12 + m13*tau13  , dummy)
       fric_p = fric_p - velx*dummy/W_xi
       call diff1y_zentral( tau12t  , dummy)
       fric_p = fric_p - velx*dummy/W_eta
       call diff1z_zentral( tau13t  , dummy)
       fric_p = fric_p - velx*dummy/W_zeta

       call diff1x_zentral( tau21t  , dummy)
       fric_p = fric_p - vely*dummy/W_xi
       call diff1y_zentral( tau22t  , dummy)
       fric_p = fric_p - vely*dummy/W_eta
       call diff1z_zentral( tau23t  , dummy)
       fric_p = fric_p - vely*dummy/W_zeta

       call diff1x_zentral( tau31t  , dummy)
       fric_p = fric_p - velz*dummy/W_xi
       call diff1y_zentral( tau32t  , dummy)
       fric_p = fric_p - velz*dummy/W_eta
       call diff1z_zentral( tau33t  , dummy)
       fric_p = fric_p - velz*dummy/W_zeta

       fric_p = (gamma-1)*fric_p/Jacobian



    else
       fric_p = nll
       fric_u = nll
       fric_v = nll
       fric_w = nll
      ! if((MB.eq.1).and.(kStep.eq.1)) then
      !    write(*,*)'-------------------'
      !    write(*,*)' Euler Calculation '
      !    write(*,*)'-------------------'
      ! end if
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! Equation of energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)

    call diff1x_zentral( velxt*u_sub(:,:,:,5) ,dummy)
    kv(:,:,:,5,kStep) = -dummy/W_xi
    call diff1y_zentral( velyt*u_sub(:,:,:,5) ,dummy)
    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep) - dummy/W_eta
    call diff1z_zentral( velzt*u_sub(:,:,:,5),dummy)
    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep) - dummy/W_zeta
    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep)*gamma

    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep) + gamma1 * ( velx*p_x + vely*p_y + velz*p_z )

    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep) + gamma1*fluxes(:,:,:,5)/(W_xi*W_eta*W_zeta)

    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep)/Jacobian

    kv(:,:,:,5,kStep) = kv(:,:,:,5,kStep) + fric_p




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Set boundary fluxes for after eq. of energy !!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(scheme_type.eq.'implicit') then
       call set_bc_flux(u,u_sub,dt,kv(:,:,:,:,kStep),kStep,fluxes,2)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! Updating u_sub!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(scheme_type.eq.'implicit') then
       call get_u_sub(u,kv,dt,scheme,kStep,u_sub)
    end if

!!!!!!!!! Recalculating p_x !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(scheme_type.eq.'implicit') then

      if (trim(Gittertypus).eq.'cylindrical') then

         call diff1x_zentral(u_sub(:,:,:,5),dummy)
         p_x = m11*dummy/W_xi
         call diff1y_zentral(u_sub(:,:,:,5),dummy)
         p_x = p_x + m21*dummy/W_eta
         call diff1z_zentral(u_sub(:,:,:,5),dummy)
         p_x = p_x + m31*dummy/W_zeta
   
         call diff1x_zentral(u_sub(:,:,:,5),dummy)
         p_y = m12*dummy/W_xi
         call diff1y_zentral(u_sub(:,:,:,5),dummy)
         p_y = p_y + m22*dummy/W_eta
         call diff1z_zentral(u_sub(:,:,:,5),dummy)
         p_y = p_y + m32*dummy/W_zeta
   
         call diff1x_zentral(u_sub(:,:,:,5),dummy)
         p_z = m13*dummy/W_xi
         call diff1y_zentral(u_sub(:,:,:,5),dummy)
         p_z = p_z + m23*dummy/W_eta
         call diff1z_zentral(u_sub(:,:,:,5),dummy)
         p_z = p_z + m33*dummy/W_zeta

      else

         call diff1x_zentral(m11*u_sub(:,:,:,5),dummy)
         p_x = dummy/W_xi
         call diff1y_zentral(m21*u_sub(:,:,:,5),dummy)
         p_x = p_x + dummy/W_eta
         call diff1z_zentral(m31*u_sub(:,:,:,5),dummy)
         p_x = p_x + dummy/W_zeta

         call diff1x_zentral(m12*u_sub(:,:,:,5),dummy)
         p_y = dummy/W_xi
         call diff1y_zentral(m22*u_sub(:,:,:,5),dummy)
         p_y = p_y + dummy/W_eta
         call diff1z_zentral(m32*u_sub(:,:,:,5),dummy)
         p_y = p_y + dummy/W_zeta

         call diff1x_zentral(m13*u_sub(:,:,:,5),dummy)
         p_z = dummy/W_xi
         call diff1y_zentral(m23*u_sub(:,:,:,5),dummy)
         p_z = p_z + dummy/W_eta
         call diff1z_zentral(m33*u_sub(:,:,:,5),dummy)
         p_z = p_z + dummy/W_zeta

      endif

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Momentum equations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp

    if(scheme_type.eq.'implicit') then
       call get_rho(u_sub,rho)
    end if


    call diff1x_zentral(velxt*rho*velx ,dummy)
    kv(:,:,:,2,kStep) = - halb*dummy/W_xi
    call diff1y_zentral(velyt*rho*velx ,dummy)
    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) - halb*dummy/W_eta
    call diff1z_zentral(velzt*rho*velx ,dummy)
    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) - halb*dummy /W_zeta

    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) - halb*rho*velxt*velx_x
    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) - halb*rho*velyt*velx_y
    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) - halb*rho*velzt*velx_z

    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) - p_x


    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) + fluxes(:,:,:,2)/(W_xi*W_eta*W_zeta)

    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep)/Jacobian/u_sub(:,:,:,1)

    kv(:,:,:,2,kStep) = kv(:,:,:,2,kStep) + fric_u


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for v


    call diff1x_zentral(velxt*rho*vely ,dummy)
    kv(:,:,:,3,kStep) = - halb*dummy/W_xi
    call diff1y_zentral(velyt*rho*vely ,dummy)
    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) - halb*dummy/W_eta
    call diff1z_zentral(velzt*rho*vely ,dummy)
    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) - halb*dummy/W_zeta

    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) - halb*rho*velxt*vely_x
    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) - halb*rho*velyt*vely_y
    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) - halb*rho*velzt*vely_z

    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) - p_y

    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) + fluxes(:,:,:,3)/(W_xi*W_eta*W_zeta)

    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep)/Jacobian/u_sub(:,:,:,1)

    kv(:,:,:,3,kStep) = kv(:,:,:,3,kStep) + fric_v


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for w


    call diff1x_zentral(velxt*rho*velz ,dummy)
    kv(:,:,:,4,kStep) = - halb*dummy/W_xi
    call diff1y_zentral(velyt*rho*velz ,dummy)
    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) - halb*dummy/W_eta
    call diff1z_zentral(velzt*rho*velz ,dummy)
    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) - halb*dummy/W_zeta

    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) - halb*rho*velxt*velz_x
    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) - halb*rho*velyt*velz_y
    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) - halb*rho*velzt*velz_z

    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) - p_z

    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) + fluxes(:,:,:,4)/(W_xi*W_eta*W_zeta)

    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep)/Jacobian/u_sub(:,:,:,1)

    kv(:,:,:,4,kStep) = kv(:,:,:,4,kStep) + fric_w




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Set boundary fluxes for after eq. of mom  !!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if(scheme_type.eq.'implicit') then
       call set_bc_flux(u,u_sub,dt,kv(:,:,:,:,kStep),kStep,fluxes,3)
    elseif(scheme_type.eq.'explicit') then
       call set_bc_flux(u,u_sub,dt,kv(:,:,:,:,kStep),kStep,fluxes,1)
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!! SPONGE TEST

    call set_sponge(kv(:,:,:,:,kStep),u)

!!!!!!!!!!!!!!!!!!! SPONGE TEST




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!  Statistics and Stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$    call get_param('av_start',avs,default=1)
!!$    if( STEP.eq.1.and.(kStep==1).and.(nt.ge.avs).and.is_true('statistik')) then 
!!$       if (is_true('statistik').or.is_true('statistik_taylor'))  then
!!$          call mittelwerte(u)
!!$       endif
!!$       if (is_true('statistik_kor')) call mittelwerte_kor(u)
!!$       if (is_true('statistik_strich')) call fluctuations(u)
!!$    endif
!!$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine iteration_step_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


  subroutine iteration_step_cmplx(u,ki,dt,kv,u_sub,kStep)
    complex(kind=rk), dimension(:,:,:,:),intent(in)             :: u
    complex(kind=rk), dimension(:,:,:,:),intent(inout)          :: ki
    real(kind=rk),optional                                      :: dt
    complex(kind=rk), dimension(:,:,:,:),optional   :: kv
    complex(kind=rk), dimension(:,:,:,:),optional   :: u_sub
    integer,optional                                            :: kStep


  end subroutine iteration_step_cmplx



  subroutine get_u_sub(u,kv,dt,scheme,kStep,ret)
    real(kind=rk), dimension(:,:,:,:),intent(in)             :: u
    real(kind=rk), dimension(:,:,:,:,:),intent(in)           :: kv
    real(kind=rk)                                            :: dt
    character(len=80)                                        :: scheme
    integer,optional                                         :: kStep
    integer                                                  :: stages,j
    real(kind=rk),dimension(:,:,:,:),intent(out)             :: ret

    stages = size(kv,5)

    ret = u

    do j=1,stages
       ret = ret + dt*a_rk(kStep,j)*kv(:,:,:,:,j)
    end do



  end subroutine get_u_sub





end module iterationstep_skew
