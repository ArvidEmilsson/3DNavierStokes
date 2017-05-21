function rhs = rhs_skew_NCP(q,~) 

    global params
    global m11 ;
    global m12 ;
    global m13 ;
    global m21 ;
    global m22 ;
    global m23 ;
    global m31 ;
    global m32 ;
    global m33 ;

    %call get_param('Gittertypus', Gittertypus,default='kartesisch')

    %call get_param('nt',nt)


    %[ n1, n2, n3, meqn ] = size(q) ;
    
    rhs = zeros(size(q)) ; 
    
    gamma = params.equation.gamma ; 
    gamma1 = gamma-1 ; 
    %Rs    = params.matetrial.Rs    ;  was gasc 
    
    %cv = params.matetrial.cv  ; 
    %cp = params.matetrial.cp  ; 
    %Pr = params.matetrial.Pr  ; 

    W_xi   = params.deriv.W_xi     ; 
    W_eta  = params.deriv.W_eta    ; 
    W_zeta = params.deriv.W_zeta   ; 
    
    %Jacobian = 1 ; 

    rho = q(:,:,:,1).^2 ;    % !call get_rho(u_sub,rho)
    
    %vel  = zeros(n1,n2,n3,3) ; 
    
    %vel(:,:,:,1) = q(:,:,:,2)./q(:,:,:,1);  %!call get_uvw(u_sub,velocities) %  %velocities(:,:,:,1)
    velx = q(:,:,:,2)./q(:,:,:,1);  % vel(:,:,:,1) ; 
    vely = q(:,:,:,3)./q(:,:,:,1);
    velz = q(:,:,:,4)./q(:,:,:,1);

    %vel_x = D_xi_multi(vel) ;  
    %velx_x = vel_x(:,:,:,1)./W_xi ;  
    
    velx_x = D_xi(velx)./W_xi ;  
    vely_x = D_xi(vely)./W_xi ;  
    velz_x = D_xi(velz)./W_xi ;  

    velx_y = D_eta(velx)./W_eta ;  
    vely_y = D_eta(vely)./W_eta ;  
    velz_y = D_eta(velz)./W_eta ;  

    velx_z = D_zeta(velx)./W_zeta ;  
    vely_z = D_zeta(vely)./W_zeta ;  
    velz_z = D_zeta(velz)./W_zeta ;  

    %m = params.geom.m ;
    
%!!!!!!!!!!!!!!!!!!!!!!! Computing projected velocities components

    velxt = m11.*velx + m12.*vely + m13.*velz ;
    velyt = m21.*velx + m22.*vely + m23.*velz ;
    velzt = m31.*velx + m32.*vely + m33.*velz ;

%!!!!!!! Computing derivatives of p with respect to physical coordinates
    p = q(:,:,:,5) ; 
    p_xi   = D_xi(  p )/W_xi ;  
    p_eta  = D_eta( p )/W_eta ;
    p_zeta = D_zeta(p )/W_zeta ;    
    %if (trim(Gittertypus).eq.'cylindrical') then
    %call diff1x_zentral(u_sub(:,:,:,5),dummy)
    p_x = m11.*p_xi + m21.*p_eta + m31.*p_zeta ; 
    p_y = m12.*p_xi + m22.*p_eta + m32.*p_zeta ; 
    p_z = m13.*p_xi + m23.*p_eta + m33.*p_zeta ; 

%  Evaluating the equation of mass   !!!!!!!!!!!!!!!!!!!!

%  RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    JsqRho =  params.geom.Jacobian.*q(:,:,:,1) ; 
    velxt_rho  =  velxt.*rho   ; 
    velyt_rho  =  velyt.*rho   ; 
    velzt_rho  =  velzt.*rho   ; 
    
    rhs(:,:,:,1) = -( D_xi(velxt_rho)./W_xi + D_eta(velyt_rho)./W_eta + D_zeta(velzt_rho)./W_zeta)./(2*JsqRho) ; % Jacobian*q(:,:,:,1) )   ; 
    
    
%     if (.not.is_true('Eulerrechnung')) then
% 
%        ! Compute mu
% 
%        call get_T(u_sub,T)
%        call visk(T,mu)
% 
% 
%        mu_d = nll
% 
%        call get_param('Pr',Pr,default=0.71_rk)
% 
% 
%        lambda  = Cp*mu/Pr
% 
%        ! Compute tau_ij
% 
% 
% 
% !!!!!! tau11
% 
%        tau11 = mu*zwei*(m11*velx_x +m21*velx_y + m31*velx_z)
% 
%        call diff1x_zentral( velxt , dummy )
%        divu = dummy/W_xi
%        call diff1y_zentral( velyt , dummy )
%        divu = divu + dummy/W_eta
%        call diff1z_zentral( velzt , dummy )
%        divu = divu + dummy/W_zeta
% 
%        tau11 = tau11 + (mu_d-zweidrittel*mu)*divu
% 
%        tau11 = tau11/Jacobian
% 
% !!!!! tau22
% 
%        tau22 = mu*zwei*( m12*vely_x + m22*vely_y + m32*vely_z )
% 
%        tau22 = tau22 + (mu_d-zweidrittel*mu)*divu
% 
%        tau22 = tau22/Jacobian
% 
% !!!!! tau33
% 
%        tau33 = mu*zwei*( m13*velz_x + m23*velz_y + m33*velz_z )
% 
%        tau33 = tau33 + (mu_d-zweidrittel*mu)*divu
% 
%        tau33 = tau33/Jacobian
% 
% !!!!!! tau12
% 
%        tau12 = mu*( m12*velx_x + m11*vely_x  +  m22*velx_y + m21*vely_y  +  m32*velx_z + m31*vely_z )
% 
%        tau12 = tau12/Jacobian
% 
% !!!!!! tau13
% 
%        tau13 = mu*( m13*velx_x + m11*velz_x  +  m23*velx_y + m21*velz_y  +  m33*velx_z + m31*velz_z )
% 
%        tau13 = tau13/Jacobian
% 
% 
% !!!!!! tau 23
% 
%        tau23 = mu*( m13*vely_x + m12*velz_x  +  m23*vely_y + m22*velz_y  + m33*vely_z + m32*velz_z )
% 
%        tau23 = tau23/Jacobian
% 
% 
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute projected Tau_i* components
% 
%        tau11t = m11*tau11 + m12*tau12 + m13*tau13  
%        tau12t = m21*tau11 + m22*tau12 + m23*tau13  
%        tau13t = m31*tau11 + m32*tau12 + m33*tau13  
% 
%        tau21t = m11*tau12 + m12*tau22 + m13*tau23  
%        tau22t = m21*tau12 + m22*tau22 + m23*tau23  
%        tau23t = m31*tau12 + m32*tau22 + m33*tau23  
% 
%        tau31t = m11*tau13 + m12*tau23 + m13*tau33  
%        tau32t = m21*tau13 + m22*tau23 + m23*tau33  
%        tau33t = m31*tau13 + m32*tau23 + m33*tau33  
% 
% 
% 
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% 
% !!!!!! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
% 
% !!!!!! Friction for momentum in u
%        call diff1x_zentral( tau11t , dummy )
%        fric_u = dummy/W_xi
%        call diff1y_zentral( tau12t , dummy )
%        fric_u = fric_u + dummy/W_eta
%        call diff1z_zentral( tau13t , dummy )
%        fric_u = fric_u + dummy/W_zeta
% 
% 
%        fric_u = fric_u/(Jacobian*u_sub(:,:,:,1))
% 
% 
% !!!!!! Friction for momentum in v
% 
%        call diff1x_zentral( tau21t , dummy )
%        fric_v = dummy/W_xi
%        call diff1y_zentral( tau22t , dummy )
%        fric_v = fric_v + dummy/W_eta
%        call diff1z_zentral( tau23t , dummy )
%        fric_v = fric_v + dummy/W_zeta
% 
% 
%        fric_v = fric_v/(Jacobian*u_sub(:,:,:,1))
% 
% 
% 
% 
% 
% !!!!!! Friction for momentum in w
% 
%        call diff1x_zentral( tau31t , dummy )
%        fric_w = dummy/W_xi
%        call diff1y_zentral( tau32t , dummy )
%        fric_w = fric_w + dummy/W_eta
%        call diff1z_zentral( tau33t , dummy )
%        fric_w = fric_w + dummy/W_zeta
% 
%        fric_w = fric_w/(Jacobian*u_sub(:,:,:,1))
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!! Friction terms for the energy equation
% 
% !!!!!! Heat Flux
% 
%        call grad_zentral(T,T_x,T_y,T_z)
% 
%        T_x = T_x/W_xi
%        T_y = T_y/W_eta
%        T_z = T_z/W_zeta
% 
% 
%        phi1 = lambda*(m11*T_x + m21*T_y + m31*T_z)
%        phi1 = phi1/Jacobian
% 
%        phi2 = lambda*(m12*T_x + m22*T_y + m32*T_z)
%        phi2 = phi2/Jacobian
% 
%        phi3 = lambda*(m13*T_x + m23*T_y + m33*T_z)
%        phi3 = phi3/Jacobian
% 
% 
% !!!!!! Friction Terms for energy equations
% 
%        ! All simple divergence terms for u_i*tau_ik and phi_k
%        call diff1x_zentral( ( velx*tau11t + vely*tau21t + velz*tau31t + (m11*phi1 + m12*phi2 + m13*phi3) )  ,dummy)
%        fric_p = dummy/W_xi
%        call diff1y_zentral( ( velx*tau12t + vely*tau22t + velz*tau32t + (m21*phi1 + m22*phi2 + m23*phi3) )  ,dummy)
%        fric_p = fric_p + dummy/W_eta
%        call diff1z_zentral( ( velx*tau13t + vely*tau23t + velz*tau33t + (m31*phi1 + m32*phi2 + m33*phi3) )  ,dummy)
%        fric_p = fric_p + dummy/W_zeta
% 
% 
%        ! u_i*dx_k (tau_ik) terms
%        call diff1x_zentral( tau11t  , dummy)         !( m11*tau11 + m12*tau12 + m13*tau13  , dummy)
%        fric_p = fric_p - velx*dummy/W_xi
%        call diff1y_zentral( tau12t  , dummy)
%        fric_p = fric_p - velx*dummy/W_eta
%        call diff1z_zentral( tau13t  , dummy)
%        fric_p = fric_p - velx*dummy/W_zeta
% 
%        call diff1x_zentral( tau21t  , dummy)
%        fric_p = fric_p - vely*dummy/W_xi
%        call diff1y_zentral( tau22t  , dummy)
%        fric_p = fric_p - vely*dummy/W_eta
%        call diff1z_zentral( tau23t  , dummy)
%        fric_p = fric_p - vely*dummy/W_zeta
% 
%        call diff1x_zentral( tau31t  , dummy)
%        fric_p = fric_p - velz*dummy/W_xi
%        call diff1y_zentral( tau32t  , dummy)
%        fric_p = fric_p - velz*dummy/W_eta
%        call diff1z_zentral( tau33t  , dummy)
%        fric_p = fric_p - velz*dummy/W_zeta
% 
%        fric_p = (gamma-1)*fric_p/Jacobian
% 
% 
% 
%     else
       fric_p =  0 ; 
       fric_u =  0 ; 
       fric_v =  0 ; 
       fric_w =  0 ; 

%    end if
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!! Equation of energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%!!! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    
    
    rhs(:,:,:,5) = -gamma*( D_xi(velxt.*p)./W_xi + D_eta(velyt.*p)./W_eta + D_zeta(velzt.*p)./W_zeta  ) ;   
    rhs(:,:,:,5) = rhs(:,:,:,5) +  gamma1 * ( velx.*p_x + vely.*p_y + velz.*p_z ) ; 
    
    rhs(:,:,:,5) = rhs(:,:,:,5)./params.geom.Jacobian ; 
    rhs(:,:,:,5) = rhs(:,:,:,5)  + fric_p  ;  


% !!!!!!!!!!!!!!!!!!  Momentum equations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
 

    rhs(:,:,:,2) = - ( D_xi(velxt_rho.*velx)./W_xi +  D_eta(velyt_rho.*velx)./W_eta + D_zeta(velzt_rho.*velx )./W_zeta ...
                     +  velxt_rho.*velx_x          + velyt_rho.*velx_y + velzt_rho.*velx_z )/2  ...
                     - p_x ; 
    rhs(:,:,:,2) = rhs(:,:,:,2)./(JsqRho) ; %Jacobian.*q(:,:,:,1)) ;                 
    rhs(:,:,:,2) = rhs(:,:,:,2) + fric_u ; 
    


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for v
    rhs(:,:,:,3) = - ( D_xi(velxt_rho.*vely)./W_xi +  D_eta(velyt_rho.*vely)./W_eta + D_zeta(velzt_rho.*vely )./W_zeta ...
                     +   velxt_rho.*vely_x          + velyt_rho.*vely_y + velzt_rho.*vely_z )/2  ...
                     - p_y ;
                 
    rhs(:,:,:,3) = rhs(:,:,:,3)./(JsqRho) ; %Jacobian.*q(:,:,:,1)) ;                 
    rhs(:,:,:,3) = rhs(:,:,:,3) + fric_v ; 
    

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for w
    rhs(:,:,:,4) = - ( D_xi(velxt_rho.*velz)./W_xi +  D_eta(velyt_rho.*velz)./W_eta + D_zeta(velzt_rho.*velz )./W_zeta ...
                     +  velxt_rho.*velz_x          + velyt_rho.*velz_y + velzt_rho.*velz_z )/2  ...
                     - p_z ;
                 
    rhs(:,:,:,4) = rhs(:,:,:,4)./(JsqRho) ; % Jacobian.*q(:,:,:,1)) ;                 
    rhs(:,:,:,4) = rhs(:,:,:,4) + fric_w ; 


%!!!!!!!!!!!!!!!! Set boundary fluxes for after eq. of mom  !!!!!!!!!!!!!!!!!! 
    rhs = setBoundary(rhs) ; 

end
