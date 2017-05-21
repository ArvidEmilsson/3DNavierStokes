function rhs = rhs_skew_NCP_rotSym(q,~) 

    global params 
    
    rhs = zeros(size(q)) ; 
    
    gamma = params.equation.gamma ; 
    gamma1 = gamma-1 ; 

    W_xi   = params.deriv.W_xi     ; 
    W_eta  = params.deriv.W_eta    ; 
    W_zeta = params.deriv.W_zeta   ; 
    
    m11 = params.geom.m(:,:,:,1,1) ; % y2.*z3  -  z2.*y3 ;
    m12 = params.geom.m(:,:,:,1,2) ; % z2.*x3  -  x2.*z3  ;
    m13 = params.geom.m(:,:,:,1,3) ; % x2.*y3  -  y2.*x3   ;
    
    m21 = params.geom.m(:,:,:,2,1) ; 
    m22 = params.geom.m(:,:,:,2,2) ; 
    m23 = params.geom.m(:,:,:,2,3) ; 
    
    m31 = params.geom.m(:,:,:,3,1) ; 
    m32 = params.geom.m(:,:,:,3,2) ; 
    m33 = params.geom.m(:,:,:,3,3) ; 
                
    Jacobian    = params.geom.Jacobian ;
    rho = q(:,:,:,1).^2 ;    % !call get_rho(u_sub,rho)
    
    
    velx = q(:,:,:,2)./q(:,:,:,1);  % u = cos(eta) U(xi,zeta) 
    vely = q(:,:,:,3)./q(:,:,:,1);  % v = sin(eta) U(xi,zeta) 
    velz = q(:,:,:,4)./q(:,:,:,1);  % w = w(xi,zeta) 

    
    velx_x = D_xi(velx)./W_xi ;  
    vely_x = D_xi(vely)./W_xi ;  
    velz_x = D_xi(velz)./W_xi ;  

    velx_y = D_eta(velx)./W_eta ;    % - s(e) U     (=0 @e=0)              %%%%% D_EAT 
    vely_y = D_eta(vely)./W_eta ;    % c(e) U                   %%%%% D_EAT 
    velz_y = D_eta(velz)./W_eta ;    % 0                      %%%%% D_EAT 

    velx_z = D_zeta(velx)./W_zeta ;  
    vely_z = D_zeta(vely)./W_zeta ;  
    velz_z = D_zeta(velz)./W_zeta ;  


%!!!!!!!!!!!!!!!!!!!!!!! Computing projected velocities components

    velxt = m11.*velx + m12.*vely + m13.*velz ;
    velyt = m21.*velx + m22.*vely + m23.*velz ;
    velzt = m31.*velx + m32.*vely + m33.*velz ;

%!!!!!!! Computing derivatives of p with respect to physical coordinates
    p = q(:,:,:,5) ; 
    p_xi   = D_xi(  p )/W_xi ;  
    p_eta  = D_eta( p )/W_eta ;                        %%%%% D_EAT 
    p_zeta = D_zeta(p )/W_zeta ;    

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
   
                                              %%%%% D_EAT 
                                              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if params.equation.dissipation
        %fric_p =  0 ; 
        %fric_u =  0 ; 
        %fric_v =  0 ; 
        %fric_w =  0 ; 

        Rs    = params.equation.Rs    ; % R specific, was gasc 
        Pr  = params.equation.Pr  ; 

        % Compute mu
        T   = p./(rho*Rs)     ;  %        call get_T(u_sub,T)
        T0  = 273.15        ;

        %mu  = ((T0 + 110.4)*params.equation.mu0)*((T/T0).^1.5)./(T + 110.4); %        call visk(T,mu)
        mu  = ((T0 + 110.4)*params.equation.mu0/(T0.^1.5))*((T).*sqrt(T))./(T + 110.4); %        call visk(T,mu)

        mu_d = 0            ;
        
        cv   = Rs/gamma1;
        cp   = cv*gamma; 

        lambda  = cp*mu/Pr ; 
        
        if params.equation.ZeroLambdaAtBound
            setZeroLambdaAtBound
        end 

%        ! Compute tau_ij
%       div U 

        divu = D_xi(  velxt)./W_xi  + ...
               D_eta( velyt)./W_eta + ...                           %%%%% D_EAT 
               D_zeta(velzt)./W_zeta ;  

% 
% !!!!!! tau11
% 
        tau11 = ( 2*mu.*(m11.*velx_x +m21.*velx_y + m31.*velx_z)  ... 
                 +(mu_d-(2/3)*mu).*divu)./Jacobian                      ;

% !!!!! tau22
% 
        tau22 =  ( 2*mu.*( m12.*vely_x + m22.*vely_y + m32.*vely_z ) ...
               + ( mu_d-(2/3)*mu).*divu)./Jacobian                      ;

% !!!!! tau33
% 
        tau33 =  (2*mu.*( m13.*velz_x + m23.*velz_y + m33.*velz_z ) ...
                  + (mu_d-(2/3)*mu).*divu)./Jacobian                    ;

% !!!!!! tau12
% 
        tau12 = mu.*( m12.*velx_x + m11.*vely_x + m22.*velx_y ... 
                    + m21.*vely_y + m32.*velx_z + m31.*vely_z )./Jacobian ;
 
% !!!!!! tau13

        tau13 = mu.*( m13.*velx_x + m11.*velz_x + m23.*velx_y + ...
                      m21.*velz_y + m33.*velx_z + m31.*velz_z )./Jacobian ;
 
% !!!!!! tau 23

        tau23 = mu.*( m13.*vely_x + m12.*velz_x + m23.*vely_y + ...
                      m22.*velz_y + m33.*vely_z + m32.*velz_z )./Jacobian ;

                  % 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute projected Tau_i* components
% 
        tau11t = m11.*tau11 + m12.*tau12 + m13.*tau13  ;
        tau12t = m21.*tau11 + m22.*tau12 + m23.*tau13  ;
        tau13t = m31.*tau11 + m32.*tau12 + m33.*tau13  ;

        tau21t = m11.*tau12 + m12.*tau22 + m13.*tau23  ;
        tau22t = m21.*tau12 + m22.*tau22 + m23.*tau23  ;
        tau23t = m31.*tau12 + m32.*tau22 + m33.*tau23  ;

        tau31t = m11.*tau13 + m12.*tau23 + m13.*tau33  ;
        tau32t = m21.*tau13 + m22.*tau23 + m23.*tau33  ;
        tau33t = m31.*tau13 + m32.*tau23 + m33.*tau33  ;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% !!!!!! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
% 
% !!!!!! Friction for momentum in u
        
        tau11t_xi   = D_xi(  tau11t)./W_xi ;   
        tau12t_eta  = D_eta( tau12t)./W_eta ;                          %%%%% D_EAT    
        tau13t_zeta = D_zeta(tau13t)./W_zeta ;   

        fric_u =( tau11t_xi + ...
                  tau12t_eta + ...
                  tau13t_zeta ) ;

% 
% !!!!!! Friction for momentum in v
% 
        fric_v =( D_xi(  tau21t)./W_xi + ...
                 D_eta( tau22t)./W_eta + ...                           %%%%% D_EAT 
                D_zeta(tau23t)./W_zeta ); % ./(Jacobian.*q(:,:,:,1) )  ;  


% 
% !!!!!! Friction for momentum in w
% 
        fric_w =( D_xi(  tau31t)./W_xi + ...
                 D_eta( tau32t)./W_eta + ...                                                  %%%%% D_EAT 
                D_zeta(tau33t)./W_zeta ); %./(Jacobian.*q(:,:,:,1) )  ;  

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!! Friction terms for the energy equation

% !!!!!! Heat Flux
% 
%        call grad_zentral(T,T_x,T_y,T_z)

        calcTFlux_rotSym                                                              %%%%% D_EAT 
                                                            
        % !!!!!! Friction Terms for energy equations
% 
%        ! All simple divergence terms for u_i*tau_ik and phi_k
         fric_p = D_xi(  velx.*tau11t + vely.*tau21t + velz.*tau31t + phi(:,:,:,1) )./W_xi ;    
         fric_p = fric_p + ...
                   D_eta( velx.*tau12t + vely.*tau22t + velz.*tau32t +  phi(:,:,:,2) )./W_eta ;                        %%%%% D_EAT 
                 
         fric_p = fric_p + ...
                   D_zeta( velx.*tau13t + vely.*tau23t + velz.*tau33t + phi(:,:,:,3) )./W_zeta ;   
         fric_p = fric_p - velx.*(fric_u ) ;     
        fric_p = fric_p - vely.*(fric_v) ;  
        fric_p = fric_p - velz.*( fric_w ) ; 
% 
        fric_p = (gamma-1)*fric_p./Jacobian ;
    else
       fric_p =  0 ; 
       fric_u =  0 ; 
       fric_v =  0 ; 
       fric_w =  0 ; 

    end 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!! Equation of energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%!!! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    
    
    rhs(:,:,:,5) = -gamma*( D_xi(velxt.*p)./W_xi + D_eta(velyt.*p)./W_eta + D_zeta(velzt.*p)./W_zeta  ) ;            %%%%% D_EAT    
    rhs(:,:,:,5) = rhs(:,:,:,5) +  gamma1 * ( velx.*p_x + vely.*p_y + velz.*p_z ) ; 
    
    rhs(:,:,:,5) = rhs(:,:,:,5)./params.geom.Jacobian ; 
    rhs(:,:,:,5) = rhs(:,:,:,5)  + fric_p  ;  


% !!!!!!!!!!!!!!!!!!  Momentum equations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
 

    rhs(:,:,:,2) = - ( D_xi(velxt_rho.*velx)./W_xi +  D_eta(velyt_rho.*velx)./W_eta + D_zeta(velzt_rho.*velx )./W_zeta ...       %%%%% D_EAT 
                     +  velxt_rho.*velx_x          + velyt_rho.*velx_y + velzt_rho.*velx_z )/2  ...
                     - p_x ; 
    rhs(:,:,:,2) = ( rhs(:,:,:,2) + fric_u)./(JsqRho) ; %Jacobian.*q(:,:,:,1)) ;                 
    


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for v
    rhs(:,:,:,3) = - ( D_xi(velxt_rho.*vely)./W_xi +  D_eta(velyt_rho.*vely)./W_eta + D_zeta(velzt_rho.*vely )./W_zeta ...               %%%%% D_EAT 
                     +   velxt_rho.*vely_x          + velyt_rho.*vely_y + velzt_rho.*vely_z )/2  ...
                     - p_y ;
                 
    rhs(:,:,:,3) = ( rhs(:,:,:,3) + fric_v )./(JsqRho) ; %Jacobian.*q(:,:,:,1)) ;                 
    

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RHS of  momentum equation for w
    rhs(:,:,:,4) = - ( D_xi(velxt_rho.*velz)./W_xi +  D_eta(velyt_rho.*velz)./W_eta + D_zeta(velzt_rho.*velz )./W_zeta ...                %%%%% D_EAT 
                     +  velxt_rho.*velz_x          + velyt_rho.*velz_y + velzt_rho.*velz_z )/2  ...
                     - p_z ;
                 
    rhs(:,:,:,4) = (rhs(:,:,:,4) + fric_w)./(JsqRho) ; % Jacobian.*q(:,:,:,1)) ;                 


%!!!!!!!!!!!!!!!! Set boundary fluxes for after eq. of mom  !!!!!!!!!!!!!!!!!! 
    
    if params.equation.polarFilter 
        kNumMin = 3 ; 
        rhs = sharpFilterFieldDynProtoNew(rhs,params.geom.r, kNumMin) ;                                     
    end 
    
    setBoundary

    
end
