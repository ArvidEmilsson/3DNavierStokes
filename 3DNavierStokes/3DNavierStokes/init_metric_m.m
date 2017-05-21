function params = init_metric(params) 
    global m11 ;
    global m12 ;
    global m13 ;
    global m21 ;
    global m22 ;
    global m23 ;
    global m31 ;
    global m32 ;
    global m33 ;
    disp('init metric..')
    
    
    XI  = params.geom.XI ; 
    % remove periodicity if artifical periodic to avoid jump at end 
    Lx = params.geom.L(1)*params.geom.forcePeriodic(1) ;  
    Ly = params.geom.L(2)*params.geom.forcePeriodic(2) ;   
    Lz = params.geom.L(3)*params.geom.forcePeriodic(3) ;   

    W_xi   = params.deriv.W_xi     ; 
    W_eta  = params.deriv.W_eta    ; 
    W_zeta = params.deriv.W_zeta   ; 

    X  = params.geom.X(:,:,:,1) - Lx*XI(:,:,:,1);     
    Y  = params.geom.X(:,:,:,2) - Ly*XI(:,:,:,2);     
    Z  = params.geom.X(:,:,:,3) - Lz*XI(:,:,:,3);     

    
    x1 = D_xi(X  )./W_xi ;  
    y1 = D_xi(Y                  )./W_xi ;  
    z1 = D_xi(Z                  )./W_xi ;  

    x2 = D_eta(X                 )./W_eta ;  
    y2 = D_eta(Y )./W_eta ;  
    z2 = D_eta(Z                 )./W_eta ;  

    x3 = D_zeta(X                )./W_zeta ;  
    y3 = D_zeta(Y                )./W_zeta ;  
    z3 = D_zeta(Z )./W_zeta ;  

%    params.geom.m = zeros( [params.geom.n,3,3] ) ;
    
    %warning('geom factors calculated inconsitent') 
    %warning('artifical periodic missing') 
    
    %params.geom.m11      = y2.*z3  -  z2.*y3 ;
    %params.geom.m11      = y2.*z3 +Ly*z3+Lz*y2 +Ly*Lz -  z2.*y3 ;
    m11      =  D_eta(Y.*z3)./W_eta  +Ly*z3+Lz*y2 +Ly*Lz -  D_zeta(Y.*z2)./W_zeta ;
    %params.geom.m12      = z2.*x3  -  x2.*z3  ;
    %params.geom.m12      = z2.*x3  -  x2.*z3  - x2*Lz;
    m12      =  D_eta(Z.*x3)./W_eta   -  D_zeta(Z.*x2)./W_zeta - x2*Lz ;
    %params.geom.m13      = x2.*y3  -  y2.*x3   ;
    %params.geom.m13      = x2.*y3  -  y2.*x3  - Ly x3;
    m13      =  D_zeta(x2.*Y)./W_zeta -  D_eta(x3.*Y)./W_eta   - Ly*x3 ;

    %params.geom.m21      = y3.*z1  -  z3.*y1 ;
    %params.geom.m21      = y3.*z1  -  z3.*y1 -Lz y1;
    m21      = D_zeta(Y.*z1)./W_zeta  -  D_xi(z3.*Y)./W_xi     -Lz*y1  ;
    %params.geom.m22      = z3.*x1  -  x3.*z1 ;
    %params.geom.m22      = z3.*x1  + z3 *Lx  + Lz x1 + Lx Lz -  x3.*z1  ;
    m22      = D_zeta(Z.*x1)./W_zeta  + z3*Lx  + Lz*x1 + Lx*Lz -  D_xi(x3.*Z)./W_xi     ;
    %params.geom.m23      = x3.*y1  -  y3.*x1 ;
    %params.geom.m23      = x3.*y1  -  y3.*x1 -y3 Lx ;    
    m23      = D_xi(x3.*Y)./W_xi      -  D_zeta(Y.*x1)./W_zeta -y3*Lx ;
 
    %params.geom.m31      = y1.*z2  -  y2.*z1 ;
    m31      = D_xi(Y.*z2)./W_xi      -  D_eta(Y.*z1)./W_eta  -Ly*z1  ;
    %params.geom.m32      = z1.*x2  -  z2.*x1 ;
    m32      = D_xi(Z.*x2)./W_xi      -  D_eta(Z.*x1)./W_eta  -z2*Lx  ;
    %params.geom.m33      = x1.*y2  -  x2.*y1 ;
    m33      = D_eta(x1.*Y)./W_eta  + Lx*y2 + x1*Ly + Lx*Ly -  D_xi(x2.*Y)./W_xi      ;

    %keyboard
    
    params.geom.Jacobian = m31.*x3 + ...
                          +m32.*y3 + ... 
                          +m33.*(z3 +Lz) ; 

    disp('  ..done.')
    
end
