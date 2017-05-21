function [params, q ] = start_vals_euler(params) 

    X = params.geom.X(:,:,:,1) ;
    Y = params.geom.X(:,:,:,2) ;
    Z = params.geom.X(:,:,:,3) ;
    
    q = zeros([params.geom.n params.equation.nVar]) ; 
    gamma = params.equation.gamma   ; 
    Ly  = params.geom.L(2); 
            
    switch params.start.type 
        case 'pulse'

            x0  = params.start.x0 ;
            y0  = params.start.y0 ;
            z0  = params.start.z0 ;
            
            
            sigX2  = params.start.sig(1)^2 ; %  (Ly/6)^2 ; 
            sigY2  = params.start.sig(2)^2 ; % (Ly/6)^2 ; 
            sigZ2  = params.start.sig(3)^2 ; % (Ly/12)^2 ;

            a = 100 ;
            p0      = 1e5 ;
            %rho0    = 1   ; 
            
            q(:,:,:, 5) = p0+ a*exp( -(X-x0).^2/sigX2 -(Y-y0).^2/sigY2 -(Z-z0).^2/sigZ2  ) ; 
            q(:,:,:, 1) = sqrt( (q(:,:,:, 5)/p0).^(1/gamma)  )   ;
        case 'Tpulse'
            x0  = params.start.x0 ;
            y0  = params.start.y0 ;
            z0  = params.start.z0 ;
            
            sigX2  = (Ly/6)^2 ; 
            sigY2  = (Ly/6)^2 ; 
            sigZ2  = (Ly/12)^2 ;

            a = .1 ;
            p0      = 1e5 ;
            %rho0    = 1   ; 
            
            q(:,:,:, 5) = p0 ; 
            q(:,:,:, 1) = sqrt( 1- a*exp( -(X-x0).^2/sigX2 -(Y-y0).^2/sigY2 -(Z-z0).^2/sigZ2  )   )   ;
        case 'linearZ'
            p0      = 1e5 ;
            gamma = params.equation.gamma ; 
            pS = (2/(gamma+1))^(gamma/(gamma-1))*p0 ; 
            
            q(:,:,:,1)  = 1 ;  %sqrt(1)
            q(:,:,:,5)  = p0 + params.geom.XI(:,:,:,3) *(pS -p0) ;
            
            r_z =   params.geom.r_xi(:,:,:,:,3) ;
            
            
            n = (sqrt( r_z(:,:,:,1).^2 + r_z(:,:,:,2).^2 + r_z(:,:,:,3).^2   )  ) ;
            r_z(:,:,:,1) = r_z(:,:,:,1)./n ; 
            r_z(:,:,:,2) = r_z(:,:,:,2)./n ; 
            r_z(:,:,:,3) = r_z(:,:,:,3)./n ; 
            
            ZETA = params.geom.XI(:,:,:,3)  ; 
            
            q(:,:,:,2) = r_z(:,:,:,1).*( ZETA*300 )  ; % SQRT RHO MISSING (1)
            q(:,:,:,3) = r_z(:,:,:,2).*( ZETA*300 )  ;
            q(:,:,:,4) = r_z(:,:,:,3).*( ZETA*300 )  ;

        case 'linearZWall'
            p0      = 1e5 ;
            gamma = params.equation.gamma ; 
            pS = (2/(gamma+1))^(gamma/(gamma-1))*p0 ; 
            
            q(:,:,:,1)  = 1 ;  %sqrt(1)
            q(:,:,:,5)  = p0 + params.geom.XI(:,:,:,3) *(pS -p0) ;
            
            r_z =   params.geom.r_xi(:,:,:,:,3) ;
            
            
            n = (sqrt( r_z(:,:,:,1).^2 + r_z(:,:,:,2).^2 + r_z(:,:,:,3).^2   )  ) ;
            r_z(:,:,:,1) = r_z(:,:,:,1)./n ; 
            r_z(:,:,:,2) = r_z(:,:,:,2)./n ; 
            r_z(:,:,:,3) = r_z(:,:,:,3)./n ; 
                    
            XI      = params.geom.XI(:,:,:,1)   ; 
            ZETA    = params.geom.XI(:,:,:,3)   ; 
            lam     = params.start.lambda       ;  
            bounLay = (tanh(XI/lam) +  tanh(-(XI-1 ) /lam) )-1; 
            q(:,:,:,2) = r_z(:,:,:,1).*( ZETA*300 ).*bounLay  ; % SQRT RHO MISSING since one here (1)
            q(:,:,:,3) = r_z(:,:,:,2).*( ZETA*300 ).*bounLay  ;
            q(:,:,:,4) = r_z(:,:,:,3).*( ZETA*300 ).*bounLay  ;
        case 'testRotSym' 
            q(:,1,:, 1) = 1 + 0.2*exp( -((X(:,1,:)-0.1e-3).^2+(Z(:,1,:)).^2)/(1e-3)^2 ) ; 
            q(:,1,:, 2) =  0.2*exp( -((X(:,1,:)+0.1e-3).^2+(Z(:,1,:)+0.004).^2)/(1e-3)^2 ) ; 
            q(:,1,:, 4) = 1 + 0.2*exp( -((X(:,1,:)+0.2e-3).^2+(Z(:,1,:)-0.006).^2)/(1e-3)^2 ) ; 
            q(:,1,:, 5) = 1 + 0.2*exp( -((X(:,1,:)+0.1e-3).^2+(Z(:,1,:)+0.005).^2)/(1e-3)^2 ) ; 
            q = setRotSym(params,q) ; 
        otherwise 
            error('no so start type')
    end
    
end