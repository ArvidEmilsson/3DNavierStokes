function params = init_geometry(params)

if params.geom.periodic(1),
    xi = (0:params.geom.n(1)-1)/ params.geom.n(1)      ;
else
    xi = (0:params.geom.n(1)-1)/(params.geom.n(1) -1)  ;    
end

if params.geom.periodic(2),
    eta = (0:params.geom.n(2)-1)/ params.geom.n(2)     ;
else
    eta = (0:params.geom.n(2)-1)/(params.geom.n(2)-1)  ;    
end

if params.geom.periodic(3),
    zeta = (0:params.geom.n(3)-1)/ params.geom.n(3)    ; 
else
    zeta = (0:params.geom.n(3)-1)/(params.geom.n(3)-1) ; 
end 
  
[XI,ETA,ZETA] = meshgrid(xi,eta,zeta) ; 

params.geom.XI = zeros([params.geom.n 3])  ; 

XI   = permute(  XI, [2,1,3]) ;  
ETA  = permute( ETA, [2,1,3]) ;  
ZETA = permute(ZETA, [2,1,3]) ; 

params.geom.XI(:,:,:,1) =   XI ;  
params.geom.XI(:,:,:,2) =  ETA ;  
params.geom.XI(:,:,:,3) = ZETA ; 

                    
params.geom.dxi(1)      =   xi(2)-  xi(1)   ;  

if params.geom.n(2)>1
    params.geom.dxi(2)      =  eta(2)- eta(1)   ;  
else
    params.geom.dxi(2)      =  1                ;   
end
                                                
if params.geom.n(3)>1
    params.geom.dxi(3)  = zeta(2)-zeta(1)   ;
else
    params.geom.dxi(3)  = 1                 ;
end

Lx = params.geom.L(1) ;  
Ly = params.geom.L(2) ;   
Lz = params.geom.L(3) ;   

switch params.geom.type
    case 'euclid'
        X =   XI*Lx ;  
        Y =  ETA*Ly ;  
        Z = ZETA*Lz ;  
    case 'test'
        alpha = params.geom.alpha ; 
        X =   XI*Lx  + alpha*sin( 2*pi*(XI+ETA+ZETA)) ;  
        Y =  ETA*Ly  + alpha*sin( 2*pi*(XI+ETA+ZETA)) ;    
        Z = ZETA*Lz  + alpha*sin( 2*pi*(XI+ETA+ZETA)) ;  

    case 'pipe'
        R   =  (ETA*2 -1)*Ly               ;
        PHI =  XI*Lx                       ;

        X = R.*cos(PHI) ;  
        Y = R.*sin(PHI) ; 
        Z = ZETA*Lz     ;
        
    case 'pipe_rXI_strech'
        tau = tanh(params.geom.sigma)   ; 
        
        R   =  (XI*2 -1)               ;
        R   = Lx*tanh(R*params.geom.sigma)/tau ; 

        PHI =  ETA*Ly                     ;

        X = R.*cos(PHI) ;  
        Y = R.*sin(PHI) ; 
        Z = ZETA*Lz     ;
        
        params.geom.r = R(:,1,1) ;   
    case 'pipe_rXI'
        
        R   =  (XI*2 -1)*Lx               ;
        PHI =  ETA*Ly                     ;

        X = R.*cos(PHI) ;  
        Y = R.*sin(PHI) ; 
        Z = ZETA*Lz     ;
        
        params.geom.r = R(:,1,1) ;   
        
    case 'capillary'
        r_nozzle  =@(z) min(z,0).^2/0.01^2*0.0035  + 0.0005 ;

        Z   = ZETA*Lz   + params.geom.z0   ;
                
        R   =  (XI*2 -1)*Lx               ;
        PHI =  ETA*Ly                       ;
        
        z   = squeeze( Z(1,1,:)  ) ; 
        r   = r_nozzle(z) ; 
                            
        for k=1:params.geom.n(3)
            R(:,:,k) = R(:,:,k)*r(k)  ; 
        end
        
        X   = R.*cos(PHI)   ;
        Y   = R.*sin(PHI)   ; 
        
        params.geom.r = R(:,1,1) ; 
    case 'capillary_strech'
        tau = tanh(1*params.geom.sigma)   ; 

        r_nozzle  =@(z) min(z,0).^2/0.01^2*0.0035  + 0.0005 ;

        Z   = ZETA*Lz   + params.geom.z0   ;
                
        R   =  (XI*2 -1)               ;
        R   = Lx*tanh(R*params.geom.sigma)/tau ; 

        PHI =  ETA*Ly                       ;
        
        z   = squeeze( Z(1,1,:)  ) ; 
        r   = r_nozzle(z) ; 
                            
        for k=1:params.geom.n(3)
            R(:,:,k) = R(:,:,k)*r(k)  ; 
        end
        
        X   = R.*cos(PHI)   ;
        Y   = R.*sin(PHI)   ; 
        
        params.geom.r = R(:,1,1) ; 
        
    case 'capillaryRoundInlet'
        r_nozzle  =@(z) min(z,0).^2/0.01^2*0.0035  + 0.0005 ;

        Z   = ZETA*Lz   + params.geom.z0   ;
                
        R   =  (XI*2 -1)*Lx               ;
        PHI =  ETA*Ly                       ;
        
        z   = squeeze( Z(1,1,:)  ) ; 
        r   = r_nozzle(z) ; 
        %keyboard                    
        for k=1:params.geom.n(3)
            R(:,:,k) = R(:,:,k)*r(k)  ;
            Z(:,:,k) = Z(:,:,k) + 30000*(r(k)-0.0005)*(R(:,:,k).^2-r(k)^2)   ;  

        end
        
        X   = R.*cos(PHI)   ;
        Y   = R.*sin(PHI)   ; 
        
        params.geom.r = R(:,1,1) ; 
        
    otherwise 
        error('no such geom type')
end

params.geom.X = zeros([params.geom.n, 3]) ; 

params.geom.X(:,:,:,1) =  X ; 
params.geom.X(:,:,:,2) =  Y ; 
params.geom.X(:,:,:,3) =  Z ; 

end