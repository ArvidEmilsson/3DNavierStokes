function params = geometry(params)

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

                    
params.geom.dxi     =   xi(2)-  xi(1) ;  
params.geom.deta    =  eta(2)- eta(1) ;  
if params.geom.n(3)>1
    params.geom.dzeta   = zeta(2)-zeta(1) ;
else
    params.geom.dzeta   = inf ;
end

Lx = params.geom.L(1) ;  
Ly = params.geom.L(2) ;   
Lz = params.geom.L(3) ;   

switch params.geom.type
    case 'euclid'
        X =   XI*Lx ;  
        Y =  ETA*Ly ;  
        Z = ZETA*Lz ;  
    case 'pipe'
        R   =  (ETA*2 -1)*Ly               ;
        PHI =  XI*Lx                       ;

        X = R.*cos(PHI) ;  
        Y = R.*sin(PHI) ; 
        Z = ZETA*Lz     ;
        
    case 'capillary'
        r_nozzle  =@(z) min(z,0).^2/0.01^2*0.0035  + 0.0005 ;

        Z   = ZETA*Lz   + params.geom.z0   ;
                
        R   =  (ETA*2 -1)*Ly               ;
        PHI =  XI*Lx                       ;
        
        z   = squeeze( Z(1,1,:)  ) ; 
        r   = r_nozzle(z) ; 
                            
        for k=1:params.geom.n(3)
            R(:,:,k) = R(:,:,k)*r(k)  ; 
        end
        
        X   = R.*cos(PHI)   ;
        Y   = R.*sin(PHI)   ; 

        

            
    otherwise 
        error('no such geom type')
end

params.geom.X =  X ; 
params.geom.Y =  Y ; 
params.geom.Z =  Z ;  

%params.geom.m = zeros(1,1,1,3,3) ; 
%params.geom.m(1,1,1,:,:) = [1,0,0;0,1,0;0,0,1] ; 
end