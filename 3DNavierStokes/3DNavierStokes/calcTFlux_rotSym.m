
        T_xi   = D_xi(  T)./W_xi ;  
        T_eta  = D_eta( T)./W_eta ;   %%%%% D_EAT 
        T_zeta = D_zeta(T)./W_zeta ;  

 
% 
        T_x  = lambda.*(m11.*T_xi + m21.*T_eta + m31.*T_zeta)./Jacobian ;
        T_y  = lambda.*(m12.*T_xi + m22.*T_eta + m32.*T_zeta)./Jacobian ;
        T_z  = lambda.*(m13.*T_xi + m23.*T_eta + m33.*T_zeta)./Jacobian ;
    
        phi = zeros([ params.geom.n 3]) ;  
        phi(:,:,:,1)    = (m11.*T_x + m12.*T_y + m13.*T_z)  ; 
        phi(:,:,:,2)    = (m21.*T_x + m22.*T_y + m23.*T_z)  ; 
        phi(:,:,:,3)    = (m31.*T_x + m32.*T_y + m33.*T_z)  ;
