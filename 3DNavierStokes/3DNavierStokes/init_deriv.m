function params = init_deriv(params) 

    n  =params.geom.n ; 
    
    D_periodic = @D_periodicFourthOrd           ; 
    D_open     = @D_openFourthOrdWeightsStrand  ; 
    
    if params.geom.periodic(1)
        params.deriv.Dxi    = D_periodic(n(1),params.geom.dxi(1)  )  ; 
        params.deriv.W_xi   = 1 ;
        
        params.deriv.W_Xi   = ones(n(1),1) ;
    else
        [w, D, fac ]        = D_open(n(1),params.geom.dxi(1) ) ; 
        params.deriv.Dxi    = sparse(diag(1./w))*D/fac         ; 
        params.deriv.W_xi   = 1 ; 
        
        params.deriv.W_Xi   = w ;             % i guess intended for flux calc, is it used? 
    end
    
    if params.geom.n(2)>1
        if params.geom.periodic(2)
            params.deriv.Deta   = D_periodic(n(2),params.geom.dxi(2) )  ; 
            params.deriv.W_eta   = 1 ; 

            params.deriv.W_Eta   = ones(n(2),1) ;
        else
            [w, D, fac ]        = D_open(n(2),params.geom.dxi(2) ) ; 
            params.deriv.Deta    = sparse(diag(1./w))*D/fac         ; 
            params.deriv.W_eta   = 1 ; 

            params.deriv.W_Eta   = w ;            

        end
    else
        params.deriv.Deta = 0 ; 
        params.deriv.W_eta = 1 ;        
        params.deriv.W_Eta = 1 ;                
    end
    
    if params.geom.n(3)>1
        if params.geom.periodic(3)
            params.deriv.Dzeta  = D_periodic(n(3),params.geom.dxi(3))  ; 
            params.deriv.W_zeta = 1 ;
            
            params.deriv.W_Zeta   = ones(n(3),1) ;
        else
            [w, D, fac ]        = D_open(n(3),params.geom.dxi(3) ) ; 
            params.deriv.Dzeta    = sparse(diag(1./w))*D/fac         ; 
            params.deriv.W_zeta   = 1 ;

            params.deriv.W_Zeta   = w ;            
        end
    else
        params.deriv.Dzeta = 0 ; 
        params.deriv.W_zeta = 1 ;        
        params.deriv.W_Zeta = 1 ;        
    end
    
end 