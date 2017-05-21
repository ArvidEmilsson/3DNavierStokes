function params = init_filter(params) 

    params.filter.A = cell(3,1) ;   
    params.filter.B = cell(3,1) ;   
    for l=1:3
        if params.geom.n(l) >1 
            if params.geom.periodic(l), 
                type = 'per' ; 
            else
                type = 'non-per' ;                 
            end
            order = 4 ; 
            alpha = 0.35;
            
            N = params.geom.n(l); 
            [params.filter.A{l},params.filter.B{l}] = pade_filter_matrices(type , order , alpha ,N) ; 
            
        end
        
        
    end



end 