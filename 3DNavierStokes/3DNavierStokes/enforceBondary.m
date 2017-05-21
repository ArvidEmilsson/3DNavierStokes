function q = enforceBondary(params,q)

    for i =1:length(params.boundaries)
        switch params.boundaries{i}.name
            case 'nonSlipFlux'
                n0      = params.boundaries{i}.n0 ;   
                n1      = params.boundaries{i}.n1 ;   
                
                q(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),2) = 0 ; 
                q(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),3) = 0 ; 
                q(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),4) = 0 ;
            
        end
    end

end 