function params = init_bound(params, q) 


for i =1:length(params.boundaries)
    %disp(params.boundaries{i}.name)
    spaceIndex = params.boundaries{i}.surface(1) ;     
    if params.geom.periodic(spaceIndex)
        error('boundary for periodic direction?')
    end
    
    n0 = [1,1,1,]       ; 
    n1 = params.geom.n  ;
    
    if params.boundaries{i}.surface(2) == 0,
        error('obsolet should be -1')
    end
    
    if  params.boundaries{i}.surface(2)==-1
        index = 1 ; 
    else
        index = params.geom.n(spaceIndex) ;
    end
    
    n0(spaceIndex) = index ;  
    n1(spaceIndex) = index ;  
    
    %%%%%%%%%%%%%%%%%% to remove commen point between bound, if needed     
    if isfield(params.boundaries{i},'removeEdge') 
        removeEdge = params.boundaries{i}.removeEdge ;  
        for l = 1:size(removeEdge,1)
            spaceIndex =  removeEdge(l,1);    
            sign       =  removeEdge(l,2);
            if sign==1
               n1(spaceIndex) = n1(spaceIndex)-1 ; 
            else
                if sign==-1
                    n0(spaceIndex) = n0(spaceIndex)+1 ; 
                else
                      error('sign is -1 or 1') 
                end
           end
               
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    params.boundaries{i}.n0 = n0 ; 
    params.boundaries{i}.n1 = n1 ; 

    
    %%% 
    switch params.boundaries{i}.name
        case {'slipEuler','nonSlipFlux'}
            sign = params.boundaries{i}.surface(2)  ;      
            
            e       = params.geom.m(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),spaceIndex,:)*sign ;
            normE   = sqrt(e(:,:,:,1,1).^2 + e(:,:,:,1,2).^2 + e(:,:,:,1,3).^2 ); 
           
            sE      = size(e) ;  
            params.boundaries{i}.e = reshape(e, [sE(1:3) 3]); 
            params.boundaries{i}.e(:,:,:,1)  = params.boundaries{i}.e(:,:,:,1)./reshape(normE, [sE(1:3)]);
            params.boundaries{i}.e(:,:,:,2)  = params.boundaries{i}.e(:,:,:,2)./reshape(normE, [sE(1:3)]);
            params.boundaries{i}.e(:,:,:,3)  = params.boundaries{i}.e(:,:,:,3)./reshape(normE, [sE(1:3)]);

        case {'nonreflecting'}
            sign = params.boundaries{i}.surface(2)  ;   
            
            e       = params.geom.r_xi(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),:,spaceIndex)*sign ; 
            normE   = sqrt(e(:,:,:,1).^2 + e(:,:,:,2).^2 + e(:,:,:,3).^2 ); 
            params.boundaries{i}.e = zeros( [n1-n0+1, 3]  ) ; 
            params.boundaries{i}.e(:,:,:,1)  = e(:,:,:,1)./normE ; %reshape(  , n1-n0+1) ; 
            params.boundaries{i}.e(:,:,:,2)  = e(:,:,:,2)./normE ; 
            params.boundaries{i}.e(:,:,:,3)  = e(:,:,:,3)./normE ; 
           
            qNow = q(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),:) ; 
            
            qRef = zeros( size(qNow)  ) ;
            
            qRef(:,:,:,1) = qNow(:,:,:,1).^2               ;
            qRef(:,:,:,2) = qNow(:,:,:,2)./qNow(:,:,:,1)   ;
            qRef(:,:,:,3) = qNow(:,:,:,3)./qNow(:,:,:,1)   ;
            qRef(:,:,:,4) = qNow(:,:,:,4)./qNow(:,:,:,1)   ;
            qRef(:,:,:,5) = qNow(:,:,:,5)                  ;
           
            
            params.boundaries{i}.qRef   = qRef  ;   
        otherwise
            error('no such bound') 

    end
    
    
end
   


end