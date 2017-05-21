%function rhs = setBoundary(rhs)
% not a dunction to grep e.g the heat fluxes 

global params 

for i =1:length(params.boundaries)
    %i
    %disp(params.boundaries{i}.name) 
    n0      = params.boundaries{i}.n0 ;   
    n1      = params.boundaries{i}.n1 ;   
    surf    = params.boundaries{i}.surface ;   
    spIndex = surf(1);
    signFac = params.boundaries{i}.surface(2)  ;      

    switch params.boundaries{i}.name
        case 'slipEuler'
            e = params.boundaries{i}.e ; 


            a = -( e(:,:,:,1).*rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),2) ...
                  +e(:,:,:,2).*rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),3) ...
                  +e(:,:,:,3).*rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),4) )      ; 
            
            rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),2) = ...
                 rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),2) +a.*e(:,:,:,1) ;  
            rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),3) = ...
                 rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),3) +a.*e(:,:,:,2) ;  
            rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),4) = ...
                 rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),4) +a.*e(:,:,:,3) ; 

           %  x = params.geom.X(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3));
           %  y = params.geom.Y(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3));
           % hold on 
           % quiver(x,y,squeeze(e(:,:,:,1)),squeeze(e(:,:,:,2)))
                       

           %  keyboard
        case 'nonSlipFlux'
           
            rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),2) = 0 ; 
            rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),3) = 0 ; 
            rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),4) = 0 ;
            
            switch spIndex
                case 1
                   
                    WW = params.deriv.W_Eta(n0(2):n1(2))*params.deriv.W_Zeta(n0(3):n1(3))';
                    
                    signFac = signFac*params.geom.dxi(2)*params.geom.dxi(3); 
                case 2
                    WW = params.deriv.W_Xi*params.deriv.W_Zeta'  ;
                    signFac = signFac*params.geom.dxi(1)*params.geom.dxi(3); 
                case 3
                    WW = params.deriv.W_Xi*params.deriv.W_Eta'  ;
                    signFac = signFac*params.geom.dxi(1)*params.geom.dxi(2); 
            end
                WW=reshape(WW,n1-n0+1); 
                %keyboard
                   rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),5) = ...
                   rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),5) -...
                  signFac*WW.*phi(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3), spIndex) ; 
        case 'nonreflecting'
                nonreflecting
                
        otherwise
            error('no such bound') 
    end

end

