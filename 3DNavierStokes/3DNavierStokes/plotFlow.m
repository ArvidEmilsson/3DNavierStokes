function plotFlow(params, q)
    
    pShape = params.plot.shape ; 

    figure(1)
    
    for l = 1:length(params.plot.vars) 
            
        var = params.plot.vars{l}; 
        switch  var 
            case 'rho'
                plotVar = q(:,:,:,1).^2  ; 
            case 'u'
                plotVar = q(:,:,:,2)./q(:,:,:,1)  ; 
            case 'v'
                plotVar = q(:,:,:,3)./q(:,:,:,1) ; 
            case 'w'
                plotVar = q(:,:,:,4)./q(:,:,:,1) ; 
            case 'p'
                plotVar = q(:,:,:,5)  ; 
            case 'Ma'
                
                gamma   = params.equation.gamma ; 
                plotVar =  sqrt(  ( q(:,:,:,2).^2 + q(:,:,:,3).^2  + q(:,:,:,4).^2)./(gamma*q(:,:,:,5)) )   ; 
            case 'T'
                Rs      = params.equation.Rs ; 
                plotVar = q(:,:,:,5)./(q(:,:,:,1).^2*Rs)     ;

            otherwise
                error(['do not know how to plot ' var])
        end

        X1 =params.geom.X(:,:,:, params.plot.spaceIndex(l,1) ) ; 
        X2 =params.geom.X(:,:,:, params.plot.spaceIndex(l,2) ) ; 
        
        subset = params.plot.subset{l} ;   
        
        subplot(pShape(1),pShape(2),l) 
%        keyboard
%subset{:}
        pcolorFull(squeeze(X1( subset{:})),squeeze(X2(subset{:})),squeeze(plotVar(subset{:}))); 
        %surf(squeeze(X1( subset{:})),squeeze(X2(subset{:})),squeeze(plotVar(subset{:})))
        shading flat %interp 
        axis equal 
        axis tight 
        colorbar
        title([var ' n=' num2str(params.time.n) ] ) 
        
    end    

    drawnow 


end
