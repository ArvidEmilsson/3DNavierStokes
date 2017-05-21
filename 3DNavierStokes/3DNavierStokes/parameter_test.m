function params = parameter_test(params) 

    % equation    %%%%%%%%%%%
    
    params.equation.dissipation  = true  ;  % to switch of expensive calculation if not needed
    params.equation.Pr      = 0.71 ; 
    params.equation.mu0     = 1.716e-5 ; 

    params.equation.Rs      = 8.3144621/0.0270  ;  
    params.equation.polarFilter     = false ;  
    % time     %%%%%%%%%%%
    
    params.time.steps           = 100 ; 
    params.time.CFL             = 0.3 ; 
    
    % geometry %%%%%%%%%%%
    %params.geom.n               = [ 168,    84,   252 ]    ; 
    params.geom.n               = [ 64,    65,  1 ]    ; 
    params.geom.L               = [2*pi, 2*pi, 2*pi ]   ; %2*pi since special treatment missing 
    
    params.geom.periodic        = [false,true,true]  ;
    params.geom.forcePeriodic   = [false,true,true]  ;
    params.geom.type            = 'test'            ;   
    params.geom.alpha           =  0.0              ;  

    % boundaries %%%%%%%%%%%

    params.boundaries = {} ; 

    %params.boundaries{1}.name    = 'slipEuler' ;
    params.boundaries{1}.name    = 'nonSlipFlux' ;
    %params.boundaries{1}.name    = 'nonreflecting' ; 
    params.boundaries{1}.surface = [1,-1] ;

    params.boundaries{2}.name    = 'nonSlipFlux' ;
    %params.boundaries{2}.name    = 'nonreflecting' ;
    params.boundaries{2}.surface =  [1,1] ;

    % start %%%%%%%%%%%
    
    params.start.type          = 'pulse' ;
    params.start.x0  = (params.geom.n(1)> 1)*params.geom.L(1)/2 ;
    params.start.y0  = (params.geom.n(2)> 1)*params.geom.L(2)/2 ;
    params.start.z0  = (params.geom.n(3)> 1)*params.geom.L(3)/2 ;

    params.start.sig    = params.geom.L/10  ;
    params.start.sig(2) = 1e6 ;  

    
    % plot %%%%%%%%%%%

    params.plot.freq    = 25 ;    
    params.plot.shape   = [2,2] ; 
    params.plot.vars    = {'u','p','v','T'} ; 
    params.plot.subset  = {{':',':',1} ,{':',':',1} ,{':',':',1} ,{':',':',1}  } ;  
    params.plot.spaceIndex = [1,2; 1,2; 1,2; 1,2]; 
    
    % io %%%%%%%%%%%%%%%
    params.io.freq = 1e10 ; 

end 