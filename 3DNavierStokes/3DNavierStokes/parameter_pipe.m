function params = parameter_pipe(params) 

    params.equation.dissipation = false ; 
    
    % time     %%%%%%%%%%%
    params.time.steps           = 1 ; 
    params.time.CFL             = 0.3 ; 
     % geometry %%%%%%%%%%%
    
    params.geom.n               = [ 168, 84,   1 ]    ; 
    params.geom.L               = [2*pi,  1,   5 ]    ; %2*pi since special treatment missing 

    params.geom.periodic        = [true,false,true]  ;
    params.geom.forcePeriodic   = [false,false,true] ;
    params.geom.type            = 'pipe'            ; 

    
    params.boundaries = {} ; 

    params.plot.freq  = 100 ; 

    params.boundaries{1}.name    = 'slipEuler' ;
    params.boundaries{1}.surface = [2,-1] ;
 
    params.boundaries{2}.name    = 'slipEuler' ;
    params.boundaries{2}.surface =  [2, 1] ;
     
    % start %%%%%%%%%%%
  
    params.start.type   = 'pulse' ;

    params.start.x0  = (params.geom.n(1)> 1)*0.25 ;
    params.start.y0  = (params.geom.n(2)> 1)*params.geom.L(2)/5 ;
    params.start.z0  = (params.geom.n(3)> 1)*params.geom.L(3)/2 ;

    params.start.sig    = [ 0.2 0.2 10]  ;

     
     params.io.freq = 100 ; 

end 