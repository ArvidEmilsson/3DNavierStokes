function params = parameter_capillary3D(params) 

     % equation    %%%%%%%%%%%
    
    params.equation.dissipation  = true  ;  % to switch of expensive calculation if not needed
    params.equation.Pr      = 0.71 ; 
    params.equation.mu0     = 1.716e-5 ; 

    params.equation.Rs    = 8.3144621/0.0270   ; 
    
    params.equation.polarFilter = true ; 
    
    params.equation.ZeroLambdaAtBound  = true ; 
    % time     %%%%%%%%%%%
    
    params.time.steps           = 10000 ; 
    params.time.CFL             = .1   ; 
    params.time.checkNaN        = 1     ; 
    
    % init_filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    params.filter.freq          = 50    ; 
    params.filter.eta           = true  ; 
    params.filter.zeta          = false ;  
    
    % geometry %%%%%%%%%%%
    % init_geometry
    params.geom.n               = [ 64, 64 ,    500  ] ; 
    params.geom.L               = [ 1  , 2*pi,   0.07 ] ; %2*pi since special treatment missing 
                                                        
    params.geom.periodic        = [false,true, false]    ;
    params.geom.forcePeriodic   = [false,false,false]   ;
    params.geom.type            = 'capillary_strech'    ; 
                                                        
    params.geom.z0              = -0.01 ;  
    params.geom.sigma           = 1.2 ;                              
    
    % boundaries %%%%%%%%%%%   
    % init_b
    params.boundaries = {} ; 

    params.boundaries{1}.name    = 'nonreflecting' ;
    params.boundaries{1}.surface = [3,-1] ;

    params.boundaries{2}.name    = 'nonreflecting' ;
    params.boundaries{2}.surface =  [3,1] ;

    %params.boundaries{1}.removeEdge = [1 1;1 -1 ]; 
    %params.boundaries{2}.removeEdge = [1 1;1 -1 ]; 

    params.boundaries{3}.name    = 'nonSlipFlux' ; %  slipEuler' ;
    params.boundaries{3}.surface = [1,-1] ;
    params.boundaries{3}.removeEdge = [3 1 ; 3 -1] ;  
    
    params.boundaries{4}.name    = 'nonSlipFlux' ;     % 'slipEuler' ;
    params.boundaries{4}.surface =  [1,1] ;
    params.boundaries{4}.removeEdge = [3 1 ; 3 -1 ] ;  

    % start %%%%%%%%%%%
    % start_vals_euler
    params.start.type   = 'pulse'   ; 
    params.start.x0     = 0         ;
    params.start.y0     = 0         ;
    params.start.z0     = 0         ;
                    
    params.start.type   = 'linearZWall' ; % 'linearZ' ;
    params.start.lambda = 0.07 ;  
%     params.start.x0  = 0 ; %(params.geom.n(1)> 1)*params.geom.L(1)/2 ;
%     params.start.y0  = 0 ; %(params.geom.n(2)> 1)*params.geom.L(2)/2 ;
%     params.start.z0  = -0.005 ; %(params.geom.n(3)> 1)*params.geom.L(3)/2 ;
%     
%     params.start.sig    = [0.001 0.001 0.001 ] ; %params.geom.L/10  ;

    % io %%%%%%%%%%%%%%
    params.io.dir       = ['/local/reiss' '/3DNavierStokes/'] ; 
    params.io.subFolder = 'capillary3Dstrech' ; 
    params.io.freq      = 20           ; 
                                        
    % plot %%%%%%%%%%%
    % plotFlow 
    params.plot.freq    = 1e10 ;    
    params.plot.shape   = [1,2] ; 
    params.plot.vars    = {'u','p'} ; %,'p','w','Ma'} ; 
    params.plot.subset  = {{':',':',1} ,{':',1,':'} ,{':',1,':'} , {':',1,':'} } ;  
    params.plot.spaceIndex = [1,2; 3,1; 1,3; 1,3]; 

end 
