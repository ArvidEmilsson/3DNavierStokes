function params = parameter_pipe_rXi(params) 

    params.equation.dissipation = false ; 
    params.equation.polarFilter = true  ; 
    
    % time     %%%%%%%%%%%
    params.time.steps           = 1000 ; 
    params.time.CFL             = 1   ; 
    params.time.checkNaN        = 1     ; 
    
    % geometry %%%%%%%%%%%
    % init_geometry
    params.filter.freq          = 1e50    ; 
     
    params.geom.n               = [ 64,  64,  1 ]    ; 
    params.geom.L               = [  .5e-3, 2*pi,  2 ]    ; %2*pi since special treatment missing 

    params.geom.periodic        = [false,true,true]  ;
    params.geom.forcePeriodic   = [false,false,true] ;
    params.geom.type            = 'pipe_rXI_strech'            ; 
    
    params.geom.sigma           = .6  ;    
    
    params.boundaries = {} ; 


    params.boundaries{1}.name    = 'slipEuler' ;
    params.boundaries{1}.surface = [1,-1] ;
 
    params.boundaries{2}.name    = 'slipEuler' ;
    params.boundaries{2}.surface =  [1, 1] ;
     
    % start %%%%%%%%%%%
  
    params.start.type   = 'pulse' ;

    params.start.x0  = (params.geom.n(1)> 1)*params.geom.L(1)/4 ;
    params.start.y0  = (params.geom.n(2)> 1)*params.geom.L(1)/4 ;
    params.start.z0  = (params.geom.n(3)> 1)*params.geom.L(3)/2 ;

    params.start.sig    = [ 0.2 0.2 10]*params.geom.L(1)   ;

     
    % io %%%%%%%%%%%%%%
    
    params.io.subFolder = 'pipe_rXi'; 
    params.io.freq      = 25      ; 
    params.io.dir       = ['/work/reiss' '/3DNavierStokes/'] ;  
    % plot %%%%%%%%%%%

    params.plot.freq    = 25    ; 
    params.plot.shape   = [1,2] ;                         
    params.plot.vars    = {'p'} ;
    params.plot.subset  = {{':',':',1} , {':',1,':'} } ; 
    params.plot.spaceIndex = [1,2 ; 3,1]     ;      
                                    
end 