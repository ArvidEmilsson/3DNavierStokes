function params = parameter_euclid(params) 

    params.geom.n               = [ 32, 32,   32 ]    ; 
    params.geom.L               = [2*pi, 2*pi,   2*pi ]    ; %2*pi since special treatment missing 

    params.geom.periodic        = [true,true,true]  ;
    params.geom.forcePeriodic   = [true,true,true] ;
    params.geom.type            = 'euclid'            ; 

    params.boundaries = {} ; 

%     params.boundaries{1}.name    = 'slipEuler' ;
%     params.boundaries{1}.surface = [1,0] ;
% 
%     params.boundaries{2}.name    = 'slipEuler' ;
%     params.boundaries{2}.surface =  [1,1] ;

end 