clear all ;  close all;  clc;

global params

path(path, genpath('./LIB/'))

% TODO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  reibungsterme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nicht ref rb 
% wall rb, mit def wärme fluss 
% filter pipe
% optimiren
% rot sym 
% gl lösen

% equation 
params.equation.nVar = 5 ; 
params.equation.gamma = 1.4 ; 

% params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = parameter_pipe(params) ;
params = parameter_pipe_rXi(params) ;
params = parameter_euclid(params) ;
params = parameter_test(params) ;
params = parameter_capillary(params) ;
params = parameter_capillary3D(params) ;
params = parameter_capillaryRotSym(params) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params      = init_geometry(params)     ; 
params      = init_deriv(params)        ; 
params      = init_metric(params)       ; 
[params, q] = start_vals_euler(params)  ;   
params      = init_bound(params,q)      ; 
params      = init_filter(params)       ; 

% time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = min(params.geom.dxi.*params.geom.L  ); % more a guess (a dirty one)

dt = params.time.CFL*dx/sqrt(params.equation.gamma*1e5/1);  
params.time.dt = dt ; 

% plotten %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  (params.plot.freq<params.time.steps || params.time.steps==0) , 
    params.time.n = 0 ; 
    plotFlow(params,q) ;  
end 

saveData(params, q, 0 ) ;  

% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for n =1:params.time.steps
    
    %%%%
    params.time.n = n ;
    %%%%%%%%%%%%%% time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    q1 = RK4(q,dt, 0, @rhs_skew_NCP ) ;
    t = toc ;
    if (  ~mod(n,params.time.checkNaN)) , 
        if any( reshape(q1(:,:,:,1)<0,[],1) )
            disp('sqrho neg')
            break 
        end
        if any(~isreal( reshape(q1(:,:,:,1),[],1)) )
            disp('imag vals')
            break 
        end
            
        if hasInfNaN(reshape(q1,[],1))
            disp(['has inf nan @ ' num2str(n) ])
            break 
        end
    end     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    q = q1 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (  ~mod(n,1)) , 
        disp(['step time ' num2str(t) ' n= ' num2str(n) ])
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  (  ~mod(n,params.plot.freq)) , 
        plotFlow(params,q) ;
        %pause
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (~mod(n,params.filter.freq))
        %disp('filter') 
        q = applyFilter(params,q)       ; 
        q = enforceBondary(params,q)    ;  
    end
    
    % save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    saveData(params, q, n ) ;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
