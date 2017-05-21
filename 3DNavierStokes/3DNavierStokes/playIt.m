clear all 
clc

params.equation.nVar = 5 ; 
params.equation.gamma = 1.4 ; 

params = parameter_capillary3D(params) ;

[ q, params]   = loadData(params,0) ;   
params.io.freq      = 20           ; 
%params.io.freq =100 ; 
%params.plot.vars{3} = 'Ma' ; 
params.plot.vars= {'rho','rho'}; 
for n =10:params.time.steps
    n
    params.time.n = n ; 
    q = loadData(params,n);     
    
    if ~isempty(q) 
        plotFlow(params,q) ; 
        %pause 
    end
    
end 