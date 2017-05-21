function [ q, params]  = loadData(params,n)


    q = [] ; 
    
    if  (params.io.freq<=params.time.steps ),
        if  (~mod(n,params.io.freq)),
            workDir     = [ params.io.dir  params.io.subFolder '/'] ;
            
            
            if n == 0 
                load([workDir 'params.mat' ] )  

            end
            
            fName  = sprintf(['/dump_'  '%010d.mat'],n )  
            
            load([workDir fName ] ) ;   
            
            
        end 
    else 
        warning('params.io.freq<=params.time.steps, any data?') 
        
    end 
 

end