function saveData(params,q,n)

persistent wasCalled
    %q_single  = single(q) ;
    
    if  (params.io.freq<=params.time.steps ),
        if  (~mod(n,params.io.freq)),
            
            workDir     = [ params.io.dir  params.io.subFolder '/'] ;
            
            if isempty(wasCalled)
                wasCalled = true ; 
                [~,~] = mkdir(workDir)  ; 
        
            end

%            fNameTemplate   = [sprintf('/dump_%d', caseNumber)  '_%010d.mat' ] ; 
            fName  = sprintf(['/dump_'  '%010d.mat'],n ) ; 
            %tic
            save([workDir fName ],'q' ,'-v6')  
            %toc
            if n == 0 
                save([workDir 'params.mat' ],'params' )  
                disp(workDir)
            end
            
        end    
    end 
 

end