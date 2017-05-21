
for l =1:length(params.boundaries)
   
    n0      = params.boundaries{l}.n0 ;   
    n1      = params.boundaries{l}.n1 ;   

    lambda(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3))  = 0 ; 
    
end
