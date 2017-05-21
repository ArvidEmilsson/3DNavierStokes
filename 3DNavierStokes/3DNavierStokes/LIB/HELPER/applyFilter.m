function [q] = applyFilter(params,q)


[n1,n2,n3, varNum] = size(q) ; 


if params.geom.n(1)> 1     
    q =  reshape(  params.filter.B{1}\(params.filter.A{1}*reshape(q, n1,n2*n3*varNum)),  n1,n2,n3,varNum) ;  
end


if (params.geom.n(2)> 1  && params.filter.eta)
    for l =1:varNum

        for k =1:n3
             q(:,:,k,l)  = (q(:,:,k,l)*params.filter.A{2}')/params.filter.B{2}  ;  
        end     
    end

end



if (params.geom.n(3)> 1  && params.filter.zeta ) 
    for l =1:varNum
        q(:,:,:,l) =  reshape( (reshape(q(:,:,:,l), n1*n2,n3)*params.filter.A{3}')/params.filter.B{3},  n1,n2,n3) ;  
    end
    %u_zeta    = reshape( reshape(u,n1*n2,n3)*Dzeta' , n1,n2,n3) ;  

end


end 