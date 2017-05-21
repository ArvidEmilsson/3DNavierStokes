function [ u_zeta ] = D_zeta( u )


global params ; 
Dzeta     = params.deriv.Dzeta;

[n1,n2,n3]  = size(u) ;
%u_zeta        = zeros(n1,n2,n3) ; 

%for i =1:n1
%    u_zeta(i,:,:)    = squeeze(u(i,:,:))*Dzeta';  
%end 
u_zeta    = reshape( reshape(u,n1*n2,n3)*Dzeta' , n1,n2,n3) ;  

end
