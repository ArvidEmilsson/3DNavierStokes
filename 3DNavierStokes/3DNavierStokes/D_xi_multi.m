function [ u_xi ] = D_xi_multi( u )


global params ; 
Dxi     = params.deriv.Dxi;

[n1,n2,n3,nVar]  = size(u) ;
%u_xi        = zeros(n1,n2,n3) ; 

%for k =1:n3
%    u_xi(:,:,k)    = Dxi*u(:,:,k);  
%end 

u_xi =  reshape( Dxi*reshape(u, n1,n2*n3*nVar),  n1,n2,n3,nVar) ;  

end
