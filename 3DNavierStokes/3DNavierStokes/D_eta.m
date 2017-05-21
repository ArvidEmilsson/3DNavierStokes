function [ u_eta ] = D_eta( u )

global params ; 

[n1,n2,n3]  = size(u) ;
if n2==1 
   u_eta = 0 ; 
   return 
end
Deta    = params.deriv.Deta ;


u_eta        = zeros(n1,n2,n3) ; 

for k =1:n3
     u_eta(:,:,k)    = u(:,:,k)*Deta';  
end 

%temp        = permute(u,[2,1,3] ); 
%temp_eta    = reshape( Deta* reshape(temp ,n2,n1*n3)  ,  n2,n1,n3  )   ; 
%u_eta       = permute(temp_eta, [2,1,3] ) ;   


%u_eta   = Deta*u;  

end
