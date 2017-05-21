function q = setRotSym(params,q)

[n1,n2,n3, varNum] = size(q) ; 


%setRotsSym sqRho, p 

for j=1:n2
    q(1:n1/2,j,:,1) = q(1:n1/2,1,:,1); 
    q(n1/2+1:n1,j,:,1) = q(n1/2:-1:1,1,:,1); 
    
    q(1:n1/2,j,:,5) = q(1:n1/2,1,:,5); 
    q(n1/2+1:n1,j,:,5) = q(n1/2:-1:1,1,:,5); 
    
    q(1:n1/2,j,:,4) = q(1:n1/2,1,:,4); 
    q(n1/2+1:n1,j,:,4) = q(n1/2:-1:1,1,:,4); 
    
end

rhoUrad = q(1:n1/2,1,:,2); 

%X = params.geom.X ;

for j=1:n2
    %j
    e_r_1 = params.geom.X(1:n1/2,j,:,1) ; 
    e_r_2 = params.geom.X(1:n1/2,j,:,2) ;
    
    n = sqrt(e_r_1.^2  + e_r_2.^2) ;
    
    e_r_1 = e_r_1./n ;  
    e_r_2 = e_r_2./n ;

   %keyboard

    q(1:n1/2,j,:,2) = -e_r_1.*rhoUrad; 
    %keyboard
    q(n1/2+1:n1,j,:,2) = e_r_1(n1/2:-1:1,:).*rhoUrad(n1/2:-1:1,:); 

    q(1:n1/2,j,:,3) = -e_r_2.*rhoUrad; 
    q(n1/2+1:n1,j,:,3) = e_r_2(n1/2:-1:1,:).*rhoUrad(n1/2:-1:1,:); 
    
    
end


end 