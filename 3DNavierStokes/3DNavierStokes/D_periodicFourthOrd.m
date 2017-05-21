function D = D_periodicFourthOrd(N,h)
% Expicit Derivate M        atrix with improved wave numbers of Tam and Webb
% Periodic Version
% D = D_TW(N,h)
% N number of Points
% h = Delta x   

alpha2= -1/12 ; %-1/(h)*0.18941; %*(2*pi)/h;
alpha3= 0 ; %1/(h)*0.02652; % *(2*pi)/h;
alpha1= 8/(12) ; %)-2*alpha2-3*alpha3   ; 

 
D= sparse(alpha1*(diag(ones(N-1,1),1)-  diag(ones(N-1,1),-1)  + diag(ones(1,1),-N+1) -  diag(ones(1,1),N-1))    +...  
   alpha2*(diag(ones(N-2,1),2)-  diag(ones(N-2,1),-2) +  diag(ones(2,1),-(N-2))- diag(ones(2,1),N-2))  +...        
   alpha3*(diag(ones(N-3,1),3)-  diag(ones(N-3,1),-3) +  diag(ones(3,1),-(N-3))- diag(ones(3,1),N-3) ));


D= 1/(h)*D ; 
end
