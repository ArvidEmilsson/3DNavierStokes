function h = pcolorFull( X,Y,u )
%PCOLORFUNCTION h = pcolorFull( X,Y,u )
% as pcolor with plotting the boundary points. 
% Boundary cells are of half size
% X, Y are to be provided

%nargin

[Nxi, Neta] = size(u);   

Xmid = zeros(Nxi+1,Neta+1) ;
Ymid = zeros(Nxi+1,Neta+1) ;

coreX = (X(1:end-1,1:end-1) + X(1:end-1,2:end) + X(2:end,1:end-1) + X(2:end,2:end) )/4  ; 
coreY = (Y(1:end-1,1:end-1) + Y(1:end-1,2:end) + Y(2:end,1:end-1) + Y(2:end,2:end) )/4  ; 

Xmid(2:end-1,2:end-1)  = coreX ;  
Ymid(2:end-1,2:end-1)  = coreY ;  

Xmid(1,2:end-1) = (X(1,1:end-1) + X(1,2:end) )/2  ;  
Ymid(1,2:end-1) = (Y(1,1:end-1) + Y(1,2:end) )/2  ;  

Xmid(end,2:end-1) = (X(end,1:end-1) + X(end,2:end) )/2  ;  
Ymid(end,2:end-1) = (Y(end,1:end-1) + Y(end,2:end) )/2  ;  

Xmid(2:end-1,1) = (X(1:end-1,1) + X(2:end,1) )/2  ;  
Ymid(2:end-1,1) = (Y(1:end-1,1) + Y(2:end,1) )/2  ;  

Xmid(2:end-1,end) = (X(1:end-1,end) + X(2:end,end) )/2  ;  
Ymid(2:end-1,end) = (Y(1:end-1,end) + Y(2:end,end) )/2  ;  

Xmid(1,1) = (X(1,1) ) ;  
Ymid(1,1) = (Y(1,1) ) ;  

Xmid(1,end) = (X(1,end) ) ;  
Ymid(1,end) = (Y(1,end) ) ;  

Xmid(end,1) = (X(end,1) ) ;  
Ymid(end,1) = (Y(end,1) ) ;  

Xmid(end,end) = (X(end,end) ) ;  
Ymid(end,end) = (Y(end,end) ) ;  

%plot(Xmid,Ymid,'r+') 


Umid = zeros(Nxi+1,Neta+1) + u(1,1)  ; % to avoid messing up the colorbar  
Umid(1:Nxi,1:Neta)  = u ; 

h = pcolor(Xmid,Ymid,Umid) ; 


end

