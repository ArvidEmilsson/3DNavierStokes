function [q ] =  sharpFilterFieldDynProtoNew(q,r ,kNumMin)
% was  sharpFilterFieldDynProtoNew in 3D acoustic code 

[n1,n2,n3, varNum ] = size(q) ;   
%global c  ;  
%global dt ; 
%global kNumMin ; 
Neta = n2 ; 

%

d = 2*pi*abs(r) ; 

DX = min(  d(1)/((kNumMin-1)*2+1)  , d(end)/Neta ) ;

kNum =  floor((d/DX-1)/2 +1) ;  % which is the max k to have in inner ring same resolution as in outer ring 

kNum = max(kNum,kNumMin )  ; 

kCutOffStart    =  1 +  kNum   ;

n1Start = find(kNum<Neta/2,1,'first') ; 
n1End   = find(kNum<Neta/2,1,'last') ; 


qhFull  = fft(q,[],2) ;   

%if 0
%c = sqrt(c2); 
%dt = params.dt ; 
%keyboard
%n2
%cflAz   = dt./max( d/Neta  , (r*2*pi)./((kNum-1)*2+1) )*c  
%cflRad  = dt./(r(2:end)-r(1:end-1))*c ; 
%keyboard
%end 

%keyboard
for n = n1Start:n1End                    
    qhFull(n, kCutOffStart(n):Neta/2 +1 ,:,: ) = 0 ; 
    qhFull(n,Neta/2 +1:end-kCutOffStart(n)+2,:,: ) = 0 ;  % not need in fortran as fftw for real returns less modes          
end

q = ifft(qhFull,[],2) ;  

end