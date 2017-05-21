function [w, D, fac ] = D_openFourthOrdWeightsStrand(N,h)
% Expicit Derivate Matrix 4th order, divisor returned 
% 
% N number of Points
% h = Delta x   

a1 = 2/3; 
a2 =-1/12;

alpha1= 2/3 ; %
alpha2= -1/12 ; %
alpha3= 0 ;  %

 
D= sparse(alpha1*(diag(ones(N-1,1),1)-  diag(ones(N-1,1),-1))      +...  
   alpha2*(diag(ones(N-2,1),2)-  diag(ones(N-2,1),-2)  )  +...        
   alpha3*(diag(ones(N-3,1),3)-  diag(ones(N-3,1),-3))) ;


ww = [17 59 43 49]/48; 
%ww = [13649/43200 12013/8640 2711/4320 5359/4320 7877/8640 43801/43200 ]; 

boundVals =  ...
    [-24/17 59/34 -4/17 -3/34  0     0 ...
    ; -1/2  0     1/2    0     0     0 ...
    ; 0     0     0      0     0     0 ... 
    ; 3/98  0    -59/98  0     32/49 -4/49 ] ;

 
WW = diag(ww) ;
boundVals = WW* boundVals ;


%keyboard   
[l, c]=size(boundVals); 

D(1:l,1:c) = boundVals; 
D(end-l+1:end,end-c+1:end) = -flipud(fliplr( boundVals )); 

b = ones(1,N)*D; 

b(1)  = b(1)+1;   
b(end)  = b(end)-1;   

D (3,1:c) = - b(1:c); 
D (end-2,end-c+1:end) = - b(end-c+1:end); 

fac=h;     
w      = ones(N,1) ; 

w(1:l)   = ww; 
w(end-l+1:end) = fliplr(ww); 


%D= 1/(h)*D ; 

end
