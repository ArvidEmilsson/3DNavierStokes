function [w, D, fac ] = D_openSixOrdWeightsStrand(N,h)
% Expicit Derivate Matrix 4th order, divisor returned 
% 
% N number of Points
% h = Delta x   

alpha1= 3/4 ; %
alpha2= -3/20 ; %
alpha3= 1/60 ;  %

 
D= sparse(alpha1*(diag(ones(N-1,1),1)-  diag(ones(N-1,1),-1))      +...  
   alpha2*(diag(ones(N-2,1),2)-  diag(ones(N-2,1),-2)  )  +...        
   alpha3*(diag(ones(N-3,1),3)-  diag(ones(N-3,1),-3))) ;


%D(1:2,1:5) = [-7 8 -1  0  0 ; -6 -1 8 -1 0 ]; 
%D(end-1:end,end-4:end) = -flipud(fliplr([-7 8 -1  0 0 ; -6 -1 8 -1 0 ])); 

ww = [13649/43200 12013/8640 2711/4320 5359/4320 7877/8640 43801/43200 ]; 

boundVals =  ...
[ -21600/13649   81763/40947 131/27298  -9143/13649 20539/81894 0 0 0 0 ...
; -81763/180195  0           7357/36039 30637/72078 -2328/12013 6611/360390  0 0 0 ...
; -131/54220    -7357/16266  0          645/2711    11237/32532 -3487/27110  0 0 0 ...
; 9143/53590    -30637/64308 -645/5359  0           13733/32154 -67/4660     72/5359 0 0 ...
; -20539/236310 2328/7877    -11237/47262 -13733/23631 0        89387/118155 -1296/7877 144/7877 0 ...
; 0             -6611/262806 3487/43801  1541/87602 -89387/131403 0          32400/43801 -6480/43801 720/43801 ];


WW = diag(ww) ;
boundVals = WW* boundVals ;


%keyboard   
[l, c]=size(boundVals); 

D(1:l,1:c) = boundVals; 
D(end-l+1:end,end-c+1:end) = -flipud(fliplr( boundVals )); 


fac=h;     
w      = ones(N,1) ; 

w(1:l)   = ww; 
w(end-l+1:end) = fliplr(ww); 


%D= 1/(h)*D ; 

end
