
qRef  = params.boundaries{i}.qRef ; 
kx    = params.boundaries{i}.e(:,:,:,1) ; 
ky    = params.boundaries{i}.e(:,:,:,2) ; 
kz    = params.boundaries{i}.e(:,:,:,3) ; 
 
qNow = q(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),:) + params.time.dt*rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),:);

primeU  = zeros(size(qNow)) ; 
    
primeU(:,:,:,1) = qNow(:,:,:,1).^2               ;
primeU(:,:,:,2) = qNow(:,:,:,2)./qNow(:,:,:,1)   ;
primeU(:,:,:,3) = qNow(:,:,:,3)./qNow(:,:,:,1)   ;
primeU(:,:,:,4) = qNow(:,:,:,4)./qNow(:,:,:,1)   ;
primeU(:,:,:,5) = qNow(:,:,:,5)                  ;

c  = sqrt(gamma*(primeU(:,:,:,5)./primeU(:,:,:,1))) ;


du = primeU  - qRef  ;

%keyboard                                                                                                                
R1 = kx.*du(:,:,:,1)                     +  kz.*du(:,:,:,3) -  ky.*du(:,:,:,4)  -            kx.*du(:,:,:,5)./(c.^2)    ;
R2 = ky.*du(:,:,:,1) -  kz.*du(:,:,:,2)                     +  kx.*du(:,:,:,4)  -            ky.*du(:,:,:,5)./(c.^2)    ;
R3 = kz.*du(:,:,:,1) +  ky.*du(:,:,:,2)  -  kx.*du(:,:,:,3)                     -            kz.*du(:,:,:,5)./(c.^2)    ;
R4 =                    kx.*du(:,:,:,2)  +  ky.*du(:,:,:,3) +  kz.*du(:,:,:,4)  +  du(:,:,:,5)./(primeU(:,:,:,1).*c)   ;
R5 =                 -  kx.*du(:,:,:,2)  -  ky.*du(:,:,:,3) -  kz.*du(:,:,:,4)  +  du(:,:,:,5)./(primeU(:,:,:,1).*c)   ;


lambdaS = kx.*primeU(:,:,:,2) + ky.*primeU(:,:,:,3) + kz.*primeU(:,:,:,4)      ;
lambdaP = kx.*primeU(:,:,:,2) + ky.*primeU(:,:,:,3) + kz.*primeU(:,:,:,4)  + c ;
lambdaM = kx.*primeU(:,:,:,2) + ky.*primeU(:,:,:,3) + kz.*primeU(:,:,:,4)  - c ;

projectS= (lambdaS>=0) ;
projectP= (lambdaP>=0); 
projectM= (lambdaM>=0); 

R1  = R1.*projectS ;
R2  = R2.*projectS ;
R3  = R3.*projectS ;
R4  = R4.*projectP ;
R5  = R5.*projectM ;
    
% 
% if(relaxation_param.ne.nll) then                                          ! relaxation_param = dx/L * sigma
% 
% !!$       R5 = zehn/(64.0_rk)/(primeU(:,:,:,2)-c) * 0.58_rk*(1-0.57_rk**2)/(zehn) * (primeU(:,:,:,5)-100000.0_rk)/(primeU(:,:,:,1))
% 
%    M  = sqrt( primeU(:,:,:,2).^2 + primeU(:,:,:,3).^2 + primeU(:,:,:,4).^2 )/c
%    call get_param('po',po)
% 
%    R5 = relaxation_param/(primeU(:,:,:,2)-c)*(eins-maxval(M).^2) * (primeU(:,:,:,5)-po)/(primeU(:,:,:,1))
% 
% end if
half = 0.5 ; 
duN     = zeros(size(qNow)) ;
                                                                                                                
duN(:,:,:,1) =   kx.*R1  +  ky.*R2  +  kz.*R3  +  primeU(:,:,:,1)./(2*c).*R4  +    primeU(:,:,:,1)./(2*c).*R5   ;
duN(:,:,:,2) =           -  kz.*R2  +  ky.*R3  +                 half*kx.*R4  -                   half*kx.*R5   ;
duN(:,:,:,3) =   kz.*R1             -  kx.*R3  +                 half*ky.*R4  -                   half*ky.*R5   ;
duN(:,:,:,4) =  -ky.*R1  +  kx.*R2             +                 half*kz.*R4  -                   half*kz.*R5   ;
duN(:,:,:,5) =                                 + half*primeU(:,:,:,1).*c.*R4  +   half*primeU(:,:,:,1).*c.*R5   ;


targetVals  = qRef + duN    ;            


targetVals(:,:,:,1) = sqrt(targetVals(:,:,:,1))                 ;
targetVals(:,:,:,2) = targetVals(:,:,:,2).*targetVals(:,:,:,1)   ;
targetVals(:,:,:,3) = targetVals(:,:,:,3).*targetVals(:,:,:,1)   ;
targetVals(:,:,:,4) = targetVals(:,:,:,4).*targetVals(:,:,:,1)   ;

rhs(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),:) = (targetVals - q(n0(1):n1(1),n0(2):n1(2),n0(3):n1(3),:))/params.time.dt ;
