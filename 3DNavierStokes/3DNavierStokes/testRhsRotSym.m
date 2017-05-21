

rhs = rhs_skew_NCP(q) ; 

rhsTilde =  rhs_skew_NCP_rotSym(q) ;

rhsRotSym = setRotSym(params, rhsTilde) ; 


diff = rhs - rhsTi; 

err = norm(reshape( diff,[],1)) ; 

disp(['err: ' num2str(err)] ) 