function [params, q ] = start_vals(params) 

    X = params.geom.X ;
    Y = params.geom.Y ;
    Z = params.geom.Z ;
    
    N = size(X) ;
    q = zeros([N params.equation.nVar]) ; 
    x0 =pi ; 
    y0 =pi ; 
    z0 =pi ; 
    sigX2  = (pi/6)^2 ; 
    sigY2  = (pi/6)^2 ; 
    sigZ2  = (pi/6)^2 ;
    a= .1 ; 
    q(:,:,:, 1) = a*exp( -(X-x0).^2/sigX2 -(Y-y0).^2/sigY2 -(Z-z0).^2/sigZ2  ) ; 
    
end