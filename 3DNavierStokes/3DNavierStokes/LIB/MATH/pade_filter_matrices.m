function [ A,B ] = pade_filter_matrices( type , order , alpha , N)
%  [ A,B ] = pade_filter_matrices( type , order , alpha , N)
% information
% reference         : GaitondeVisbal2000
%
% choice of alpha   : alpha = 0     --> explicit filter
%                   : alpha = 0.5   --> no filtering
%
% Bu' = Au --> u' = B\Au
%
% 2013-06-19 (ml), 2013-XX-XX (??), ...
%
% TODO: testen, erweiterung auf höhere ordnungen, ...
%% parameter
% set default
default.type  = 'non-per';
default.order = 4   ;
default.alpha = 0.35;

if nargin == 0
    type  = default.type;
    order = default.order;
    alpha = default.alpha;
end

if nargin == 1
    order = default.order;
    alpha = default.alpha;
end

if nargin == 2
    alpha = default.alpha;
end

% check alpha  values
if alpha < -0.5 || alpha > 0.5
    disp(['error in pade_filter.m: non-suitable choice of alpha: ',num2str(alpha)])
    disp( '                        possible alpha -0.5 < alpha < 0.5')
    return
end

% determine size and reshape to vector

%% filter
% B == matrix =============================================================
switch type
    case 'non-per'
        d = ones(N+0,1);
        s = ones(N-1,1)*alpha;
        
        % build matrix (sparse)
        B = gallery('tridiag',s,d,s);
        
        % do not modify bounds
        B(1,2)          = 0;
        B(end,end-1)    = 0;
        
    case 'per'
        e = ones(N,1);
        
        % build matrix (sparse)
        B = sparse(1:N,[2:N 1],alpha*e,N,N);
        B = B + B' + sparse(1:N,1:N,e,N,N);
        
    otherwise
        disp(['error in pade_filter.m: type not implemented: ',type])
        disp( '                        possible: per and non-per')
        return
end

% A == matrix =============================================================
% coeffients (interior)
a = zeros(1,6 );
b = zeros(5,11);

switch order
    case 2 %...............................................................
        a(1) = +1/2     + alpha;
        a(2) = +1/2     + alpha;        
        as   = 1; % stencil width
        
    case 4 %...............................................................
        a(1) = +5/8     + 3*alpha/4;
        a(2) = +1/2     + alpha;
        a(3) = -1/8     + alpha/4;        
        as   = 2; % stencil width
        
        % point 2 (non-per case)
        b(1,1) = +1/16  + 7*alpha/8;
        b(1,2) = +3/4   + alpha/2;
        b(1,3) = +3/8   + alpha/4;
        b(1,4) = -1/4   + alpha/2;
        b(1,5) = +1/16  - alpha/8;
        
    case 6 %...............................................................
        a(1) = +11/16   + 5*alpha/8;
        a(2) = +15/32   + 17*alpha/16;
        a(3) = -3/16    + 3*alpha/8;        
        a(4) = +1/32    - alpha/16;
        as   = 3; % stencil width
        
        % point 2 (non-per case)
        b(1,1) = 1/64  + 31*alpha/32;
        b(1,2) = 29/32 + 3*alpha/16;
        b(1,3) = 15/64 + 17*alpha/32;
        b(1,4) = -5/16 + 5*alpha/8;
        b(1,5) = 15/64 - 15*alpha/32;
        b(1,6) = -3/32 + 3*alpha/16;
        b(1,7) = 1/64  - 1*alpha/32;
        
        % point 3 (non-per case)
        b(2,1) = -1/64  + 1*alpha/32;
        b(2,2) = 3/32   + 13*alpha/16;
        b(2,3) = 49/64  + 15*alpha/32;
        b(2,4) = 5/16   + 3*alpha/8;
        b(2,5) = -15/64 + 15*alpha/32;
        b(2,6) = 3/32   - 3*alpha/16;
        b(2,7) = -1/64  + 1*alpha/32;
    case 10
        % coefs from VisbalGaitonde2002; only periodic case in this paper 
        % boundary terms are in GaitondeVisbal2000

        a(1) = (193 + 126*alpha)/256 ;
        a(2) = (105 + 302*alpha)/256 ;
        a(3) = (15*(-1+2*alpha))/64  ;    
        a(4) = (45*(1-2*alpha))/512  ; 
        a(5) = (5*(-1+2*alpha))/256  ; 
        a(6) =      (1-2*alpha)/512  ; 
        as   = 5 ; % stencil width, as defined by Matze 
        
        % non per case missing!
        if strcmp(type,'non-per'), error('10th order filter only non-per')  ; end 
            
    otherwise
        disp(['error in pade_filter.m: order not implemented: ',num2str(order)])
        disp( '                        possible: 2, 4, 6, 10(non-per)')
        return
end

% A -- matrix -------------------------------------------------------------
a = 0.5*a;
switch type
    case 'non-per'
        % build matrix (sparse)
        A = sparse(1:N,1:N,a(1)*ones(N,1),N,N);
        for i = 1:5
            A = spdiags(a(1+i)*ones(N-i,1),-i,A);
        end
        A = A+A';
        
        % do not modify bounds
        A(1,:) = 0;
        A(N,:) = 0;
        A(1,1) = 1;
        A(N,N) = 1;
        
        for j = 1:(as-1)
            A(j+1, 1       :(as*2+1)) =        b(j,1:(as*2+1));
            A(N-j, N-(as*2):N       ) = fliplr(b(j,1:(as*2+1)));
        end        
        
    case 'per'
        e = ones(N,1);
        
        % build matrix (sparse)
        A =   sparse(1:N,[2:N 1]        ,a(2)*e,N,N) ...
            + sparse(1:N,[3:N 1 2]      ,a(3)*e,N,N) ...
            + sparse(1:N,[4:N 1 2 3]    ,a(4)*e,N,N) ...
            + sparse(1:N,[5:N 1 2 3 4]  ,a(5)*e,N,N) ...
            + sparse(1:N,[6:N 1 2 3 4 5],a(6)*e,N,N);
        
        A = A+A'+sparse(1:N,1:N,2*a(1)*e,N,N);
        
    otherwise
        disp(['error in pade_filter.m: type not implemented: ',type])
        disp( '                        possible: per and non-per')
        return
end

%% filter

end