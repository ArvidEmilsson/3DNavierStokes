function u1 =RK4(u0,dt, t,rhs)
% Runge Kutta 4. Ordnung

% die vier Aufrufe der rechten Seiten; 
% egal, ob rhs zahl, vektor  oder vektoren=matrix zurueckgibt: 
% solange Dimension u entspricht geht alles glatt.  
 %disp('calc k1')
k1=rhs(u0          , t     );
 %disp('calc k2')
k2=rhs(u0+(dt/2)* k1 , t+dt/2);
 %disp('calc k3')
k3=rhs(u0+(dt/2)* k2 , t+dt/2);
 %disp('calc k4')
k4=rhs(u0+dt* k3   , t+dt  );



u1 = u0+(dt/6) * (k1 + 2*(k2+k3) +k4) ;

