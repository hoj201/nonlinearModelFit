function [dxdt, dfdx, dfdp] = testDubins(x,u,p)
    dxdt = [p(1)*sin(x(3))
            p(2)*cos(x(3));
            u(1)] ;
        
    dfdx = [0 0 p(1)*cos(x(3));
            0 0 -p(2)*sin(x(3));
            0 0 0] ;
        
    dfdp = [sin(x(3))   0 ;
            0           cos(x(3)) ;
            0           0] ;
end