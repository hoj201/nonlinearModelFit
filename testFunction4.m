function [dxdt, dfdx, dfdp] = testFunction4(x,u,p)
    dxdt = [-p(1)*u(1)*x(2); -p(1)*u(1)*x(1)] ;
    dfdx = [0 -p(1)*u(1) ; -p(1)*u(1) 0] ;
    dfdp = [-u(1)*x(2) ; -u(1)*x(1)] ;
end