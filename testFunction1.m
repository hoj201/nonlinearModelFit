function [dxdt, dfdx, dfdp] = testFunction1(x,u,p)
    dxdt = [-p(1)*x(1)*u(1); -p(2)*x(2)*u(1)] ;
    dfdx = [-p(1)*u(1) 0; 0 -p(2)*u(1)] ;
    dfdp = [-x(1)*u(1) 0; 0 -x(2)*u(1)] ;
end