function [dxdt, dfdx, dfdp] = testFunction2(x,u,p)
    dxdt = -(p(1)+x(1)*u(1))^2 ;
    dfdx = -2*(p(1)+x(1)*u(1))*u(1) ;
    dfdp = -2*(p(1)+x(1)*u(1)) ;
end