function [dxdt, dfdx, dfdp] = testFunction3(x,u,p)
    dxdt = -exp(p(1)*x(1)/u(1)) ;
    dfdx = -(p(1)/u(1))*exp(p(1)*x(1)/u(1)) ;
    dfdp = -(x(1)/u(1))*exp(p(1)*x(1)/u(1)) ;
end