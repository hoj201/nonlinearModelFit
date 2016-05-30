function [dxdt, dfdx, dfdp] = testNumJacobFunction(x,u,p)
    dxdt = [x(1)*u(1)*p(1);
            x(2)*u(1)*p(2);
            x(3)*u(1)*p(3);
            x(4)*u(1)*p(1)] ;
        
    dfdx = [u(1)*p(1) 0 0 0;
            0 u(1)*p(2) 0 0;
            0 0 u(1)*p(3) 0;
            0 0 0 u(1)*p(1)];
        
    dfdp = [x(1)*u(1) 0 0;
            0 x(2)*u(1) 0;
            0 0 x(3)*u(1);
            x(4)*u(1) 0 0];
end