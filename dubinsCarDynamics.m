function [dxdt, dfdx, dfdp] = dubinsCarDynamics(x,u,p)
    dxdt = [p(3)*u(2)*sin(p(1)*x(3)) ;
            p(3)*u(2)*cos(p(2)*x(3)) ;
            u(1)                    ];
        
    dfdx = [0, 0,  p(3)*u(2)*p(1)*cos(p(1)*x(3)) ;
            0, 0, -p(3)*u(2)*p(2)*sin(p(2)*x(3)) ;
            0, 0,  0                       ];
        
    dfdp = [p(3)*u(2)*x(3)*cos(p(1)*x(3)), 0, u(2)*sin(p(1)*x(3));
            0, -p(3)*u(2)*x(3)*sin(p(2)*x(3)), u(2)*cos(p(2)*x(3));
            0, 0, 0];
end