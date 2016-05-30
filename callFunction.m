function out = callFunction(fun,in)
    in = num2cell(in) ;
    out = fun(in{:}) ;
end