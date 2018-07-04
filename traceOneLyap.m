function [ Q ] = traceOneLyap( X, A )
% This function finds the trace one solution and the corresponding c 
%of the following lyapunov equation XQ + Qx + A = cI and A

    tr0 = trace(lyap(X, A));
    tr1 = trace(lyap(X, A - eye(size(A, 1))));
    cVal = 1 - (0 - 1) / (tr0 - tr1) *(tr1 - 1);
    Q = lyap(X, A - cVal*eye(size(A, 1)));
end

