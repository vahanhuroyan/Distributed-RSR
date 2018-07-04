function [ outMat ] = consGMSF( X, C )  
% This function soles the gms problem with the linear part added with
% coefficient C
    reg_param = 1e-10;
    [D, n] = size(X);
    iterNumb = 20;
    Q_0 = eye(D) / D;
    for iter = 1:iterNumb
        temp = sum(( X' * Q_0).^2, 2);
        covMat = (X.*repmat(0.5*min(temp.^-0.5,1 / reg_param), 1, D)' * X');
        
        if(any(any(covMat == Inf)) || any(any(isnan(covMat))))
            break;
        end
        Q_1 = traceOneLyap(covMat, C);
        Q_0 = Q_1;
    end
    outMat = Q_1;
    outMat = outMat / trace(outMat);
end

