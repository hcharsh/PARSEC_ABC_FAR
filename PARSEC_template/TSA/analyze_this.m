function QT_interest = analyze_this(y, var_interest)
    % y = max(y, 0) + 1e-8;
    % nphase = size(y, 1);
    % var_n = size(var_interest, 1);
    % for var_ind = 1:var_n
    %     var_index = var_interest(var_ind, 1);
    %     QT_interest(((var_ind - 1)*nphase) + 1: (var_ind*nphase) , 1) = y(:, var_index);
    % end



    
    y = max(y, 0) + 1e-8;
    nphase = size(y, 1);    
%     A = y(:, 1);
%     B = y(:, 2);
%     C = y(:, 3);
%     QT_interest(((1 - 1)*nphase) + 1: (1*nphase) , 1)   = log10(A);
%     QT_interest(((2 - 1)*nphase) + 1: (2*nphase) , 1)   = log10(B);
%     QT_interest(((3 - 1)*nphase) + 1: (3*nphase) , 1)   = log10(C);
    QT_interest(((1 - 1)*nphase) + 1: (1*nphase) , 1)   = y(:, 2);
    QT_interest(((2 - 1)*nphase) + 1: (2*nphase) , 1)   = y(:, 3);
%     QT_interest(((3 - 1)*nphase) + 1: (3*nphase) , 1)   = y(:, 3);

end