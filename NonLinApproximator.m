function [P, sgP] = NonLinApproximator (y,r,fun, P_0)
% P_0 is a line-array with starting guess for values of parameters
% p_number - число параметров P
N = size(y, 2);
p_number = length(P_0);

for iterations = 1 : 420
    
    F = zeros(1, N);
    ii = 1 : N;
    F(ii) = fun(r(:, ii), P_0); 


    % якобиан
    J = jacobi(fun, r, p_number, P_0);

    % из уравнения G*(P - P_0) = J*(y - fun(r, P_0))
    Y = J'*(y' - F');
    G = J'*J;
    P = P_0' + G\Y;

    if isnan(P) || (iterations > 1 && abs(P_0 - P_last) < 1e-10)
        break
    end

    P_last = P;
    P_0 = P';
end

sgP = 'have no idea how(';
end

function [number] = find_arg_number(fun, N)
ii = 1;
exception_identifier = 'MATLAB:badsubscript';
while strmatch(exception_identifier, 'MATLAB:badsubscript')
    exception_identifier = ' ';
    try fun(zeros(1, N), zeros(1, ii))
    catch E
        exception_identifier = E.identifier;
        ii = ii + 1;
    end
end

number = ii;
end

function [J] = jacobi(fun, r, p_number, P_0)
% h - шаг вдоль оси
h = 1e-6; 
N = length(r(1, :));
J = zeros(N, p_number);

for istr = 1 : N
    for icol = 1 : p_number
        delta = zeros(1, p_number);
            delta(icol) = h*abs(P_0(icol));
            J(istr, icol) = (fun(r(:, istr),P_0 + delta) - fun(r(:, istr),P_0 - delta))/(2*h);
    end   
end

end