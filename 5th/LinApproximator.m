function [P,sgP] = LinApproximator(y,r,funcs)
N = size(y, 2);
% K = size(r, 1);
M  = size(funcs, 1);
% y - вектор из N эксперементальных точек
% r - вектор из N векторов размера k
% funcs - вектор из M функций степени k


auxilary = zeros(N, M);
for istr = 1 : N
    for icol = 1 : M
        f = cell2mat(funcs(icol));
        vec = num2cell(r(:, istr));
        auxilary(istr, icol) =  f(vec{:});
    end
end

Y = (y*auxilary)';
G = auxilary'*auxilary;
a = G\Y;
P = a;

Err = sumsqr(y - auxilary*a)/N;
sgP = sqrt(Err/(2*trace(G)));

end
