function [Q,D] = ElectroStaticDipoles(XYZ,R,F)
XYZ = XYZ';
N = length(R);
F = [F; zeros(3 * N, 1)];

Elem = 1:N;
Elem = Elem * ones(1, N);
Distances = norm(XYZ(:, Elem) - XYZ(:, Elem'));
SumR = R(Elem) + R(Elem');

if Distances < SumR
    error('Data is not possible, check the distances between balls')
end


M = zeros(4 * N);
for istr = 1 : N
    for icol = 1 : 4 * N
        j = mod(icol, N);
        if(j == 0)
            j = N;
        end
        co = (icol - j) / N;
        r = XYZ(1:3,istr)-XYZ(1:3,j);
        if(istr == icol)
            M(istr, icol) = 1 / R(istr);
        elseif(co == 0)
            M(istr, icol) = 1 / norm(r);
        elseif(j == istr)
            M(istr, icol) = 0;
        else
            M(istr, icol) = r(co) / norm(r) ^ 3; 
        end   
    end
end


for istr = (N + 1) : (4 * N)
    for icol = 1 : (4 * N)
        j = mod(icol, N);
        i = mod(istr, N);
        if(j == 0)
            j = N;
        end
        if(i == 0)
            i = N;
        end
        co_p = (icol - j) / N;
        co_E = (istr - i) / N;
        r = XYZ(1:3,i)-XYZ(1:3,j);
        if(i == icol)
            M(istr, icol) = 0;
        elseif(istr == icol)
            M(istr, icol) = 1 / R(i) ^ 3;
        elseif(co_p == 0)
            M(istr, icol) = r(co_E) / norm(r) ^ 3;
        elseif(j == i)
            M(istr, icol) = 0;
        elseif(co_E == co_p)
            M(istr, icol) = (3 * r(co_E) ^ 2 - norm(r)^2 ) / norm(r)^5;
        else
            M(istr, icol) = 3 * r(co_E) * r(co_p)/ norm(r)^5;
        end 
    end
end

Pr = M \ F;
Q = Pr(1:N, 1);
Dx = Pr(N + 1:2 * N,1)';
Dy = Pr(2 * N + 1:3 * N,1)';
Dz = Pr(3 * N + 1:4 * N,1)';
D = [Dx; Dy; Dz]';
D = -D;
end



