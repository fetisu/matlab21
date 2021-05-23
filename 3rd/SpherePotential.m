function [F,X,Y,P] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)
XYZ = XYZ';
b = b - a .* dot(a,b) ./ dot(a,a);

P = zeros(2,3);
P(1,:) = a;
P(2,:) = b;
dx = (Dx(2) - Dx(1)) / (Nxy(2) - 1);
dy = (Dy(2) - Dy(1)) / (Nxy(1) - 1);
X = Dx(1):dx:Dx(2);
Y = Dy(1):dy:Dy(2);


F = zeros(Nxy(1),Nxy(2));
for icol = 1:Nxy(2)
    for istr = 1:Nxy(1)
        x0 = r0(1) + a(1) * X(icol) + b(1) * Y(istr);
        y0 = r0(2) + a(2) * X(icol) + b(2) * Y(istr);
        z0 = r0(3) + a(3) * X(icol) + b(3) * Y(istr);
        ro = [x0; y0; z0];
        
        k = 1:N;
        k_dist = norm(XYZ(:, k) - ro);
        k_dist = k_dist - R(k);
        k_dist = R(k) + max(k_dist, 0);
        F(istr, icol) = F(istr, icol) + Q(k) / k_dist;
    end
end      


X = repmat(X, Nxy(1), 1);
Y = repmat(Y', 1, Nxy(2));
end

