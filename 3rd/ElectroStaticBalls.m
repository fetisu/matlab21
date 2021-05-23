function Q = ElectroStaticBalls(XYZ,R,F)
XYZ = XYZ';
N = length(R);

Elem = 1:N;
Elem = Elem * ones(1, N);
Distances = norm(XYZ(:, Elem) - XYZ(:, Elem'));
SumR = R(Elem) + R(Elem');

if Distances < SumR
    error('Data is not possible, check the distances between balls')
end


Distances = Distances + diag(ones(1, N));
SLE = 1./Distances;
sle = 1./R(1:N);
sle = diag(sle);
SLE = SLE - diag(ones(1, N)) + sle;
            
Q = SLE \ F;
