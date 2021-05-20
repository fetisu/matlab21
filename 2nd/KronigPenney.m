function [E] = KronigPenney(kk, mmin, aain, bbin, U0in, Emax)
global mm;
global aa;
global bb;
global U0;
global int_k;
int_k = kk(1);
eV_to_J = 1.602176634e-19;
aa = aain*10^(-9); %from nm to m%
bb = bbin*10^(-9); %from nm to m%
U0 = U0in*eV_to_J;
mm = mmin;

Eees = U0in+0.01:0.001:Emax;
ff = Eees;

for i = 1:(length(Eees))
    Ee = Eees(i)*eV_to_J;
    ff(i) = funcf(mm, aa, bb, U0, Ee);
end
figure, hold on, grid on, plot(ff, Eees);


extremums = U0in+0.01;
dif = diff(ff);
for j = 1:(length(dif)-1)
    if (sign(dif(j)) * sign(dif(j+1))) <= 0
        extremums = [extremums, Eees(j)];
    end
end


fk = @funck;
k_len= length(kk);
for q = 1:(length(extremums)-1)
    E_int = [extremums(q), extremums(q+1)];
    for qk = 1:k_len
        int_k = kk(qk);
        if (sign(funck(E_int(1))*sign(funck(E_int(end))))) < 0
            E(q, qk) = fzero(fk, E_int);
        end
    end
end

end



function f = funcf(mm, aa, bb, U0, E)
h = 1.054571817e-34;

mu_q = 2*mm*E/(h*h);
lambda_q = 2*mm*(E-U0)/(h*h);
mu = sqrt(mu_q);
lambda = sqrt(lambda_q);

f = cos(mu*aa)*cos(lambda*bb) - (lambda_q + mu_q)/(2*mu*lambda)*sin(mu*aa)*sin(lambda*bb);
end



function fk = funck(Ee)
eV_to_J = 1.602176634e-19;
global int_k;
global mm;
global aa;
global bb;
global U0;
h = 1.054571817e-34;
E = Ee*eV_to_J;

mu_q = 2*mm*E/(h*h);
lambda_q = 2*mm*(E-U0)/(h*h);
mu = sqrt(mu_q);
lambda = sqrt(lambda_q);

fk = cos(mu*aa)*cos(lambda*bb) - (lambda_q + mu_q)/(2*mu*lambda)*sin(mu*aa)*sin(lambda*bb) - cos(int_k*(aa+bb));
end

