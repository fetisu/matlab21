t = 0:0.21:20^(2);
a = 0:0.03:2.83;
z = (t.^(1/2)).*-1;
z = [z.*-1, z, a, a.*-1];
z = sort(z);
x = fresnelc(z);
y = fresnels(z);
figure; hold on; grid on; plot(x, y);
