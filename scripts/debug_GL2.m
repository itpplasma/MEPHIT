old = load('../run/33353_2325/mephit_old.h5');
new = load('../run/33353_2325/mephit_new.h5');
n = double(old.mesh.n);

kf = 30;
k_min = old.mesh.kt_low(kf) + 1;
k_max = old.mesh.kt_low(kf) + old.mesh.kt_max(kf);
k_off = old.mesh.npoint - 1;
ndim = old.mesh.kt_max(kf);
area = old.mesh.area(k_min:k_max);
orient = double(old.mesh.orient(k_min:k_max)) * 2 - 1;
e_f = old.mesh.tri_edge(1, k_min:k_max);

x_old = old.debug_currn.x(k_min:k_max);
d_old = old.debug_currn.d(k_min:k_max);
du_old = old.debug_currn.du(k_min:k_max);
I_f_old = old.debug_currn.I(e_f) .* orient;
rem_old = x_old + I_f_old;
I_old = old.debug_currn.I((k_off+k_min):(k_off+k_max));
grad_pn_old = old.debug_currn.grad_pn(:, k_min:k_max);
lorentz_old = old.debug_currn.lorentz(:, k_min:k_max);
terms_old = old.debug_currn.terms(:, :, k_min:k_max);
x_new = new.debug_currn.x(k_min:k_max);
d_new = new.debug_currn.d(k_min:k_max);
du_new = new.debug_currn.du(k_min:k_max);
I_f_new = new.debug_currn.I(e_f) .* orient;
coeff_f = new.debug_currn.coeff_f(k_min:k_max);
rem_new = x_new + I_f_new .* coeff_f;
I_new = new.debug_currn.I((k_off+k_min):(k_off+k_max));
grad_pn_new = new.debug_currn.grad_pn(:, k_min:k_max);
lorentz_new = new.debug_currn.lorentz(:, k_min:k_max);
terms_new = new.debug_currn.terms(:, :, k_min:k_max);

icol = repmat([1:ndim], [2, 1]);
irow = icol;
irow(2, :) = shift(irow(2, :), 1);
K_old = vertcat(d_old, shift(du_old, 1));
K_old = sparse(irow(:), icol(:), K_old(:), ndim, ndim);
K_new = vertcat(d_new, shift(du_new, 1));
K_new = sparse(irow(:), icol(:), K_new(:), ndim, ndim);

figure;
subplot(1, 2, 1);
hold on;
semilogy(abs(x_old), '-k');
semilogy(abs(x_new), ':r');
semilogy(abs(I_f_new .* coeff_f), '--b');
hold off;
legend('old', 'new', 'I_f');
ylabel('abs inhom / statA');
subplot(1, 2, 2);
hold on;
plot(arg(x_old), '-k');
plot(arg(x_new), ':r');
plot(arg(I_f_new .* coeff_f), '--b');
hold off;
legend('old', 'new', 'I_f');
ylabel('arg inhom / statA');

figure;
hold on;
subplot(1, 2, 1);
hold on;
plot(imag(d_new), ':r');
plot(imag(d_old), '-k');
hold off;
legend('old', 'new');
ylabel('Im d');
subplot(1, 2, 2);
hold on;
plot(imag(du_new), ':r');
plot(imag(du_old), '-k');
hold off;
legend('old', 'new');
ylabel('Im du');

figure;
subplot(1, 2, 1);
hold on;
semilogy(abs(I_old), '-k');
semilogy(abs(I_new), ':r');
semilogy(abs(K_old \ x_new.'), '--b');
semilogy(abs(K_new \ x_old.'), '-.g');
hold off;
legend('old', 'new', 'old K', 'new K');
ylabel('abs I_{pol} / statA');
subplot(1, 2, 2);
hold on;
plot(arg(I_old), '-k');
plot(arg(I_new), ':r');
plot(arg(K_old \ x_new.'), '--b');
plot(arg(K_new \ x_old.'), '-.g');
hold off;
legend('old', 'new', 'old K', 'new K');
ylabel('arg I_{pol} / statA');

sum_terms_old = sum(terms_old, 2);
sum_terms_new = sum(terms_new, 2);
figure;
for k = 1:3
  subplot(1, 3, k);
  hold on;
  semilogy(abs(sum_terms_old(k, :)), '-k');
  semilogy(abs(sum_terms_new(k, :)), ':r');
  hold off;
  legend('old', 'new');
  ylabel(['abs term ', num2str(k)]);
end
figure;
for k = 1:3
  subplot(1, 3, k);
  hold on;
  plot(arg(sum_terms_old(k, :)), '-k');
  plot(arg(sum_terms_new(k, :)), ':r');
  hold off;
  legend('old', 'new');
  ylabel(['arg term ', num2str(k)]);
end

coords = {'R', 'phi', 'Z'};
for k = 1:3
  figure;
  subplot(1, 2, 1);
  hold on;
  semilogy(abs(lorentz_old(k, :)), '-k');
  semilogy(abs(lorentz_new(k, :)), ':r');
  semilogy(abs(grad_pn_old(k, :)), '--b');
  hold off;
  legend('old', 'new', 'grad pn');
  ylabel(['abs f_{', coords{k}, '} / dyn cm^{-2}']);
  subplot(1, 2, 2);
  hold on;
  plot(arg(lorentz_old(k, :)), '-k');
  plot(arg(lorentz_new(k, :)), ':r');
  plot(arg(grad_pn_old(k, :)), '--b');
  hold off;
  legend('old', 'new', 'grad pn');
  ylabel(['arg f_{', coords{k}, '} / dyn cm^{-2}']);
end
