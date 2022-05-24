old = load('../run/33353_2325/mephit_old.h5');
new = load('../run/33353_2325/mephit_new.h5');
mde = load('../run/33353_2325/mephit.h5');
n = double(old.mesh.n);

kf = 30;
k_min = old.mesh.kt_low(kf) + 1;
k_max = old.mesh.kt_low(kf) + old.mesh.kt_max(kf);
k_off = old.mesh.npoint - 1;
area = old.mesh.area(k_min:k_max);
orient = double(old.mesh.orient(k_min:k_max)) * 2 - 1;
e_f = old.mesh.tri_edge(1, k_min:k_max);

I_f_old = old.debug_currn.I(e_f) .* orient;
I_old = old.debug_currn.I((k_off+k_min):(k_off+k_max));
grad_pn_old = old.debug_currn.grad_pn(:, k_min:k_max);
lorentz_old = old.debug_currn.lorentz(:, k_min:k_max);
I_f_new = new.debug_currn.I(e_f) .* orient;
I_new = new.debug_currn.I((k_off+k_min):(k_off+k_max));
grad_pn_new = new.debug_currn.grad_pn(:, k_min:k_max);
lorentz_new = new.debug_currn.lorentz(:, k_min:k_max);
I_f_mde = mde.debug_currn.I(e_f) .* orient;
I_mde = mde.debug_currn.I((k_off+k_min):(k_off+k_max));
grad_pn_mde = mde.debug_currn.grad_pn(:, k_min:k_max);
lorentz_mde = mde.debug_currn.lorentz(:, k_min:k_max);

figure;
subplot(1, 2, 1);
hold on;
semilogy(abs(I_f_old), '-k');
semilogy(abs(I_f_new), ':r');
semilogy(abs(I_f_mde), '--b');
hold off;
legend('old', 'new', 'MDE');
ylabel('abs I_{f} / statA');
subplot(1, 2, 2);
hold on;
plot(arg(I_f_old), '-k');
plot(arg(I_f_new), ':r');
plot(arg(I_f_mde), '--b');
hold off;
legend('old', 'new', 'MDE');
ylabel('arg I_{f} / statA');

figure;
subplot(1, 2, 1);
hold on;
semilogy(abs(I_old), '-k');
semilogy(abs(I_new), ':r');
semilogy(abs(I_mde), '--b');
hold off;
legend('old', 'new', 'MDE');
ylabel('abs I_{pol} / statA');
subplot(1, 2, 2);
hold on;
plot(arg(I_old), '-k');
plot(arg(I_new), ':r');
plot(arg(I_mde), '--b');
hold off;
legend('old', 'new', 'MDE');
ylabel('arg I_{pol} / statA');

coords = {'R', 'phi', 'Z'};
for k = 1:3
  figure;
  subplot(1, 2, 1);
  hold on;
  semilogy(abs(lorentz_old(k, :)), '-k');
  semilogy(abs(lorentz_new(k, :)), ':r');
  semilogy(abs(lorentz_mde(k, :)), '--b');
  semilogy(abs(grad_pn_mde(k, :)), '-.g');
  hold off;
  legend('old', 'new', 'MDE', 'grad pn');
  ylabel(['abs f_{', coords{k}, '} / dyn cm^{-2}']);
  subplot(1, 2, 2);
  hold on;
  plot(arg(lorentz_old(k, :)), '-k');
  plot(arg(lorentz_new(k, :)), ':r');
  plot(arg(lorentz_mde(k, :)), '--b');
  plot(arg(grad_pn_mde(k, :)), '-.g');
  hold off;
  legend('old', 'new', 'MDE', 'grad pn');
  ylabel(['arg f_{', coords{k}, '} / dyn cm^{-2}']);
end