old = load('../run/33353_2325/mephit_old.h5');
mde = load('../run/33353_2325/mephit.h5');
n = double(old.mesh.n);

kf = 30;
k_min = old.mesh.kt_low(kf) + 1;
k_max = old.mesh.kt_low(kf) + old.mesh.kt_max(kf);
k_off = old.mesh.npoint - 1;
area = old.mesh.area(k_min:k_max);
orient = double(old.mesh.orient(k_min:k_max)) * 2 - 1;
e_f = old.mesh.tri_edge(1, k_min:k_max);
% old HDF5 file had no prefix
I_f_old = old.debug_currn.I(e_f) .* orient;
I_old = old.debug_currn.I((k_off+k_min):(k_off+k_max));
grad_pn_old = old.debug_currn.grad_pn(:, k_min:k_max);
lorentz_old = old.debug_currn.lorentz(:, k_min:k_max);
I_f_mde = mde.debug_currn_000.I(e_f) .* orient;
I_mde = mde.debug_currn_000.I((k_off+k_min):(k_off+k_max));
grad_pn_mde = mde.debug_currn_000.grad_pn(:, k_min:k_max);
lorentz_mde = mde.debug_currn_000.lorentz(:, k_min:k_max);

figure;
subplot(1, 2, 1);
semilogy(abs(I_f_old), '-k');
hold on;
semilogy(abs(I_f_mde), ':r');
hold off;
legend('old', 'MDE');
ylabel('abs I_{f} / statA');
subplot(1, 2, 2);
plot(arg(I_f_old), '-k');
hold on;
plot(arg(I_f_mde), ':r');
hold off;
legend('old', 'MDE');
ylabel('arg I_{f} / rad');

figure;
subplot(1, 2, 1);
semilogy(abs(I_old), '-k');
hold on;
semilogy(abs(I_mde), ':r');
hold off;
legend('old', 'MDE');
ylabel('abs I_{pol} / statA');
subplot(1, 2, 2);
plot(arg(I_old), '-k');
hold on;
plot(arg(I_mde), ':r');
hold off;
legend('old', 'MDE');
ylabel('arg I_{pol} / rad');

coords = {'R', 'phi', 'Z'};
for k = 1:3
  figure;
  subplot(1, 2, 1);
  semilogy(abs(lorentz_old(k, :)), '-k');
  hold on;
  semilogy(abs(lorentz_mde(k, :)), ':r');
  semilogy(abs(grad_pn_mde(k, :)), '--b');
  hold off;
  legend('old', 'MDE', 'grad p_{n}');
  ylabel(['abs f_{', coords{k}, '} / dyn cm^{-2}']);
  subplot(1, 2, 2);
  plot(arg(lorentz_old(k, :)), '-k');
  hold on;
  plot(arg(lorentz_mde(k, :)), ':r');
  plot(arg(grad_pn_mde(k, :)), '--b');
  hold off;
  legend('old', 'MDE', 'grad p_{n}');
  ylabel(['arg f_{', coords{k}, '} / rad']);
end

k_min_in = mde.mesh.kp_low(kf - 1) + 1;
k_max_in = mde.mesh.kp_low(kf - 1) + mde.mesh.kp_max(kf - 1);
k_min_out = mde.mesh.kp_low(kf) + 1;
k_max_out = mde.mesh.kp_low(kf) + mde.mesh.kp_max(kf);
theta_in = mde.mesh.node_theta_geom(k_min_in:k_max_in);
pn_in = mde.iter.pn_000.L1_DOF(k_min_in:k_max_in);
theta_out = mde.mesh.node_theta_geom(k_min_out:k_max_out);
pn_out = mde.iter.pn_000.L1_DOF(k_min_out:k_max_out);
figure;
subplot(1, 2, 1);
plot(theta_in, real(pn_in), '-k');
hold on;
plot(theta_out, real(pn_out), ':r');
hold off;
legend('inner', 'outer');
xlabel('\theta / rad');
ylabel('Re p_{n} / dyn cm^{-2}');
subplot(1, 2, 2);
plot(theta_in, imag(pn_in), '-k');
hold on;
plot(theta_out, imag(pn_out), ':r');
hold off;
legend('inner', 'outer');
xlabel('\theta / rad');
ylabel('Im p_{n} / dyn cm^{-2}');

diff_re = interp1(theta_out, real(pn_out), theta_in, 'cubic') - real(pn_in);
diff_im = interp1(theta_out, imag(pn_out), theta_in, 'cubic') - imag(pn_in);
figure;
subplot(1, 2, 1);
plot(theta_in, hypot(diff_im, diff_re), '-k');
xlabel('\theta / rad');
ylabel('abs \Delta p_{n} / dyn cm^{-2}');
subplot(1, 2, 2);
plot(theta_in, atan2(diff_im, diff_re), '-k');
xlabel('\theta / rad');
ylabel('arg \Delta p_{n} / rad');
