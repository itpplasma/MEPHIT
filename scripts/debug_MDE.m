mde = load('../run/33353_2325/mephit.h5');
n = double(mde.mesh.n);

kf = 151;
% poloidal edge index = point index - 1
k_min = mde.mesh.kp_low(kf);
k_max = mde.mesh.kp_low(kf) + mde.mesh.kp_max(kf) - 1;

B0 = mde.cache.mid_fields.B0(k_min:k_max);
h = [mde.cache.mid_fields.B0_R(k_min:k_max); ...
     mde.cache.mid_fields.B0_phi(k_min:k_max); ...
     mde.cache.mid_fields.B0_Z(k_min:k_max)];
dp0_dpsi = mde.cache.fs.dp_dpsi(kf);
Bn_psi = mde.debug_presn.Bn_psi_contravar(k_min+1:k_max+1);
grad_pn_mde = mde.debug_presn.grad_pn(:, k_min+1:k_max+1);
lhs = dot(h, grad_pn_mde) ./ B0;
rhs = -Bn_psi ./ B0 * dp0_dpsi;

figure;
subplot(1, 2, 1);
hold on;
semilogy(abs(lhs), '-k');
semilogy(abs(rhs), '--r');
hold off;
legend('LHS', 'RHS');
ylabel('abs MDE p_{n} / dyn cm^{-3}');
subplot(1, 2, 2);
hold on;
plot(arg(lhs), '-k');
plot(arg(rhs), '--r');
hold off;
legend('LHS', 'RHS');
ylabel('arg MDE p_{n} / dyn cm^{-3}');

figure;
subplot(1, 2, 1);
hold on;
semilogy(abs(grad_pn_mde(1, :)), '-k');
semilogy(abs(grad_pn_mde(2, :)), ':b');
semilogy(abs(grad_pn_mde(3, :)), '--r');
hold off;
legend('R', 'phi', 'Z');
ylabel('abs grad p_{n} / dyn cm^{-3}');
subplot(1, 2, 2);
hold on;
plot(arg(grad_pn_mde(1, :)), '-k');
plot(arg(grad_pn_mde(2, :)), ':b');
plot(arg(grad_pn_mde(3, :)), '--r');
hold off;
legend('R', 'phi', 'Z');
ylabel('arg grad p_{n} / dyn cm^{-3}');

