data = load('../run/33353_2325/mephit.h5');
n = double(data.mesh.n);

kf = 30;
% poloidal edge index = point index - 1
k_min = data.mesh.kp_low(kf);
k_max = data.mesh.kp_low(kf) + data.mesh.kp_max(kf) - 1;

B0 = data.cache.mid_fields.B0(k_min:k_max);
h = [data.cache.mid_fields.B0_R(k_min:k_max) ./ B0; ...
     data.cache.mid_fields.B0_phi(k_min:k_max) ./ B0; ...
     data.cache.mid_fields.B0_Z(k_min:k_max) ./ B0];
dp0_dpsi = data.cache.fs.dp_dpsi(kf);
Bn_psi = data.debug_MDE.Bn_psi_contravar(k_min:k_max);
grad_pn = data.debug_MDE.grad_pn(:, k_min:k_max);
R = data.mesh.cntr_R;
Z = data.mesh.cntr_Z;
B_0 = [data.cache.cntr_fields.B0_R; ...
       data.cache.cntr_fields.B0_phi; ...
       data.cache.cntr_fields.B0_Z];
j_0 = [data.cache.cntr_fields.j0_R; ...
       data.cache.cntr_fields.j0_phi; ...
       data.cache.cntr_fields.j0_Z];
B_n = [data.vac.Bn.comp_R; data.vac.Bn.RT0_comp_phi; data.vac.Bn.comp_Z];
jn_perp = (2.99792458e10 * cross(B_0, data.debug_currn.grad_pn) ...
           + repmat(dot(j_0, B_0), [3, 1]) .* B_n ...
           - repmat(dot(B_n, B_0), [3, 1]) .* j_0) ...
           ./ repmat(dot(B_0, B_0), [3, 1]);
grid_R = linspace(min(data.mesh.node_R), max(data.mesh.node_R), 200);
grid_Z = linspace(min(data.mesh.node_Z), max(data.mesh.node_Z), 350);
[RR, ZZ] = meshgrid(grid_R, grid_Z);
jnperp_R = griddata(R, Z, jn_perp(1, :), RR, ZZ, 'linear');
jnperp_phi = griddata(R, Z, jn_perp(2, :), RR, ZZ, 'linear');
jnperp_Z = griddata(R, Z, jn_perp(3, :), RR, ZZ, 'linear');
[djnperpR_dR, ~] = gradient(jnperp_R, grid_R, grid_Z);
[~, djnperpZ_dZ] = gradient(jnperp_Z, grid_R, grid_Z);
div_jnperp_oct = interp2(RR, ZZ, djnperpR_dR + jnperp_R ./ RR + djnperpZ_dZ + ...
                 1i * n * jnperp_phi, data.mesh.mid_R, data.mesh.mid_Z);
div_jnperp = data.debug_MDE.div_jnperp(k_min:k_max);
div_jnperp_RT0 = data.debug_MDE.div_jnperp_RT0(k_min:k_max);
grad_jnpar = data.debug_MDE.grad_jnpar(:, k_min:k_max);
inhom = data.debug_currn.x(k_min+1:k_max+1);

lhs = dot(h, grad_pn);
rhs = -Bn_psi ./ B0 * dp0_dpsi;

figure;
subplot(1, 2, 1);
semilogy(abs(lhs), '-k');
hold on;
semilogy(abs(rhs), '--r');
hold off;
legend('LHS', 'RHS');
ylabel('abs MDE p_{n} / dyn cm^{-3}');
subplot(1, 2, 2);
plot(angle(lhs), '-k');
hold on;
plot(angle(rhs), '--r');
hold off;
legend('LHS', 'RHS');
ylabel('arg MDE p_{n} / dyn cm^{-3}');

figure;
subplot(1, 2, 1);
semilogy(abs(grad_pn(1, :)), '-k');
hold on;
semilogy(abs(grad_pn(2, :)), ':b');
semilogy(abs(grad_pn(3, :)), '--r');
hold off;
legend('R', 'phi', 'Z');
ylabel('abs grad p_{n} / dyn cm^{-3}');
subplot(1, 2, 2);
plot(arg(grad_pn(1, :)), '-k');
hold on;
plot(arg(grad_pn(2, :)), ':b');
plot(arg(grad_pn(3, :)), '--r');
hold off;
legend('R', 'phi', 'Z');
ylabel('arg grad p_{n} / dyn cm^{-3}');

lhs = B0 .* dot(h, grad_jnpar);
rhs = inhom;
rhs_var = -div_jnperp;
rhs_oct = -div_jnperp_oct(k_min+1:k_max+1);
rhs_RT0 = -div_jnperp_RT0;

figure;
subplot(1, 2, 1);
semilogy(abs(lhs), '-k');
hold on;
semilogy(abs(rhs), '--r');
semilogy(abs(rhs_var), ':b');
semilogy(abs(rhs_RT0), '-.g');
semilogy(abs(rhs_oct), '--y');
hold off;
legend('LHS', 'RHS', 'RHS (stencil)', 'RHS (RT0)', 'RHS (oct)');
ylabel('abs MDE j_{n}^{||} / statA cm^{-3}');
subplot(1, 2, 2);
plot(angle(lhs), '-k');
hold on;
plot(angle(rhs), '--r');
plot(angle(rhs_var), ':b');
plot(angle(rhs_RT0), '-.g');
plot(angle(rhs_oct), '--y');
hold off;
legend('LHS', 'RHS', 'RHS (stencil)', 'RHS (RT0)', 'RHS (oct)');
ylabel('arg MDE j_{n}^{||} / statA cm^{-3}');

figure;
subplot(1, 2, 1);
semilogy(abs(grad_jnpar(1, :)), '-k');
hold on;
semilogy(abs(grad_jnpar(2, :)), ':b');
semilogy(abs(grad_jnpar(3, :)), '--r');
hold off;
legend('R', 'phi', 'Z');
ylabel('abs grad j_{n}^{||} / statA cm^{-3}');
subplot(1, 2, 2);
plot(angle(grad_jnpar(1, :)), '-k');
hold on;
plot(angle(grad_jnpar(2, :)), ':b');
plot(angle(grad_jnpar(3, :)), '--r');
hold off;
legend('R', 'phi', 'Z');
ylabel('arg grad j_{n}^{||} / statA cm^{-3}');

