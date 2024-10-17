% fig 4; f-plane analysis

close all; clear

Figs_path = 'Figs/';

x = ncread('Test_data/Case_1L0_all.nc', 'x');
y = ncread('Test_data/Case_1L0_all.nc', 'y');
t = ncread('Test_data/Case_1L0_all.nc', 't');

eval_A = @(f) (f + f(:, [1 end:-1:2], :))/2;
eval_S = @(f) (f - f(:, [1 end:-1:2], :))/2;
max_A = @(f) squeeze(max(max(abs(eval_A(f)))));
d_dy = @(f) d_dx(f, y, 2, 2);

lw = 2;
it = 501;
window_1 = [-2 2 -2 2];
n = 7;

Filt = @(f) filter2(ones(n)/n^2, squeeze(f), 'same');

[ix, iy] = cropped_index(x, y, window_1);

options = optimset('TolFun',1e-10, 'TolX', 1e-10);

Q_1 = ncread('Test_data/Case_1L0_all.nc','Q'); Q_1_1_A = max_A(Q_1);
[Q_c_1, x_c_1, y_c_1] = interp_max2(Q_1, x, y);
Q_S_1 = eval_S(Q_1(:,:,1)); dQ_SdY_1 = max(max(d_dy(Q_S_1)));
DeltaY_1 = Q_1_1_A / dQ_SdY_1;
Max_res_1 = t;
for i = 1:length(t)
    Q1_A = eval_A(Q_1(:, :, i)); dQ1_SdY = d_dy(eval_S(Q_1(:,:,i)));
    f1 = @(x) sum(sum((Q1_A + x * dQ1_SdY).^2));
    d1 = fminsearch(f1, 2e-7, options);
    Max_res_1(i) = max(max(abs(Q1_A + d1 * dQ1_SdY)));
end

Q_1 = ncread('Test_data/Case_2L0:0_all.nc','Q_1'); Q_1_3_A = max_A(Q_1);
[Q_c_3, x_c_3, y_c_3] = interp_max2(Q_1, x, y);
Q_S_3 = eval_S(Q_1(:,:,1)); dQ_SdY_3 = max(max(d_dy(Q_S_3)));
DeltaY_3 = Q_1_3_A / dQ_SdY_3;
Max_res_3 = t;
for i = 1:length(t)
    Q1_A = eval_A(Q_1(:, :, i)); dQ1_SdY = d_dy(eval_S(Q_1(:,:,i)));
    f1 = @(x) sum(sum((Q1_A + x * dQ1_SdY).^2));
    d1 = fminsearch(f1, 2e-7, options);
    Max_res_3(i) = max(max(abs(Q1_A + d1 * dQ1_SdY)));
end

Q_1 = ncread('Test_data/Case_2L0:1_all.nc','Q_1'); Q_1_5_A = max_A(Q_1);
Q_2 = ncread('Test_data/Case_2L0:1_all.nc','Q_2'); Q_2_5_A = max_A(Q_2);
[Q1_c_5, x1_c_5, y1_c_5] = interp_max2(Q_1, x, y);
[Q2_c_5, x2_c_5, y2_c_5] = interp_max2(Q_2, x, y);
Q1_S_5 = eval_S(Q_1(:,:,1)); dQ1_SdY_5 = max(max(d_dy(Q1_S_5)));
Q2_S_5 = eval_S(Q_2(:,:,1)); dQ2_SdY_5 = max(max(d_dy(Q2_S_5)));
DeltaY1_5 = Q_1_5_A / dQ1_SdY_5;
DeltaY2_5 = Q_2_5_A / dQ2_SdY_5;
Max_res_51 = t; Max_res_52 = t;
for i = 1:length(t)
    Q1_A = eval_A(Q_1(:, :, i)); dQ1_SdY = d_dy(eval_S(Q_1(:,:,i)));
    Q2_A = eval_A(Q_2(:, :, i)); dQ2_SdY = d_dy(eval_S(Q_2(:,:,i)));
    f1 = @(x) sum(sum((Q1_A + x * dQ1_SdY).^2));
    f2 = @(x) sum(sum((Q2_A + x * dQ2_SdY).^2));
    d1 = fminsearch(f1, 2e-7, options);
    d2 = fminsearch(f2, 2e-7, options);
    Max_res_51(i) = max(max(abs(Q1_A + d1 * dQ1_SdY)));
    Max_res_52(i) = max(max(abs(Q2_A + d2 * dQ2_SdY)));
end

Q1_A = eval_A(Q_1(:, :, it)); dQ1_SdY = d_dy(eval_S(Q_1(:,:,it)));
Q2_A = eval_A(Q_2(:, :, it)); dQ2_SdY = d_dy(eval_S(Q_2(:,:,it)));

f1 = @(x) sum(sum((Q1_A + x * dQ1_SdY).^2));
f2 = @(x) sum(sum((Q2_A + x * dQ2_SdY).^2));

DeltaY1 = fminsearch(f1, 2e-7, options);
DeltaY2 = fminsearch(f2, 2e-7, options);

%semilogy(t, Q_1_1_A, t, Q_1_3_A, t, Q_1_5_A, t, Q_2_5_A, 'LineWidth', lw); grid
%xlabel('T'); ylabel('Max|Q_A|'); xlim([t(1) t(end)])
%legend('1L0', '2L0:0', '2L0:1 (layer 1)', '2L0:1 (layer 2)', 'location', 'NorthWest')
%set(gca,'FontSize',12,'linewidth',0.7);
%save_and_close([Figs_path 'Fig_4_a'])

figure;
semilogy(t, DeltaY_1, t, DeltaY_3, t, DeltaY1_5, t, DeltaY2_5, 'LineWidth', lw); grid
xlabel('T'); ylabel('\DeltaY_j'); xlim([t(1) t(end)])
legend('1L0', '2L0:0', '2L0:1 (layer 1)', '2L0:1 (layer 2)', 'location', 'NorthWest')
set(gca,'FontSize',12,'linewidth',0.7);
save_and_close([Figs_path 'Fig_4_a'])

figure;
semilogy(t, Max_res_1./DeltaY_1, t, Max_res_3./DeltaY_3, t, Max_res_51./DeltaY1_5, t, Max_res_52./DeltaY2_5, 'LineWidth', lw); grid
xlabel('T'); ylabel('Max|\DeltaQ_j|/\DeltaY_j'); xlim([t(1) t(end)])
legend('1L0', '2L0:0', '2L0:1 (layer 1)', '2L0:1 (layer 2)', 'location', 'NorthWest')
set(gca,'FontSize',12,'linewidth',0.7);
save_and_close([Figs_path 'Fig_4_b'])

f = Filt(Q1_A(ix-67,iy));
splot2(f, x(ix), y(iy), 'X', 'Y')
hold on
contour(x(ix), y(iy), f', [-5e-5 -3e-5 -2e-6 2e-6 3e-5 5e-5], 'k')
colormap(cmap2(0, 0))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_4_c'])

f = Filt(Q1_A(ix-67,iy) + DeltaY1 * dQ1_SdY(ix-67,iy));
splot2(f, x(ix), y(iy), 'X', 'Y')
hold on
contour(x(ix), y(iy), f', [-3e-8 -2e-8 -3e-9 3e-9 2e-8 3e-8], 'k')
colormap(cmap2(0, 0))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_4_d'])

% f = Q2_A(ix,iy) + DeltaY2 * dQ2_SdY(ix,iy);
% splot2(f, x(ix), y(iy), 'X', 'Y')
% hold on
% contour(x(ix), y(iy), f', [-3e-9 -2e-9 -1e-9 0 1e-9], 'k')
% colormap(cmap2(0, 0, 0, 0, 0))
% plot_window(window_1)
% hold off
% save_and_close([Figs_path 'Fig_4_d'])


function [ix, iy] = cropped_index(x, y, window)

    Nx = length(x); Ny = length(y);

    ix1 = sum(x < window(1)); ix2 = 1 + Nx - sum(x > window(2));
    iy1 = sum(y < window(3)); iy2 = 1 + Ny - sum(y > window(4));

    ix = ix1:ix2;
    iy = iy1:iy2;

end

function plot_window(lims)

    axis equal

    xlim(lims(1:2))
    ylim(lims(3:4))

    line([lims(1) lims(1)], [lims(3) lims(4)], 'Color', 'black')
    line([lims(2) lims(2)], [lims(3) lims(4)], 'Color', 'black')
    line([lims(1) lims(2)], [lims(3) lims(3)], 'Color', 'black')
    line([lims(1) lims(2)], [lims(4) lims(4)], 'Color', 'black')

    set(gca,'FontSize',12,'linewidth',0.7);

end


function save_and_close(savename)

    %save_figure(savename, 'eps', 1)
    save_figure(savename, 'png', 1)
    
    %savefig(savename)

    close(1)

end