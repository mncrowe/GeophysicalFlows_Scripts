% Creates Fig 2 for paper with Georgi

close all; clear

Figs_path = 'Figs/';

Lx = 20.48;
lw = 2;
T_range = [0, 500];
log_range = [1e-9, 1e1];

x = ncread('Test_data/Case_1L0_all.nc', 'x');
y = ncread('Test_data/Case_1L0_all.nc', 'y');

t = ncread('Test_data/Case_1L0_all.nc', 't'); dt = t(2) - t(1);
t2 = (t(2:end) + t(1:end-1))/2;

% Calculate max Q_A for each case:

eval_A = @(f) (f + f(:, [1 end:-1:2], :))/2;
eval_A2 = @(f) (f + f(:, end:-1:1, :))/2;
max_A = @(f) squeeze(max(max(abs(eval_A(f)))));

Q_1 = ncread('Test_data/Case_1L1_all.nc','Q'); Q_1_2_A = max_A(Q_1);
[Q_c_2, x_c_2, y_c_2] = interp_max2(Q_1, x, y);
x_c_2 = unloop(x_c_2, Lx);
x_c_2 = smooth(x_c_2); y_c_2 = smooth(y_c_2);
dXdt_2 = smooth(diff(x_c_2)) + 1;

Q_1 = ncread('Test_data/Case_2L1:0_all.nc','Q_1'); Q_1_4_A = max_A(Q_1);
[Q_c_4, x_c_4, y_c_4] = interp_max2(Q_1, x, y);
x_c_4 = unloop(x_c_4, Lx);
x_c_4 = smooth(x_c_4); y_c_4 = smooth(y_c_4);
dXdt_4 = smooth(diff(x_c_4)) + 1;

Q_1 = ncread('Test_data/Case_2L1:1_all.nc','Q_1'); Q_1_6_A = max_A(Q_1);
Q_2 = ncread('Test_data/Case_2L1:1_all.nc','Q_2'); Q_2_6_A = max_A(Q_2);
[Q_c_6, x_c_6, y_c_6] = interp_max2(Q_1, x, y);
x_c_6 = unloop(x_c_6, Lx);
x_c_6 = smooth(x_c_6); y_c_6 = smooth(y_c_6);
dXdt_6 = smooth(diff(x_c_6)) + 1;

x2 = linspace(-6, 2, 801)';
y2 = linspace(-2, 2, 401)';
%mask = (abs(atan2(y2', x2)) > 7*pi/8) .* ((x2.^2 + y2'.^2) > 4);

Q_1 = ncread('Test_data/Case_1L1_window.nc', 'Q');
Q_1A = eval_A2(Q_1);
Q_2_w = smooth(max(squeeze(Q_1A(1, :, 1:5:end))));
[~, i_m] = max(Q_2_w); t_c_2 = t(i_m);

Q_1 = ncread('Test_data/Case_2L1:0_window.nc', 'Q_1');
Q_1A = eval_A2(Q_1);
Q_4_w = smooth(max(squeeze(Q_1A(1, :, 1:5:end))));
[~, i_m] = max(Q_4_w); t_c_4 = t(i_m);

Q_1 = ncread('Test_data/Case_2L1:1_window.nc', 'Q_1');
Q_1A = eval_A2(Q_1);
Q_6_w = smooth(max(squeeze(Q_1A(1, :, 1:5:end))));
[~, i_m] = max(Q_6_w); t_c_6 = t(i_m);


% fig 2a

figure;
semilogy(t, Q_1_2_A, t, Q_1_4_A, t, Q_1_6_A, t, Q_2_6_A, 'LineWidth', lw); grid
xlabel('T'); ylabel('Max|Q_A|'); xlim(T_range); ylim(log_range);
xline([t_c_2, t_c_4, t_c_6])
legend('1L1', '2L1:0', '2L1:1 (layer 1)', '2L1:1 (layer 2)', 'location', 'NorthWest')
set(gca,'FontSize',12,'linewidth',0.7);
save_and_close([Figs_path 'Fig_2_a'])

% fig 2b

figure;
semilogy(t, Q_2_w, t, Q_4_w, t, Q_6_w, 'LineWidth', lw); grid
xlabel('T'); ylabel('Max|Q_A(X = X_c - 6)|'); xlim(T_range); ylim(log_range);
xline([t_c_2, t_c_4, t_c_6])
legend('1L1', '2L1:0', '2L1:1', 'location', 'NorthWest')
set(gca,'FontSize',12,'linewidth',0.7);
save_and_close([Figs_path 'Fig_2_b'])

% fig 2c

figure;
plot(t, y_c_2, t, y_c_4, t, y_c_6, 'LineWidth', lw); grid
xlabel('T'); ylabel('Y_c'); xlim(T_range)
xline([t_c_2, t_c_4, t_c_6])
legend('1L1', '2L1:0', '2L1:1', 'location', 'NorthWest')
set(gca,'FontSize',12,'linewidth',0.7);
save_and_close([Figs_path 'Fig_2_c'])

% fig 2d

figure;
plot(t2, dXdt_2, t2, dXdt_4, t2, dXdt_6, 'LineWidth', lw); grid
xlabel('T'); ylabel('dX_c/dT'); xlim(T_range)
xline([t_c_2, t_c_4, t_c_6])
legend('1L1', '2L1:0', '2L1:1', 'location', 'SouthWest')
set(gca,'FontSize',12,'linewidth',0.7);
save_and_close([Figs_path 'Fig_2_d'])


function save_and_close(savename)

    %save_figure(savename, 'eps', 1)
    save_figure(savename, 'png', 1)
    
    %savefig(savename)

    close(1)

end

function x = unloop(x, Lx)

    % un-loops a periodic x

    N = length(x);

    for i = 2:N

        while abs(x(i) - x(i-1)) > Lx / 2

            if x(i) > x(i-1)
                x(i) = x(i) - Lx;
            else
                x(i) = x(i) + Lx;
            end

        end

    end

end
