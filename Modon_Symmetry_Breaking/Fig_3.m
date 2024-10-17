% creates plot of maximum Q and position of maximum Q

close all; clear

Figs_path = 'Figs/';

window_1 = [-6 2 -2 2];
sz = [600 300];
n = 7;
f_exp = 0.6;
Nc = 256;

it1 = 1367;
it2 = 2496;

Q_1 = ncread('Test_data/Case_2L1:1_window.nc', 'Q_1');
Q_2 = ncread('Test_data/Case_2L1:1_window.nc', 'Q_2');

Asymm = @(f) (f + f(:, end:-1:1, :)) / 2;
Symm = @(f) (f - f(:, end:-1:1, :)) / 2;
Filt = @(f) filter2(ones(n)/n^2, squeeze(f), 'same');

x = linspace(-6, 2, 801)';
y = linspace(-2, 2, 401)';

% [-15 -10 -0.01 0.01 10 15]
f = Filt(Q_1(:, :, it1) + y');
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-20 20])
contour(x, y, f', [-10 -1 -0.1 0.1 1 10], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_a'])

f = Filt(Q_1(:, :, it2) + y');
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-20 20])
contour(x, y, f', [-10 -1 -0.1 0.1 1 10], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_b'])

f = Filt(Symm(Q_1(:, :, it1)));
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-20 20])
contour(x, y, f', [-15 -2 2 15], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_c'])

f = Filt(Symm(Q_1(:, :, it2)));
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-20 20])
contour(x, y, f', [-15 -2 2 15], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_d'])

f = Filt(Asymm(Q_1(:, :, it1)));
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-4.3 4.3])
contour(x, y, f', [-1 1], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_e'])

f = Filt(Asymm(Q_1(:, :, it2)));
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-2.5 2.5])
contour(x, y, f', [-1 1], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_f'])

f = Q_2(:, :, it1);
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-0.3 0.3])
contour(x, y, f', [-0.2 -0.1 0.1 0.2], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_g'])

f = Q_2(:, :, it2);
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-0.7 0.7])
contour(x, y, f', [-0.2 -0.1 0.1 0.2], 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_h'])

f = Filt(Asymm(Q_2(:, :, it1)));
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-0.33 0.33])
%contour(x, y, f', 3, 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_i'])

f = Filt(Asymm(Q_2(:, :, it2)));
splot2(f, x, y, 'X', 'Y', 0, 1, sz)
hold on
clim([-3e-3 3e-3])
%contour(x, y, f', 3, 'k')
colormap(cmap2(0, 0, @(x) x.^f_exp, Nc))
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_3_j'])


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