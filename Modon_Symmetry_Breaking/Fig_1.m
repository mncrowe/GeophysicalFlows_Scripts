% Creates Fig 1 according to Georgi's suggestions

clear; close all

Q_1 = ncread('Test_data/Case_2L1:1_all.nc', 'Q_1');
Q_2 = ncread('Test_data/Case_2L1:1_all.nc', 'Q_2');

psi_1 = ncread('Test_data/Case_2L1:1_all.nc', 'psi_1');
psi_2 = ncread('Test_data/Case_2L1:1_all.nc', 'psi_2');

x = ncread('Test_data/Case_2L1:1_all.nc', 'x');
y = ncread('Test_data/Case_2L1:1_all.nc', 'y');

PV_1 = Q_1 + y';
PV_2 = Q_2 + y';

window_1 = [-3 3 -3 3];

Figs_path = 'Figs/';

d_dy = @(f, y) d_dx(f, y, 2, 2);

[ix, iy] = cropped_index(x, y, window_1);

u0_1 = -d_dy(psi_1(:,:,1), y);
u0_2 = -d_dy(psi_2(:,:,1), y);

splot2(u0_1(ix,iy), x(ix), y(iy), 'X', 'Y')
hold on
contour(x(ix), y(iy), u0_1(ix,iy)', 5, 'k')
colormap(cmap2(0, 0))
plot_circle(1)
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_1_a'])

splot2(PV_1(ix,iy,1), x(ix), y(iy), 'X', 'Y')
hold on
contour(x(ix), y(iy), PV_1(ix,iy,1)', 5, 'k')
colormap(cmap2(0, 0))
plot_circle(1)
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_1_b'])

splot2(u0_2(ix,iy), x(ix), y(iy), 'X', 'Y')
hold on
contour(x(ix), y(iy), u0_2(ix,iy)', 5, 'k')
colormap(cmap2(0, 0))
plot_circle(1)
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_1_c'])

splot2(PV_2(ix,iy,1), x(ix), y(iy), 'X', 'Y')
hold on
contour(x(ix), y(iy), PV_2(ix,iy,1)', 5, 'k')
colormap(cmap2(0, 0))
plot_circle(1)
plot_window(window_1)
hold off
save_and_close([Figs_path 'Fig_1_d'])


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

function plot_circle(a)

    fplot(@(t) a*sin(t), @(t) a*cos(t), [0,2*pi], 'k--')

end

function save_and_close(savename)

    %save_figure(savename, 'eps', 1)
    save_figure(savename, 'png', 1)
    
    %savefig(savename)

    close(1)

end