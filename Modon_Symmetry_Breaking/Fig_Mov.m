% Creates frames for movie files

clear; close all

system('rm frames/*');

n = 5;
T_start = 200;
%ext = 'png';
%N = 256;
f_exp = 0.6;

Filt = @(f) filter2(ones(n)/n^2, squeeze(f), 'same');
eval_A = @(f) (f + f(:, end:-1:1, :))/2;
eval_S = @(f) (f - f(:, end:-1:1, :))/2;
f_scale = @(f) sign(f) .* abs(f).^f_exp;

filename = 'Test_data/Case_2L1:1_window.nc';
framename = 'frames/frame_';

T = ncread(filename, 't');
i0 = sum(T<T_start) + 1;

Q1 = ncread(filename, 'Q_1');
Q2 = ncread(filename, 'Q_2');
Y = ncread(filename, 'y')';

% Q1_A

F = f_scale(eval_A(Q1));
for i = 1:length(T)
    F(:, :, i) = Filt(F(:, :, i));
end
%Fm = max(max(max(abs(F)))); scale = [-Fm Fm];
save_frames(F(2:end, 2:end, i0:end), [framename 'Q1_A']); %, ext, cmap2([-1 1], 0, @(x) x.^f_exp, N));
%add_text2frame([framename 'Q1_A'], 'T = ', T(i0:end), [10 10], [framename 'Q1_A'], 4)

% Q1_S

F = f_scale(eval_S(Q1));
for i = 1:length(T)
    F(:, :, i) = Filt(F(:, :, i));
end
%Fm = max(max(max(abs(F)))); scale = [-Fm Fm];
save_frames(F(2:end, 2:end, i0:end), [framename 'Q1_S']); %, ext, cmap2([-1 1], 0, @(x) x.^f_exp, N));
%add_text2frame([framename 'Q1_S'], 'T = ', T(i0:end), [10 10], [framename 'Q1_S'], 4)

% Q1 + Y

F = f_scale(Q1 + Y);
for i = 1:length(T)
    F(:, :, i) = Filt(F(:, :, i));
end
%Fm = max(max(max(abs(F)))); scale = [-Fm Fm];
save_frames(F(2:end, 2:end, i0:end), [framename 'Q1_Y']); %, ext, cmap2([-1 1], 0, @(x) x.^f_exp, N));
%add_text2frame([framename 'Q1_Y'], 'T = ', T(i0:end), [10 10], [framename 'Q1_Y'], 4)

% Q2

F = f_scale(Q2);
for i = 1:length(T)
    F(:, :, i) = Filt(F(:, :, i));
end
%Fm = max(max(max(abs(F)))); scale = [-Fm Fm];
save_frames(F(2:end, 2:end, i0:end), [framename 'Q2_T']); %, ext, cmap2([-1 1], 0, @(x) x.^f_exp, N));
%add_text2frame([framename 'Q2_T'], 'T = ', T(i0:end), [10 10], [framename 'Q2_T'], 4)

% Q2_A

F = f_scale(eval_A(Q2));
for i = 1:length(T)
    F(:, :, i) = Filt(F(:, :, i));
end
%Fm = max(max(max(abs(F)))); scale = [-Fm Fm];
save_frames(F(2:end, 2:end, i0:end), [framename 'Q2_A']); %, ext, cmap2([-1 1], 0, @(x) x.^f_exp, N));
%add_text2frame([framename 'Q2_A'], 'T = ', T(i0:end), [10 10], [framename 'Q2_A'], 4)

close all
