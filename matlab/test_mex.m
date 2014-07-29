clearvars; clc; close all;
load 'bunny_data.mat'

% % Use the same matrix as the c++ code
M_PC2 = [ +8.38535964e-001  -7.75283277e-002  -8.10404643e-002  +7.44949058e-002 ; ...
          +7.93413445e-002  +8.78419638e-001  -1.93956606e-002  -2.78396234e-002 ; ...
          +1.03379093e-001  +1.39857326e-002  +1.05629706e+000  -3.56172025e-003 ; ...
          +0.00000000e+000  +0.00000000e+000  +0.00000000e+000  +1.00000000e+000 ];
% Just translation and rotation
% M_PC2 = [ +9.91173983e-001  -9.16407481e-002  -9.57921892e-002  +7.44949058e-002 ; ...
%           +8.99348930e-002  +9.95705009e-001  -2.19853427e-002  -2.78396234e-002 ; ...
%           +9.73955095e-002  +1.31762382e-002  +9.95158553e-001  -3.56172025e-003 ; ...
%           +0.00000000e+000  +0.00000000e+000  +0.00000000e+000  +1.00000000e+000 ];
M_PC2_normal_mat = inv(M_PC2)';

% Affine Transform position
pc2 = (M_PC2(1:3,1:3) * bunny_data.vert' + repmat(M_PC2(1:3,4), 1, size(bunny_data.vert,1)))';
% Affine Transform normal
npc2 = M_PC2(1:3,1:3) * bunny_data.norm';
len = sqrt(sum(abs(npc2).^2,1));
npc2 = (npc2 ./ repmat(len, 3, 1))';

icp_data = struct;
icp_data.num_iterations = 29;
icp_data.pc1 = bunny_data.vert';
icp_data.npc1 = bunny_data.norm';  % Optional parameter
icp_data.pc2 = pc2';
icp_data.npc2 = npc2';  % Optional parameter
icp_data.m_pc2_initial = eye(4);
icp_data.min_distance_sq = 1e-6;  % Optional parameter
icp_data.max_distance_sq = 1e+3;  % Optional parameter
icp_data.cos_normal_threshold = 0.9;  % Optional parameter
icp_data.method = 1;  % Optional parameter

M_PC2_icp = icp(icp_data);

% Affine Transform position
pc2_icp = (M_PC2_icp(1:3,1:3) * pc2' + repmat(M_PC2_icp(1:3,4), 1, size(bunny_data.vert,1)))';

figure;
scatter3(bunny_data.vert(:,1), bunny_data.vert(:,3), bunny_data.vert(:,2),'r.'); hold on;
scatter3(pc2(:,1), pc2(:,3), pc2(:,2),'b.');
scatter3(pc2_icp(:,1), pc2_icp(:,3), pc2_icp(:,2),'g.');
% view(0,90);
set(gcf,'renderer','opengl'); axis vis3d;

mean_err_before = mean(sum((pc2 - bunny_data.vert).^2, 2));
mean_err_after = mean(sum((pc2_icp - bunny_data.vert).^2, 2));

display(['Mean error before ICP: ', num2str(mean_err_before)]);
display(['Mean error after ICP: ', num2str(mean_err_after)]);