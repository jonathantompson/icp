clearvars; close all; clc;

dim = 3;
npts = 1;
pc1 = sym('pc1_', [dim+1, npts]);
pc2 = sym('pc2_', [dim+1, npts]);
weights = sym('weight', [1,npts]);
pc1(dim+1,:) = 1;
pc2(dim+1,:) = 1;
scale = sym('scale');
cur_translation_scale = sym('cur_translation_scale_');
PI = sym('M_PI');

c = sym('c', [1,9]);

mat = euler2RotMat(2*PI*c(4), 2*PI*c(5), 2*PI*c(6));
mat = rightMultScale(mat, c(7), c(8), c(9));
mat = leftMultTranslation(mat, c(1)*cur_translation_scale, ...
  c(2)*cur_translation_scale, c(3)*cur_translation_scale);

delta_x = pc1 - mat * pc2;
err_per_point = scale * weights .* sum(delta_x.^2,1);  % Dot product 1st dim
err = sum(sum(err_per_point));

for i = 1:9
  dy_dx = diff(err, c(i));
  code = ccode(dy_dx);
  code = strrep(code, 't0 = ', ['j[',num2str(i-1),'] += ']);
  disp(code);
end