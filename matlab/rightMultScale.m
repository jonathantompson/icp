function [ m ] = rightMultScale( m, x, y, z )
% % ROW COLUMN_MAJOR
% m(1) = m(1) * x;
% m(2) = m(2) * y;
% m(3) = m(3) * z;
% m(4) = m(4) * x;
% m(5) = m(5) * y;
% m(6) = m(6) * z;
% m(7) = m(7) * x;
% m(8) = m(8) * y;
% m(9) = m(9) * z;
% m(10) = m(10) * x;
% m(11) = m(11) * y;
% m(12) = m(12) * z;
% COLUMN_MAJOR
m(1) = m(1) * x;
m(2) = m(2) * x;
m(3) = m(3) * x;
m(4) = m(4) * x;
m(5) = m(5) * y;
m(6) = m(6) * y;
m(7) = m(7) * y;
m(8) = m(8) * y;
m(9) = m(9) * z;
m(10) = m(10) * z;
m(11) = m(11) * z;
m(12) = m(12) * z;

end

