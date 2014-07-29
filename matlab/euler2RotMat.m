function [ m ] = euler2RotMat( x_angle, y_angle, z_angle )

% http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToMatrix/index.htm
c1 = cos(x_angle);  % Heading
s1 = sin(x_angle);
c2 = cos(y_angle);  % Attitude
s2 = sin(y_angle);
c3 = cos(z_angle);  % bank
s3 = sin(z_angle);

m = sym('m_',[4,4]);
% COLUMN_MAJOR
m(1) = c1*c2;
m(5) = -c1*s2*c3 + s1*s3;
m(9) = c1*s2*s3 + s1*c3;
m(13) = 0;
m(2) = s2;
m(6) = c2*c3;
m(10) = -c2*s3;
m(14) = 0;
m(3) = -s1*c2;
m(7) = s1*s2*c3 + c1*s3;
m(11) = -s1*s2*s3 + c1*c3;
m(15) = 0;
m(4) = 0;
m(8) = 0;
m(12) = 0;
m(16) = 1;
% % ROW_MAJOR
% m(1) = c1*c2;
% m(2) = -c1*s2*c3 + s1*s3;
% m(3) = c1*s2*s3 + s1*c3;
% m(4) = 0;
% m(5) = s2;
% m(6) = c2*c3;
% m(7) = -c2*s3;
% m(8) = 0;
% m(9) = -s1*c2;
% m(10) = s1*s2*c3 + c1*s3;
% m(11) = -s1*s2*s3 + c1*c3;
% m(12) = 0;
% m(13) = 0;
% m(14) = 0;
% m(15) = 0;
% m(16) = 1;

end

