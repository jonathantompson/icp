function [ m ] = leftMultTranslation( m, x, y, z )
% % ROW_MAJOR
% m(4) = m(4) + x;
% m(8) = m(8) + y;
% m(12) = m(12) + z;
% COLUMN_MAJOR
m(13) = m(13) + x;
m(14) = m(14) + y;
m(15) = m(15) + z;
end

