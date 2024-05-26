% surf visualization function
function visualization(X, Y, Z, xlabel_, ylabel_, zlabel_, xlim_, ylim_, title_, zlim_)
% plot to show concentration in a specific plain
% X, Y: X-Y plain coordinate
% Z: concentration

surf(X, Y, Z);
colorbar;
shading flat;
xlabel(xlabel_);
ylabel(ylabel_);
zlabel(zlabel_);
xlim(xlim_);
ylim(ylim_);
if (nargin==10)
    % normally we do not constrain zlim except for comparing directly
    zlim(zlim_);
end
title(title_);
drawnow;

end