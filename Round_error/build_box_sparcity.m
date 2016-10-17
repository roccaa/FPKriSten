function [I,J] = build_box_sparcity( nvar,nparam )

var = 1:nvar;
param = nvar+1:nvar+nparam;
I = repmat(var,nparam,1);
I = [I param'];
J = I;
end

