%% Build the [0;1]^n hyper box (...)
function [G] = create_unitBox(n)

G = eye(n,n);
G = [G ones(n,1)];

end

