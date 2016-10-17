function [ostr,lvarstr] = build_sdpvar( nvars,nparams )

ostr = 'sdpvar ';
lvarstr = '[';
for i=1:(nvars+nparams)
    ostr = [ostr 'x' num2str(i) ' '];
    lvarstr = [lvarstr 'x' num2str(i) ' '];
end
ostr = [ostr ';'];
lvarstr = [lvarstr '];'];


end

