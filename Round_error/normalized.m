function [ normalized_poly ] = normalized( q,vars,interval )
    a=0;
    b=0;
    n=length(vars);
    normalized_poly=q;
    for i=1:n 
        a=interval(i,1);
        b=interval(i,2);
        normalized_poly = subs(normalized_poly,vars(i),vars(i)*(b-a)+a);
    end
end