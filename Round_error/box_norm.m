function [ pn ] = box_norm( p,v,int )
    n=length(v);
    pn=p;
    for i=1:n 
        a=int(i,1);
        b=int(i,2);
        pn = replace(pn,v(i),v(i)*(b-a)+a);
    end
end
