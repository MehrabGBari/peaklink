function k = kfun(u,v,p1,p2)
k = tanh(p1*(u*v')+p2);