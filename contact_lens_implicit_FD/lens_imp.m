#! /usr/bin/octave

#This should be run with octave

#Parameters
#Contact Length
L = 1.0;

#Diffusion Constant
D = 0.1;


#n_vars is the number of variables *not including* the endpoints
#so the variables should be indexed from 1.
n_vars = 500;
dx = L/(n_vars-1);
dt = 0.01;

#Construct matrix for implicit finite difference method
alpha = D*dt/dx^2;
A = full(spdiags(ones(n_vars,1)*[-alpha,1+2*alpha,-alpha],[-1,0,1],n_vars,n_vars));

#Setup main loop
C_post = 0;
C_pre  = 0;

c_old = ones(n_vars,1);
save "data.dat" c_old


#Each iteration
for i = 1:100
    b = c_old;
    b(1)     += alpha*C_post;
    b(end-1) += alpha*C_pre;

    #c_new = bicgstab(A, b, 1e-06, 500);
    c_new = A\ b;

    save "data.dat" -append c_new;
    c_old = c_new;
endfor
