syms x z

brunt=0.5;

rhonot=1000;

g=9.81;

rho=brunt^2*rhonot*z/g+rhonot;

%rhobar=1000+z;
%rhoprime=


u=1/2*sin(x*pi/10)^2*rho;

uxint=int(u,x,0,20);

disp(uxint)


KE=int(uxint,z,0,5);

disp(KE)

w=1/2*sin(z*2*pi/5)^2*rho;

wxint=int(w,x,0,20);

disp(wxint)

KEw=int(wxint,z,0,5);

disp(KEw)



PEintx=int(rho*z*9.81,x,0,20);

disp(PEintx)

PE=int(PEintx,z,0,5);

disp(PE)