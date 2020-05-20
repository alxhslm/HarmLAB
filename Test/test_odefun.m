function Mydot = test_odefun(t,y,w0,U,hbm,problem)
NInput = size(U,2);
w = hbm.harm.kHarm*w0(:);
Wu = repmat(1i*w,1,NInput);

ph = repmat(permute(exp(1i*w*t), [3 2 1]),NInput,1);
amp     = repmat(permute(U,         [2 3 1]),1,length(t));
ampdot  = repmat(permute(U.*Wu,     [2 3 1]),1,length(t));
ampddot = repmat(permute(U.*(Wu.^2),[2 3 1]),1,length(t));
u     = sum(real(amp     .* ph) ,3);
udot  = sum(real(ampdot  .* ph) ,3);
uddot = sum(real(ampddot .* ph) ,3);

%extract x and xdot
NDof = problem.NDof;
x = y(1:NDof,:); 
xdot = y(NDof+1:end,:);

Fe = problem.Ku*u + problem.Cu*udot + problem.Mu*uddot;
Fl = -problem.C*xdot - problem.K*x;

States = struct('t',t,'x',x,'xdot',xdot,'u',u,'udot',udot,'uddot',uddot,'w0',w0);
Fnl  = -test_model('nl' ,States,hbm,problem);

Ftot = Fl + Fnl + Fe;

%finally assemble the xdot vector
Mydot = [xdot;
         Ftot];