function Mydot = test_odefun(t,y,w0,U,hbm,problem)
NInput = size(U,2);
if length(w0)>1
    w = hbm.harm.kHarm(:,1)*w0(1) + hbm.harm.kHarm(:,2)*w0(2);
else
    w = (0:hbm.harm.NHarm)'*w0;
end
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
xdot = y(NDof+1:2*NDof,:);
xalg = y(2*NDof+1:end,:);

Fe = problem.Ku*u + problem.Cu*udot + problem.Mu*uddot;
Fl = -problem.C*xdot - problem.K*x;
Fnl  = -test_model('nl' ,t,x,xdot,0*xdot,u,udot,uddot,xalg,hbm,problem,w0);
Falg =  test_model('alg',t,x,xdot,0*xdot,u,udot,uddot,xalg,hbm,problem,w0);

Ftot = Fl + Fnl + Fe;

%finally assemble the xdot vector
Mydot = [xdot;
         Ftot;
         Falg];