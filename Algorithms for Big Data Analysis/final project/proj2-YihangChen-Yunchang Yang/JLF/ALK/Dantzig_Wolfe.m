function [x, fval, bound, exit_flag] = Dantzig_Wolfe(mast, sub, K)
%Dantzig Wolfe Decomp solves a linear programming with a special structure:
% min c'*x
% s.t. A*x <= b
% x >= 0
%where A is an m*n matrix which can be write as block angular form:
% __ __ __ __ __ __ __ __
% | L1 L2 ... Lk | | x1 | | c1 | | b0 |
% | A1 | | x2 | | c2 | | b1 |
% A = | A2 |, x = | : |, c =| : |, b = | b2 |
% | ... | | : | | : | | : |
% |__ Ak __| |__ xk __| |__ ck __| |__ bk __|
%so the LP can be decomposed into a master problem(MP) and k subproblems(SPk),
%we can rewrite the MP as a restricted master problem by Resolution Theorem
%
% Inputs:
% mast is a struct includes MP's coefficients, i.e. L1,...,Lk, b0
% sub is a struct includes coefficients of SPks, i.e. c1,...,ck, A1,...,Ak,
% b1,...,bk and the initial extreme points v1, ..., vk.
% K is the number of subproblems
%
% Outputs:
% x = n*1 vector, final solution
% fval is a scalar, final objective value
% bound is a matrix includes all LBs and UBs for each iteration
% exit_flag describes the exit condition of the problem as follows:
% 0 - optimal solution
% 1 - LP is unboundeds
x=[]; fval=[]; bound=[]; exit_flag=[];
%% Step 0: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate an initial basis B for the master problem.
% Let x_B be the basic variables and B_bar the index set of the basic variables
% Set all other variables to zero to get the restricted master problem.
% Go to STEP 1.
s = mast.b; %slack variables for inequality constraints
x_B = [s; ones(K,1)]; %initial basic variables
%x_Bflag is an index of the basic variables in the restricted master problem
%associated with linking constraints and subproblem k, i.e. slack variables of
%linking constraints are initially basic and other basic variables
%associated with subproblems are set to 1.
x_Bflag = [zeros(length(s),1); [1:K]'];
f_sub=[];
for k=1:K
%obtain initial extreme points from struct sub for subproblems these are
%zero vectors, so initial objective function values will be zero, v(k)
%is initial extreme point of subproblem k
v_sub{k}=sub.v{k};
f_sub=cat(1, f_sub, sub.c{k}'*v_sub{k});
for a=1:length(s)
%generating initial extreme point for linking constraint
%v_L{a,k}= initial extreme point of linking constraint of Lk
%b(k) is the RHS vector in A_k*x^k=b^k
v_L{a,k}=zeros(length(sub.b{k}),1);
end
end
f_s=zeros(length(s),1);
f_B=[f_s; f_sub]; %initial f_B i.e. the objective coefficient of
%the restricted MP
B=eye(length(x_B)); %initial basis for master problem
iter_num=0; %counter
options=optimset('LargeScale', 'on');%choose largescale LP
%solver in linprog
while iter_num>=0
%% Step 1: Simplex Multiplier Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for duals by solving the system B^T*pie=f_B, then Go to STEP2.
pie=B'\f_B; %solve B^T*pie=f_B
pie_sub=pie(end-K+1:end);%duals of kth convexity constraints, pie_k_2
pie(end-K+1:end)=[]; %duals of linking constraints, pie_1
%% Step 2: Optimality Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each k=1,...,K, using the revised simplex method to solve SP_k i.e.
% min sig_k=[(c^k)^T-(pie^1)^T*L_k]*x^k
% s.t.A_k*x^k <= b^k
% x^k >= 0
% If SP_k is unbounded, then Go to STEP3, else let x^k=(v_i*)^k denote the
% optimal basic feasible solution, compute(r^k)_* = (sig_k)^* - pie_k^2.
% If r_min={(r^k)_*}>=0, then the current basis B is optimal,
% else Go to STEP3.
for k=1:K
c_sub=[sub.c{k}'-pie'*mast.L{k}]'; %update the objective coefficient
[x_sub{iter_num+1, k}, sig(k) exitflag(k)] = ... %call linprog solver
linprog(c_sub, sub.A{k}, sub.b{k},[],[],zeros(length(c_sub),1),[],[],options);
sig(k)=sig(k)-pie_sub(k); %computer (r^k)_*
end
r_min=min(sig); %minimum ratio test to obtain r_min
r_minflag=find(sig==r_min);
if isempty(find(r_min < 0)) || abs(r_min) <= 1e-8 %reduced cost>=0, optimal
disp('problem solved')
fval=0;
x_Bflag_s=x_Bflag(1:length(s));
for k=1:K %convert to optimal solution for original LP
x{k,1}=x_B(length(s)+k)*v_sub{k};
for a=1:length(x_Bflag_s)
if x_Bflag_s(a)==k
    x{k}=x{k}+x_B(a)*v_L{a, k};
end
end
%convert to optimal obj val for original LP
fval=fval+sub.c{k}'*x{k};
end
x=cell2mat(x);
exit_flag=0;
break
else
%% Step 3: Column Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If all subproblems are bounded and r_min={(r^k)_*}<0, then let t be
% the index of k in SP_k such that r_min=(r^t)_*.
% Let alpha_bar=[q_i*^t e_t]^T =[L_t*v_i*^t e_t]^T where v_i*^t is the
% optimal extreme points of SP_t and Go to STEP4. Else there is a
% subproblem SP_s that is unbounded and so an extreme direction d_j*^s
% will be generated such that [(c^s)^T-(pie^1)^T*L_s]*d_j*^s < 0 and
% so let alpha_bar = [(q_bar)_j*^s 0]^T = [L_s*d_j*^s 0]^T and
% go to STEP4.
if length(find(exitflag==1)) == K %if subproblems bounded and r_min<0
t=r_minflag(1); %subproblem t such that r_min=r_*_t
q_t=mast.L{t}*x_sub{iter_num+1,t}; %generate q_i*_t
e=zeros(length(x_B)-length(q_t),1);
e(t)=1; %generate e_t
alpha_bar=[q_t; e]; %generate alpha_bar
else %if any subproblems is unbounded
disp('unbouded subproblem exist')
unboundflag=find(exitflag==-3);
t=unboundflag(1); %subproblem s with extreme direction d_j*^s
q_t=mast.L{t}*x_sub{iter_num+1,t}; %generate (q_bar)_j*^s
%generate alpha_bar
alpha_bar=[q_t; zeros(length(x_B)-length(q_t),1)];
end
end
%% Step 4: Descent Direction Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for d in the linear system Bd = -alpha_bar.
% If d>=0, then the LP is unbounded STOP, else go to STEP5.
d=-B\alpha_bar; %solve Bd=-alpha_bar
d_flag=find(d<0);
if isempty(d_flag) %if d>=0, unbounded
disp('unbounded LP')
exit_flag=1;
return
else %else Go to STEP 5
%% Step 5: Step Length Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer the step length alpha by minimum ratio test. Let l* be
% the index of the basic variables then attains the minimum ratio
% alpha. Go to STEP 6.
alpha=min(x_B(d_flag)./abs(d(d_flag))); %minimum ratio test
%% Step 6: Update Basic Variables and Basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let x_B = x_B + alpha*d. Go to STEP 7.
x_B=x_B+alpha*d; %get new basis variables
delta=1e-30; %computation error tolerance
leave=find(abs(x_B)<=delta); %index of leave variable
while isempty(leave)
delta=10*delta;
leave=find(abs(x_B)<=delta);
end
x_B(leave(1))=alpha;
x_Bflag(leave(1))=t;
if leave(1) <= length(s) %update f_s and extreme point
f_s(leave(1))=sub.c{t}'*x_sub{iter_num+1,t};
v_L{leave(1),t}=x_sub{iter_num+1,t};
else
f_sub(leave(1)-length(s))=sub.c{t}'*x_sub{iter_num+1,t};
v_sub{leave(1)-length(s)}=x_sub{iter_num+1,t};
end
%% Step 7: Basis Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let B_l* be the column associated with the leaving variable x_l*.
% Update the basis matrix B by removing B_l* and adding the column
% alpha_bar, and update B_set.
B(:,leave(1))=alpha_bar; %update the basis B
end
iter_num=iter_num+1;
f_B=[f_s; f_sub]; %update f_B for next iteration
bound(:,iter_num)=[f_B'*x_B + sum(sig); f_B'*x_B];%new lower/upper bound
end %Go to STEP 1
end
