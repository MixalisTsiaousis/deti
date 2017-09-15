function [Ai,Bi,Ci,Di] = deti(u,y)
format long
if length(u) ~= length(y)
        fprintf('\n\n Identification could not be completed. Length of input signal u and output signal y must be the same.')
        Ai = [];
        Bi = [];
        Ci = [];
        Di = [];
        return
end
J = [];
sys_ranks = [];
matfit = [];
sub = [];
s = size(u,1);
m = size(u,2);
l = size(y,2);
i = 2;
j = s-2*i+1;
t_steps = 0:(length(u)-1);
n = 1;
a = 2*m*i+n;
inside_while = false;
while j >= a
    H1 = zeros(i*m+i*l,j);
    H2 = zeros(i*m+i*l,j);
    U_1 = cell(1,m);
    U_2 = cell(1,m);
    A = zeros(i,j,m);
    B = zeros(i,j,m);
    for t = 1:m
        u_1(1:i,t:t) = u(1:i,t:t);  
        u_2(1:j,t:t) = u(i:i+j-1,t:t);   
        u_3(1:i,t:t) = u(i+1:2*i,t:t);   
        u_4(1:j,t:t) = u(2*i:2*i+j-1,t:t);
        U_1{t} = zeros(i,j);
        U_2{t} = zeros(i,j);
        U_1{t} = hankel(u_1(1:i,t:t),u_2(1:j,t:t));  
        U_2{t} = hankel(u_3(1:i,t:t),u_4(1:j,t:t));
        A(:,:,t) = cell2mat(U_1(t));  
        B(:,:,t) = cell2mat(U_2(t));
        H1(t:m+l:end,:) = A(:,:,t);
        H2(t:m+l:end,:) = B(:,:,t);
    end
%-------------------------------------------
    Y_1 = cell(1,l);
    Y_2 = cell(1,l);
    C = zeros(i,j,l);
    D = zeros(i,j,l);
    for t = 1:l
        y_1(1:i,t:t) = y(1:i,t:t);    
        y_2(1:j,t:t) = y(i:i+j-1,t:t);             
        y_3(1:i,t:t) = y(i+1:2*i,t:t);
        y_4(1:j,t:t) = y(2*i:2*i+j-1,t:t);
        Y_1{t} = hankel(y_1(1:i,t:t),y_2(1:j,t:t));
        Y_2{t} = hankel(y_3(1:i,t:t),y_4(1:j,t:t));
        C(:,:,t) = cell2mat(Y_1(t));
        D(:,:,t) = cell2mat(Y_2(t));
        H1(m+t:m+l:end,:) = C(:,:,t);
        H2(m+t:m+l:end,:) = D(:,:,t); 
    end
%--------------------------------------------   
    H = [H1;H2];
    rank_H1 = rank(H1);
    n = rank_H1-m*i;
    a = 2*m*i+n;
    if a > j
        Ai = [];
        Bi = [];
        Ci = [];
        Di = [];
        fprintf('\n\n Identification could not be completed. Insert more input/output measurements.\n')
        return
    end
    [U,S] = svd(H);  
    U11 = U(1:m*i+l*i,1:2*m*i+n);
    U12 = U(1:m*i+l*i,2*m*i+n+1:end);
    S11 = S(1:2*m*i+n,1:2*m*i+n);
    U12_t = transpose(U12);
    K = U12_t*U11*S11;
    [U2,~] = svd(K);   
    Uq = U2(1:2*l*i-n,1:n);
    Uq_t = transpose(Uq);
    O = Uq_t*U12_t*U(m+l+1:(i+1)*(m+l),:)*S;  
    R = U((m*i+l*i+m+1):(m+l)*(i+1),:)*S;
    N = Uq_t*U12_t*U(1:m*i+l*i,:)*S;
    T = U(m*i+l*i+1:m*i+l*i+m,:)*S;
    A_least_square = [N;T];
    B_least_square = [O;R];
    sol = B_least_square/A_least_square;
    Ai = sol(1:n,1:n);
    Bi = sol(1:n,n+1:end);
    Ci = sol(n+1:end,1:n);
    Di = sol(n+1:end,n+1:end);
    idcell = {{Ai Bi Ci Di}};
    J = [J idcell];
    sysi = ss(Ai,Bi,Ci,Di,[]);
    yi = lsim(sysi,u,t_steps);
    fit = goodnessOfFit(yi,y,'NRMSE');
    matfit = [matfit fit];
    sys_ranks = [sys_ranks n];
    sub = [sub i];
    i = i+1;
    j = s-2*i+1;
    a = 2*m*i+n;
    inside_while = true;
end
if inside_while == false
    Ai = [];
    Bi = [];
    Ci = [];
    Di = [];
    fprintf('\n\n Identification could not be completed. Insert more input/output measurements.\n')
    return
end  
if l > 1
    means = mean(matfit);
    [~,I] = max(means);
    bestfit = matfit(:,I);
else
    [bestfit,I] = max(matfit);
end      
if bestfit < 0.9999*ones(size(bestfit))
    Ai = [];
    Bi = [];
    Ci = [];
    Di = [];
    fprintf('\n\n Identification could not be completed. Insert more input/output measurements.\n')
    return
end   
Ai = J{I}{1};
Bi = J{I}{2};
Ci = J{I}{3};
Di = J{I}{4};
n = sys_ranks(I);
i = sub(I);
if isempty(Ai) == 0
    fprintf('\n\nThe rank of the system is %d',n)
    fprintf('\n\nGoodness of fit:%.3f',bestfit);
    fprintf('\n\nSubscript i is: %d\n',i);
end
end    