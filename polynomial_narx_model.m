clc 
clear 
close all

data = load("iddata-05.mat");
id_data  = data.id;
val_data = data.val;

y_id  = id_data.y(:);   u_id  = id_data.u(:);
y_val = val_data.y(:);  u_val = val_data.u(:);


na = 1;
nb = 2;
nk = 1;

maxm = 20;        

[X_id,  Y_id]  = xy_fcn(y_id,  u_id,  na, nb);
[X_val, Y_val] = xy_fcn(y_val, u_val, na, nb);

n = max(na,nb);

%% 

MSE_pred_val_all = zeros(maxm,1);
MSE_sim_val_all  = zeros(maxm,1);
theta_all = cell(maxm,1); 

for m = 1:maxm
    
    PHI_id = phi_fcn(X_id, m);
    theta  = PHI_id \ Y_id;
    theta_all{m} = theta; 

    PHI_val = phi_fcn(X_val, m);
    y_pred_val = PHI_val * theta;
    MSE_pred_val_all(m) = mean((Y_val - y_pred_val).^2); 

    y_sim_val = sim_fcn(y_val, u_val, na, nb, m, theta);
    MSE_sim_val_all(m) = mean((y_val(n+1:end) - y_sim_val(n+1:end)).^2);
end

[mini, m_opt] = min(MSE_pred_val_all); 

figure;
plot(MSE_pred_val_all,'-o'); 
hold on;
plot(MSE_sim_val_all,'-o');
hold on;
plot(m_opt, MSE_pred_val_all(m_opt),'x','MarkerSize',18,'LineWidth',2);
grid on;
legend("MSE pred VAL", "MSE sim VAL", "m optim");
xlabel('Polynomial degree m');
ylabel('Mean Squared Error (MSE)');
title(sprintf('MSE on validation (na=%d, nb=%d)', na, nb));

%% m optim

m = 3;
theta = theta_all{m};

PHI_id  = phi_fcn(X_id,  m);
PHI_val = phi_fcn(X_val, m);
 
y_pred_id  = PHI_id  * theta;
y_pred_val = PHI_val * theta;

MSE_pred_id  = mean((Y_id  - y_pred_id ).^2);
MSE_pred_val = mean((Y_val - y_pred_val).^2);

y_sim_id  = sim_fcn(y_id,  u_id,  na, nb, m, theta);
y_sim_val = sim_fcn(y_val, u_val, na, nb, m, theta);

MSE_sim_id  = mean((y_id(n+1:end)  - y_sim_id(n+1:end)).^2); 
MSE_sim_val = mean((y_val(n+1:end) - y_sim_val(n+1:end)).^2);


figure;
plot(Y_id,'b');
hold on; 
plot(y_pred_id,'r');
grid on;
title(sprintf('PREDICTION ID \n(na=%d, nb=%d, m=%d)', na, nb, m));
xlabel('k (sample index)'); 
ylabel('y (output)');
legend("Y_{id} (real)","y_{pred,id}");

figure;
plot(Y_val,'b'); 
hold on; 
plot(y_pred_val,'r');
grid on;
title(sprintf('PREDICTION VAL \n(na=%d, nb=%d, m=%d)', na, nb, m));
xlabel('k (sample index)'); 
ylabel('y (output)');
legend("Y_{val} (real)","y_{pred,val}");

figure;
plot(y_id,'b');
hold on;
plot(y_sim_id,'r');
grid on;
title(sprintf('SIMULATION ID \n(na=%d, nb=%d, m=%d)', na, nb, m));
xlabel('k (sample index)');
ylabel('y (output)');
legend("y_{id} (real)","y_{sim,id}");

figure;
plot(y_val,'b');
hold on; 
plot(y_sim_val,'r');
grid on;
title(sprintf('SIMULATION VAL \n(na=%d, nb=%d, m=%d)', na, nb, m));
xlabel('k (sample index)');
ylabel('y (output)');
legend("y_{val} (real)","y_{sim,val}");

%% MSE for m optim

disp(["MSE prediction ID  = " num2str(MSE_pred_id)]);
disp(["MSE prediction VAL = " num2str(MSE_pred_val)]);
disp(["MSE simulation ID  = " num2str(MSE_sim_id)]);
disp(["MSE simulation VAL = " num2str(MSE_sim_val)]);

%% functions

function [X, Y] = xy_fcn(y, u, na, nb)

    n = max(na,nb);
    N = length(y);
    rows = N - n;

    X = zeros(rows, na+nb);
    Y = zeros(rows, 1);

    for k = n+1:N
        r = k-n;

        for i = 1:na
            X(r,i) = -y(k-i);
        end

        for j = 1:nb
            X(r,na+j) = u(k-j);
        end

        Y(r) = y(k);
    end
end

function phi_vect = phi_row_fcn(x, m)

    d = length(x);
    phi_vect = [];

    exp = zeros(1,d);
    gata = 0; 

    while gata == 0
        if sum(exp) <= m
            term = 1;
            for j = 1:d
                if exp(j) > 0
                    term = term * (x(j)^exp(j));
                end
            end
            phi_vect = [phi_vect, term];
        end

        index = 1;
        while index <= d
            exp(index) = exp(index) + 1; 

            if exp(index) <= m
                break;
            else
                exp(index) = 0;
                index = index + 1;
            end
        end

        if index > d
            gata = 1;
        end
    end
end


function PHI = phi_fcn(X, m)

    rows = size(X,1); 
    PHI = [];

    for k = 1:rows 
        x = X(k,:);
        phi_vect = phi_row_fcn(x, m);
        PHI = [PHI; phi_vect]; 
    end
end

function y_sim = sim_fcn(y, u, na, nb, m, theta)

    n = max(na,nb);
    N = length(y);
    d = na + nb;

    y_sim = zeros(N,1);
    y_sim(1:n) = y(1:n);

    for k = n+1:N
        x = zeros(1,d);

        for i = 1:na
            x(i) = -y_sim(k-i);
        end
        for j = 1:nb
            x(na+j) = u(k-j); 
        end

        phi_vect = phi_row_fcn(x, m);
        y_sim(k) = phi_vect * theta;
    end
end
