% load_quasar_data
%
% Loads the data in the quasar data files
%
% Upon completion of this script, the matrices and data are as follows:
%
% lambdas - A length n = 450 vector of wavelengths {1150, ..., 1599}
% train_qso - A size m-by-n matrix, where m = 200 and n = 450, of noisy
%      observed quasar spectra for training.
% test_qso - A size m-by-n matrix, where m = 200 and n = 450, of noisy observed
%       quasar spectra for testing.

load quasar_train.csv;
lambdas = quasar_train(1, :)';
train_qso = quasar_train(2:end, :);
load quasar_test.csv;
test_qso = quasar_test(2:end, :);

%Theta Implementation
y = train_qso(1,:)';

x0 = [ones(rows(lambdas),1),lambdas];

minw = min(x0(:,2));
maxw = max(x0(:,2));

mx = 2/(maxw-minw);
bx = 1-mx*maxw;
NM = [1 bx; 0 mx];

X=x0*NM;


%Nomralize Y
minRes = min(y);
maxRes = max(y);
my = 2/(maxRes - minRes);
by = 1 - my*maxRes;

y = my*y + by;


theta = inv(X'*X)* X' * y

figure(1);
  hold off;
  plot(lambdas,y,"xk","linewidth",1);
  hold on;
  plot(lambdas,X*theta,"r","linewidth",2);
  xlabel('Lambdas (Angstrom (Å))');
  ylabel('Flujo de Espectral');


#Tau = 5  
figure(2);
hold off;

xpon = [ones(rows(lambdas),1),lambdas(:,1)];
ypon = train_qso(1,:)';
plot(lambdas,ypon,"+r");
hold on;
tau=5;
yar=[];  
    for xite = xpon'      
        w = exp(-(xpon(:,2)-xite(2,1)).^2/(2*tau^2));
        W = diag(w, 0);
        theta = inv(xpon' * W * xpon) * xpon' * W * ypon;
        yar = [yar theta'*xite];
    end
yp2 = yar';
plot(lambdas, yp2,'b', 'linewidth', 3);
hold on;
xlabel('Lambdas (Angstrom (Å))');
ylabel('Flujo de Espectral');
title ("Tau = 5");


#Tau = 1  
figure(3);
hold off;

xpon = [ones(rows(lambdas),1),lambdas(:,1)];
ypon = train_qso(1,:)';
plot(lambdas,ypon,"+r");
hold on;
tau=1;
yar=[];  
    for xite = xpon'      
        w = exp(-(xpon(:,2)-xite(2,1)).^2/(2*tau^2));
        W = diag(w, 0);
        theta = inv(xpon' * W * xpon) * xpon' * W * ypon;
        yar = [yar theta'*xite];
    end
yp2 = yar';
plot(lambdas, yp2,'b', 'linewidth', 3);
hold on;
xlabel('Lambdas (Angstrom (Å))');
ylabel('Flujo de Espectral');
title ("Tau = 1");

#Tau = 10  
figure(4);
hold off;

xpon = [ones(rows(lambdas),1),lambdas(:,1)];
ypon = train_qso(1,:)';
plot(lambdas,ypon,"+r");
hold on;
tau=10;
yar=[];  
    for xite = xpon'      
        w = exp(-(xpon(:,2)-xite(2,1)).^2/(2*tau^2));
        W = diag(w, 0);
        theta = inv(xpon' * W * xpon) * xpon' * W * ypon;
        yar = [yar theta'*xite];
    end
yp2 = yar';
plot(lambdas, yp2,'b', 'linewidth', 3);
hold on;
xlabel('Lambdas (Angstrom (Å))');
ylabel('Flujo de Espectral');
title ("Tau = 10");

#Tau = 100  
figure(5);
hold off;

xpon = [ones(rows(lambdas),1),lambdas(:,1)];
ypon = train_qso(1,:)';
plot(lambdas,ypon,"+r");
hold on;
tau=100;
yar=[];  
    for xite = xpon'      
        w = exp(-(xpon(:,2)-xite(2,1)).^2/(2*tau^2));
        W = diag(w, 0);
        theta = inv(xpon' * W * xpon) * xpon' * W * ypon;
        yar = [yar theta'*xite];
    end
yp2 = yar';
plot(lambdas, yp2,'b', 'linewidth', 3);
hold on;
xlabel('Lambdas (Angstrom (Å))');
ylabel('Flujo de Espectral');
title ("Tau = 100");


#Tau = 1000 
figure(6);
hold off;

xpon = [ones(rows(lambdas),1),lambdas(:,1)];
ypon = train_qso(1,:)';
plot(lambdas,ypon,"+r");
hold on;
tau=1000;
yar=[];  
    for xite = xpon'      
        w = exp(-(xpon(:,2)-xite(2,1)).^2/(2*tau^2));
        W = diag(w, 0);
        theta = inv(xpon' * W * xpon) * xpon' * W * ypon;
        yar = [yar theta'*xite];
    end
yp2 = yar';
plot(lambdas, yp2,'b', 'linewidth', 3);
hold on;
xlabel('Lambdas (Angstrom (Å))');
ylabel('Flujo de Espectral');
title ("Tau = 1000");
