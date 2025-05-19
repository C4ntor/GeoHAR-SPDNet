close all; clear all; clc; pause(0.01);
confPath;
rng('default');
rng(0) ;
format long;
cd('C:\Users\andre\OneDrive - Università degli Studi di Macerata\Articolo Zhang Palma\GeoHAR-SPDNet_v2\')

%--ARGS--%
opts.loss_function= "loge"; % values: mse, loge, frob
n_lags=5;
compute_geohar=true;
if compute_geohar==true
    n_lags=3;  %as the diag block will always contain 3 matrices of size nxn where n is the number of stocks
end
method = 'log-euclidean';
data_filename = "RCOV20.csv"; %RCOVReal.csv
opts.dataDir = fullfile('./data') ;
opts.imdbPathtrain = fullfile(opts.dataDir, data_filename);
opts.batchSize = 1; 
[X,Y] = dataset_builder(n_lags, opts.imdbPathtrain, compute_geohar, method);
opts.data = struct('X', X, 'Y', Y);

opts.training_index= 1983;  % - determines the number of (test) predictions to be made as length(X) - training_index
opts.numEpochs = opts.training_index+1;
%opts.numEpochs = 100;
opts.gpus = [] ;
opts.learningRate = 0.01*ones(1,opts.training_index);
opts.weightDecay = 0.0005 ;
opts.continue = 0;

%spdnet initialization (new network)
net = spdnet_init_afew(opts) ;
[net, info, train_predictions, val_predictions] = spdnet_train_afew(net, opts, X, Y);


n = length(opts.data(1).Y);
n= n*(n+1)/2;
headers = cell(1, n); 
for i = 1:n
    headers{i} = sprintf('y%d', i);
end
ntest = length(val_predictions);
predictions = array2table(zeros(ntest, n), 'VariableNames', headers');
for i=1 : length(val_predictions)
    row= vech(val_predictions{i})';
    predictions{i,:} = row;
end
writetable(predictions, '5.SPDNetLEHARLE20.csv');