function [net, info, train_predictions, predictions] = neuralnet_train_afew(net, opts, X, Y)
opts.errorLabels = {'top1e'};
c=1;
epoch=1;
predictions={}
for test_index=opts.numEpochs:length(X)
    epoch=c;
    learningRate = opts.learningRate(epoch);
    Xtrain=X(c:test_index-1);
    Xtest = X(test_index:test_index);
    Ytrain =Y(c:test_index-1);
    Ytest = Y(test_index:test_index);
    train= struct('X', Xtrain, 'Y', Ytrain, 'B', Xtrain);
    val=struct('X', Xtest, 'Y', Ytest, 'B', Xtest);
    c=c+1;
     % fast-forward to last checkpoint
     modelPath = @(ep) fullfile(opts.dataDir, sprintf('net-epoch-%d.mat', ep));
     modelFigPath = fullfile(opts.dataDir, 'net-train.pdf') ;
     if opts.continue
         if exist(modelPath(epoch),'file')
             if epoch == opts.numEpochs
                 load(modelPath(epoch), 'net', 'info') ;
             end
             continue ;
         end
         if epoch > 1
             fprintf('resuming by loading epoch %d\n', epoch-1) ;
             load(modelPath(epoch-1), 'net', 'info') ;
         end
     end
    
    
    [net,stats.train, train_predictions] = process_epoch(opts, train, learningRate, net) ;
    [net,stats.val, val_predictions] = process_epoch(opts, val, 0, net) ;
   
    predictions{end+1} = val_predictions{1};
    % save
    evaluateMode = 0;
    
    if evaluateMode, sets = {'train'} ; else sets = {'train', 'val'} ; end

    for f = sets
        f = char(f) ;
        n = numel(eval(f)) ; %
        info.(f).objective(epoch) = stats.(f)(2) / n ;
        info.(f).error(:,epoch) = stats.(f)(3:end) / n ;
    end
    %if ~evaluateMode, save(modelPath(epoch), 'net', 'info') ; end

    figure(1) ; clf ;
    hasError = 1 ;
    subplot(1,1+hasError,1) ;
    if ~evaluateMode
        semilogy(1:epoch, info.train.objective, '.-', 'linewidth', 2) ;
        hold on ;
    end
    semilogy(1:epoch, info.val.objective, '.--') ;
    xlabel('training epoch') ; ylabel('energy') ;
    grid on ;
    h=legend(sets) ;
    set(h,'color','none');
    title('objective') ;
    if hasError
        subplot(1,2,2) ; leg = {} ;
        if ~evaluateMode
            plot(1:epoch, info.train.error', '.-', 'linewidth', 2) ;
            hold on ;
            leg = horzcat(leg, strcat('train ', opts.errorLabels)) ;
        end
        plot(1:epoch, info.val.error', '.--') ;
        leg = horzcat(leg, strcat('val ', opts.errorLabels)) ;
        set(legend(leg{:}),'color','none') ;
        grid on ;
        xlabel('training epoch') ; ylabel('error') ;
        title('error') ;
    end
    drawnow ;
    %print(1, modelFigPath, '-dpdf') ;
    disp('---')
    fprintf('EPOCH: %d',epoch);
    fprintf(' TRAIN-MSE: %.3f', stats.train(3)) ;
    fprintf(' TRAIN-loss: %.3f', stats.train(2)) ;
    fprintf(' TEST-MSE: %.3f', stats.val(3)) ;
    fprintf(' TEST-loss: %.3f', stats.val(2)) ;
end
end
    
    
function [net,stats, out_array] = process_epoch(opts, dataset, learningRate, net)

training = learningRate > 0 ;
if training, mode = 'training' ; else mode = 'validation' ; end

stats = [0 ; 0 ; 0] ; % [totalTime; totalLoss; totalError]
numGpus = numel(opts.gpus) ;
if numGpus >= 1
    one = gpuArray(single(1)) ;
else
    one = single(1) ;
end

batchSize = opts.batchSize;
if ~exist('batchSize','var') || isempty(batchSize), batchSize = 1; end
N_total = length(dataset);
totalSamples=0;
totalGrad = cell(size(net.layers));
isGradInit = false;

errors = 0;
numDone = 0 ;
out_array=cell(1,N_total);
for ib = 1 : batchSize : N_total
    %fprintf('%s: epoch %02d: batch %3d/%3d:', mode, epoch, ib,length(dataset)) ;
    batchTime = tic ;
    res = [];
    if (ib+batchSize> N_total)
        batchSize_r = N_total-ib+1;
    else
        batchSize_r = batchSize;
    end

    totalSamples = totalSamples + batchSize_r;
    batch_data = cell(batchSize_r,1);
    batch_target = cell(batchSize_r,1);
    batch_baseline = cell(batchSize_r,1);
    
    for j = 1 : batchSize_r
        dataset_index = ib + j - 1; 
        batch_data{j} = dataset(dataset_index).X;  
        batch_target{j} = dataset(dataset_index).Y;  
        batch_baseline{j} = dataset(dataset_index).B;
    end
    net.layers{end}.class = batch_target;

    %forward/backward spdnet
    if training, dzdy = one; else dzdy = [] ; end
    res = vl_myforbackward(net, batch_data, dzdy, res) ; 


    %accumulating gradients
    for l = 1:numel(net.layers)
        if ~isempty(res(l).dzdw)
            if ~isGradInit
                for ll = 1:numel(net.layers)
                    if ~isempty(res(ll).dzdw)
                        totalGrad{ll} = zeros(size(res(ll).dzdw), 'like', res(ll).dzdw);
                    else
                        totalGrad{ll}=[];
                    end
                end
                isGradInit=true;
            end
            totalGrad{l} = totalGrad{l} + (res(l).dzdw );
        end
    end


    % accumulate training errors
    %
    predictions = gather(res(end-1).x) ;
    out_array{ib} = predictions{1};
    error = 0;
    for ixx=1 : length(predictions)
        n=length(predictions{ixx});
        diff = predictions{ixx}-batch_target{ixx};
        error = error + sum(diff(:).^2)/ (n^2); 
    end

    numDone = numDone + batchSize_r ;
    errors = error/length(predictions);
    batchTime = toc(batchTime) ;
    speed = batchSize/batchTime ;
    stats = stats+[batchTime ; res(end).x ; error]; % works even when stats=[]
    
end
    
    
    %fprintf(' %.2f s (%.1f data/s)', batchTime, speed) ;

    %fprintf(' MSE: %.5f', stats(3)) ;
    %fprintf(' loss: %.5f', stats(2)) ;
    %fprintf(' [%d/%d]', numDone, batchSize_r);
    %fprintf('\n') ;
    
if training && isGradInit
    for l = 1:numel(net.layers)
            if ~isempty(totalGrad{l})
                avgGrad = totalGrad{l} / totalSamples; % average per-sample gradient
                
                if ~isfield(net.layers{l}, 'learningRate') 
                    net.layers{l}.learningRate = 1 ; 
                end
                if ~isfield(net.layers{l}, 'weightDecay')
                    net.layers{l}.weightDecay = 0 ;
                end
                thisLR = learningRate * net.layers{l}.learningRate ;
                if isfield(net.layers{l}, 'weight') && net.layers{l}.weightDecay ~= 0 
                    avgGrad = avgGrad + net.layers{l}.weightDecay * net.layers{l}.weight; 
                end
                if isfield(net.layers{l}, 'weight') && ~isempty(net.layers{l}.weight)
                    net.layers{l}.weight = net.layers{l}.weight - thisLR * avgGrad ;
                end
            end
    end
end
    
end

function y=gather(x)
    try 
        if isa(x, 'gpuArray')
            y=builtin('gather',x);
        else
            y=x;
        end
    catch 
        y=x;
    end
end 
