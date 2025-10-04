function net = neuralnet_init_afew(opts)
%initializes an MLP
rng('default');
rng(0) ;
opts.layernum = 3;
Winit = cell(opts.layernum,1);
opts.datadim = [length(opts.data(1).X), 80, 60, length(opts.data(1).Y)];
% for RCOV20 use [length(opts.data(1).X), 40, 30, length(opts.data(1).Y)];

for iw = 1 : opts.layernum
    curDim = opts.datadim(iw);
    nextDim = opts.datadim(iw+1);
    %fprintf('layer %d: curDim = %d, nextDim = %d\n', iw, curDim, nextDim);
    Winit{iw} = 0.1*randn(opts.datadim(iw), opts.datadim(iw+1));  %rescale weights, othw gradient might explode
end


net.layers = {} ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{1}) ;
net.layers{end+1} = struct('type', 'rec') ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{2}) ;
net.layers{end+1} = struct('type', 'rec') ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{3}) ;
net.layers{end+1} = struct('type', opts.loss_function) ;






