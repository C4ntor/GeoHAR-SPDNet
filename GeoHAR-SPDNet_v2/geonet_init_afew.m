function net = geonet_init_afew(opts)
% spdnet_init initializes a spdnet
rng('default');
rng(0) ;
opts.layernum = 3;
Winit = cell(opts.layernum,1);
opts.datadim = [length(opts.data(1).X), 80, 60, length(opts.data(1).Y)];
% for RCOV20 use [length(opts.data(1).X), 40, 30, length(opts.data(1).Y)];


for iw = 1 : opts.layernum
    A = rand(opts.datadim(iw));
    [U1, S1, V1] = svd(A * A');
    curDim = opts.datadim(iw);
    nextDim = opts.datadim(iw+1);
    fprintf('layer %d: curDim = %d, nextDim = %d\n', iw, curDim, nextDim);
    Winit{iw} = U1(:,1:opts.datadim(iw+1));
end


net.layers = {} ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{1}) ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{2}) ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{3}) ;
net.layers{end+1} = struct('type', opts.loss_function) ;






