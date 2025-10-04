function [X, Y] = dataset_builder(N_LAGS, DATA_PATH, GEOHAR, method, stride)
    % DATASET_BUILDER Prepares the dataset with SPD matrices.
    % method - 'procrustes' or 'log-euclidean' for Fréchet mean computation

    if nargin < 4
        method = 'log-euclidean'; % Default method
    end

    data = readtable(DATA_PATH);
    dataset_length = height(data);

    if N_LAGS<=0 || N_LAGS>= dataset_length
       error('number of lags must be greater or equal than 1 and smaller than'+(dataset_length))
    end

    if stride<1 || stride >= dataset_length-N_LAGS
       error('stride must be greater or equal than 1 and smallereq than'+(dataset_length-N_LAGS))
    end


    for k = 1:dataset_length
        time_series(:,:,k) = invech(data(k,:)); % Build time series of RCOV matrices
    end

    % Build Train and Test sets
    Y = {};
    X = {};

    if GEOHAR
        for i = 23:dataset_length-stride-1
            obs = time_series(:,:,i-22:i-1);
            obs_squeezed = squeeze(num2cell(obs, [1 2]));
            Y{end+1} = time_series(:,:,i-1+stride);
            X{end+1} = diagblockhar(obs_squeezed{:}, method);
        end
    else
        for i = N_LAGS+1:dataset_length-stride-1
            obs = time_series(:,:,i-N_LAGS:i-1);
            obs_squeezed = squeeze(num2cell(obs, [1 2]));
            Y{end+1} = time_series(:,:,i-1+stride);
            X{end+1} = diagblock(obs_squeezed{:});
        end
    end
end


