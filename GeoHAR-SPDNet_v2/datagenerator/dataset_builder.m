function [X, Y] = dataset_builder(N_LAGS, DATA_PATH, GEOHAR, method)
    % DATASET_BUILDER Prepares the dataset with SPD matrices.
    % method - 'procrustes' or 'log-euclidean' for Fréchet mean computation

    if nargin < 4
        method = 'log-euclidean'; % Default method
    end

    data = readtable(DATA_PATH);
    dataset_length = height(data);

    for k = 1:dataset_length
        time_series(:,:,k) = invech(data(k,:)); % Build time series of RCOV matrices
    end

    % Build Train and Test sets
    Y = {};
    X = {};

    if GEOHAR
        for i = 23:dataset_length
            obs = time_series(:,:,i-22:i-1);
            obs_squeezed = squeeze(num2cell(obs, [1 2]));
            Y{end+1} = time_series(:,:,i);
            X{end+1} = diagblockhar(obs_squeezed{:}, method);
        end
    else
        for i = N_LAGS+1:dataset_length
            obs = time_series(:,:,i-N_LAGS:i-1);
            obs_squeezed = squeeze(num2cell(obs, [1 2]));
            Y{end+1} = time_series(:,:,i);
            X{end+1} = diagblock(obs_squeezed{:});
        end
    end
end


