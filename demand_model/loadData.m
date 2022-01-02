function [datTrain,varargout] = loadData(file,sheet,varargin)

[num,~,raw] = xlsread(file,sheet,'','basic');

timestamps = datetime(num(:,1),'ConvertFrom','Excel');
if nargin > 2 && ~isempty(varargin{1}) % training and testing data are different
    endTime = varargin{1};
    includeTraining = varargin{2};
    
    indTrain = find(timestamps<=endTime);
    indTest = find(timestamps>endTime);
    if includeTraining
        indTest = [indTrain;indTest];
    end
    
    % training data
    datTrain.timestamps = timestamps(indTrain);
    datTrain.kW = num(indTrain,2);
    
    % testing data
    datTest.timestamps = timestamps(indTest);
    datTest.kW = num(indTest,2);
    varargout{1} = datTest;
        
else
    datTrain.timestamps = timestamps;
    datTrain.kW = num(:,2);
end