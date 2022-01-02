function hf = formatPlots(hf,sz,xl,yl,varargin)

N = length(hf.Children);


linkaxes(hf.Children)
for ii = 1:N
   hf.Children(ii).FontSize = sz;
   hf.Children(ii).XLabel.String = xl;
   hf.Children(ii).YLabel.String =yl;
   hf = gcf;
   
   % rotation
    if nargin > 4 && ~isempty(varargin{1})
        hf.Children(ii).XTickLabelRotation = varargin{1};
    end
    
    % Ctick labels
    if nargin > 5 && ~isempty(varargin{2})
        hf.Children(ii).XTickLabel = varargin{2};
    end
    
    if nargin > 6 && ~isempty(varargin{3})
        hf.Children(ii).YTickLabel = varargin{3};
    end
   
end

