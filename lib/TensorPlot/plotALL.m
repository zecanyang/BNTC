function [] = plotALL(varargin)

nTensor = 4;
% scrsz = get(0,'ScreenSize');


strName = cell(1,4);
strName{1} = 'Ground Truth';
strName{2} = 'Recovery';
strName{3} = 'Noise';
strName{4} = 'Observation';

parThresh = zeros(1,nTensor);
parThresh(3) =0.8;
parThresh(4) =0.7;

if length(varargin{8})==1
    for i = 1:nTensor
        subplot(2,nTensor,i);
        voxel3(abs(varargin{i}),'thresh',parThresh(i), 'degree',5);
        title(strName{i},'Fontsize',15,'FontWeight','bold');
        xlabel(''); ylabel(''); zlabel('');
        set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
        box on;
    end
else
    subplot(2,nTensor,2);
    voxel3(abs(varargin{2}),'thresh',parThresh(2), 'degree',5);
    title(strName{2},'Fontsize',15,'FontWeight','bold');
    xlabel(''); ylabel(''); zlabel('');
    set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
    box on;
end



subplot(2,4,5); imagesc(varargin{5}), colorbar; title('Factor-1','Fontsize',15); 
subplot(2,4,6); imagesc(varargin{6}), colorbar; title('Factor-2','Fontsize',15); 
subplot(2,4,7); imagesc(varargin{7}), colorbar; title('Factor-3','Fontsize',15);  
subplot(2,4,8); plot(varargin{8}, '-*R','LineWidth',1);title('Model Error','Fontsize',15); xlabel('Iteration'); ylabel('RMSE');

