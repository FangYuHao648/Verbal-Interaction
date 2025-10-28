% savefolder={'B2B_WF','B2B_NF','F2F_WF','F2F_NF','TO_WF','TO_NF','TWE_WF','TWE_NF','VoC_WF','VoC_NF','ViC_WF','ViC_NF'};
savefolder={'F2F_WF','F2F_NF'};
locfile = 'D:\桌面\Matlab Working path\data&code&experiment\1.maincal&plot_code\32loc.loc';
correct=0;
erro=cell(30,1);
ISC_delta_all =zeros(12,30);
ISC_theta_all =zeros(12,30);
ISC_alpha_all =zeros(12,30);
ISC_beta_all  =zeros(12,30);
ISC_gamma_all =zeros(12,30);
% isc_mutisub_mutiband=zeros(5,12);

% 总文件夹数
total_folders = length(savefolder);
total_files = 0;

% 先统计总文件数
for h = 1:total_folders
    path = ['D:\桌面\Matlab Working path\data&code&experiment\8.vi_scalp_segment_data\' savefolder{h} '\*.set'];
    namelist = dir(path);
    total_files = total_files + length(namelist);
end

% 创建waitbar
hWaitbar = waitbar(0, 'Processing EEG data...');
processed_files = 0;
for h = 1:total_folders
    path = ['D:\桌面\Matlab Working path\data&code&experiment\8.vi_scalp_segment_data\' savefolder{h} '\*.set'];
    namelist = dir(path);
    len = length(namelist);
    filename=cell(1,len);
    for k = 1:len
        tempfilepath = [ path(1:end-5) namelist(k).name];
        filename{k}= tempfilepath;
    end
    folder = dir(fullfile(path));
    folder = {folder.name};

     %% 计算两个人ISC
     count=1;
    for i=1:2:len
         % if i==100
         %     count=count-2;
         %     continue;
         % else
            count=i;
            % 导入被试A/被试B
            EEG_A = pop_loadset(filename{i});
            EEG_B = pop_loadset(filename{i+1});
            % 带通滤波
            EEG_A_delta = pop_eegfilt(EEG_A,1,3,[],[],1);
            EEG_A_theta = pop_eegfilt(EEG_A,4,7,[],[],1);
            EEG_A_alpha = pop_eegfilt(EEG_A,8,12,[],[],1);
            EEG_A_beta = pop_eegfilt(EEG_A,13,30,[],[],1);
            % EEG_A_gamma = pop_eegfilt(EEG_A,31,49,[],[],1);
            
            EEG_B_delta = pop_eegfilt(EEG_B,1,3,[],[],1);
            EEG_B_theta = pop_eegfilt(EEG_B,4,7,[],[],1);
            EEG_B_alpha = pop_eegfilt(EEG_B,8,12,[],[],1);
            EEG_B_beta = pop_eegfilt(EEG_B,13,30,[],[],1);
            % EEG_B_gamma = pop_eegfilt(EEG_B,31,49,[],[],1);
            
            %% delta频段
            DATA_A_delta = double(EEG_A_delta.data');
            DATA_B_delta = double(EEG_B_delta.data');
            s1 = size(DATA_A_delta,1);
            s2 = size(DATA_B_delta,1);
            if s1~=s2
                min_s = min(s1,s2);
                DATA_A_delta = DATA_A_delta(s1-min_s+1:s1,:);
                DATA_B_delta = DATA_B_delta(s2-min_s+1:s2,:);
                outputline=['delta频带',savefolder{h},'条件下','sub',num2str((i+1)/2),'A/B长度不匹配'];
                erro(correct+1)=cellstr(outputline);
                correct=correct+1;
            end
            DATA_delta = cat(3,DATA_A_delta,DATA_B_delta);
            EEGData_delta.X = DATA_delta;
            EEGData_delta.fs = 250;
            EEGData_delta.badchannels = {};
            EEGData_delta.eogchannels = [];            
            save('D:\桌面\Matlab Working path\data&code&experiment\EEGData_delta','-struct','EEGData_delta')
            clear EEG_A_delta EEG_B_delta DATA_A_delta DATA_B_delta DATA_delta EEGData_delta 
            
            % datafile = 'D:\桌面\Matlab Working path\EEGData_delta.mat';
            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\EEGData_delta.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig =0;
            [ISC_delta,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
            ISC_delta_all(h,(i+1)/2) = ISC_delta(1)+ISC_delta(2)+ISC_delta(3);
            %ISC_delta_all(h,(count+1)/2) = ISC_delta(1);
            %% theta频段
            DATA_A_theta = double(EEG_A_theta.data');
            DATA_B_theta = double(EEG_B_theta.data');
            s1 = size(DATA_A_theta,1);
            s2 = size(DATA_B_theta,1);
            if s1~=s2
                min_s = min(s1,s2);
                DATA_A_theta = DATA_A_theta(s1-min_s+1:s1,:);
                DATA_B_theta = DATA_B_theta(s2-min_s+1:s2,:);
                outputline=['theta频带',savefolder{h},'条件下','sub',num2str((i+1)/2),'A/B长度不匹配'];
                erro(correct+1)=cellstr(outputline);
                correct=correct+1;
            end
            DATA_theta = cat(3,DATA_A_theta,DATA_B_theta);
            EEGData_theta.X = DATA_theta;
            EEGData_theta.fs = 250;
            EEGData_theta.badchannels = {};
            EEGData_theta.eogchannels = [];            
            save('D:\桌面\Matlab Working path\data&code&experiment\EEGData_theta','-struct','EEGData_theta')
            clear EEG_A_theta EEG_B_theta DATA_A_theta DATA_B_theta DATA_theta EEGData_theta 

            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\EEGData_theta.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 0;
            [ISC_theta,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);  
            ISC_theta_all(h,(i+1)/2) = ISC_theta(1)+ISC_theta(2)+ISC_theta(3);
            %ISC_theta_all(h,(count+1)/2) = ISC_theta(1);
            %% alpha频段
            DATA_A_alpha = double(EEG_A_alpha.data');
            DATA_B_alpha = double(EEG_B_alpha.data');
            s1 = size(DATA_A_alpha,1);
            s2 = size(DATA_B_alpha,1);
            if s1~=s2
                min_s = min(s1,s2);
                DATA_A_alpha = DATA_A_alpha(s1-min_s+1:s1,:);
                DATA_B_alpha = DATA_B_alpha(s2-min_s+1:s2,:);
                outputline=['alpha频带',savefolder{h},'条件下','sub',num2str((i+1)/2),'A/B长度不匹配'];
                erro(correct+1)=cellstr(outputline);
                correct=correct+1;
            end
            DATA_alpha = cat(3,DATA_A_alpha,DATA_B_alpha);
            EEGData_alpha.X = DATA_alpha;
            EEGData_alpha.fs = 250;
            EEGData_alpha.badchannels = {};
            EEGData_alpha.eogchannels = [];            
            save('D:\桌面\Matlab Working path\data&code&experiment\EEGData_alpha','-struct','EEGData_alpha')
            clear EEG_A_alpha EEG_B_alpha DATA_A_alpha DATA_B_alpha DATA_alpha EEGData_alpha 

            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\EEGData_alpha.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 1;
            [ISC_alpha,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
            ISC_alpha_all(h,(i+1)/2) = ISC_alpha(1)+ISC_alpha(2)+ISC_alpha(3);
            %ISC_alpha_all(h,(count+1)/2) = ISC_alpha(1);
            % colormap(nclCM(141));
            % set(gcf,'color',[235,241,223]./255);
            %% beta频段
            DATA_A_beta = double(EEG_A_beta.data');
            DATA_B_beta = double(EEG_B_beta.data');
            s1 = size(DATA_A_beta,1);
            s2 = size(DATA_B_beta,1);
            if s1~=s2
                min_s = min(s1,s2);
                DATA_A_beta = DATA_A_beta(s1-min_s+1:s1,:);
                DATA_B_beta = DATA_B_beta(s2-min_s+1:s2,:);
                outputline=['beta频带',savefolder{h},'条件下','sub',num2str((i+1)/2),'A/B长度不匹配'];
                erro(correct+1)=cellstr(outputline);
                correct=correct+1;
            end
            DATA_beta = cat(3,DATA_A_beta,DATA_B_beta);
            EEGData_beta.X = DATA_beta;
            EEGData_beta.fs = 250;
            EEGData_beta.badchannels = {};
            EEGData_beta.eogchannels = [];            
            save('D:\桌面\Matlab Working path\data&code&experiment\EEGData_beta','-struct','EEGData_beta')
            clear EEG_A_beta EEG_B_beta DATA_A_beta DATA_B_beta DATA_beta EEGData_beta 

            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\EEGData_beta.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 0;
            [ISC_beta,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
            ISC_beta_all(h,(i+1)/2) = ISC_beta(1)+ISC_beta(2)+ISC_beta(3);
            %ISC_beta_all(h,(count+1)/2) = ISC_beta(1);
            % %% gamma频段
            % DATA_A_gamma = double(EEG_A_gamma.data');
            % DATA_B_gamma = double(EEG_B_gamma.data');
            % s1 = size(DATA_A_gamma,1);
            % s2 = size(DATA_B_gamma,1);
            % if s1~=s2
            %     min_s = min(s1,s2);
            %     DATA_A_gamma = DATA_A_gamma(s1-min_s+1:s1,:);
            %     DATA_B_gamma = DATA_B_gamma(s2-min_s+1:s2,:);
            %     outputline=['gamma频带',savefolder{h},'条件下','sub',num2str((i+1)/2),'A/B长度不匹配'];
            %     erro(correct+1)=cellstr(outputline);
            %     correct=correct+1;
            % end
            % DATA_gamma = cat(3,DATA_A_gamma,DATA_B_gamma);
            % EEGData_gamma.X = DATA_gamma;
            % EEGData_gamma.fs = 250;
            % EEGData_gamma.badchannels = {};
            % EEGData_gamma.eogchannels = [];            
            % 
            % clear EEG_A_gamma EEG_B_gamma DATA_A_gamma DATA_B_gamma DATA_gamma EEGData_gamma 
            % 
            % datafile = 'D:\桌面\Matlab Working path\EEGData_gamma.mat';
            % gamma = 0.4;
            % Nsec = 1;
            % plotfig = 0;
            % [ISC_gamma,ISC_persubject,ISC_persecond,W,A] = isceeg2(datafile,locfile,gamma,Nsec,plotfig);
            % ISC_gamma_all(h,(i+1)/2) = ISC_gamma(1)+ISC_gamma(2)+ISC_gamma(3);
            % ISC_gamma_all(h,(count+1)/2) = ISC_gamma(1);
            
            % 更新waitbar
            processed_files = processed_files + 2; % 每次处理两个文件
            waitbar(processed_files / total_files, hWaitbar, sprintf('Processing %d of %d files', processed_files, total_files));
             count=count+2;
         end
    end
%
% 关闭waitbar
close(hWaitbar);