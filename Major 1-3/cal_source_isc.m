savename={'B2B_F','B2B_N','F2F_F','F2F_N','TO_F','TO_N','TWE_F','TWE_N','VoC_F','VoC_N','ViC_F','ViC_N'};
%savename={'ViC_N'};

locfile = 'D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data';
locfile1 = 'D:\桌面\Matlab Working path\data&code&experiment\1.maincal&plot_code\32loc.loc';
for h=1:12   
  % 从子文件夹收集所有 .mat 文件，并组装成 DATA.(Var###).Value
L = dir(fullfile(locfile, savename{h}, '*.mat'));
if isempty(L)
    error('目录不存在或无 .mat 文件：%s', fullfile(locfile, savename{h}));
end
% 按文件名排序（你的文件名有前导0，字母序即自然序）
[~,ix] = sort({L.name}); 
L = L(ix);

DATA = struct();
for k = 1:numel(L)
    fpath = fullfile(L(k).folder, L(k).name);
    S = load(fpath);

    % 优先取 struct.Value；否则取第一个二维数值矩阵
    picked = [];
    fn = fieldnames(S);
    for t = 1:numel(fn)
        v = S.(fn{t});
        if isstruct(v) && isfield(v,'Value') && isnumeric(v.Value) && ismatrix(v.Value)
            picked = v.Value; 
            break;
        end
    end
    if isempty(picked)
        for t = 1:numel(fn)
            v = S.(fn{t});
            if isnumeric(v) && ismatrix(v)
                picked = v; 
                break;
            end
        end
    end
    if isempty(picked)
        error('文件中未找到二维数值矩阵：%s', fpath);
    end

    % 填到 DATA 中，字段名 Var001、Var002...
    DATA.(sprintf('Var%03d', k)).Value = picked;
end

filename = fieldnames(DATA);
l = numel(filename);
    for i=1:2:60
            
            DATA_A = DATA.(filename{i}).Value;
            DATA_B = DATA.(filename{i+1}).Value;
            DATA_A=DATA_A';
            DATA_B=DATA_B';
            fs = 250;
            %%delta
            DATA_A_delta = bp_band(DATA_A, fs, 1, 3);
            DATA_B_delta = bp_band(DATA_B, fs, 1, 3);
            DATA_delta = cat(3,DATA_A_delta,DATA_B_delta);
            EEGData_delta.X = DATA_delta;
            EEGData_delta.eogchannels = [];
            EEGData_delta.badchannels = {};
            EEGData_delta.fs = 250;
            save('D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_delta','-struct','EEGData_delta')
            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_delta.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 0;
            [ISC_delta,W,A] = isceeg2(datafile,locfile1,gamma,Nsec,plotfig);
            ISC_delta_sum3((i+1)/2,h) = ISC_delta(1)+ISC_delta(2)+ISC_delta(3);
            ISC_delta_1st((i+1)/2,h) = ISC_delta(1);
            
            %%theta
            DATA_A_theta = bp_band(DATA_A, fs, 4, 7);
            DATA_B_theta = bp_band(DATA_B, fs, 4, 7);
            DATA_theta = cat(3,DATA_A_theta,DATA_B_theta);
            EEGData_theta.X = DATA_theta;
            EEGData_theta.eogchannels = [];
            EEGData_theta.badchannels = {};
            EEGData_theta.fs = 250;
            save('D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_theta','-struct','EEGData_theta')
            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_theta.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 0;
            [ISC_theta,W,A] = isceeg2(datafile,locfile1,gamma,Nsec,plotfig);
            ISC_theta_sum3((i+1)/2,h) = ISC_theta(1)+ISC_theta(2)+ISC_theta(3);
            ISC_theta_1st((i+1)/2,h) = ISC_theta(1);

            %%alpha
            DATA_A_alpha = bp_band(DATA_A, fs, 8, 12);
            DATA_B_alpha = bp_band(DATA_B, fs, 8, 12);
            DATA_alpha = cat(3,DATA_A_alpha,DATA_B_alpha);
            EEGData_alpha.X = DATA_alpha;
            EEGData_alpha.eogchannels = [];
            EEGData_alpha.badchannels = {};
            EEGData_alpha.fs = 250;
            save('D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_alpha','-struct','EEGData_alpha')
            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_alpha.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 0;
            [ISC_alpha,W,A] = isceeg2(datafile,locfile1,gamma,Nsec,plotfig);
            ISC_alpha_sum3((i+1)/2,h) = ISC_alpha(1)+ISC_alpha(2)+ISC_alpha(3);
            ISC_alpha_1st((i+1)/2,h) = ISC_alpha(1);

            %%beta
            DATA_A_beta = bp_band(DATA_A, fs, 13, 30);
            DATA_B_beta = bp_band(DATA_B, fs,13, 30);
            DATA_beta = cat(3,DATA_A_beta,DATA_B_beta);
            EEGData_beta.X = DATA_beta;
            EEGData_beta.eogchannels = [];
            EEGData_beta.badchannels = {};
            EEGData_beta.fs = 250;
            save('D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_beta','-struct','EEGData_beta')
            datafile = 'D:\桌面\Matlab Working path\data&code&experiment\9.brainstrom_vi_source_data\EEGData_beta.mat';
            gamma = 0.4;
            Nsec = 1;
            plotfig = 0;
            [ISC_beta,W,A] = isceeg2(datafile,locfile1,gamma,Nsec,plotfig);
            ISC_beta_sum3((i+1)/2,h) = ISC_beta(1)+ISC_beta(2)+ISC_beta(3);
            ISC_beta_1st((i+1)/2,h) = ISC_beta(1);
           
    end
end
% function Xf = bp_band(X, fs, f1, f2)
% % X: [time x channels]
% % fs: 采样率
% % f1,f2: 带通边界
%     [b,a] = butter(4, [f1 f2]/(fs/2), 'bandpass');  % 4阶足够稳健
%     Xf = filtfilt(b, a, double(X));                 % 零相位
% end
function Xf = bp_band(X, fs, f1, f2)
% X: [time x channels]  fs: 采样率  f1,f2: 带通边界(Hz)
    EEG = struct();
    EEG.data   = double(X)';                 % [chan x time]
    EEG.srate  = fs;
    EEG.nbchan = size(EEG.data,1);
    EEG.pnts   = size(EEG.data,2);
    EEG.trials = 1;
    EEG.event  = [];                         % 必要占位
    EEG.urevent= [];

    if exist('pop_eegfiltnew','file')
        EEG = pop_eegfiltnew(EEG,'locutoff',f1,'hicutoff',f2,'plotfreqz',0);
    else
        EEG = pop_eegfilt(EEG, f1, f2, [], [], 1);   % 1 = 零相位
    end

    Xf = double(EEG.data)';                  % 还原为 [time x channels]
end
