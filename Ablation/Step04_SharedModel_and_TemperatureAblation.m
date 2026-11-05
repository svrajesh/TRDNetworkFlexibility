%% Step04 Compute Effective Temperature
thresholds = 50; %round([25,30,35,37.5,40,45,50,55,60,62.5,65,70,75]);
regs = 1:8;
bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};
% pmem feature folder and binary data dir
% Loop thresh and band
for thresh = 1:length(thresholds)
    for ablation = 1:length(regs)
        for b = 1:length(bands)
            binarydatadir = dir(betterLoadSave([USER_DIR,'/WeightedReconBinarizedData/Ablation/',task,'/',pt,'/','Reg_',num2str(ablation),'/binarized_',bands{b},'/']));
            binarydatadir = binarydatadir(~startsWith({binarydatadir.name},'.'));
            load(betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/pMEM/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_result.mat']));
            energyfolder = betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/Energy/reg_',num2str(ablation),'_Energy_features_',num2str(thresholds(thresh)),'/',bands{b},'/']);
            nodeNum = length(results(1).parameters.h);
            allhvals = zeros(length(results),nodeNum);
            allJvals = zeros(length(results),nodeNum,nodeNum);
            allrecs = {};
            %% Gather the data
            for rec = 1:length(results)
                fn = results(rec).filename;
                [fp,fn] = fileparts(fn);
                binarydatafile = binarydatadir(find(contains({binarydatadir.name},fn)));
                load(betterLoadSave(fullfile(binarydatafile.folder,binarydatafile.name)));
                allrecs{rec} = double(bitpack_helper(binarydata,'decode'));
                % Process the gathered data
                allhvals(rec, :) = results(rec).parameters.h; % Store h values
                allJvals(rec, :, :) = results(rec).parameters.J; % Store J values
            end
            
            %%
            % h_star = mean(allhvals,1);
            % J_star = squeeze(mean(allJvals,1));
            pairwiseStructureCorrs = structureCorrelations(allhvals,allJvals);
            if ~isfolder(betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/Figures/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/']))
                mkdir(betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/Figures/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/']));
            end
            figname = betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/Figures/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_StructureCorrelations.png']);
            
            figure;
            set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]); % Set figure background to white and size to 800x600
            
            mask = triu(true(size(pairwiseStructureCorrs)),1);
            
            
            % Use the mask to extract upper triangle elements
            upperTri = pairwiseStructureCorrs(mask);
            
            histogram(upperTri(:), 'Normalization', 'probability'); % Plot histogram
            xlabel('Correlation Coefficient');
            ylabel('Probability');
            title('Histogram of Pairwise Structure Correlations');
            saveas(gcf, figname); % Save figure as figname
            
            
            
            %%
            revisedparams = pfunc_02_Inferrer_ML_MultiRecording(allrecs,ablation);
            beta = revisedparams.beta;
            temperature = 1./revisedparams.beta;
            for rec = 1:length(results)
                results(rec).parameters.h_shared = revisedparams.h_star;
                results(rec).parameters.J_shared = revisedparams.J_star;
                results(rec).parameters.effectiveTemperature = temperature(rec);
                % model eval
                h_eff = beta(rec)*revisedparams.h_star;
                J_eff = beta(rec)*revisedparams.J_star;
                [probN, prob1, prob2, rD, r, D1, D2] = pfunc_03_Accuracy(h_eff, J_eff, allrecs{rec});
                mets.d = rD;
                mets.KLD1 = D1;
                mets.KLD2 = D2;
                mets.r = r;
                R2imem = corrcoef(probN(:), prob1(:));
                R2imem = R2imem(1,2)^2;
                R2pmem = corrcoef(probN(:), prob2(:));
                R2pmem = R2pmem(1,2)^2;
                mets.R2imem = R2imem;
                mets.R2pmem = R2pmem;
                results(rec).metrics.sharedMetrics = mets;
            end
            
            %% shared energy features
            h = revisedparams.h_star;
            J = revisedparams.J_star;
            E = mfunc_Energy(h,J);
            nodeNumber = length(h);
            vectorList = mfunc_VectorList(nodeNumber);
            [LocalMinIndex, BasinGraph, AdjacentList] = mfunc_LocalMin(nodeNumber, E);
            basicSharedFeatures = struct;
            basicSharedFeatures.basinInds = LocalMinIndex;
            Nmin = length(LocalMinIndex);
            basicSharedFeatures.basinCount = Nmin;
            close all
            basins = {};
            basinSizes = [];
            
            basinEnergies = E(LocalMinIndex);
            GMenergy = min(basinEnergies);
            E_diff = basinEnergies-GMenergy;
            if length(E_diff) >1
                E_diff = mean(E_diff(E_diff~=0));
            else
                E_diff = NaN;
            end
            basicSharedFeatures.E_diff = E_diff;
            
            
            for m = 1:Nmin
                % Basin
                bM = vectorList(:,LocalMinIndex(m));
                bM = (bM+1)/2;
                basins{m} = sprintf('%d', bM);
                
                % Basin Size
                basinSizes(m) = sum(BasinGraph(:,3) == LocalMinIndex(m));
                
            end
            basicSharedFeatures.basins = basins';
            basicSharedFeatures.basinSizes = basinSizes';
            basicSharedFeatures.basin_sd = std(basinSizes);
            
            
            
            %% Calculate Disconnectivity by Dijkstra
            % [~, ~] = mfunc_DisconnectivityGraph(E, LocalMinIndex, AdjacentList);
            % close all %% next iter will save these
            
            for rec = 1:length(results)
                fn = results(rec).filename;
                [fp,fn] = fileparts(fn);
                enfn = fullfile(energyfolder,['Energy_',fn,'.mat']);
                load(enfn);
                
                % State Number
                StateNumber = features.StateNumber;
                % Basin Number
                BasinNumber = mfunc_GetBasinNumber(StateNumber, BasinGraph, LocalMinIndex);
                % Dynamics
                Dynamics = mfunc_GetDynamics(BasinNumber, Nmin);
                features.sharedFeatures = basicSharedFeatures;
                features.sharedFeatures.StateNumber = StateNumber;
                features.sharedFeatures.BasinNumber = BasinNumber;
                features.sharedFeatures.Dynamics = Dynamics;
                save(enfn,"features",'-v7.3');
                
            end
            save(betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/pMEM/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_result.mat']),'results','-v7.3');
        end
    end
end



