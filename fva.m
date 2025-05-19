models = [];

thresholds = [0.85, 0.90, 0.95];
base_directory = "C:\Users\A\Documents\MSc_Thesis";
ecGEM_directory = fullfile(base_directory, "ecGEMs");

function fva_result = run_fva(Model, threshold, inputDir, outputDir, uniquename)
    sol_opt = solveLP(Model);
    objFlux = sol_opt.x(Model.c>0);
    bioModel = setParam(Model,'lb',Model.rxns(Model.c>0),threshold*objFlux);
    bioModel = setParam(bioModel,'ub',Model.rxns(Model.c>0),objFlux);
    model = bioModel;
    pos_cols = find(all(model.S >= 0, 1) & any(model.S ~= 0, 1));
    neg_cols = find(all(model.S <= 0, 1) & any(model.S ~= 0, 1));
    exc_rxns = [pos_cols, neg_cols];
    exc_rxns = sort(exc_rxns);

    % Initialize arrays with correct size
    n_exc = length(exc_rxns);
    LB_exc = zeros(n_exc, 1);
    UB_exc = zeros(n_exc, 1);
    rxn_names = cell(n_exc, 1);
    direction = cell(n_exc, 1);
    excluded_counter = 0;

    % Perform FVA only on exchange reactions
    for idx = 1:n_exc
        i = exc_rxns(idx);
        if startsWith(model.rxns{i},'sink_') || startsWith(model.rxns{i},'prot_pool_exchange') 
            excluded_counter = excluded_counter + 1;
            continue
        end
        rxn_names{idx} = model.rxns{i};  % Store reaction names in order
        iRxn = model.rxns(i);
        if ismember(i, pos_cols)
            direction{idx} = 'excreted';
        else
            direction{idx} = 'absorbed';
        end
        iModel = setParam(model,'obj',iRxn,1);
        iSol = solveLP(iModel);
        UB_exc(idx) = iSol.x(ismember(iModel.rxns,iRxn));
        iModel = setParam(model,'obj',iRxn,-1);
        iSol = solveLP(iModel);
        LB_exc(idx) = iSol.x(ismember(iModel.rxns,iRxn));
    end

    % Trim arrays to actual size
    rxn_names = rxn_names(1:n_exc-excluded_counter);
    LB_exc = LB_exc(1:n_exc-excluded_counter);
    UB_exc = UB_exc(1:n_exc-excluded_counter);
    direction = direction(1:n_exc-excluded_counter);

    % Create table with matched dimensions
    fva_result = table(rxn_names, LB_exc, UB_exc, direction, ...
        'VariableNames', {'Reaction', 'Min_Flux', 'Max_Flux', 'Direction'});
    cd (outputDir)
    writetable(fva_result, sprintf('%s_fva_results%.2f.csv', uniquename, threshold));
    fprintf('FVA results saved to: %s/%s\n', outputDir, sprintf('%s_fva_results%.2f.csv', uniquename, threshold));
    cd (inputDir)
end

function filtered_fva = filter_fva(fva_file, keep_zeros)
    % Read the original FVA results into a table
    T = readtable(fva_file);

    % Preallocate a logical mask
    to_keep = false(height(T),1);

    for i = 1:height(T)
        LB = T.Min_Flux(i);
        UB = T.Max_Flux(i);

        if keep_zeros
            % Keep lines where both LB & UB are ≥0 or both ≤0
            if ((LB >= 0 && UB >= 0) || (LB <= 0 && UB <= 0)) && ~(LB == 0 && UB == 0)
                to_keep(i) = true;
            end
        else
            % Exclude any zero; keep only strictly >0 or strictly <0
            if (LB > 0 && UB > 0) || (LB < 0 && UB < 0)
                to_keep(i) = true;
            end
        end
    end

    % Create a filtered table
    filtered_fva = T(to_keep, :);

    % Build an output filename to indicate zero behavior
    if keep_zeros
        outFile = strrep(fva_file, '.csv', '_filtered_zeros.csv');
    else
        outFile = strrep(fva_file, '.csv', '_filtered_nozeros.csv');
    end

    % Write the filtered table
    writetable(filtered_fva, outFile);
    fprintf('Filtered results saved to: %s\n', outFile);
end

[samplePath, sampleName, ~] = fileparts(pwd);
[~, cohort, ~] = fileparts(samplePath);
inputDir = fullfile(samplePath, sampleName);
outputDir = fullfile(base_directory, 'fva_results', cohort, sampleName);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end


% Loop through all elements in the models list and run FVA for each
for k = 1:length(models)
    current_model = strcat(models(k), '.yml');
    % Load the ecModel based on the current model (without preset names)
    for i = 1:length(thresholds)
        out_file_name = strcat(current_model, '_fva_results', sprintf('%.2f', thresholds(i)), '.csv');
        out_file_path = fullfile(outputDir, out_file_name);
        if ~isfile(out_file_path)
            fprintf('Creating files for %s, threshold %.2f\n', current_model, thresholds(i))
            if ~exist('ecModel', 'var')
                % If model has not been loaded yet, load it
                ecModel = readYAMLmodel(current_model);    
            end
            run_fva(ecModel, thresholds(i), inputDir, outputDir, current_model);
        else
            fprintf('Files already created for %s, threshold %.2f\n', current_model, thresholds(i))
        end
    clearvars ecModel
    end
end

cd (outputDir)

file_list = dir('*Ready_ecGEM.yml_fva_results0.??.csv');
for i = 1:length(file_list)
    filter_fva(file_list(i).name, false);
    filter_fva(file_list(i).name, true);
end

disp('Results saved successfully to:')
disp(pwd)
%exit