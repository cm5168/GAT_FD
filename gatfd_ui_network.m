classdef gatfd_ui_network < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GATFD_network_UIFigure         matlab.ui.Figure
        Button_loadfiles               matlab.ui.control.Button
        LoadedFIlesTextAreaLabel       matlab.ui.control.Label
        TextArea_LoadedFiles           matlab.ui.control.TextArea
        UIAxes                         matlab.ui.control.UIAxes
        LoadTemporalMaskButton         matlab.ui.control.Button
        CalculateNetworkPropertiesButton  matlab.ui.control.Button
        NetworkAverageGlobalMeasuresPanel  matlab.ui.container.Panel
        GlobalEfficiencyCheckBox       matlab.ui.control.CheckBox
        LocalEfficiencyCheckBox        matlab.ui.control.CheckBox
        ClusteringCoefficientCheckBox  matlab.ui.control.CheckBox
        DegreeCheckBox                 matlab.ui.control.CheckBox
        SmallWorldCoefficientCheckBox  matlab.ui.control.CheckBox
        ModularityCoefficientCheckBox  matlab.ui.control.CheckBox
        NormalizedClusteringCoefficientCheckBox  matlab.ui.control.CheckBox
        NormalizedPathLengthCheckBox   matlab.ui.control.CheckBox
        TransitivityCoefficientCheckBox  matlab.ui.control.CheckBox
        AssortativityCoefficientCheckBox  matlab.ui.control.CheckBox
        CharacteristicPathLengthCheckBox  matlab.ui.control.CheckBox
        NodalMeasuresPanel             matlab.ui.container.Panel
        ListBox                        matlab.ui.control.ListBox
        SelectNodesButton              matlab.ui.control.Button
        NodalGlobalEfficiencyCheckBox  matlab.ui.control.CheckBox
        NodalLocalEfficiencyCheckBox   matlab.ui.control.CheckBox
        NodalClusteringCoefficientCheckBox  matlab.ui.control.CheckBox
        NodalDegreeCheckBox            matlab.ui.control.CheckBox
        BetweennessCentralityCheckBox  matlab.ui.control.CheckBox
        ThresholdingPanel              matlab.ui.container.Panel
        ThresholdingMethodDropDownLabel  matlab.ui.control.Label
        ThresholdingMethodDropDown     matlab.ui.control.DropDown
        ThresholdRangeLabel            matlab.ui.control.Label
        ThresholdUpperEditField        matlab.ui.control.NumericEditField
        ThresholdLowerEditField        matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel_2      matlab.ui.control.Label
        ThresholdStepLabel             matlab.ui.control.Label
        ThresholdStepEditField         matlab.ui.control.NumericEditField
        AbsoluteCheckBox               matlab.ui.control.CheckBox
        ParallelCheckBox               matlab.ui.control.CheckBox
    end

    
    properties (Access = private)
        gatn_setting % Description
        gatn_measure
    end
    
    methods (Access = private)
        % Dialog for choosing nodes
        function results = choose_node(app)
            d=uifigure("Position",[100 100 292 714],'Name','Select nodes for calculation');
            
            % Create SelectNodesforCalculationListBoxLabel
            SelectNodesforCalculationListBoxLabel = uilabel(d);
            SelectNodesforCalculationListBoxLabel.HorizontalAlignment = 'right';
            SelectNodesforCalculationListBoxLabel.Position = [69 666 156 22];
            SelectNodesforCalculationListBoxLabel.Text = 'Select Nodes for Calculation';

            % Create SelectNodesforCalculationListBox
            SelectNodesforCalculationListBox = uilistbox(d);
            SelectNodesforCalculationListBox.Multiselect = 'on';
            SelectNodesforCalculationListBox.Position = [31 79 231 588];
            SelectNodesforCalculationListBox.Items = app.gatn_setting.node_list;
            SelectNodesforCalculationListBox.Value = {};

            % Create DoneButton
            DoneButton = uibutton(d, 'push',...
                                  'Position',[162 34 100 22],...
                                  'Text','Done',...
                                  'ButtonPushedFcn',@fnc_done);
   
            function fnc_done(src,event)
                app.gatn_setting.node_list_selected=SelectNodesforCalculationListBox.Value;
                app.ListBox.Items = app.gatn_setting.node_list_selected;
                close(d);
            end
        end

        % Network Global Efficiency
        function results = calc_global_efficiency(~,temp_net)
            results = efficiency_bin(temp_net,0);
            return
        end
        
        % Network Local Efficiency
        function results = calc_network_local_efficiency(~,temp_net)
            results = mean(efficiency_bin(temp_net,1));
            return
        end
        
        % Network Clustering Coefficient
        function results = calc_network_clustering_coefficient(~,temp_net)
            results = mean(clustering_coef_bu(temp_net));
            return
        end
        
        % Network Averaged Degree
        function results = calc_network_average_degree(~,temp_net)
            results = mean(degrees_und(temp_net));
            return
        end
                
        % Network characteristic path length
        function results = calc_network_characteristic_path(~,temp_net)
            results = charpath(distance_bin(temp_net),0,0);
            return
        end
        
        % Network SW Normalized CC
        function results = calc_sw_norm_cc(~,temp_net)
            num_edges=round(sum(temp_net(:))/2);
            num_nodes=length(temp_net);
            cc=0;
            for i=1:20
                temp_rand=makerandCIJ_und(num_nodes,num_edges);
                cc=cc+mean(clustering_coef_bu(temp_rand));
            end
            cc=cc/20;
            results = mean(clustering_coef_bu(temp_net))/cc;
            return
        end
        
        % Network SW Normalized Path Length
        function results = calc_sw_norm_pl(~,temp_net)
            num_edges=round(sum(temp_net(:))/2);
            num_nodes=length(temp_net);
            pl=0;
            for i=1:20
                temp_rand=makerandCIJ_und(num_nodes,num_edges);
                pl=pl+charpath(distance_bin(temp_rand),0,0);
            end
            pl=pl/20;
            results = charpath(distance_bin(temp_net),0,0)/pl;
            return
        end
        
        % Network SW Coefficient
        function results = calc_sw_coefficient(~,temp_net)
            num_edges=round(sum(temp_net(:))/2);
            num_nodes=length(temp_net);
            pl=0;
            cc=0;
            for i=1:20
                temp_rand=makerandCIJ_und(num_nodes,num_edges);
                pl=pl+charpath(distance_bin(temp_rand),0,0);
                cc=cc+mean(clustering_coef_bu(temp_rand));
            end
            pl=pl/20;
            cc=cc/20;
            results = (mean(clustering_coef_bu(temp_net))/cc)/(charpath(distance_bin(temp_net),0,0)/pl);
            return
        end
        
        % Network Modularity
        function results = calc_modularity(~,temp_net)
            [~,results]=modularity_und(temp_net);
            return
        end
        
        % Network Transitivity
        function results = calc_transitivity(~,temp_net)
            results = transitivity_bu(temp_net);
            return
        end
        
        % Network Assortativity
        function results = calc_assortativity(~,temp_net)
            results = assortativity_bin(temp_net,0);
            return
        end
        
        % Nodal Global Efficiency / Nodal Efficiency
        function results = calc_nodal_efficiency(~,temp_net)    % From Brain Connectivity Toolbox
            n=length(temp_net);     %number of nodes
            temp_net(1:n+1:end)=0;      %clear diagonal
            temp_net=double(temp_net~=0);       %enforce double precision
            l=1;        %path length
            Lpath=temp_net;     %matrix of paths l
            D=temp_net;     %distance matrix
            n_=length(temp_net);
            disp(['[DEBUG] Nodal efficiency calculation for ', num2str(n_), ' nodes']);
            
            Idx=true;
            while any(Idx(:))
                l=l+1;
                Lpath=Lpath*temp_net;
                Idx=(Lpath~=0)&(D==0);
                D(Idx)=l;
            end
            
            D(~D | eye(n_))=inf;        %assign inf to disconnected nodes and to diagonal
            D=1./D;     %invert distance
            results=sum(D)./n;
            return
        end
        
        % Nodal Local Efficiency
        function results = calc_nodal_local_efficiency(~,temp_net)
            results = efficiency_bin(temp_net,1);
            return
        end
        
        % Nodal Clustering Coefficient
        function results = calc_nodal_clustering_coefficient(~,temp_net)
            results = clustering_coef_bu(temp_net);
            return
        end
        
        function results = calc_nodal_degree(~,temp_net)
            results = degrees_und(temp_net);
            return
        end
        
        function results = calc_nodal_betweenness(~,temp_net)
            results = betweenness_bin(temp_net);
            return
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Update UI
            movegui(app.GATFD_network_UIFigure,'southwest');
            app.ThresholdingMethodDropDown.ItemsData = {1,2,3};
            app.NormalizedClusteringCoefficientCheckBox.Text = ['Normalized Clustering Coefficient (',char(947),')'];
            app.NormalizedPathLengthCheckBox.Text = ['Normalized Path Length (',char(955),')'];
            
            % Check if brain connectivity toolbox installed
            prec_bct=exist('degrees_und.m','file');
            if ~(prec_bct==2)
                errordlg('Please check if Brain Connectivity Toolbox is correctly installed','Prerequisite not satisfied');
            end

            % Initialize settings
            app.gatn_setting.threshold_lower=0.1;
            app.gatn_setting.threshold_upper=0.4;
            app.gatn_setting.threshold_step=0.02;
            app.gatn_setting.threshold_method=1;  % 0:cost, 1:absolute, 2:relative
            app.gatn_setting.threshold_absolute=0;
            
            app.gatn_setting.parallel=0;
            app.gatn_setting.file_list={};
            app.gatn_setting.file_path_list = [];
            app.gatn_setting.condition = [];
            app.gatn_setting.node_list={};
            app.gatn_setting.node_list_unselected={};
            app.gatn_setting.node_list_selected={};
            app.gatn_setting.iscondition=0;           
            
            % Predefined global/network measures
            app.gatn_measure.glob_name={'global_efficiency','network_local_efficiency','network_clustering_coefficient',...
                                        'network_average_degree', 'characteristic_path_length','small_world_coefficient',...
                                        'normalized_clustering_coefficient','normalized_path_length','transitivity',...
                                        'assortativity','modularity'};
            app.gatn_measure.glob_count=length(app.gatn_measure.glob_name);
            app.gatn_measure.glob_measure=zeros(1,app.gatn_measure.glob_count);         % Measurement flag
            app.gatn_measure.glob_label={'glo_eff','glo_lef','glo_clc','glo_deg','glo_cpl','sw_coef','norm_cc','norm_pl','glo_tra','glo_ast','glo_mod'};
            app.gatn_measure.glob_func={@app.calc_global_efficiency,...
                                        @app.calc_network_local_efficiency,...
                                        @app.calc_network_clustering_coefficient,...
                                        @app.calc_network_average_degree,...
                                        @app.calc_network_characteristic_path,...
                                        @app.calc_sw_coefficient,...
                                        @app.calc_sw_norm_cc,...
                                        @app.calc_sw_norm_pl,...
                                        @app.calc_transitivity,...
                                        @app.calc_assortativity,...
                                        @app.calc_modularity};
            % Predefined nodal measures
            app.gatn_measure.nod_name={'nodal_efficiency','nodal_local_efficiency','nodal_clustering_coefficient','nodal_degree','nodal_betweenness'};
            app.gatn_measure.nod_label={'nod_eff','nod_lef','nod_clc','nod_deg','nod_bet'};
            app.gatn_measure.nod_count=length(app.gatn_measure.nod_name);
            app.gatn_measure.nod_measure=zeros(1,app.gatn_measure.nod_count);
            app.gatn_measure.nod_func={@app.calc_nodal_efficiency,...
                                       @app.calc_nodal_local_efficiency,...
                                       @app.calc_nodal_clustering_coefficient,...
                                       @app.calc_nodal_degree,...
                                       @app.calc_nodal_betweenness};
            
        end

        % Button pushed function: Button_loadfiles
        function gatn_loadfile(app, event)
            [fnc_temp_file,fnc_temp_path] = uigetfile('*.mat','Select One or More Files', 'MultiSelect','on');
            if isequal(fnc_temp_file,0)
                disp('Selection Canceled')
            else
                if ~iscell(fnc_temp_file)
                    fnc_temp_file={fnc_temp_file};
                end
                app.gatn_setting.file_list = fnc_temp_file;
                app.gatn_setting.file_path_list = fnc_temp_path;
                app.TextArea_LoadedFiles.Value = app.gatn_setting.file_list;
                temp_file_path=fullfile(fnc_temp_path,app.gatn_setting.file_list{1});
                temp_feature=load(temp_file_path);
                %disp("subj_data d_corr:"); disp(temp_feature.subj_data.d_corr(:,:,1));
                if isfield(temp_feature.subj_data,'d_atlas_list')
                    app.gatn_setting.node_list=temp_feature.subj_data.d_atlas_list;
                else
                    for i=1:length(temp_feature.subj_data.d_corr(1,:,1))
                        app.gatn_setting.node_list(i)={num2str(i)};
                    end
                end
                app.gatn_setting.node_list_selected=app.gatn_setting.node_list;
                app.ListBox.Items = app.gatn_setting.node_list_selected;
            end
        end

        % Button pushed function: LoadTemporalMaskButton
        function gatn_loadconditions(app, event)
            % This function upload design matrix from file that generated
            % by TrFNC-Design
            [fnc_temp_file,fnc_temp_path,~] = uigetfile('*.mat','Select Conditions');
            if isequal(fnc_temp_file,0)
                app.gatn_setting.iscondition=0;
                return
            else
                file_fullpath=fullfile(fnc_temp_path,fnc_temp_file);
                temp=load(file_fullpath,'window_condition');
                try
                    app.gatn_setting.condition=temp.window_condition.dfnc_window_condi;
                    plot(app.UIAxes,temp.window_condition.dfnc_window_condi);
                    app.gatn_setting.iscondition=1;
                catch
                    errordlg("Please check task design matrix is correct","Error");
                end
            end
        end

        % Button pushed function: CalculateNetworkPropertiesButton
        function gatn_run(app, event)
            % Check if loaded file
            if isempty(app.gatn_setting.file_list)
                errordlg("Please load the input file first","Error");
                return
            end

            % Update Settings
            % Global Measures
            app.gatn_measure.glob_measure(1)=app.GlobalEfficiencyCheckBox.Value;
            app.gatn_measure.glob_measure(2)=app.LocalEfficiencyCheckBox.Value;
            app.gatn_measure.glob_measure(3)=app.ClusteringCoefficientCheckBox.Value;
            app.gatn_measure.glob_measure(4)=app.DegreeCheckBox.Value;
            app.gatn_measure.glob_measure(5)=app.CharacteristicPathLengthCheckBox.Value;
            app.gatn_measure.glob_measure(6)=app.SmallWorldCoefficientCheckBox.Value;
            app.gatn_measure.glob_measure(7)=app.NormalizedClusteringCoefficientCheckBox.Value;
            app.gatn_measure.glob_measure(8)=app.NormalizedPathLengthCheckBox.Value;
            app.gatn_measure.glob_measure(9)=app.TransitivityCoefficientCheckBox.Value;
            app.gatn_measure.glob_measure(10)=app.AssortativityCoefficientCheckBox.Value;
            app.gatn_measure.glob_measure(11)=app.ModularityCoefficientCheckBox.Value;
            meas_glob_list=find(app.gatn_measure.glob_measure);
            run_if_glob=not(isempty(meas_glob_list));
            
            % Nodal Measures
            app.gatn_measure.nod_measure(1)=app.NodalGlobalEfficiencyCheckBox.Value;
            app.gatn_measure.nod_measure(2)=app.NodalLocalEfficiencyCheckBox.Value;
            app.gatn_measure.nod_measure(3)=app.NodalClusteringCoefficientCheckBox.Value;
            app.gatn_measure.nod_measure(4)=app.NodalDegreeCheckBox.Value;
            app.gatn_measure.nod_measure(5)=app.BetweennessCentralityCheckBox.Value;
            meas_nod_list=find(app.gatn_measure.nod_measure);
            run_if_nod=not(isempty(meas_nod_list) || isempty(app.gatn_setting.node_list_selected));
           
            if not(run_if_glob || run_if_nod)
                errordlg("Please select network properties","Error");
                return
            end
            
            [fnc_out_file,fnc_out_path] = uiputfile('*.csv','Select directory to save settings','network_properties.csv');
            if isequal(fnc_out_file,0)
                return
            end
            
            fnc_filelength=length(app.gatn_setting.file_list);
            
            % Initialize process bar
            process_fig = uifigure;
            process_d = uiprogressdlg(process_fig,'Title','Calculating Network');
            
            % ################################
            % # Check if design matches data #
            % ################################

            % Other Settings
            app.gatn_setting.threshold_lower=app.ThresholdLowerEditField.Value;
            app.gatn_setting.threshold_upper=app.ThresholdUpperEditField.Value;
            app.gatn_setting.threshold_step=app.ThresholdStepEditField.Value;
            app.gatn_setting.threshold_method=app.ThresholdingMethodDropDown.Value;  % 0:cost, 1:absolute, 2:relative
            app.gatn_setting.threshold_absolute=app.AbsoluteCheckBox.Value;
            app.gatn_setting.parallel=app.ParallelCheckBox.Value;
            
            % Load number of nodes
            if run_if_nod
                data_num_nodes=length(app.gatn_setting.node_list_selected);
                nodal_idx=zeros(1,data_num_nodes);
                for idx=1:data_num_nodes
                    nodal_idx(idx)=find(contains(app.gatn_setting.node_list,app.gatn_setting.node_list_selected{idx}));
                end
            end
            
            % Load Settings
            fnc_in_file=fullfile(app.gatn_setting.file_path_list,app.gatn_setting.file_list{1});
            temp_sub_data=load(fnc_in_file);
            
            % Check data
            dnet_data_threshold_list=app.gatn_setting.threshold_lower:app.gatn_setting.threshold_step:app.gatn_setting.threshold_upper;
            data_num_thres_steps=length(dnet_data_threshold_list);
            data_window_size=size(temp_sub_data.subj_data.d_corr,1);
            
            if run_if_glob
                dnet_data_measures_glob=app.gatn_measure.glob_name(meas_glob_list);
            else
                dnet_data_measures_glob=[];
            end
            
            if run_if_nod
                dnet_data_nodes_list=app.gatn_setting.node_list_selected;
                dnet_data_measures_nod=app.gatn_measure.nod_name(meas_nod_list);
            else
                dnet_data_nodes_list=[];
                dnet_data_measures_nod=[];
            end
            
            dnet_data_files=app.gatn_setting.file_list;

            if app.gatn_setting.parallel==1 % Parallel Processing
                process_d.Value=0.3;
                process_d.Message = {'Running in parallel'};
                
                % Create Local variable
                local_threshold_absolute=app.gatn_setting.threshold_absolute;
                local_threshold_method=app.gatn_setting.threshold_method;
                
                local_file_path=app.gatn_setting.file_path_list;
                local_file_list=app.gatn_setting.file_list;
                
                if run_if_glob
                    local_glob_func=app.gatn_measure.glob_func;
                    local_glob_count=length(meas_glob_list);
                    local_glob_name=app.gatn_measure.glob_name(meas_glob_list);
                end
                
                if run_if_nod
                    local_nod_func=app.gatn_measure.nod_func;
                    local_nod_count=length(meas_nod_list);
                    local_nod_name=app.gatn_measure.nod_name(meas_nod_list);
                end

                parfor idx_file=1:fnc_filelength              
                    % Load data file
                    fnc_in_file=fullfile(local_file_path,local_file_list{idx_file});
                    temp_mat=matfile(fnc_in_file,'writable',true);
                    temp_subj_data=temp_mat.subj_data;
                    
                    % Data type:  frames * measures * threshold steps * nodes * subjects
                    if run_if_glob
                        network_data_mat_global_para=zeros(data_window_size,local_glob_count,data_num_thres_steps);
                    else
                        network_data_mat_global_para=[];
                    end
                    if run_if_nod
                        network_data_mat_nodal_para=zeros(data_window_size,local_nod_count,data_num_thres_steps,data_num_nodes);
                    else
                        network_data_mat_nodal_para=[];
                    end
                
                    for idx_thr=1:numel(dnet_data_threshold_list)
                        % Update Progress
 
                        for idx_frame=1:data_window_size
                                                    
                            % Thresholding
                            temp_net=squeeze(temp_subj_data.d_corr(idx_frame,:,:));
                            temp_net(isnan(temp_net))=0;
                            
                            %% Absolute
                            if local_threshold_absolute==1
                                temp_net=abs(temp_net);
                            end
                            
                            %% Methods
                            if local_threshold_method==3
                                temp_net_thresh_para=prctile(temp_net(:),(1-dnet_data_threshold_list(idx_thr))*100);
                            elseif  local_threshold_method==2
                                temp_net_thresh_para=max(temp_net(:))*dnet_data_threshold_list(idx_thr);                    
                            else
                                temp_net_thresh_para=dnet_data_threshold_list(idx_thr);
                            end

                            temp_net=threshold_absolute(temp_net,temp_net_thresh_para);
                            temp_net=weight_conversion(temp_net,'binarize');
                            
                            % Calculate Properties
                            % Global
                            if run_if_glob
                                for idx_meas=1:local_glob_count
                                    network_data_mat_global_para(idx_frame,idx_meas,idx_thr)=local_glob_func{meas_glob_list(idx_meas)}(temp_net);
                                end
                            end
                            % Nodal
                            if run_if_nod
                                for idx_meas=1:local_nod_count
                                    disp(['[DEBUG] File: ', app.gatn_setting.file_list{idx_file}, ', Window: ', num2str(idx_frame), ', Threshold: ', num2str(dnet_data_threshold_list(idx_thr))]);
                                    disp(['[DEBUG] Nodal measure: ', local_nod_name{idx_meas}]);
                                    disp(['[DEBUG] temp_net size: ', mat2str(size(temp_net))]);
                                    disp('[DEBUG] temp_net (binarized):');
                                    disp(temp_net);
                                    %disp(['[DEBUG] Selected nodes: ', strjoin(app.gatn_setting.node_list_selected, ', ')]);
                                    temp_prop=app.gatn_measure.nod_func{meas_nod_list(idx_meas)}(temp_net);
                                    disp(['[DEBUG] Raw nodal measure output: ', mat2str(temp_prop)]);
                                    disp(['[DEBUG] Node indices: ', mat2str(nodal_idx)]);
                                    disp(['[DEBUG] Selected node values: ', mat2str(temp_prop(nodal_idx))]);
                                    network_data_mat_nodal_para(idx_frame,idx_meas,idx_thr,:)=temp_prop(nodal_idx);
                                end
                            end
                        end    % End of for-loop (Sliding window)
                    end    % End of for-loop (Threshold list)
                    if run_if_glob
                        temp_mat.network_data_mat_global=network_data_mat_global_para;
                        temp_mat.network_measure_global=local_glob_name;
                    end
                    if run_if_nod
                        temp_mat.network_data_mat_nodal=network_data_mat_nodal_para;
                        temp_mat.network_measure_nodal=local_nod_name;
                    end
                end
                    
            else % Non-parallel
                % Generate Settings
                if run_if_glob
                    local_glob_count=length(meas_glob_list);
                    local_glob_name=app.gatn_measure.glob_name(meas_glob_list);
                end
                if run_if_nod
                    local_nod_count=length(meas_nod_list);
                    local_nod_name=app.gatn_measure.nod_name(meas_nod_list);
                end
                
                for idx_file=1:fnc_filelength
                    % Update progress info
                    process_d.Value=(idx_file-1)/fnc_filelength*0.7;
                    process_d.Message = {[num2str(idx_file),'/',num2str(fnc_filelength),': Calculating Network Properties']};
                    
                    % Load data file
                    fnc_in_file=fullfile(app.gatn_setting.file_path_list,app.gatn_setting.file_list{idx_file});
                    
                    temp_mat=matfile(fnc_in_file,'writable',true);
                    temp_subj_data=temp_mat.subj_data;
                    
                    if run_if_glob
                        network_data_mat_global=zeros(data_window_size,local_glob_count,data_num_thres_steps);
                    end
                    if run_if_nod
                        network_data_mat_nodal=zeros(data_window_size,local_nod_count,data_num_thres_steps,data_num_nodes);
                    end
                        
                    for idx_thr=1:data_num_thres_steps
                        % Update Progress
                        process_d.Value=(idx_file-0.7*(1-(dnet_data_threshold_list(idx_thr)-app.gatn_setting.threshold_lower)/(app.gatn_setting.threshold_upper-app.gatn_setting.threshold_lower)))/fnc_filelength*(7/10);
                        process_d.Message = {[num2str(idx_file),'/',num2str(fnc_filelength),': Calculating Network Properties - Cost:',num2str(dnet_data_threshold_list(idx_thr))]};
                       
                        for idx_frame=1:data_window_size
                            % Thresholding
                            temp_net=squeeze(temp_subj_data.d_corr(idx_frame,:,:)); 
                            %disp("temp_net"); disp(temp_net); %good so far
                            temp_net(isnan(temp_net))=0;
                            
                            %% Absolutecnn
                            if app.gatn_setting.threshold_absolute==1
                                temp_net=abs(temp_net);
                            end
                            %% Methods
                            %disp(['[DEBUG] Processing file: ', app.gatn_setting.file_list{idx_file}, ', Window: ', num2str(idx_frame), ', Threshold: ', num2str(dnet_data_threshold_list(idx_thr))]);
                            %disp(['[DEBUG] temp_net (before thresholding):']); disp(temp_net);
                            if app.gatn_setting.threshold_method==3
                                temp_net_thresh=prctile(temp_net(:),(1-dnet_data_threshold_list(idx_thr))*100);
                            elseif  app.gatn_setting.threshold_method==2
                                temp_net_thresh=max(temp_net(:))*dnet_data_threshold_list(idx_thr);                    
                            elseif app.gatn_setting.threshold_method==1
                                temp_net_thresh=dnet_data_threshold_list(idx_thr);
                            end
                            
                            temp_net=threshold_absolute(temp_net,temp_net_thresh);
                            temp_net=weight_conversion(temp_net,'binarize');
                            
                            % Calculate Properties
                            % Global
                            if run_if_glob
                                for idx_meas=1:local_glob_count
                                    network_data_mat_global(idx_frame,idx_meas,idx_thr)=app.gatn_measure.glob_func{meas_glob_list(idx_meas)}(temp_net);
                                end
                            end
                            % Nodal
                            if run_if_nod
                                for idx_meas=1:local_nod_count
                                    %disp(['[DEBUG] File: ', app.gatn_setting.file_list{idx_file}, ', Window: ', num2str(idx_frame), ', Threshold: ', num2str(dnet_data_threshold_list(idx_thr))]);
                                    %disp(['[DEBUG] Nodal measure: ', local_nod_name{idx_meas}]);
                                    temp_prop=app.gatn_measure.nod_func{meas_nod_list(idx_meas)}(temp_net);
                                    %disp("temp_prop (raw nodal measure output):");
                                    %disp(temp_prop);
                                    %disp(['[DEBUG] Raw nodal measure output: ', mat2str(temp_prop)]);
                                    % disp(['[DEBUG] Node indices: ', mat2str(nodal_idx)]);
                                    %disp(['[DEBUG] Selected node values: ', mat2str(temp_prop(nodal_idx))]);
                                    network_data_mat_nodal(idx_frame,idx_meas,idx_thr,:)=temp_prop(nodal_idx);
                                end
                            end
                        end
                    end
                    if run_if_glob
                        temp_mat.network_data_mat_global=network_data_mat_global;
                        temp_mat.network_measure_global=local_glob_name;
                    end
                    if run_if_nod
                        temp_mat.network_data_mat_nodal=network_data_mat_nodal;
                        temp_mat.network_measure_nodal=local_nod_name;
                    end
                end
            end
            
            
            process_d.Value=0.7;
            process_d.Message = {'Update Data'};
            
            % Write output
            out_path=fullfile(fnc_out_path,fnc_out_file);
            save([out_path,'.mat'],"dnet_data_threshold_list","dnet_data_nodes_list","dnet_data_measures_glob","dnet_data_measures_nod","dnet_data_files",'-v7.3');
            out_mat=matfile([out_path,'.mat'],'writable',true);
            
            if run_if_glob
                out_mat.dnet_data_data_mat_global=zeros(data_window_size,local_glob_count,data_num_thres_steps,fnc_filelength);
            end
            if run_if_nod
                out_mat.dnet_data_data_mat_nodal=zeros(data_window_size,local_nod_count,data_num_thres_steps,data_num_nodes,fnc_filelength);
            end

            % Create output data
            if app.gatn_setting.iscondition==1
                out_data={};
                out_data(1)={'filename'};
                % Thresholding
                if run_if_glob
                    for idx_meas=1:local_glob_count
                        out_data(end+1)={[app.gatn_measure.glob_label{meas_glob_list(idx_meas)},'_var']};
                        out_data(end+1)={[app.gatn_measure.glob_label{meas_glob_list(idx_meas)},'_mean']};
                    end
                end
                if run_if_nod
                    for idx_node=1:data_num_nodes
                        for idx_meas=1:local_nod_count
                            out_data(end+1)={[app.gatn_measure.nod_label{meas_nod_list(idx_meas)},'_',app.gatn_setting.node_list_selected{idx_node},'_var']};
                            out_data(end+1)={[app.gatn_measure.nod_label{meas_nod_list(idx_meas)},'_',app.gatn_setting.node_list_selected{idx_node},'_mean']};
                        end
                    end
                end
                out_data=strrep(out_data,'.','_');
            end    % End of if (condition is loaded)
            
            for idx_file=1:fnc_filelength            % Load data file
                process_d.Value=0.7+0.3*idx_file/fnc_filelength;
                process_d.Message = {[num2str(idx_file),'/',num2str(fnc_filelength),': Update Data']};
                
                fnc_in_file=fullfile(app.gatn_setting.file_path_list,app.gatn_setting.file_list{idx_file});
                temp_mat=matfile(fnc_in_file);
                
                if run_if_glob
                    network_data_mat_global=temp_mat.network_data_mat_global;
                end
                if run_if_nod
                    network_data_mat_nodal=temp_mat.network_data_mat_nodal;
                end
                if fnc_filelength==1
                    if run_if_glob
                        out_mat.dnet_data_data_mat_global(:,:,:)=network_data_mat_global;
                    end
                    if run_if_nod
                        out_mat.dnet_data_data_mat_nodal(:,:,:,:)=network_data_mat_nodal;
                    end
                else
                    if run_if_glob
                        out_mat.dnet_data_data_mat_global(:,:,:,idx_file)=network_data_mat_global;
                    end
                    if run_if_nod
                        out_mat.dnet_data_data_mat_nodal(:,:,:,:,idx_file)=network_data_mat_nodal;
                    end
                end
                
                % Write Data
                if app.gatn_setting.iscondition==1      % If condition is loaded
                    out_counter=1;
                    out_data(idx_file+1,out_counter)=app.gatn_setting.file_list(idx_file);
                    out_counter=out_counter+1;
                
                    condit_index=app.gatn_setting.condition==1;
                    if run_if_glob
                        for idx_meas=1:local_glob_count
                            out_data(idx_file+1,out_counter)={mean(var(squeeze(network_data_mat_global(condit_index,idx_meas,:)),0,1))};
                            out_counter=out_counter+1;
                            out_data(idx_file+1,out_counter)={mean(mean(squeeze(network_data_mat_global(condit_index,idx_meas,:)),1))};
                            out_counter=out_counter+1;
                        end
                    end
                    
                    if run_if_nod
                        for idx_node=1:data_num_nodes
                        % Calculate Properties
                            for idx_meas=1:local_nod_count
                                out_data(idx_file+1,out_counter)={mean(var(squeeze(network_data_mat_nodal(condit_index,idx_meas,:,idx_node)),0,1))};
                                out_counter=out_counter+1;
                                out_data(idx_file+1,out_counter)={mean(mean(squeeze(network_data_mat_nodal(condit_index,idx_meas,:,idx_node)),1))};
                                out_counter=out_counter+1;
                            end     % End of for-loop (measure)
                        end    % End of for-loop (nodes)
                    end    % End of if (nod)
                end    % End of if (if condition is loaded)
            end    % End of for-loop (Load files)
            
            if app.gatn_setting.iscondition==1
                data_cell=cell2table(out_data(2:end,:));
                data_cell.Properties.VariableNames=out_data(1,:);
                % Output data
                % Issue: other function if older version
                writetable(data_cell,out_path,'Delimiter','comma');
            end
            
            close(process_d);
            close(process_fig);
            msgbox("Finished");
        end

        % Button pushed function: SelectNodesButton
        function gatn_select_nodes(app, event)
            app.choose_node();
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GATFD_network_UIFigure and hide until all components are created
            app.GATFD_network_UIFigure = uifigure('Visible', 'off');
            app.GATFD_network_UIFigure.Position = [100 100 650 788];
            app.GATFD_network_UIFigure.Name = 'GAT-FD - Network Property Calculation v0.2a';

            % Create Button_loadfiles
            app.Button_loadfiles = uibutton(app.GATFD_network_UIFigure, 'push');
            app.Button_loadfiles.ButtonPushedFcn = createCallbackFcn(app, @gatn_loadfile, true);
            app.Button_loadfiles.Position = [131 522 100 22];
            app.Button_loadfiles.Text = 'Load Files';

            % Create LoadedFIlesTextAreaLabel
            app.LoadedFIlesTextAreaLabel = uilabel(app.GATFD_network_UIFigure);
            app.LoadedFIlesTextAreaLabel.HorizontalAlignment = 'right';
            app.LoadedFIlesTextAreaLabel.Position = [60 740 75 22];
            app.LoadedFIlesTextAreaLabel.Text = 'Loaded FIles';

            % Create TextArea_LoadedFiles
            app.TextArea_LoadedFiles = uitextarea(app.GATFD_network_UIFigure);
            app.TextArea_LoadedFiles.Editable = 'off';
            app.TextArea_LoadedFiles.Position = [60 554 171 183];
            app.TextArea_LoadedFiles.Value = {'Loaded Files'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GATFD_network_UIFigure);
            title(app.UIAxes, 'State Conditions')
            xlabel(app.UIAxes, 'Time')
            ylabel(app.UIAxes, 'Condition')
            app.UIAxes.PlotBoxAspectRatio = [2.09126984126984 1 1];
            app.UIAxes.Position = [238 554 349 200];

            % Create LoadTemporalMaskButton
            app.LoadTemporalMaskButton = uibutton(app.GATFD_network_UIFigure, 'push');
            app.LoadTemporalMaskButton.ButtonPushedFcn = createCallbackFcn(app, @gatn_loadconditions, true);
            app.LoadTemporalMaskButton.Position = [450 522 126 22];
            app.LoadTemporalMaskButton.Text = 'Load Temporal Mask';

            % Create CalculateNetworkPropertiesButton
            app.CalculateNetworkPropertiesButton = uibutton(app.GATFD_network_UIFigure, 'push');
            app.CalculateNetworkPropertiesButton.ButtonPushedFcn = createCallbackFcn(app, @gatn_run, true);
            app.CalculateNetworkPropertiesButton.Position = [414 27 172 22];
            app.CalculateNetworkPropertiesButton.Text = 'Calculate Network Properties';

            % Create NetworkAverageGlobalMeasuresPanel
            app.NetworkAverageGlobalMeasuresPanel = uipanel(app.GATFD_network_UIFigure);
            app.NetworkAverageGlobalMeasuresPanel.Title = 'Network Average (Global) Measures';
            app.NetworkAverageGlobalMeasuresPanel.Position = [60 239 530 128];

            % Create GlobalEfficiencyCheckBox
            app.GlobalEfficiencyCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.GlobalEfficiencyCheckBox.Text = 'Global Efficiency';
            app.GlobalEfficiencyCheckBox.Position = [10 82 112 22];

            % Create LocalEfficiencyCheckBox
            app.LocalEfficiencyCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.LocalEfficiencyCheckBox.Text = 'Local Efficiency';
            app.LocalEfficiencyCheckBox.Position = [10 57 106 22];

            % Create ClusteringCoefficientCheckBox
            app.ClusteringCoefficientCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.ClusteringCoefficientCheckBox.Text = 'Clustering Coefficient';
            app.ClusteringCoefficientCheckBox.Position = [10 32 137 22];

            % Create DegreeCheckBox
            app.DegreeCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.DegreeCheckBox.Text = 'Degree';
            app.DegreeCheckBox.Position = [10 7 60 22];

            % Create SmallWorldCoefficientCheckBox
            app.SmallWorldCoefficientCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.SmallWorldCoefficientCheckBox.Text = 'Small World Coefficient';
            app.SmallWorldCoefficientCheckBox.Position = [157 57 147 22];

            % Create ModularityCoefficientCheckBox
            app.ModularityCoefficientCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.ModularityCoefficientCheckBox.Text = 'Modularity Coefficient';
            app.ModularityCoefficientCheckBox.Position = [382 32 137 22];

            % Create NormalizedClusteringCoefficientCheckBox
            app.NormalizedClusteringCoefficientCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.NormalizedClusteringCoefficientCheckBox.Text = 'Normalized Clustering Coefficient ';
            app.NormalizedClusteringCoefficientCheckBox.Position = [157 32 229 22];

            % Create NormalizedPathLengthCheckBox
            app.NormalizedPathLengthCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.NormalizedPathLengthCheckBox.Text = 'Normalized Path Length';
            app.NormalizedPathLengthCheckBox.Position = [157 7 209 22];

            % Create TransitivityCoefficientCheckBox
            app.TransitivityCoefficientCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.TransitivityCoefficientCheckBox.Text = 'Transitivity Coefficient';
            app.TransitivityCoefficientCheckBox.Position = [383 82 139 22];

            % Create AssortativityCoefficientCheckBox
            app.AssortativityCoefficientCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.AssortativityCoefficientCheckBox.Text = 'Assortativity Coefficient';
            app.AssortativityCoefficientCheckBox.Position = [383 57 147 22];

            % Create CharacteristicPathLengthCheckBox
            app.CharacteristicPathLengthCheckBox = uicheckbox(app.NetworkAverageGlobalMeasuresPanel);
            app.CharacteristicPathLengthCheckBox.Text = 'Characteristic Path Length';
            app.CharacteristicPathLengthCheckBox.Position = [157 82 165 22];

            % Create NodalMeasuresPanel
            app.NodalMeasuresPanel = uipanel(app.GATFD_network_UIFigure);
            app.NodalMeasuresPanel.Title = 'Nodal Measures';
            app.NodalMeasuresPanel.Position = [60 54 530 177];

            % Create ListBox
            app.ListBox = uilistbox(app.NodalMeasuresPanel);
            app.ListBox.Items = {'Selected Nodes'};
            app.ListBox.Position = [8 8 231 118];
            app.ListBox.Value = 'Selected Nodes';

            % Create SelectNodesButton
            app.SelectNodesButton = uibutton(app.NodalMeasuresPanel, 'push');
            app.SelectNodesButton.ButtonPushedFcn = createCallbackFcn(app, @gatn_select_nodes, true);
            app.SelectNodesButton.Tooltip = {'Select Nodes for nodal network properties calculation.'};
            app.SelectNodesButton.Position = [8 131 100 22];
            app.SelectNodesButton.Text = 'Select Nodes';

            % Create NodalGlobalEfficiencyCheckBox
            app.NodalGlobalEfficiencyCheckBox = uicheckbox(app.NodalMeasuresPanel);
            app.NodalGlobalEfficiencyCheckBox.Text = 'Global Efficiency';
            app.NodalGlobalEfficiencyCheckBox.Position = [254 125 112 22];

            % Create NodalLocalEfficiencyCheckBox
            app.NodalLocalEfficiencyCheckBox = uicheckbox(app.NodalMeasuresPanel);
            app.NodalLocalEfficiencyCheckBox.Text = 'Local Efficiency';
            app.NodalLocalEfficiencyCheckBox.Position = [254 95 106 22];

            % Create NodalClusteringCoefficientCheckBox
            app.NodalClusteringCoefficientCheckBox = uicheckbox(app.NodalMeasuresPanel);
            app.NodalClusteringCoefficientCheckBox.Text = 'Clustering Coefficient';
            app.NodalClusteringCoefficientCheckBox.Position = [254 65 137 22];

            % Create NodalDegreeCheckBox
            app.NodalDegreeCheckBox = uicheckbox(app.NodalMeasuresPanel);
            app.NodalDegreeCheckBox.Text = 'Degree';
            app.NodalDegreeCheckBox.Position = [254 35 60 22];

            % Create BetweennessCentralityCheckBox
            app.BetweennessCentralityCheckBox = uicheckbox(app.NodalMeasuresPanel);
            app.BetweennessCentralityCheckBox.Text = 'Betweenness Centrality';
            app.BetweennessCentralityCheckBox.Position = [254 5 148 22];

            % Create ThresholdingPanel
            app.ThresholdingPanel = uipanel(app.GATFD_network_UIFigure);
            app.ThresholdingPanel.Title = 'Thresholding';
            app.ThresholdingPanel.Position = [60 372 530 145];

            % Create ThresholdingMethodDropDownLabel
            app.ThresholdingMethodDropDownLabel = uilabel(app.ThresholdingPanel);
            app.ThresholdingMethodDropDownLabel.Position = [11 95 118 22];
            app.ThresholdingMethodDropDownLabel.Text = 'Thresholding Method';

            % Create ThresholdingMethodDropDown
            app.ThresholdingMethodDropDown = uidropdown(app.ThresholdingPanel);
            app.ThresholdingMethodDropDown.Items = {'Absolute', 'Proportional', 'Cost'};
            app.ThresholdingMethodDropDown.Position = [353 95 160 22];
            app.ThresholdingMethodDropDown.Value = 'Cost';

            % Create ThresholdRangeLabel
            app.ThresholdRangeLabel = uilabel(app.ThresholdingPanel);
            app.ThresholdRangeLabel.Position = [11 65 97 22];
            app.ThresholdRangeLabel.Text = 'Threshold Range';

            % Create ThresholdUpperEditField
            app.ThresholdUpperEditField = uieditfield(app.ThresholdingPanel, 'numeric');
            app.ThresholdUpperEditField.Position = [442 65 71 22];
            app.ThresholdUpperEditField.Value = 0.4;

            % Create ThresholdLowerEditField
            app.ThresholdLowerEditField = uieditfield(app.ThresholdingPanel, 'numeric');
            app.ThresholdLowerEditField.Position = [339 65 71 22];
            app.ThresholdLowerEditField.Value = 0.1;

            % Create ThresholdEditFieldLabel_2
            app.ThresholdEditFieldLabel_2 = uilabel(app.ThresholdingPanel);
            app.ThresholdEditFieldLabel_2.HorizontalAlignment = 'center';
            app.ThresholdEditFieldLabel_2.Position = [417 65 19 22];
            app.ThresholdEditFieldLabel_2.Text = '-';

            % Create ThresholdStepLabel
            app.ThresholdStepLabel = uilabel(app.ThresholdingPanel);
            app.ThresholdStepLabel.Position = [11 35 87 22];
            app.ThresholdStepLabel.Text = 'Threshold Step';

            % Create ThresholdStepEditField
            app.ThresholdStepEditField = uieditfield(app.ThresholdingPanel, 'numeric');
            app.ThresholdStepEditField.Position = [442 35 71 22];
            app.ThresholdStepEditField.Value = 0.02;

            % Create AbsoluteCheckBox
            app.AbsoluteCheckBox = uicheckbox(app.ThresholdingPanel);
            app.AbsoluteCheckBox.Text = 'Use absolute value of correlation coefficient';
            app.AbsoluteCheckBox.Position = [11 5 260 22];

            % Create ParallelCheckBox
            app.ParallelCheckBox = uicheckbox(app.GATFD_network_UIFigure);
            app.ParallelCheckBox.Text = 'Run in parallel';
            app.ParallelCheckBox.Position = [65 27 98 22];

            % Show the figure after all components are created
            app.GATFD_network_UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = gatfd_ui_network

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.GATFD_network_UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GATFD_network_UIFigure)
        end
    end
end