classdef gatfd_ui_design < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GATFD_design_UIFigure         matlab.ui.Figure
        DesignDurationSequencesLabel  matlab.ui.control.Label
        DurationSequenceEditField     matlab.ui.control.EditField
        DesignConditionSequence0forrestEditFieldLabel  matlab.ui.control.Label
        DesignConditionSequenceEditField  matlab.ui.control.EditField
        UIAxes                        matlab.ui.control.UIAxes
        TRsEditFieldLabel             matlab.ui.control.Label
        TRsEditField                  matlab.ui.control.NumericEditField
        UpdateDesignButton            matlab.ui.control.Button
        SaveDesignMatrixButton        matlab.ui.control.Button
        WindowSizeTREditFieldLabel    matlab.ui.control.Label
        WindowSizeField               matlab.ui.control.NumericEditField
        StepSizeTREditFieldLabel      matlab.ui.control.Label
        StepSizeField                 matlab.ui.control.NumericEditField
        DefaulSettingButton           matlab.ui.control.Button
        PlotOptionPanel               matlab.ui.container.Panel
        TaskDesignCheckBox            matlab.ui.control.CheckBox
        HemodynamicResponseCheckBox   matlab.ui.control.CheckBox
        TemporalMaskCheckBox          matlab.ui.control.CheckBox
        TabGroup                      matlab.ui.container.TabGroup
        AutomaticTab                  matlab.ui.container.Tab
        CheckBoxActivation            matlab.ui.control.CheckBox
        EstimatedActivationLevelThreshold01EditFieldLabel  matlab.ui.control.Label
        EstimatedActivationLevelThreshold01EditField  matlab.ui.control.NumericEditField
        EstimatedActivationCoveragePercentageThresholdEditField  matlab.ui.control.NumericEditField
        EstimatedActivationCoveragePercentageThresholdEditFieldLabel  matlab.ui.control.Label
        CheckBoxCondition             matlab.ui.control.CheckBox
        ConditionCoveragePercentageThresholdEditFieldLabel  matlab.ui.control.Label
        ConditionCoveragePercentageThresholdEditField  matlab.ui.control.NumericEditField
        ManualTab                     matlab.ui.container.Tab
        StageConditionSequenceEditFieldLabel  matlab.ui.control.Label
        StageConditionSequenceEditField  matlab.ui.control.EditField
        StageDurationSequenceTRsEditFieldLabel  matlab.ui.control.Label
        StageDurationSequenceTRsEditField  matlab.ui.control.EditField
        WiththecurrentsettingthetotalnumberofframesisLabel  matlab.ui.control.Label
        totalstepsEditField           matlab.ui.control.NumericEditField
        CalculateButton               matlab.ui.control.Button
        TaskConditionSpecificationLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        gatd_setting        % App Parameters
        gatd_plot           % Plot info holder
    end
    
    
    methods (Access = private)
        
        function parameters_default(app)
            % Slding Window Setting
            app.gatd_setting.tr=1;
            app.gatd_setting.window_size=20;
            app.gatd_setting.step_size=1;
            
            % Design
            app.gatd_setting.design_condition_list='';
            app.gatd_setting.design_duration_list='';
            
            % Design Auto/Manual
            app.gatd_setting.stage_method=1;    % 1: auto; 2: manual
            
            % Condition
            app.gatd_setting.activation_enable=1;
            app.gatd_setting.activation_level=0.8;
            app.gatd_setting.activation_percent=80;
            app.gatd_setting.condition_enable=1;
            app.gatd_setting.condition_percent=80;
            
            app.gatd_setting.stage_condition_list='';
            app.gatd_setting.stage_duration_list='';
        end
        
        function parameters_update_field(app)
            % Slding Window Setting
            app.TRsEditField.Value = app.gatd_setting.tr;
            app.WindowSizeField.Value = app.gatd_setting.window_size;
            app.StepSizeField.Value = app.gatd_setting.step_size;
           
            % Design
            app.DesignConditionSequenceEditField.Value = app.gatd_setting.design_condition_list;
            app.DurationSequenceEditField.Value = app.gatd_setting.design_duration_list;
            
            % Condition
            app.CheckBoxActivation.Value = app.gatd_setting.activation_enable;
            app.EstimatedActivationLevelThreshold01EditField.Value = app.gatd_setting.activation_level;
            app.EstimatedActivationCoveragePercentageThresholdEditField.Value = app.gatd_setting.activation_percent;
            app.CheckBoxCondition.Value = app.gatd_setting.condition_enable;
            app.ConditionCoveragePercentageThresholdEditField.Value = app.gatd_setting.condition_percent;
            
            app.StageConditionSequenceEditField.Value = app.gatd_setting.stage_condition_list;
            app.StageDurationSequenceTRsEditField.Value = app.gatd_setting.stage_duration_list;

        end
        
        function parameters_update_setting(app)
            % Slding Window Setting
            app.gatd_setting.tr=app.TRsEditField.Value;
            app.gatd_setting.window_size=app.WindowSizeField.Value;
            app.gatd_setting.step_size=app.StepSizeField.Value;
            
            % Design
            app.gatd_setting.design_condition_list=app.DesignConditionSequenceEditField.Value;
            app.gatd_setting.design_duration_list=app.DurationSequenceEditField.Value;
            
            % Design Auto/Manual
            selectedTab = app.TabGroup.SelectedTab;
            switch selectedTab.Title
                case 'Automatic'
                    app.gatd_setting.stage_method=1;
                    
                    % Condition
                    app.gatd_setting.activation_enable=app.CheckBoxActivation.Value;
                    app.gatd_setting.activation_level=app.EstimatedActivationLevelThreshold01EditField.Value;
                    app.gatd_setting.activation_percent=app.EstimatedActivationCoveragePercentageThresholdEditField.Value;
                    app.gatd_setting.condition_enable=app.CheckBoxCondition.Value;
                    app.gatd_setting.condition_percent=app.ConditionCoveragePercentageThresholdEditField.Value;
                case 'Manual'
                    app.gatd_setting.stage_method=2;
                    
                    % Condition
                    app.gatd_setting.stage_condition_list=app.StageConditionSequenceEditField.Value;
                    app.gatd_setting.stage_duration_list=app.StageDurationSequenceTRsEditField.Value;
            end
            

        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            movegui(app.GATFD_design_UIFigure,'northeast');
            
            % Default Setting
            app.parameters_default();
            app.parameters_update_field();
        end

        % Button pushed function: DefaulSettingButton
        function gatd_default(app, event)
            app.parameters_default();
            app.parameters_update_field();
        end

        % Value changed function: CheckBoxCondition
        function gatd_s_condition(app, event)
            app.gatd_setting.condition_enable = app.CheckBoxCondition.Value;
            if app.gatd_setting.condition_enable==1
                app.ConditionCoveragePercentageThresholdEditField.Enable = 'on';
                app.ConditionCoveragePercentageThresholdEditFieldLabel.Enable = 'on';
            else
                app.ConditionCoveragePercentageThresholdEditField.Enable = 'off';
                app.ConditionCoveragePercentageThresholdEditFieldLabel.Enable = 'off';
            end
        end

        % Value changed function: CheckBoxActivation
        function gatd_s_activation(app, event)
            app.gatd_setting.activation_enable = app.CheckBoxActivation.Value;
            if app.gatd_setting.activation_enable==1
                app.EstimatedActivationLevelThreshold01EditFieldLabel.Enable = 'on';
                app.EstimatedActivationLevelThreshold01EditField.Enable = 'on';
                app.EstimatedActivationCoveragePercentageThresholdEditField.Enable = 'on';
                app.EstimatedActivationCoveragePercentageThresholdEditFieldLabel.Enable = 'on';
            else
                app.EstimatedActivationLevelThreshold01EditFieldLabel.Enable = 'off';
                app.EstimatedActivationLevelThreshold01EditField.Enable = 'off';
                app.EstimatedActivationCoveragePercentageThresholdEditField.Enable = 'off';
                app.EstimatedActivationCoveragePercentageThresholdEditFieldLabel.Enable = 'off';
            end
        end

        % Button pushed function: CalculateButton
        function gatd_update_frames(app, event)
            app.parameters_update_setting();
            try
                design_duration_list_in_s=str2num(app.gatd_setting.design_duration_list);
                dfnc_length=floor(sum(design_duration_list_in_s)/app.gatd_setting.tr);
                app.totalstepsEditField.Value=dfnc_length-app.gatd_setting.window_size+1;
            catch
                app.totalstepsEditField.Value=0;
            end
        end

        % Button pushed function: UpdateDesignButton
        function gatd_update_design(app, event)
            % Update settings
            app.parameters_update_setting();
            
            % Calculate hrf
            dfnc_hrf=spm_hrf(app.gatd_setting.tr);
                        
            % Update Condition
            if isempty(app.gatd_setting.design_condition_list)
                errordlg("Please enter condition sequence");
                return
            end
            
            try
                run_design_condition_list=str2num(app.gatd_setting.design_condition_list);
            catch
                errordlg("Please enter the condition sequence with only one space between numbers","Error");
            end
            
            % Update Duration       
            if isempty(app.gatd_setting.design_duration_list)
                errordlg("Please enter duration sequence");
                return
            end
            
            try
                run_design_duration_list_in_s=str2num(app.gatd_setting.design_duration_list);
            catch
                errordlg("Please enter the duration sequence with only one space between numbers","Error");
            end

            % Check if two sequences match
            if length(run_design_condition_list)~=length(run_design_duration_list_in_s)
                errordlg("The condition sequence and duration sequence do not match","Error");
                return
            end
            
            % Load and check if manual
            if app.gatd_setting.stage_method==2
                if isempty(app.gatd_setting.stage_condition_list)
                    errordlg("Please enter condition sequence");
                    return
                end
                
                try
                    run_stage_condition_list=str2num(app.gatd_setting.stage_condition_list);
                catch
                    errordlg("Please enter the condition sequence with only one space between numbers","Error");
                end
                
                % Update Duration
                if isempty(app.gatd_setting.stage_duration_list)
                    errordlg("Please enter duration sequence");
                    return
                end
                
                try
                    run_stage_duration_list=str2num(app.gatd_setting.stage_duration_list);
                catch
                    errordlg("Please enter the duration sequence with only one space between numbers","Error");
                end
    
                % Check if two sequences match
                if length(run_stage_condition_list)~=length(run_stage_duration_list)
                    errordlg("The condition sequence and duration sequence do not match","Error");
                    return
                end
            end
            
            % Update design
            dfnc_length=floor(sum(run_design_duration_list_in_s)/app.gatd_setting.tr);    % Get full design length
            temp_duration_list_s=cumsum(run_design_duration_list_in_s);
            temp_duration_list=zeros(size(run_design_duration_list_in_s));
            for i=1:(length(temp_duration_list)-1)
                temp_duration_list(i)=round(temp_duration_list_s(i)/app.gatd_setting.tr);
            end
            temp_duration_list(end)=floor(temp_duration_list_s(end)/app.gatd_setting.tr);
            
%             for i=1:(length(temp_duration_list)-1)
%                 temp_duration_list(length(temp_duration_list)-i+1)=temp_duration_list(length(temp_duration_list)-i+1)-temp_duration_list(length(temp_duration_list)-i);
%             end
            run_design_duration_list=[temp_duration_list(1) diff(temp_duration_list)];
            
            dfnc_design=zeros(1,dfnc_length);
            duration_list_cs=cumsum(run_design_duration_list);    % Get transition position
            dfnc_design(1:duration_list_cs(1))=run_design_condition_list(1);  % Load first condition
            for idx=2:length(run_design_condition_list)
                dfnc_design(duration_list_cs(idx-1)+1:duration_list_cs(idx))=run_design_condition_list(idx);
            end
            app.gatd_setting.dfnc_design=dfnc_design;   % Save design to settings
            
            % Update HRF
            dfnc_reponse=conv(app.gatd_setting.dfnc_design,dfnc_hrf);
            app.gatd_setting.dfnc_reponse=dfnc_reponse(1:dfnc_length);      % Cut data in the end
            
            % Calculate Condition
            switch app.gatd_setting.stage_method
                case 1
                    dfnc_condi=app.gatd_setting.dfnc_reponse>app.gatd_setting.activation_level;
                    dfnc_condi=double(dfnc_condi);
                    dfnc_window_condi=ones(1,dfnc_length+1-app.gatd_setting.window_size);
                    
                    if app.CheckBoxActivation.Value==1
                        for idx=1:length(dfnc_window_condi)
                            if sum(dfnc_condi(idx:idx+app.gatd_setting.window_size-1))<(app.gatd_setting.window_size*app.gatd_setting.condition_percent/100)
                                dfnc_window_condi(idx)=0;
                            end
                        end
                    end
                    
                    if app.CheckBoxCondition.Value==1
                        for idx=1:length(dfnc_window_condi)
                            if sum(dfnc_design(idx:idx+app.gatd_setting.window_size-1))<(app.gatd_setting.window_size*app.gatd_setting.condition_percent/100)
                                dfnc_window_condi(idx)=0;
                            end
                        end
                    end

                case 2
                    stage_length=sum(run_stage_duration_list); % Get full stage length
                    if stage_length==(dfnc_length+1-app.gatd_setting.window_size)
                        stage_design=zeros(1,stage_length);
                        duration_list_cs=cumsum(run_stage_duration_list);    % Get transition position
                        stage_design(1:duration_list_cs(1))=run_stage_condition_list(1);  % Load first condition
                        for idx=2:length(run_stage_condition_list)
                            stage_design(duration_list_cs(idx-1)+1:duration_list_cs(idx))=run_stage_condition_list(idx);
                        end
                        dfnc_window_condi=stage_design;
                    else
                        errordlg("The specified stage does not match design and window size","Error");
                        return
                    end
            end
            app.gatd_setting.dfnc_window_condi=dfnc_window_condi;
            
            dfnc_condi_with_window_n=conv(ones(1,app.gatd_setting.window_size),dfnc_window_condi);
            dfnc_condi_with_window=dfnc_condi_with_window_n>0;
            dfnc_condi_with_window=double(dfnc_condi_with_window);
            
            % Plot data
            % Calculate design changing timing point
            design_switch_point=find(diff(app.gatd_setting.dfnc_design)~=0);
            design_time_point=[1,repelem(design_switch_point+0.5,2),dfnc_length];
            design_value_point=[app.gatd_setting.dfnc_design(1),app.gatd_setting.dfnc_design(1),repelem(app.gatd_setting.dfnc_design(design_switch_point+1),2)];
            design_value_point(end)=app.gatd_setting.dfnc_design(end);%design_value_point(1:end-1);
            

            
            % Calculate y limits
            dfnc_reponse_max=max(app.gatd_setting.dfnc_reponse);            % Used for plotting y-axis
            dfnc_reponse_min=min(app.gatd_setting.dfnc_reponse);            % Used for plotting y-axis
            
            app.UIAxes.cla;
            axis(app.UIAxes,[0,dfnc_length,dfnc_reponse_min-0.1,dfnc_reponse_max+0.1]);
            hold(app.UIAxes,'on');
            %app.gatd_plot.pholder_window=area(app.UIAxes,dfnc_condi_with_window*dfnc_reponse_max*100,'LineStyle',':');
            %app.gatd_plot.pholder_window=patch(app.UIAxes,[1:dfnc_length,dfnc_length:-1:1],[dfnc_condi_with_window*(dfnc_reponse_max+0.1),zeros(1,dfnc_length)],[0.52 1 0.36],'LineStyle',':');
            app.gatd_plot.pholder_window={};
            app.gatd_plot.max_overlap_window=max(dfnc_condi_with_window_n);
            
            for i=1:app.gatd_plot.max_overlap_window
                % Calculate design changing timing point
                temp_condi=double(dfnc_condi_with_window_n==i);
                window_switch_point=find(diff(temp_condi)~=0);
                window_time_point=[1,repelem(window_switch_point+0.5,2),dfnc_length];
                window_value_point=[temp_condi(1),temp_condi(1),repelem(temp_condi(window_switch_point+1),2)];
                window_value_point(end)=temp_condi(end);%=window_value_point(1:end-1);
                
                % Plot
                app.gatd_plot.pholder_window{i}=patch(app.UIAxes,[window_time_point,fliplr(window_time_point)],[window_value_point*(dfnc_reponse_max-dfnc_reponse_min+0.2)+dfnc_reponse_min-0.1,ones(1,length(window_time_point))*(dfnc_reponse_min-0.1)],[0.52 1 0.36],'LineStyle','none');
                app.gatd_plot.pholder_window{i}.FaceAlpha=0.1+i*0.8/app.gatd_plot.max_overlap_window;
                app.gatd_plot.pholder_window{i}.FaceColor=[0.4 1 0.3];
            end
            
            app.gatd_plot.pholder_hrf=plot(app.UIAxes,1:dfnc_length,app.gatd_setting.dfnc_reponse,'--r','LineWidth',2);
            app.gatd_plot.pholder_hrf.Color=[1 0.32 0.16];
            %app.gatd_plot.pholder_design=plot(app.UIAxes,1:dfnc_length,app.gatd_setting.dfnc_design,'k','LineWidth',1);
            app.gatd_plot.pholder_design=plot(app.UIAxes,design_time_point,design_value_point,'k','LineWidth',1);
            hold(app.UIAxes,'off');
            
            app.TaskDesignCheckBox.Enable = 'on';
            app.TemporalMaskCheckBox.Enable = 'on';
            app.HemodynamicResponseCheckBox.Enable = 'on';
        end

        % Value changed function: TaskDesignCheckBox
        function gatd_p_task(app, event)
            value = app.TaskDesignCheckBox.Value;
            switch value
                case 1
                    app.gatd_plot.pholder_design.Visible='on';
                case 0
                    app.gatd_plot.pholder_design.Visible='off';
            end
        end

        % Value changed function: HemodynamicResponseCheckBox
        function gatd_p_hrf(app, event)
            value = app.HemodynamicResponseCheckBox.Value;
            switch value
                case 1
                    app.gatd_plot.pholder_hrf.Visible='on';
                case 0
                    app.gatd_plot.pholder_hrf.Visible='off';
            end
        end

        % Value changed function: TemporalMaskCheckBox
        function gatd_p_window(app, event)
            value = app.TemporalMaskCheckBox.Value;
            switch value 
                case 1
                    for i=1:app.gatd_plot.max_overlap_window
                        app.gatd_plot.pholder_window{i}.Visible='on';
                    end
                case 0
                    for i=1:app.gatd_plot.max_overlap_window
                        app.gatd_plot.pholder_window{i}.Visible='off';
                    end
            end
        end

        % Button pushed function: SaveDesignMatrixButton
        function gatd_output_mat(app, event)
            [fnc_out_file,fnc_out_path] = uiputfile('*.mat','Select directory to save settings','dynamic_condition.mat');
            if isequal(fnc_out_file,0)
                return
            else
                window_condition=app.gatd_setting;
                temp_file_path=fullfile(fnc_out_path,fnc_out_file);
                save(temp_file_path,'window_condition');
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GATFD_design_UIFigure and hide until all components are created
            app.GATFD_design_UIFigure = uifigure('Visible', 'off');
            app.GATFD_design_UIFigure.Position = [100 100 600 680];
            app.GATFD_design_UIFigure.Name = 'GAT-FD - Task Design v0.2a';

            % Create DesignDurationSequencesLabel
            app.DesignDurationSequencesLabel = uilabel(app.GATFD_design_UIFigure);
            app.DesignDurationSequencesLabel.Position = [60 496 166 22];
            app.DesignDurationSequencesLabel.Text = 'Design Duration Sequence (s)';

            % Create DurationSequenceEditField
            app.DurationSequenceEditField = uieditfield(app.GATFD_design_UIFigure, 'text');
            app.DurationSequenceEditField.Position = [300 496 240 22];

            % Create DesignConditionSequence0forrestEditFieldLabel
            app.DesignConditionSequence0forrestEditFieldLabel = uilabel(app.GATFD_design_UIFigure);
            app.DesignConditionSequence0forrestEditFieldLabel.Position = [60 526 214 22];
            app.DesignConditionSequence0forrestEditFieldLabel.Text = 'Design Condition Sequence (0 for rest)';

            % Create DesignConditionSequenceEditField
            app.DesignConditionSequenceEditField = uieditfield(app.GATFD_design_UIFigure, 'text');
            app.DesignConditionSequenceEditField.Position = [300 526 240 22];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GATFD_design_UIFigure);
            title(app.UIAxes, 'Design Block')
            xlabel(app.UIAxes, 'Frame (TR)')
            ylabel(app.UIAxes, 'Condition')
            app.UIAxes.PlotBoxAspectRatio = [2.82295081967213 1 1];
            app.UIAxes.Position = [30 71 540 230];

            % Create TRsEditFieldLabel
            app.TRsEditFieldLabel = uilabel(app.GATFD_design_UIFigure);
            app.TRsEditFieldLabel.Position = [60 616 39 22];
            app.TRsEditFieldLabel.Text = 'TR(s) ';

            % Create TRsEditField
            app.TRsEditField = uieditfield(app.GATFD_design_UIFigure, 'numeric');
            app.TRsEditField.Position = [450 616 90 22];
            app.TRsEditField.Value = 1;

            % Create UpdateDesignButton
            app.UpdateDesignButton = uibutton(app.GATFD_design_UIFigure, 'push');
            app.UpdateDesignButton.ButtonPushedFcn = createCallbackFcn(app, @gatd_update_design, true);
            app.UpdateDesignButton.Position = [250 311 100 22];
            app.UpdateDesignButton.Text = 'Update Design';

            % Create SaveDesignMatrixButton
            app.SaveDesignMatrixButton = uibutton(app.GATFD_design_UIFigure, 'push');
            app.SaveDesignMatrixButton.ButtonPushedFcn = createCallbackFcn(app, @gatd_output_mat, true);
            app.SaveDesignMatrixButton.Position = [400 312 120 22];
            app.SaveDesignMatrixButton.Text = 'Save Design Matrix';

            % Create WindowSizeTREditFieldLabel
            app.WindowSizeTREditFieldLabel = uilabel(app.GATFD_design_UIFigure);
            app.WindowSizeTREditFieldLabel.Position = [60 586 102 22];
            app.WindowSizeTREditFieldLabel.Text = 'Window Size (TR)';

            % Create WindowSizeField
            app.WindowSizeField = uieditfield(app.GATFD_design_UIFigure, 'numeric');
            app.WindowSizeField.Position = [450 586 90 22];
            app.WindowSizeField.Value = 20;

            % Create StepSizeTREditFieldLabel
            app.StepSizeTREditFieldLabel = uilabel(app.GATFD_design_UIFigure);
            app.StepSizeTREditFieldLabel.Position = [60 556 84 22];
            app.StepSizeTREditFieldLabel.Text = 'Step Size (TR)';

            % Create StepSizeField
            app.StepSizeField = uieditfield(app.GATFD_design_UIFigure, 'numeric');
            app.StepSizeField.Position = [450 556 90 22];
            app.StepSizeField.Value = 1;

            % Create DefaulSettingButton
            app.DefaulSettingButton = uibutton(app.GATFD_design_UIFigure, 'push');
            app.DefaulSettingButton.ButtonPushedFcn = createCallbackFcn(app, @gatd_default, true);
            app.DefaulSettingButton.Position = [90 311 100 22];
            app.DefaulSettingButton.Text = 'Defaul Setting';

            % Create PlotOptionPanel
            app.PlotOptionPanel = uipanel(app.GATFD_design_UIFigure);
            app.PlotOptionPanel.Title = 'Plot Option';
            app.PlotOptionPanel.Position = [60 30 480 42];

            % Create TaskDesignCheckBox
            app.TaskDesignCheckBox = uicheckbox(app.PlotOptionPanel);
            app.TaskDesignCheckBox.ValueChangedFcn = createCallbackFcn(app, @gatd_p_task, true);
            app.TaskDesignCheckBox.Enable = 'off';
            app.TaskDesignCheckBox.Text = 'Task Design';
            app.TaskDesignCheckBox.Position = [20 0 87 22];
            app.TaskDesignCheckBox.Value = true;

            % Create HemodynamicResponseCheckBox
            app.HemodynamicResponseCheckBox = uicheckbox(app.PlotOptionPanel);
            app.HemodynamicResponseCheckBox.ValueChangedFcn = createCallbackFcn(app, @gatd_p_hrf, true);
            app.HemodynamicResponseCheckBox.Enable = 'off';
            app.HemodynamicResponseCheckBox.Text = 'Hemodynamic Response';
            app.HemodynamicResponseCheckBox.Position = [155 0 157 22];
            app.HemodynamicResponseCheckBox.Value = true;

            % Create TemporalMaskCheckBox
            app.TemporalMaskCheckBox = uicheckbox(app.PlotOptionPanel);
            app.TemporalMaskCheckBox.ValueChangedFcn = createCallbackFcn(app, @gatd_p_window, true);
            app.TemporalMaskCheckBox.Enable = 'off';
            app.TemporalMaskCheckBox.Text = 'Temporal Mask';
            app.TemporalMaskCheckBox.Position = [359 0 104 22];
            app.TemporalMaskCheckBox.Value = true;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GATFD_design_UIFigure);
            app.TabGroup.Position = [60 346 480 120];

            % Create AutomaticTab
            app.AutomaticTab = uitab(app.TabGroup);
            app.AutomaticTab.Title = 'Automatic';

            % Create CheckBoxActivation
            app.CheckBoxActivation = uicheckbox(app.AutomaticTab);
            app.CheckBoxActivation.ValueChangedFcn = createCallbackFcn(app, @gatd_s_activation, true);
            app.CheckBoxActivation.Text = '';
            app.CheckBoxActivation.Position = [15 68 25 22];
            app.CheckBoxActivation.Value = true;

            % Create EstimatedActivationLevelThreshold01EditFieldLabel
            app.EstimatedActivationLevelThreshold01EditFieldLabel = uilabel(app.AutomaticTab);
            app.EstimatedActivationLevelThreshold01EditFieldLabel.Position = [45 68 233 22];
            app.EstimatedActivationLevelThreshold01EditFieldLabel.Text = 'Estimated Activation Level Threshold (0-1)';

            % Create EstimatedActivationLevelThreshold01EditField
            app.EstimatedActivationLevelThreshold01EditField = uieditfield(app.AutomaticTab, 'numeric');
            app.EstimatedActivationLevelThreshold01EditField.Position = [390 68 80 22];
            app.EstimatedActivationLevelThreshold01EditField.Value = 0.8;

            % Create EstimatedActivationCoveragePercentageThresholdEditField
            app.EstimatedActivationCoveragePercentageThresholdEditField = uieditfield(app.AutomaticTab, 'numeric');
            app.EstimatedActivationCoveragePercentageThresholdEditField.Position = [390 38 80 22];
            app.EstimatedActivationCoveragePercentageThresholdEditField.Value = 80;

            % Create EstimatedActivationCoveragePercentageThresholdEditFieldLabel
            app.EstimatedActivationCoveragePercentageThresholdEditFieldLabel = uilabel(app.AutomaticTab);
            app.EstimatedActivationCoveragePercentageThresholdEditFieldLabel.Position = [45 39 317 22];
            app.EstimatedActivationCoveragePercentageThresholdEditFieldLabel.Text = 'Estimated Activation Coverage Percentage Threshold (%)';

            % Create CheckBoxCondition
            app.CheckBoxCondition = uicheckbox(app.AutomaticTab);
            app.CheckBoxCondition.ValueChangedFcn = createCallbackFcn(app, @gatd_s_condition, true);
            app.CheckBoxCondition.Text = '';
            app.CheckBoxCondition.Position = [15 8 25 22];
            app.CheckBoxCondition.Value = true;

            % Create ConditionCoveragePercentageThresholdEditFieldLabel
            app.ConditionCoveragePercentageThresholdEditFieldLabel = uilabel(app.AutomaticTab);
            app.ConditionCoveragePercentageThresholdEditFieldLabel.Position = [45 8 258 22];
            app.ConditionCoveragePercentageThresholdEditFieldLabel.Text = 'Condition Coverage Percentage Threshold (%)';

            % Create ConditionCoveragePercentageThresholdEditField
            app.ConditionCoveragePercentageThresholdEditField = uieditfield(app.AutomaticTab, 'numeric');
            app.ConditionCoveragePercentageThresholdEditField.Position = [390 8 80 22];
            app.ConditionCoveragePercentageThresholdEditField.Value = 80;

            % Create ManualTab
            app.ManualTab = uitab(app.TabGroup);
            app.ManualTab.Title = 'Manual';

            % Create StageConditionSequenceEditFieldLabel
            app.StageConditionSequenceEditFieldLabel = uilabel(app.ManualTab);
            app.StageConditionSequenceEditFieldLabel.Position = [15 38 148 22];
            app.StageConditionSequenceEditFieldLabel.Text = 'Stage Condition Sequence';

            % Create StageConditionSequenceEditField
            app.StageConditionSequenceEditField = uieditfield(app.ManualTab, 'text');
            app.StageConditionSequenceEditField.Tooltip = {'Use 1 for desired timepoints, 0 for unwanted timepoints)'};
            app.StageConditionSequenceEditField.Position = [220 38 240 22];

            % Create StageDurationSequenceTRsEditFieldLabel
            app.StageDurationSequenceTRsEditFieldLabel = uilabel(app.ManualTab);
            app.StageDurationSequenceTRsEditFieldLabel.Position = [15 11 176 22];
            app.StageDurationSequenceTRsEditFieldLabel.Text = 'Stage Duration Sequence (TRs)';

            % Create StageDurationSequenceTRsEditField
            app.StageDurationSequenceTRsEditField = uieditfield(app.ManualTab, 'text');
            app.StageDurationSequenceTRsEditField.Position = [220 11 240 22];

            % Create WiththecurrentsettingthetotalnumberofframesisLabel
            app.WiththecurrentsettingthetotalnumberofframesisLabel = uilabel(app.ManualTab);
            app.WiththecurrentsettingthetotalnumberofframesisLabel.Position = [15 67 295 22];
            app.WiththecurrentsettingthetotalnumberofframesisLabel.Text = 'With the current setting, the total number of frames is ';

            % Create totalstepsEditField
            app.totalstepsEditField = uieditfield(app.ManualTab, 'numeric');
            app.totalstepsEditField.Editable = 'off';
            app.totalstepsEditField.Position = [337 67 46 22];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.ManualTab, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @gatd_update_frames, true);
            app.CalculateButton.Position = [392 67 65 22];
            app.CalculateButton.Text = 'Calculate';

            % Create TaskConditionSpecificationLabel
            app.TaskConditionSpecificationLabel = uilabel(app.GATFD_design_UIFigure);
            app.TaskConditionSpecificationLabel.Position = [60 466 182 22];
            app.TaskConditionSpecificationLabel.Text = 'Task Condition Specification';

            % Show the figure after all components are created
            app.GATFD_design_UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = gatfd_ui_design

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.GATFD_design_UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GATFD_design_UIFigure)
        end
    end
end