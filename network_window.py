# network_window.py

from PyQt6.QtWidgets import (QWidget, QLabel, QPushButton,QVBoxLayout, QTextEdit, QListWidget, QMessageBox, QFileDialog, QComboBox, QCheckBox, QLineEdit, QGroupBox, QGridLayout, QApplication, QDialog, QProgressDialog)
from PyQt6.QtCore import Qt
from scipy.io import loadmat
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import sys
import numpy as np
import os
import bct  # bctpy package for Brain Connectivity Toolbox functions
import scipy.io as sio
import pandas as pd
def choose_node(node_list, parent=None):
    dialog = QDialog(parent)
    dialog.setWindowTitle("Select nodes for calculation")
    dialog.setGeometry(100, 100, 292, 714)

    layout = QVBoxLayout(dialog)

    label = QLabel("Select Nodes for Calculation", dialog)
    label.setAlignment(Qt.AlignmentFlag.AlignCenter)
    layout.addWidget(label)

    listbox = QListWidget(dialog)
    listbox.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
    listbox.addItems([str(node) for node in node_list])
    layout.addWidget(listbox)

    done_button = QPushButton("Done", dialog)
    layout.addWidget(done_button)

    selected_nodes = []

    def on_done():
        nonlocal selected_nodes
        selected_nodes = [item.text() for item in listbox.selectedItems()]
        dialog.accept()

    done_button.clicked.connect(on_done)

    if dialog.exec():
        return selected_nodes
    else:
        return []
class NetworkWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.createComponents()
        self.startupFcn()
    #  Button pushed function: Button_loadfiles
    def calc_global_efficiency(self, temp_net):
        return bct.efficiency_bin(temp_net, 0)

    def calc_network_local_efficiency(self, temp_net):
        return float(np.mean(bct.efficiency_bin(temp_net, 1)))

    def calc_network_clustering_coefficient(self, temp_net):
        return float(np.mean(bct.clustering_coef_bu(temp_net)))

    def calc_network_average_degree(self, temp_net):
        return float(np.mean(bct.degrees_und(temp_net)))

    def calc_network_characteristic_path(self, temp_net):
        D = bct.distance_bin(temp_net)
        return bct.charpath(D, 0, 0)[0]

    def calc_sw_norm_cc(self, temp_net):
        n_nodes = temp_net.shape[0]
        n_edges = int(np.round(temp_net.sum() / 2))
        rand_cc = 0
        for _ in range(20):
            r = bct.makerandCIJ_und(n_nodes, n_edges)
            rand_cc += np.mean(bct.clustering_coef_bu(r))
        rand_cc /= 20
        return float(np.mean(bct.clustering_coef_bu(temp_net)) / rand_cc)

    def calc_sw_norm_pl(self, temp_net):
        n_nodes = temp_net.shape[0]
        n_edges = int(np.round(temp_net.sum() / 2))
        rand_pl = 0
        for _ in range(20):
            r = bct.makerandCIJ_und(n_nodes, n_edges)
            D = bct.distance_bin(r)
            rand_pl += bct.charpath(D, 0, 0)[0]
        rand_pl /= 20
        D0 = bct.distance_bin(temp_net)
        return float(bct.charpath(D0, 0, 0)[0] / rand_pl)

    def calc_sw_coefficient(self, temp_net):
        n_nodes = temp_net.shape[0]
        n_edges = int(np.round(temp_net.sum() / 2))
        rand_cc = rand_pl = 0
        for _ in range(20):
            r = bct.makerandCIJ_und(n_nodes, n_edges)
            rand_cc += np.mean(bct.clustering_coef_bu(r))
            D = bct.distance_bin(r)
            rand_pl += bct.charpath(D, 0, 0)[0]
        rand_cc /= 20
        rand_pl /= 20
        num_cc = np.mean(bct.clustering_coef_bu(temp_net))
        D0 = bct.distance_bin(temp_net)
        num_pl = bct.charpath(D0, 0, 0)[0]
        return float((num_cc / rand_cc) / (num_pl / rand_pl))

    def calc_modularity(self, temp_net):
        _, mod = bct.modularity_und(temp_net)
        return mod

    def calc_transitivity(self, temp_net):
        return float(bct.transitivity_bu(temp_net))

    def calc_assortativity(self, temp_net):
        return float(bct.assortativity_bin(temp_net, 0))

    def calc_nodal_efficiency(self, temp_net):
        A = temp_net.copy()
        np.fill_diagonal(A, 0)
        A = (A != 0).astype(float)
        # use bctpy efficiency for local
        return bct.efficiency_bin(A, 1)

    def calc_nodal_local_efficiency(self, temp_net):
        return bct.efficiency_bin(temp_net, 1)

    def calc_nodal_clustering_coefficient(self, temp_net):
        return bct.clustering_coef_bu(temp_net)

    def calc_nodal_degree(self, temp_net):
        return bct.degrees_und(temp_net)

    def calc_nodal_betweenness(self, temp_net):
        return bct.betweenness_bin(temp_net)
    def gatn_loadfile(self):
        # Open file dialog with multi-selection enabled
        files, _ = QFileDialog.getOpenFileNames(self, "Select One or More Files", "", "MAT files (*.mat)")
        if not files:
            print("Selection Canceled")
            return
        # Extract filenames and directory path
        file_paths = files
        file_names = [os.path.basename(f) for f in file_paths]
        file_dir = os.path.dirname(file_paths[0])
        # Save to app settings
        self.gatn_setting['file_list'] = file_names
        self.gatn_setting['file_path_list'] = file_dir
        # Update loaded files textarea
        self.textarea_loaded_files.setPlainText('\n'.join(file_names))

        # Load first file
        first_file = file_paths[0]
        temp_feature = loadmat(first_file, squeeze_me=True, struct_as_record=False)
        subj_data = temp_feature.get('subj_data', None)
        node_list = []
        d_atlas_list = subj_data.d_atlas_list
        node_list = d_atlas_list.tolist()

        # Store selected nodes initially as full list
        self.gatn_setting['node_list'] = node_list
        self.gatn_setting['node_list_selected'] = node_list
        # Update list widget
        self.list_nodes.clear()
        self.list_nodes.addItems(node_list)
    
    def gatn_loadconditions(self):
        # Open file dialog to select .mat file
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Conditions", "", "MAT files (*.mat)")
        if not file_path:
            self.gatn_setting['iscondition'] = 0
            return
        mat_data = loadmat(file_path, squeeze_me=True, struct_as_record=False)
        window_condition = mat_data.get('window_condition', None)
        if window_condition is None:
            QMessageBox.warning(self, "Error", "'window_condition' not found in file.")
            return
        if isinstance(window_condition, np.ndarray):
            window_condition = window_condition.item()
        dfnc_window_condi = None
        if hasattr(window_condition, 'dfnc_window_condi'):
            dfnc_window_condi = window_condition.dfnc_window_condi
        else:
            # Print available fields for debugging
            print("Available fields:", dir(window_condition))
            QMessageBox.warning(self, "Error", "'dfnc_window_condi' not found in window_condition.")
            return
        self.gatn_setting['condition'] = dfnc_window_condi
        self.gatn_setting['iscondition'] = 1

        # Plot on UIAxes (matplotlib)
        ax = self.figure.axes[0]
        ax.clear()
        ax.set_title('State Conditions')
        ax.set_xlabel('Time')
        ax.set_ylabel('Condition')
        ax.plot(np.array(dfnc_window_condi).flatten())
        self.canvas.draw()
    
    def gatn_run(self):
        # Check if loaded file
        if not self.gatn_setting['file_list']:
            QMessageBox.warning(self, "Error", "No files loaded.")
            return

        # Update Settings from UI
        self.gatn_measure['glob_measure'][0] = int(self.cb_global_eff.isChecked())
        self.gatn_measure['glob_measure'][1] = int(self.cb_local_eff.isChecked())
        self.gatn_measure['glob_measure'][2] = int(self.cb_clustering.isChecked())
        self.gatn_measure['glob_measure'][3] = int(self.cb_degree.isChecked())
        self.gatn_measure['glob_measure'][4] = int(self.cb_char_path.isChecked())
        self.gatn_measure['glob_measure'][5] = int(self.cb_small_world.isChecked())
        self.gatn_measure['glob_measure'][6] = int(self.cb_norm_clustering.isChecked())
        self.gatn_measure['glob_measure'][7] = int(self.cb_norm_path.isChecked())
        self.gatn_measure['glob_measure'][8] = int(self.cb_transitivity.isChecked())
        self.gatn_measure['glob_measure'][9] = int(self.cb_assortativity.isChecked())
        self.gatn_measure['glob_measure'][10] = int(self.cb_modularity.isChecked())
        meas_glob_list = [i for i, v in enumerate(self.gatn_measure['glob_measure']) if v]
        run_if_glob = bool(meas_glob_list)

        self.gatn_measure['nod_measure'][0] = int(self.cb_nodal_global_eff.isChecked())
        self.gatn_measure['nod_measure'][1] = int(self.cb_nodal_local_eff.isChecked())
        self.gatn_measure['nod_measure'][2] = int(self.cb_nodal_clustering.isChecked())
        self.gatn_measure['nod_measure'][3] = int(self.cb_nodal_degree.isChecked())
        self.gatn_measure['nod_measure'][4] = int(self.cb_nodal_betweenness.isChecked())
        meas_nod_list = [i for i, v in enumerate(self.gatn_measure['nod_measure']) if v]
        run_if_nod = bool(meas_nod_list) and bool(self.gatn_setting['node_list_selected'])

        if not (run_if_glob or run_if_nod):
            QMessageBox.warning(self, "Error", "No measures selected.")
            return

        # Output file dialog
        out_file, _ = QFileDialog.getSaveFileName(self, "Select directory to save settings", "network_properties.csv", "CSV files (*.csv)")
        if not out_file:
            return

        fnc_filelength = len(self.gatn_setting['file_list'])

        # Threshold settings
        try:
            threshold_lower = float(self.edit_lower.text())
            threshold_upper = float(self.edit_upper.text())
            threshold_step = float(self.edit_step.text())
        except Exception:
            QMessageBox.warning(self, "Error", "Invalid threshold values.")
            return
        threshold_method = self.dropdown_method.currentIndex()  # 0:Absolute, 1:Proportional, 2:Cost
        threshold_absolute = int(self.cb_absolute.isChecked())
        parallel = int(self.cb_parallel.isChecked())
        self.gatn_setting['threshold_lower'] = threshold_lower
        self.gatn_setting['threshold_upper'] = threshold_upper
        self.gatn_setting['threshold_step'] = threshold_step
        self.gatn_setting['threshold_method'] = threshold_method
        self.gatn_setting['threshold_absolute'] = threshold_absolute
        self.gatn_setting['parallel'] = parallel

        # Load first file for node info and window size
        first_file = os.path.join(self.gatn_setting['file_path_list'], self.gatn_setting['file_list'][0])
        temp_feature = loadmat(first_file, squeeze_me=True, struct_as_record=False)
        subj_data = temp_feature.get('subj_data', None)
        if subj_data is None:
            QMessageBox.warning(self, "Error", "subj_data not found in first file.")
            return
        d_atlas_list = subj_data.d_atlas_list
        node_list = d_atlas_list.tolist()
        data_num_nodes = len(node_list)
        d_corr = subj_data.d_corr
        data_window_size = d_corr.shape[0]

        # Threshold list
        dnet_data_threshold_list = np.arange(threshold_lower, threshold_upper + threshold_step/2, threshold_step)
        data_num_thres_steps = len(dnet_data_threshold_list)

        # Output arrays
        if run_if_glob:
            local_glob_count = len(meas_glob_list)
            dnet_data_data_mat_global = np.zeros((data_window_size, local_glob_count, data_num_thres_steps, fnc_filelength))
        if run_if_nod:
            local_nod_count = len(meas_nod_list)
            dnet_data_data_mat_nodal = np.zeros((data_window_size, local_nod_count, data_num_thres_steps, len(self.gatn_setting['node_list_selected']), fnc_filelength))

        # Save meta info
        dnet_data_measures_glob = [self.gatn_measure['glob_name'][i] for i in meas_glob_list] if run_if_glob else []
        dnet_data_measures_nod = [self.gatn_measure['nod_name'][i] for i in meas_nod_list] if run_if_nod else []
        dnet_data_files = self.gatn_setting['file_list']
        dnet_data_nodes_list = self.gatn_setting['node_list_selected']

        # Progress dialog
        progress = QProgressDialog("Calculating Network", "Cancel", 0, fnc_filelength, self)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.show()
    
        # Main calculation loop
        for idx_file, file_name in enumerate(self.gatn_setting['file_list']):
            if progress.wasCanceled():
                break
            file_path = os.path.join(self.gatn_setting['file_path_list'], file_name)
            temp_feature = loadmat(file_path, squeeze_me=True, struct_as_record=False)
            subj_data = temp_feature.get('subj_data', None)
            if subj_data is None:
                continue
            d_corr = subj_data.d_corr
            # For each threshold
            for idx_thres, thres in enumerate(dnet_data_threshold_list):
                # For each window
                for idx_win in range(data_window_size):
                    temp_net = d_corr[idx_win]
                    temp_net = np.nan_to_num(temp_net, nan=0.0)
                    # Absolute thresholding
                    if threshold_absolute == 1:
                        temp_net = np.abs(temp_net)
                    # Thresholding methods
                    if threshold_method == 2:  # Cost (percentile)
                        temp_net_thresh = np.percentile(temp_net, (1-thres)*100)
                    elif threshold_method == 1:  # Proportional
                        temp_net_thresh = np.max(temp_net) * thres
                    else:  # Absolute
                        temp_net_thresh = thres
                    temp_net = np.where(temp_net >= temp_net_thresh, temp_net, 0)
                    # Binarization
                    temp_net[temp_net != 0] = 1
                    # Global measures
                    if run_if_glob:
                        for i, measure_idx in enumerate(meas_glob_list):
                            func = self.gatn_measure['glob_func'][measure_idx]
                            try:
                                val = func(temp_net)
                            except Exception:
                                val = np.nan
                            dnet_data_data_mat_global[idx_win, i, idx_thres, idx_file] = val
                    # Nodal measures
                    if run_if_nod:
                        for i, measure_idx in enumerate(meas_nod_list):
                            func = self.gatn_measure['nod_func'][measure_idx]
                            try:
                                vals = func(temp_net)
                                node_indices = [node_list.index(n) for n in self.gatn_setting['node_list_selected']]
                                vals_selected = np.array(vals)[node_indices]
                            except Exception:
                                vals_selected = np.full(len(self.gatn_setting['node_list_selected']), np.nan)
                            dnet_data_data_mat_nodal[idx_win, i, idx_thres, :, idx_file] = vals_selected
            progress.setValue(idx_file + 1)

        # Save output
        out_mat_path = out_file + '.mat'
        sio.savemat(out_mat_path, {
            'dnet_data_threshold_list': dnet_data_threshold_list,
            'dnet_data_nodes_list': dnet_data_nodes_list,
            'dnet_data_measures_glob': dnet_data_measures_glob,
            'dnet_data_measures_nod': dnet_data_measures_nod,
            'dnet_data_files': dnet_data_files,
            'dnet_data_data_mat_global': dnet_data_data_mat_global if run_if_glob else None,
            'dnet_data_data_mat_nodal': dnet_data_data_mat_nodal if run_if_nod else None
        }, do_compression=True)

        # Optionally save CSV if condition is loaded
        if self.gatn_setting['iscondition'] == 1:
            # Example: flatten and save global measures
            if run_if_glob:
                flat_global = dnet_data_data_mat_global.reshape(-1, dnet_data_data_mat_global.shape[-1])
                df_global = pd.DataFrame(flat_global)
                df_global.to_csv(out_file, index=False)

        progress.close()
        QMessageBox.information(self, "Finished", "Network calculation finished.")

    def startupFcn(self):
        self.gatn_setting = {
            'threshold_lower': 0.1,
            'threshold_upper': 0.4,
            'threshold_step': 0.02,
            'threshold_method': 1,  # 0:cost, 1:absolute, 2:relative
            'threshold_absolute': 0,
            'parallel': 0,
            'file_list': [],
            'file_path_list': [],
            'condition' : [],
            'node_list': [],
            'node_list_unselected': [],
            'node_list_selected': [],
            'iscondition': 0,   
        }    
        self.gatn_measure = {}
        self.gatn_measure['glob_name'] = [
            'global_efficiency', 'network_local_efficiency', 'network_clustering_coefficient',
            'network_average_degree', 'characteristic_path_length', 'small_world_coefficient',
            'normalized_clustering_coefficient', 'normalized_path_length', 'transitivity',
            'assortativity', 'modularity'
        ]
        self.gatn_measure['glob_label'] = [
            'glo_eff', 'glo_lef', 'glo_clc', 'glo_deg', 'glo_cpl', 'sw_coef',
            'norm_cc', 'norm_pl', 'glo_tra', 'glo_ast', 'glo_mod'
        ]
        self.gatn_measure['glob_count'] = len(self.gatn_measure['glob_name'])
        self.gatn_measure['glob_measure'] = [0] * self.gatn_measure['glob_count']
        self.gatn_measure['glob_func'] = [
            self.calc_global_efficiency,
            self.calc_network_local_efficiency,
            self.calc_network_clustering_coefficient,
            self.calc_network_average_degree,
            self.calc_network_characteristic_path,
            self.calc_sw_coefficient,
            self.calc_sw_norm_cc,
            self.calc_sw_norm_pl,
            self.calc_transitivity,
            self.calc_assortativity,
            self.calc_modularity
        ]
        # Nodal
        self.gatn_measure['nod_name'] = [ 'nodal_efficiency', 'nodal_local_efficiency',  'nodal_clustering_coefficient', 'nodal_degree', 'nodal_betweenness' ]
        self.gatn_measure['nod_label'] = [ 'nod_eff', 'nod_lef', 'nod_clc', 'nod_deg', 'nod_bet' ]
        self.gatn_measure['nod_count'] = len(self.gatn_measure['nod_name'])
        self.gatn_measure['nod_measure'] = [0] * self.gatn_measure['nod_count']
        self.gatn_measure['nod_func'] = [
            self.calc_nodal_efficiency,
            self.calc_nodal_local_efficiency,
            self.calc_nodal_clustering_coefficient,
            self.calc_nodal_degree,
            self.calc_nodal_betweenness
        ]

    #    Button pushed function: SelectNodesButton
    def gatn_select_nodes(self):
        node_list = self.gatn_setting['node_list']
        selected = choose_node(node_list, self)
        if selected:
            self.gatn_setting['node_list_selected'] = selected
            self.list_nodes.clear()
            self.list_nodes.addItems(selected)
    def createComponents(self):
        self.gatn_setting = {}
        self.setWindowTitle('GAT-FD - Network Property Calculation v0.2a')
        self.setGeometry(100, 100, 650, 788)
        # Main layout
        layout = QGridLayout(self)

        # --- Loaded Files section ---
        self.label_loaded_files = QLabel('Loaded Files:', self)
        layout.addWidget(self.label_loaded_files, 0, 0)
        self.textarea_loaded_files = QTextEdit(self)
        self.textarea_loaded_files.setReadOnly(True)
        self.textarea_loaded_files.setPlainText('Loaded Files')
        layout.addWidget(self.textarea_loaded_files, 1, 0, 3, 1)
        self.button_loadfiles = QPushButton('Load Files', self)
        self.button_loadfiles.clicked.connect(self.gatn_loadfile)
        layout.addWidget(self.button_loadfiles, 4, 0)

        # --- State Conditions plot ---
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        ax = self.figure.add_subplot(111)
        ax.set_title('State Conditions')
        ax.set_xlabel('Time')
        ax.set_ylabel('Condition')
        layout.addWidget(self.canvas, 0, 1, 3, 1)
        self.button_load_mask = QPushButton('Load Temporal Mask', self)
        self.button_load_mask.clicked.connect(self.gatn_loadconditions)
        layout.addWidget(self.button_load_mask, 4, 1)

        # --- Thresholding section ---
        self.group_threshold = QGroupBox('Thresholding', self)
        grid_thresh = QGridLayout()
        grid_thresh.addWidget(QLabel('Thresholding Method', self), 0, 0)
        self.dropdown_method = QComboBox(self)
        self.dropdown_method.addItems(['Absolute', 'Proportional', 'Cost'])
        self.dropdown_method.setCurrentText('Cost')
        grid_thresh.addWidget(self.dropdown_method, 0, 1)
        grid_thresh.addWidget(QLabel('Threshold Range', self), 1, 0)
        self.edit_lower = QLineEdit(self)
        self.edit_lower.setText('0.1')
        self.edit_upper = QLineEdit(self)
        self.edit_upper.setText('0.4')
        grid_thresh.addWidget(self.edit_lower, 1, 1)
        grid_thresh.addWidget(QLabel('-', self), 1, 2)
        grid_thresh.addWidget(self.edit_upper, 1, 3)
        grid_thresh.addWidget(QLabel('Threshold Step', self), 2, 0)
        self.edit_step = QLineEdit(self)
        self.edit_step.setText('0.02')
        grid_thresh.addWidget(self.edit_step, 2, 1)
        self.cb_absolute = QCheckBox('Use absolute value of correlation coefficient', self)
        grid_thresh.addWidget(self.cb_absolute, 3, 0, 1, 4)
        self.group_threshold.setLayout(grid_thresh)
        layout.addWidget(self.group_threshold, 5, 0, 1, 2)

        # --- Network Average (Global) Measures ---
        self.group_global = QGroupBox('Network Average (Global) Measures', self)
        grid_global = QGridLayout()
        self.cb_global_eff = QCheckBox('Global Efficiency', self)
        self.cb_local_eff = QCheckBox('Local Efficiency', self)
        self.cb_clustering = QCheckBox('Clustering Coefficient', self)
        self.cb_degree = QCheckBox('Degree', self)
        self.cb_small_world = QCheckBox('Small World Coefficient', self)
        self.cb_modularity = QCheckBox('Modularity Coefficient', self)
        self.cb_norm_clustering = QCheckBox('Normalized Clustering Coefficient', self)
        self.cb_norm_path = QCheckBox('Normalized Path Length', self)
        self.cb_transitivity = QCheckBox('Transitivity Coefficient', self)
        self.cb_assortativity = QCheckBox('Assortativity Coefficient', self)
        self.cb_char_path = QCheckBox('Characteristic Path Length', self)
        grid_global.addWidget(self.cb_global_eff, 0, 0)
        grid_global.addWidget(self.cb_local_eff, 1, 0)
        grid_global.addWidget(self.cb_clustering, 2, 0)
        grid_global.addWidget(self.cb_degree, 3, 0)
        grid_global.addWidget(self.cb_small_world, 0, 1)
        grid_global.addWidget(self.cb_modularity, 1, 1)
        grid_global.addWidget(self.cb_norm_clustering, 2, 1)
        grid_global.addWidget(self.cb_norm_path, 3, 1)
        grid_global.addWidget(self.cb_transitivity, 0, 2)
        grid_global.addWidget(self.cb_assortativity, 1, 2)
        grid_global.addWidget(self.cb_char_path, 2, 2)
        self.group_global.setLayout(grid_global)
        layout.addWidget(self.group_global, 6, 0, 1, 2)

        # --- Nodal Measures ---
        self.group_nodal = QGroupBox('Nodal Measures', self)
        grid_nodal = QGridLayout()
        self.button_select_nodes = QPushButton('Select Nodes', self)
        self.button_select_nodes.clicked.connect(self.gatn_select_nodes)
        grid_nodal.addWidget(self.button_select_nodes, 0, 0)
        self.list_nodes = QListWidget(self)
        self.list_nodes.addItem('Selected Nodes')
        grid_nodal.addWidget(self.list_nodes, 1, 0, 4, 1)
        self.cb_nodal_global_eff = QCheckBox('Global Efficiency', self)
        self.cb_nodal_local_eff = QCheckBox('Local Efficiency', self)
        self.cb_nodal_clustering = QCheckBox('Clustering Coefficient', self)
        self.cb_nodal_degree = QCheckBox('Degree', self)
        self.cb_nodal_betweenness = QCheckBox('Betweenness Centrality', self)
        grid_nodal.addWidget(self.cb_nodal_global_eff, 0, 1)
        grid_nodal.addWidget(self.cb_nodal_local_eff, 1, 1)
        grid_nodal.addWidget(self.cb_nodal_clustering, 2, 1)
        grid_nodal.addWidget(self.cb_nodal_degree, 3, 1)
        grid_nodal.addWidget(self.cb_nodal_betweenness, 4, 1)
        self.group_nodal.setLayout(grid_nodal)
        layout.addWidget(self.group_nodal, 7, 0, 1, 2)

        # --- Run in parallel and Calculate button ---
        self.cb_parallel = QCheckBox('Run in parallel', self)
        layout.addWidget(self.cb_parallel, 8, 0)
        self.button_calculate = QPushButton('Calculate Network Properties', self)
        self.button_calculate.clicked.connect(self.gatn_run)
        layout.addWidget(self.button_calculate, 8, 1)
