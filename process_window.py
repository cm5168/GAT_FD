# process_window.py
import os
import numpy as np
import nibabel as nib
from scipy.ndimage import affine_transform
from scipy.io import savemat, loadmat
from scipy.signal import butter, filtfilt
from PyQt6.QtWidgets import QWidget, QLabel,QComboBox, QFileDialog, QPushButton, QListWidget, QProgressDialog,QLineEdit, QFormLayout, QMessageBox
from PyQt6.QtGui import QIntValidator 
from nilearn.input_data import NiftiLabelsMasker
from nilearn.image import resample_to_img
from nibabel.orientations import io_orientation, axcodes2ornt, ornt_transform, apply_orientation

class ProcessWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.createComponents()
        self.startupFcn()
      
    # Reset default value. (This only update the parameters not the field)
    def parameters_default(self):
        self.settings = {
        'file_type': 0,
        'file_list': [],
        'file_path':[],
        'atlas_type': 0,
        'atlas_list': [],
        'atlas_path': [],
        'TR' : 1,
        'window_size': 20,
        'step_size': 1,
        'filter_type': 0,
        'filter_setting1': '',
        'filter_setting2': '',
        'kernel': 0,
        'kernel_setting': ''
        }

    # Update parameters from field
    def parameters_update_setting(self):
        self.settings['TR'] = int(self.TRs_input.text()) if self.TRs_input.text().isdigit() else 1
        self.settings['window_size'] = int(self.window_size_input.text()) if self.window_size_input.text().isdigit() else 20
        self.settings['step_size'] = int(self.step_size_input.text()) if self.step_size_input.text().isdigit() else 1
        self.settings['filter_type'] = self.filter_input.currentIndex()
        self.settings['filter_setting1'] = self.filter_edit_field_1.text() if self.filter_edit_field_1.text() else ""
        self.settings['filter_setting2'] = self.filter_edit_field_2.text() if self.filter_edit_field_2.text() else ""
        self.settings['kernel'] = self.window_kernel_input.currentIndex()
        self.settings['kernel_setting'] = self.sigma_edit_field.text() if self.sigma_edit_field.text() else ""

    # Update field from parameters
    def parameters_update_field(self):
        # File dropdown
        self.file_format_dropdown.setCurrentIndex(self.settings.get('file_type', 0))
        index = self.settings.get('file_list', 0)
        if isinstance(index, list):
            index = index[0] if index else 0  # get first item or default to 0
        self.file_list_widget.setCurrentRow(index)
        # Atlas dropdown
        self.atlas_input.setCurrentIndex(self.settings.get('atlas_type', 0))
        index = self.settings.get('atlas_list', 0)
        try:
            self.atlas_list_widget.setCurrentRow(int(index))
        except (ValueError, TypeError):
            self.atlas_list_widget.setCurrentRow(0)

        self.TRs_input.setText(str(self.settings.get('TR', 1)))
        self.window_size_input.setText(str(self.settings.get('window_size', 20)))
        self.step_size_input.setText(str(self.settings.get('step_size', 1)))
        self.filter_input.setCurrentIndex(self.settings.get('filter_type', 0))
        self.filter_edit_field_1.setText(str(self.settings.get('filter_setting1', 0)))
        self.filter_edit_field_2.setText(str(self.settings.get('filter_setting2', 0)))
        self.window_kernel_input.setCurrentIndex(self.settings.get('kernel', 0))
        self.sigma_edit_field.setText(str(self.settings.get('kernel_setting', 0)))
    
    # Gaussian Kernel Function
    def gauss_kernel(g_x, g_m, g_s):
        result = np.exp(-((g_x - g_m) ** 2) / (2 * g_s ** 2))
        return result
    
    # Code that executes after component creation
    def startupFcn(self):
        self.parameters_default()
        # initialize the atlas dropdown interface (doesn't affect much the function, just the entrance appearace)
        path = os.path.join(self.settings.get('path', ''), 'atlas', 'BN_atlas_1_25mm.nii.txt')
        if os.path.exists(path):
            with open(path, 'r') as f:
                roi_names = [line.strip() for line in f.readlines()]
            self.settings['atlas_list'] = roi_names # can be replace with self.atlas_list (actually self.settings['___'] = self.___ )
        else:
            roi_names = ['[File not found]']
        self.settings['atlas_path'] = []
        self.atlas_list_widget.clear()
        self.atlas_list_widget.addItems(roi_names)

    # Value changed function: FileFormatDropDown
    def gatp_file_type(self):
        self.file_type = self.file_format_dropdown.currentIndex() 
        if self.file_type == 0:
            # Enable atlas dropdown, disable atlas list and button
            self.atlas_input.setEnabled(True)
            self.atlas_list_widget.setEnabled(False)
            self.add_atlas_button.setEnabled(False)
            self.add_atlas_button.setText("Add Masks")
            self.atlas_input.setCurrentIndex(1)  # Default to Brainnetome
            # Load Brainnetome ROI names
            path = os.path.join(self.settings.get('path', ''), 'atlas', 'BN_atlas_1_25mm.nii.txt')
            if os.path.exists(path):
                with open(path, 'r') as f:
                    roi_names = [line.strip() for line in f.readlines()]
                self.settings['atlas_list'] = roi_names
            else:
                roi_names = ['[File not found]']
            self.atlas_list_widget.clear()
            self.atlas_list_widget.addItems(roi_names)
            # Clear file list
            self.settings['file_list'] = []
            self.settings['file_path_list'] = []
            self.file_list_widget.clear()
        elif self.file_type == 1:
            # Disable atlas dropdown, enable atlas list and button
            self.atlas_input.setEnabled(False)
            self.atlas_list_widget.setEnabled(True)
            self.add_atlas_button.setEnabled(True)
            self.add_atlas_button.setText("Add Labels")
            # Clear atlas and file lists
            self.settings['atlas_list'] = []
            self.atlas_list_widget.clear()
            self.settings['file_list'] = []
            self.settings['file_path_list'] = []
            self.file_list_widget.clear()
    # Open image or matrix file function  
    def gatp_open_file(self):
        self.file_type = self.file_format_dropdown.currentIndex() 
        if self.file_type == 0: #image
            filter_str = "NIfTI Files (*.nii *.nii.gz)"
        elif self.file_type == 1: #matrix
            filter_str = "MATLAB Matrix Files (*.mat)"
        else:
            return
        # Open file dialog
        files, _ = QFileDialog.getOpenFileNames(self, "Select One or More Files", "", filter_str)
        if not files:
            return
        self.settings['file_list'] = files
        self.settings['file_path'] = [f.rsplit('/', 1)[0] for f in files]
        self.file_list_widget.clear()
        self.file_list_widget.addItems([f.rsplit('/', 1)[-1] for f in files])
    
    # Value changed function: AtlasDropDown
    def gatp_atlas_type(self):
        value = self.atlas_input.currentIndex()
        if value ==0:
            self.atlas_list_widget.setEnabled(False)
            self.add_atlas_button.setEnabled(False)
            path = os.path.join(self.settings.get('path', ''), 'atlas', 'BN_atlas_1_25mm.nii.txt')
            if os.path.exists(path):
                with open(path, 'r') as f:
                    roi_names = [line.strip() for line in f.readlines()]
                self.settings['atlas_list'] = roi_names
            else:
                roi_names = ['[File not found]']
            self.settings['atlas_path'] = []
            self.atlas_list_widget.clear()
            self.atlas_list_widget.addItems(roi_names)
        elif value == 1:
            self.atlas_list_widget.setEnabled(False)
            self.add_atlas_button.setEnabled(False)
            path = os.path.join(self.settings.get('path', ''), 'atlas', 'AAL2v1_2mm.nii.txt')
            if os.path.exists(path):
                with open(path, 'r') as f:
                    roi_names = [line.strip() for line in f.readlines()]
                self.settings['atlas_list'] = roi_names
            else:
                roi_names = ['[File not found]']
            self.settings['atlas_path'] = []
            self.atlas_list_widget.clear()
            self.atlas_list_widget.addItems(roi_names)

        elif value == 2:
            # Load AAL ROI names from file
            self.atlas_list_widget.setEnabled(True)
            self.add_atlas_button.setEnabled(True)
            self.settings['atlas_list'] = []
            self.settings['atlas_path'] = []
            self.atlas_list_widget.clear()
        
    # Button pushed function: Button_AddAtlas
    def gatp_open_atlas(self):
        self.file_type = self.file_format_dropdown.currentIndex() 
        if self.file_type == 0:
            # If input is image: allow multiple .nii or .nii.gz files
            files, _ = QFileDialog.getOpenFileNames(self, "Select One or More Atlas Files", "", "NIfTI Files (*.nii *.nii.gz)")
            if not files:
                return
            self.settings['atlas_list'] = files
            self.settings['atlas_path'] = [f.rsplit('/', 1)[0] for f in files]
            self.atlas_list_widget.clear()
            self.atlas_list_widget.addItems([f.rsplit('/', 1)[-1] for f in files])

        elif self.file_type == 1:
            # If input is .txt ROI names
            file, _ = QFileDialog.getOpenFileName(self, "Select ROI Name File", "", "Text Files (*.txt)")
            if not file:
                return
            with open(file, 'r') as f:
                roi_names = [line.strip() for line in f.readlines()]
            self.settings['atlas_list'] = roi_names
            self.settings['atlas_path'] = file.rsplit('/', 1)[0]
            self.atlas_list_widget.clear()
            self.atlas_list_widget.addItems(roi_names)
    # Value changed function: WindowKernelDropDown
    def gatp_kernel(self):
        self.kernel = self.window_kernel_input.currentIndex()   
        # No kernel
        if self.kernel == 0:  
            self.sigma_edit_field.setEnabled(False)
            self.sigma_label.setEnabled(False)
        # Gaussian kernel
        if self.kernel == 1:  
            self.sigma_edit_field.setEnabled(True)
            self.sigma_label.setEnabled(True)

    # Value changed function: FilterDropDown
    def gatp_filter_type(self):
        self.filter_type = self.filter_input.currentIndex()  
        if self.filter_type == 0:  # None
            self.filter_label_2.setEnabled(False)
            self.filter_edit_field_2.setEnabled(False)
            self.filter_label_1.setEnabled(False)
            self.filter_edit_field_1.setEnabled(False)
            self.filter_label_2.setText("Filter Option")
            self.filter_label_1.setText("Filter Option")

        elif self.filter_type == 1:  # Bandpass
            self.filter_label_2.setEnabled(True)
            self.filter_edit_field_2.setEnabled(True)
            self.filter_label_1.setEnabled(True)
            self.filter_edit_field_1.setEnabled(True)
            self.filter_label_2.setText("Lower Limit (1/Hz)")
            self.filter_label_1.setText("Upper Limit (1/Hz)")

        elif self.filter_type == 2:  # Highpass
            self.filter_label_2.setEnabled(False)
            self.filter_edit_field_2.setEnabled(False)
            self.filter_label_1.setEnabled(True)
            self.filter_edit_field_1.setEnabled(True)
            self.filter_label_2.setText("Filter Option")
            self.filter_label_1.setText("Cut-off Frequency (1/Hz)")

        elif self.filter_type == 3:  # Lowpass
            self.filter_label_2.setEnabled(False)
            self.filter_edit_field_2.setEnabled(False)
            self.filter_label_1.setEnabled(True)
            self.filter_edit_field_1.setEnabled(True)
            self.filter_label_2.setText("Filter Option")
            self.filter_label_1.setText("Cut-off Frequency (1/Hz)")

        elif self.filter_type == 4:  # Wavelet
            self.filter_label_2.setEnabled(True)
            self.filter_edit_field_2.setEnabled(True)
            self.filter_label_1.setEnabled(True)
            self.filter_edit_field_1.setEnabled(True)
            self.filter_label_2.setText("Wavelet Selections")
            self.filter_label_1.setText("Wavelet Levels")
    # Button pushed function: Button_Default
    def gatp_defaultset(self):
        self.parameters_default()
        self.parameters_update_field()

    # Button pushed function: SaveSettingsButton
    def gatp_savesetting(self):
        self.parameters_update_setting() #Update settings from fields
        print("setting:", self.settings)

    # Code that executes after component creation

    # Button pushed function: Button_LoadSetting

    #  Component initialization & Create UIFigure and components
    def createComponents(self):
        # Create GATFD_process_UIFigure 
        self.setWindowTitle('GAT-FD - Sliding Window Correlation Analysis v0.2a')
        self.setFixedSize(800, 600)
        self.layout = QFormLayout()
        
        # File Format dropdown
        self.file_format_dropdown = QComboBox()
        self.file_format_dropdown.addItems(["Image (.nii/.nii.gz)", "Matrix (.mat) (Time series)"])
        self.file_format_dropdown.setCurrentIndex(0)
        self.file_format_dropdown.currentIndexChanged.connect(self.gatp_file_type)
        self.layout.addRow("File Format", self.file_format_dropdown )
        # File List Widget
        self.file_list_widget = QListWidget()
        self.layout.addRow("File List", self.file_list_widget )
        # Load file button
        self.open_file_button = QPushButton("Load File(s)")
        self.open_file_button.clicked.connect(self.gatp_open_file)
        self.layout.addRow(self.open_file_button )

        # TRs input
        self.TRs_input = QLineEdit()
        self.TRs_input.setValidator(QIntValidator())
        self.TRs_input.setText("1")
        self.layout.addRow("TR(s)", self.TRs_input )

        # Window Size input
        self.window_size_input = QLineEdit()
        self.window_size_input.setValidator(QIntValidator())
        self.window_size_input.setText("20")
        self.layout.addRow("Window Size", self.window_size_input)

        # Step Size input
        self.step_size_input = QLineEdit()
        self.step_size_input.setValidator(QIntValidator())
        self.step_size_input.setText("1")
        self.layout.addRow("Step Size", self.step_size_input)

        # filter dropdown input
        self.filter_input = QComboBox()
        self.filter_input.addItems(['None', 'Bandpass', 'Highpass', 'Lowpass', 'Wavelet'])
        self.filter_input.setCurrentIndex(0)
        self.filter_input.currentIndexChanged.connect(self.gatp_filter_type)
        self.layout.addRow("Filter", self.filter_input)

        # filter edit field 1
        self.filter_label_1 =QLabel()
        self.filter_edit_field_1 = QLineEdit()
        self.filter_label_1.setEnabled(False)
        self.filter_edit_field_1.setEnabled(False)
        self.filter_label_1.setText("Filter Option")

        self.filter_label_2 =QLabel()
        self.filter_edit_field_2 = QLineEdit()
        self.filter_label_2.setEnabled(False)
        self.filter_edit_field_2.setEnabled(False)
        self.filter_label_2.setText("Filter Option")
        
        self.layout.addRow(self.filter_label_1, self.filter_edit_field_1)
        self.layout.addRow(self.filter_label_2, self.filter_edit_field_2)
        
        # Kernel input
        self.window_kernel_input = QComboBox()
        self.window_kernel_input.addItems(['None', 'Gaussian'])
        self.window_kernel_input.setCurrentIndex(0)
        self.window_kernel_input.currentIndexChanged.connect(self.gatp_kernel)
        self.layout.addRow("Window Kernel", self.window_kernel_input)

        # Kernel edit field
        self.sigma_label =QLabel()
        self.sigma_edit_field = QLineEdit()
        self.sigma_label.setEnabled(False)
        self.sigma_edit_field.setEnabled(False)
        self.sigma_label.setText("\sigma")
        self.layout.addRow(self.sigma_label, self.sigma_edit_field)

        # Atlas dropdown input
        self.atlas_input = QComboBox()
        self.atlas_input.addItems(['Brainnetome 1.25mm', 'AAL 2mm', 'Custom'])
        self.atlas_input.setCurrentIndex(0)
        self.atlas_input.currentIndexChanged.connect(self.gatp_atlas_type)
        self.layout.addRow("Atlas", self.atlas_input)
        # Atlas List Widget
        self.atlas_list_widget = QListWidget()
        self.atlas_list_widget.setEnabled(False)
        self.layout.addRow("Atlas List", self.atlas_list_widget )
        # Add Atlas button
        self.add_atlas_button = QPushButton('Add Mask')
        self.add_atlas_button.clicked.connect(self.gatp_open_atlas)
        self.layout.addRow(self.add_atlas_button)
    
        # Default Setting Button
        self.default_settings_button = QPushButton('Default Setting')
        self.default_settings_button.clicked.connect(self.gatp_defaultset)
        self.layout.addRow(self.default_settings_button)

        # Save button 
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.gatp_savesetting)
        self.layout.addRow(self.save_button)

        # Run button 
        self.run_button = QPushButton('Run')
        self.run_button.clicked.connect(self.gatp_run_process)
        self.layout.addRow(self.run_button)

        self.setLayout(self.layout)

    # Run process
    def gatp_run_process(self):
        # ——— Update settings and parameters ———
        fnc_pro_filelength = len(self.settings['file_list'])
        self.parameters_update_setting()

        # Load filter settings based on user selection
        filter_type = self.settings['filter_type']
        if filter_type == 0:
            filter_setting1 = filter_setting2 = 0
        elif filter_type == 1:
            filter_setting1 = float(self.settings['filter_setting1'])
            filter_setting2 = float(self.settings['filter_setting2'])
            # Validate bandpass: lower < upper
            if filter_setting1 >= filter_setting2:
                QMessageBox.critical(self, "Error", "Please enter correct band pass parameters")
                return
        elif filter_type in (2, 3):
            filter_setting1 = float(self.settings['filter_setting1'])
            filter_setting2 = 0
        elif filter_type == 4:
            filter_setting1 = int(self.settings['filter_setting1'])
            filter_setting2 = int(self.settings['filter_setting2'])
        else:
            # Unexpected filter type
            print("filter type error")
            return

        # Determine atlas and window size depending on file format
        if self.file_format_dropdown.currentIndex() == 0:
            atlas_index = self.atlas_input.currentIndex()
            atlas_files = self.settings['atlas_list']
            if atlas_index == 2 and atlas_files:
                # Custom: merge multiple atlas masks into one integer-labeled mask
                masks = []
                for path in atlas_files:
                    img = nib.load(path)
                    masks.append((img.get_fdata() > 0).astype(np.int32))
                fnc_pro_atlas_masks = np.sum([m * (i+1)
                                            for i, m in enumerate(masks)], axis=0)
                atlas_img = fnc_pro_atlas_masks  # use array directly
            elif atlas_index == 0:
                # Use default Brainetome atlas
                atlas_img = os.path.join(self.settings.get('path', ''), 'atlas', 'BN_atlas_1_25mm.nii')
                fnc_pro_atlas_masks = nib.load(atlas_img).get_fdata().astype(np.int32)
            elif atlas_index == 1:
                # Use default AAL2 atlas
                atlas_img = os.path.join(self.settings.get('path', ''), 'atlas', 'AAL2v1_2mm.nii.txt')
                fnc_pro_atlas_masks = nib.load(atlas_img).get_fdata().astype(np.int32)
            else:
                QMessageBox.critical(self, "Error", "Invalid atlas selection")
                return
            atl_len = int(fnc_pro_atlas_masks.max())
            fnc_window_size = self.settings['window_size']
        else:
            # Matrix input: atlas_list length defines number of regions
            atl_len = len(self.settings['atlas_list'])
            if atl_len < 1:
                QMessageBox.critical(self, 'Error', 'Please load label')
                return
            fnc_window_size = self.settings['window_size']

        # Prompt for output folder
        fnc_out_path = QFileDialog.getExistingDirectory(self, "Select Output Folder")
        if not fnc_out_path:
            return

        # Initialize progress dialog
        progress_dialog = QProgressDialog("Processing...", None, 0, fnc_pro_filelength, self)
        progress_dialog.setWindowTitle("Processing")
        progress_dialog.show()

        for i in range(fnc_pro_filelength):
            file_path = os.path.join(self.settings['file_path'][i], self.settings['file_list'][i])

            if self.file_format_dropdown.currentIndex() == 0:
                # Load NIfTI with memory mapping
                func_img = nib.load(file_path, mmap=True)
                func_data = func_img.get_fdata()
                fnc_raw_len = func_data.shape[3]
                fnc_window_count = fnc_raw_len - fnc_window_size + 1

                # --- Resample atlas to match functional image (MATLAB-style) ---
                # --- Reorient both atlas and func images to RAS+ orientation ---
                if isinstance(atlas_img, str):
                    atlas_nib = nib.load(atlas_img)
                else:
                    atlas_nib = nib.Nifti1Image(atlas_img, func_img.affine)
                atlas_data = atlas_nib.get_fdata()
                func_data = func_img.get_fdata()
                # Print voxel size and orientation info for debug
                print("Atlas voxel size:", atlas_nib.header.get_zooms()[:3])
                print("Atlas affine:\n", atlas_nib.affine)
                print("Atlas orientation:", nib.aff2axcodes(atlas_nib.affine))
                print("Atlas array shape:", atlas_data.shape)
                print("Func voxel size:", func_img.header.get_zooms()[:3])
                print("Func affine:\n", func_img.affine)
                print("Func orientation:", nib.aff2axcodes(func_img.affine))
                print("Func array shape:", func_data.shape)
                # --- Reorient to LAS (Left-Anterior-Superior) to match MATLAB convention ---
                from nibabel.orientations import io_orientation, axcodes2ornt, ornt_transform, apply_orientation
                target_ornt = axcodes2ornt(('L', 'A', 'S'))
                # Atlas reorient
                atlas_ornt = io_orientation(atlas_nib.affine)
                atlas_transform = ornt_transform(atlas_ornt, target_ornt)
                atlas_data_las = apply_orientation(atlas_data, atlas_transform)
                print("Atlas orientation after LAS reorient:", ('L', 'A', 'S'))
                print("Atlas array shape after LAS reorient:", atlas_data_las.shape)
                # Func reorient
                func_ornt = io_orientation(func_img.affine)
                func_transform = ornt_transform(func_ornt, target_ornt)
                func_data_las = apply_orientation(func_data, func_transform)
                print("Func orientation after LAS reorient:", ('L', 'A', 'S'))
                print("Func array shape after LAS reorient:", func_data_las.shape)
                # Use func_data_las and atlas_data_las for downstream processing
                func_data = func_data_las
                atlas_data = atlas_data_las
                func_shape = func_data_las.shape[:3]  # Define func_shape for output_shape
                # --- Affine transform and cropping as before, but use reoriented data ---
                # Step 1: Resample atlas to its own shape (like imwarp in MATLAB)
                atlas_affine = atlas_nib.affine
                func_affine = func_img.affine
                # Compute the transform from atlas to func space
                affine_transform_mat = np.linalg.inv(func_affine) @ atlas_affine
                # Resample atlas to the shape it would have after imwarp (let output_shape be None)
                resampled_atlas = affine_transform(
                    atlas_data_las,
                    matrix=affine_transform_mat[:3, :3],
                    offset=affine_transform_mat[:3, 3],
                    order=0
                )
                fnc_rawdata_tt = resampled_atlas
                resampled_shape = resampled_atlas.shape
                target_shape = func_shape
                print('fnc_rawdata_tt shape:', fnc_rawdata_tt.shape)
                print('atlas (target) shape:', target_shape)
                # Step 2: Crop resampled atlas to match func_shape, using MATLAB-style indices
                def get_crop_indices(src, tgt):
                    if src > tgt:
                        shift = (src - tgt) // 2
                        a_start, a_end = 0, tgt
                        d_start, d_end = shift, shift + tgt
                    else:
                        shift = (tgt - src) // 2
                        a_start, a_end = shift, shift + src
                        d_start, d_end = 0, src
                    return a_start, a_end, d_start, d_end
                axs, axe, dxs, dxe = get_crop_indices(target_shape[0], resampled_shape[0])
                ays, aye, dys, dye = get_crop_indices(target_shape[1], resampled_shape[1])
                azs, aze, dzs, dze = get_crop_indices(target_shape[2], resampled_shape[2])
                print(f'axs: {axs+1} to {axe}')
                print(f'ays: {ays+1} to {aye}')
                print(f'azs: {azs+1} to {aze}')
                print(f'dxs: {dxs+1} to {dxe}')
                print(f'dys: {dys+1} to {dye}')
                print(f'dzs: {dzs+1} to {dze}')
                cropped_atlas = np.zeros(target_shape, dtype=resampled_atlas.dtype)
                cropped_atlas[axs:axe, ays:aye, azs:aze] = resampled_atlas[dxs:dxe, dys:dye, dzs:dze]
                fnc_pro_atlas_masks = cropped_atlas.astype(np.int32)

                # Shift and trim data to match atlas (debug)
                print('Before shift/trim:')
                print('fnc_rawdata_tt shape:', fnc_rawdata_tt.shape)
                print('atlas (target) shape:', func_data[..., 0].shape)
                # In Python, nilearn's resample_to_img already does this alignment, but for debug, print the region indices if you were to crop manually
                # Example indices (not used, for debug only):
                # axs, axe, ays, aye, azs, aze = 0, fnc_rawdata_tt.shape[0], 0, fnc_rawdata_tt.shape[1], 0, fnc_rawdata_tt.shape[2]
                print('No manual crop indices needed in Python (nilearn handles this)')

                # Use ROI order from .txt file (MATLAB style)
                n_x, n_y, n_z, n_t = func_data.shape
                flat_func = func_data.reshape(-1, n_t)  # (n_voxels, n_timepoints)
                flat_atlas = fnc_pro_atlas_masks.flatten()  # (n_voxels,)
                roi_labels = np.arange(1, len(self.settings['atlas_list']) + 1)  # 1-based labels
                atlased_data = np.zeros((n_t, len(roi_labels)), dtype=np.float32)
                for idx, label in enumerate(roi_labels):
                    mask = (flat_atlas == label)
                    if np.any(mask):
                        ts = flat_func[mask, :]
                        atlased_data[:, idx] = ts.mean(axis=0)
                    else:
                        atlased_data[:, idx] = 0
            else:
                # .mat input: load full matrix (time × regions)
                mat = loadmat(file_path)
                key = next(k for k in mat if not k.startswith('__'))
                atlased_data = mat[key].astype(np.float32)
                fnc_raw_len, _ = atlased_data.shape
                fnc_window_count = fnc_raw_len - fnc_window_size + 1
                atl_len = atlased_data.shape[1]

            progress_dialog.setValue(i+1)

            # ——— Apply temporal filter if requested ———
            tr = self.settings['TR']
            fs = 1.0 / tr
            nyq = fs / 2.0
            if filter_type == 1:
                # Bandpass filter
                low = (1.0 / filter_setting2) / nyq
                high = (1.0 / filter_setting1) / nyq
                b, a = butter(2, [low, high], btype='band')
                atlased_data = filtfilt(b, a, atlased_data, axis=0)
            elif filter_type == 2:
                # Highpass filter
                cutoff = (1.0 / filter_setting1) / nyq
                b, a = butter(2, cutoff, btype='high')
                atlased_data = filtfilt(b, a, atlased_data, axis=0)
            elif filter_type == 3:
                # Lowpass filter
                cutoff = (1.0 / filter_setting1) / nyq
                b, a = butter(2, cutoff, btype='low')
                atlased_data = filtfilt(b, a, atlased_data, axis=0)
            elif filter_type == 5:
                # Wavelet denoising per column
                for col in range(atlased_data.shape[1]):
                    coeffs = pywt.wavedec(atlased_data[:, col], 'db1', level=filter_setting1)
                    for idx in range(len(coeffs)):
                        if idx+1 != filter_setting2:
                            coeffs[idx] = np.zeros_like(coeffs[idx])
                    atlased_data[:, col] = pywt.waverec(coeffs, 'db1')[:fnc_raw_len]

            # ——— Sliding-window correlation ———
            corr_data = np.zeros((fnc_window_count, atl_len, atl_len), dtype=np.float32)
            for w in range(fnc_window_count):
                window = atlased_data[w:w+fnc_window_size, :]
                if self.settings['kernel'] == 2:
                    idx = np.arange(1, fnc_window_size+1)
                    center = fnc_window_size/2 + 1
                    gauss = np.exp(-((idx-center)**2)/(2*(self.settings['kernel_setting']**2)))
                    window = window * gauss[:, None]
                tmp = np.corrcoef(window, rowvar=False)
                tmp[np.isnan(tmp)] = 0
                corr_data[w] = tmp

            # Save results to .mat
            subj_data = {
                'd_atlas': atlased_data,
                'd_corr': corr_data,
                'd_setting': self.settings,
                'd_windowsize': self.settings['window_size'],
                'd_stepsize': self.settings['step_size'],
                'd_kernel': self.settings['kernel'],
                'd_atlas_list': self.settings['atlas_list']
            }
            fname = os.path.splitext(self.settings['file_list'][i])[0]
            savemat(os.path.join(fnc_out_path, f"{fname}.mat"), {'subj_data': subj_data})
            progress_dialog.setLabelText(f"{i+1}/{fnc_pro_filelength}: Done")

        # Close dialog and notify completion
        progress_dialog.close()
        QMessageBox.information(self, "Finished", "Finished")

