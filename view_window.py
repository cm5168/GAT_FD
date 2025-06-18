# view_window.py
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QLabel, QPushButton, QSlider, QVBoxLayout,
    QHBoxLayout, QGridLayout, QFileDialog, QMessageBox, QComboBox
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QKeyEvent
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import scipy.io

class ViewWindow(QMainWindow    ):
    def __init__(self):
        super().__init__()
        self.createComponents()

    def gatv_load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(  self, "Select Conditions", "","MAT files (*.mat)")
        if not file_path:
            return
        # Load .mat file with options to mimic MATLAB structure
        mat_data = scipy.io.loadmat(file_path, squeeze_me=True, struct_as_record=False)
        subj_data = mat_data['subj_data']
        if isinstance(subj_data, np.ndarray):
            subj_data = subj_data.item()  # unwrap 1x1 array to struct
        d_corr = subj_data.d_corr
        self.gatv_setting['networks'] = d_corr
        self.gatv_setting['frame_limit'] = [1, len(d_corr)]
        self.frame_slider.setMinimum(0)
        self.frame_slider.setMaximum(len(d_corr))
        # print("Loaded file:", file_path.split("/")[-1])
        # print("Total frames:", d_corr)
        self.gatv_update_network()
        
    def gatv_update_network(self):
        if 'networks' in self.gatv_setting:
            self.gatv_setting['frame_cur'] = int(self.frame_slider.value())
            self.frame_slider_cursor.setText(str(self.gatv_setting['frame_cur']))
            self.network_ax.clear()

            data = self.gatv_setting['networks']
            if data.ndim == 3:
                matrix = data[self.gatv_setting['frame_cur'] - 1, :, :]
            else:
                matrix = data  # just one matrix, no frame dimension

            self.network_ax.imshow(matrix, cmap='jet')
            self.network_canvas.draw()


    def createComponents(self):
        # Create GATFD_process_UIFigure 
        self.setWindowTitle("GAT-FD - Result Display v0.2a")
        self.setGeometry(100, 100, 832, 810)

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        layout = QGridLayout(self.central_widget)

        self.gatv_setting = {}
        self.gatv_data = {}
        self.gatv_plot = {}

        # UIAxes_Network
        self.network_canvas = FigureCanvas(Figure())
        self.network_ax = self.network_canvas.figure.add_subplot(111)
        self.network_ax.set_title("Network")
        layout.addWidget(self.network_canvas, 0, 1)

        # FrameSlider and Label
        self.frame_slider_label = QLabel("Frame")
        layout.addWidget(self.frame_slider_label, 1, 0)

        self.frame_slider = QSlider(Qt.Orientation.Horizontal)
        self.frame_slider.setMinimum(0)
        self.frame_slider.setMaximum(100)
        # self.frame_slider.valueChanged.connect(self.gatv_update_network)
        layout.addWidget(self.frame_slider, 1, 1)

        self.frame_slider_cursor = QLabel("0")
        layout.addWidget(self.frame_slider_cursor, 1, 2)

        # LoadProcessedFileButton
        self.load_processed_button = QPushButton("Load Processed File")
        self.load_processed_button.clicked.connect(self.gatv_load_file)
        layout.addWidget(self.load_processed_button, 2, 0)

        # UIAxes_TaskDesign
        self.task_design_canvas = FigureCanvas(Figure())
        self.task_design_ax = self.task_design_canvas.figure.add_subplot(111)
        layout.addWidget(self.task_design_canvas, 0, 0)

        # UIAxes_Timeseries
        self.timeseries_canvas = FigureCanvas(Figure())
        self.timeseries_ax = self.timeseries_canvas.figure.add_subplot(111)
        layout.addWidget(self.timeseries_canvas, 3, 0, 1, 3)

        # LoadTemporalMaskButton
        self.load_temporal_button = QPushButton("Load Temporal Mask")
        # self.load_temporal_button.clicked.connect(self.gatv_load_design)
        layout.addWidget(self.load_temporal_button, 2, 1)

        # LoadNetworkPropertyFileButton
        self.load_property_button = QPushButton("Load Network Property File")
        # self.load_property_button.clicked.connect(self.gatv_load_measure)
        layout.addWidget(self.load_property_button, 4, 0)

        # Dropdown Labels and Widgets
        self.threshold_label = QLabel("Threshold")
        layout.addWidget(self.threshold_label, 5, 0)

        self.threshold_dropdown = QComboBox()
        self.threshold_dropdown.addItem("Mean")
        layout.addWidget(self.threshold_dropdown, 5, 1)

        self.measure_label = QLabel("Measure")
        layout.addWidget(self.measure_label, 6, 0)

        self.measure_dropdown = QComboBox()
        layout.addWidget(self.measure_dropdown, 6, 1)

        self.subject_label = QLabel("Subject")
        layout.addWidget(self.subject_label, 7, 0)

        self.subject_dropdown = QComboBox()
        self.subject_dropdown.addItem("GroupAverage")
        layout.addWidget(self.subject_dropdown, 7, 1)

        self.update_button = QPushButton("Update")
        # self.update_button.clicked.connect(self.gatv_update_timeseries)
        layout.addWidget(self.update_button, 8, 0)

        self.show()