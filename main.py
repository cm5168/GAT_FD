import sys
from PyQt6.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout

from process_window import ProcessWindow
from network_window import NetworkWindow
from view_window import ViewWindow
from design_window import DesignWindow

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('GAT_FD')
        self.setFixedSize(240, 360)
        self.init_ui()
    def init_ui(self):
        layout = QVBoxLayout()

        btn_swa = QPushButton('Sliding-Window Analysis')
        btn_swa.clicked.connect(self.open_process_window)
        layout.addWidget(btn_swa)

        btn_network = QPushButton('Network Analysis')
        btn_network.clicked.connect(self.open_network_window)
        layout.addWidget(btn_network)

        btn_display = QPushButton('Display')
        btn_display.clicked.connect(self.open_view_window)
        layout.addWidget(btn_display)

        btn_design = QPushButton('Task Design')
        btn_design.clicked.connect(self.open_design_window)
        layout.addWidget(btn_design)
        
        self.setLayout(layout)
    def open_process_window(self):
        self.proc_win = ProcessWindow()
        self.proc_win.show()
    
    def open_network_window(self):
        self.net_win = NetworkWindow()
        self.net_win.show()
    
    def open_view_window(self):
        self.view_win = ViewWindow()
        self.view_win.show()
    
    def open_design_window(self):
        self.design_win = DesignWindow()
        self.design_win.show()
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_win = MainWindow()
    main_win.show()
    sys.exit(app.exec())