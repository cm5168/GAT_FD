# network_window.py
from PyQt6.QtWidgets import QWidget, QLabel, QPushButton, QVBoxLayout

class NetworkWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Network Analysis')
        self.setFixedSize(800, 600)

        layout = QVBoxLayout()

        label = QLabel('Network Analysis UI goes here')
        layout.addWidget(label)

        # Example button
        button = QPushButton('Calculate Network Properties')
        layout.addWidget(button)

        self.setLayout(layout)
