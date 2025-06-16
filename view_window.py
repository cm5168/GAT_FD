# view_window.py
from PyQt6.QtWidgets import QWidget, QLabel, QPushButton, QVBoxLayout

class ViewWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Display')
        self.setFixedSize(800, 600)

        layout = QVBoxLayout()

        label = QLabel('Display UI goes here')
        layout.addWidget(label)

        # Example button
        button = QPushButton('Load Processed File')
        layout.addWidget(button)

        self.setLayout(layout)
