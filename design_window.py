# design_window.py
from PyQt6.QtWidgets import QWidget, QLabel, QPushButton, QVBoxLayout

class DesignWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Task Design')
        self.setFixedSize(800, 600)

        layout = QVBoxLayout()

        label = QLabel('Task Design UI goes here')
        layout.addWidget(label)

        # Example button
        button = QPushButton('Save Design Matrix')
        layout.addWidget(button)

        self.setLayout(layout)
    