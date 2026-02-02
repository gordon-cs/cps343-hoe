#!/usr/bin/env python3
# Jonathan Senning <jonathan.senning@gordon.edu>
# Written: 2024-02-21

import tkinter as tk
import numpy as np
from time import sleep

def rgb2hex(rgb):
    return f'#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}'

def hex2rgb(hex):
    r = int(hex[1:3], base=16)
    g = int(hex[3:5], base=16)
    b = int(hex[5:7], base=16)
    return (r, g, b)

class Matrix:
    def __init__(self, nrows, ncols, x, y, ppe, canvas, color):
        self.nrows = nrows      # rows in matrix
        self.ncols = ncols      # columns in matrix
        self.x = x              # x position of upper-left matrix corner
        self.y = y              # y position of upper-left matrix corner
        self.ppe = ppe          # pixels per matrix entry
        self.canvas = canvas    # TK canvas
        self.X = [[canvas.create_rectangle(
            x + i * self.ppe, y + j * self.ppe,
            x + (i+1) * self.ppe, y + (j+1) * self.ppe, fill=color)
                                      for i in range(self.ncols)]
                                      for j in range(self.nrows)]

    def setColor(self, color):
        """Set color of each pixel in matrix entry"""
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.canvas.itemconfig(self.X[i][j], fill = color)

    def setEntry(self, i, j, color):
        """Set matrix entry color"""
        self.canvas.itemconfig(self.X[i][j], fill = color)

    def getEntry(self, i, j):
        """Get matrix entry color"""
        return self.canvas.itemcget(self.X[i][j], 'fill')

class MatrixProductHeatMap():
    M = 3         # Product matrix is MxN
    P = 6         # First factor matrix is MxP
    N = 4         # Second factor matrix is PxN
    PPE = 20      # Pixels per matrix entry (on each row and column) 
    MARGIN = 20      # Tk windo margin
    SPACING = 10     # Horizintal spacing between window elements
    animate = False  # Animation toggle
    mode = 'ijk'     # Default product mode

    width = 2 * MARGIN + 3 * SPACING + (2*N + P) * PPE
    height = 2 * MARGIN + max(M, P) * PPE

    i = j = k = 0

    start_rgb = np.array([240, 240, 240])   # Initial product entry color
    end_rgb = np.array([255, 0, 0])         # Final product entry color
    delta_rgb = np.fix((end_rgb - start_rgb) / P).astype(int) # color step
    
    default_color = '#f0f0f0'               # default matrix factor color
    accessed_color = '#10c010'              # color when entry is accessed
    modified_color = '#d00000'              # color when entry is modified
    start_color = rgb2hex(start_rgb)        # Initial product entry color
    end_color = rgb2hex(end_rgb)            # Final product entry color

    def __init__(self, window):
        """Initialize App window"""
        self.window = window
        self.frame = tk.Frame(self.window)

        # Frames for mode (e.g. ijk) and command (e.g. start/stop) buttons
        self.buttonFrame = tk.Frame(self.window)
        self.cmdFrame = tk.Frame(self.window)
        self.canvas = tk.Canvas(self.frame,
                                width = self.width,
                                height = self.height)
        self.canvas.pack()

        self.modeLabel = tk.Label(self.frame, text = self.mode, fg = 'blue',
                                  padx=2, pady=10)
        self.modeLabel.pack()
        
        # Mode (e.g. ijk) Buttons
        self.ijkButton = tk.Button(self.buttonFrame, text = 'ijk', width = 3,
                                   command = self.set_mode_ijk)
        self.jikButton = tk.Button(self.buttonFrame, text = 'jik', width = 3,
                                   command = self.set_mode_jik)
        self.ikjButton = tk.Button(self.buttonFrame, text = 'ikj', width = 3,
                                   command = self.set_mode_ikj)
        self.kijButton = tk.Button(self.buttonFrame, text = 'kij', width = 3,
                                   command = self.set_mode_kij)
        self.jkiButton = tk.Button(self.buttonFrame, text = 'jki', width = 3,
                                   command = self.set_mode_jki)
        self.kjiButton = tk.Button(self.buttonFrame, text = 'kji', width = 3,
                                   command = self.set_mode_kji)
        self.ijkButton.pack(side=tk.LEFT)
        self.jikButton.pack(side=tk.LEFT)
        self.ikjButton.pack(side=tk.LEFT)
        self.kijButton.pack(side=tk.LEFT)
        self.jkiButton.pack(side=tk.LEFT)
        self.kjiButton.pack(side=tk.LEFT)

        # Command Buttons
        self.startButton = tk.Button(self.cmdFrame, text = 'Start', width = 5,
                                    command = self.start_multiply)
        self.stopButton = tk.Button(self.cmdFrame, text = 'Stop', width = 5,
                                    command = self.stop_multiply)
        self.quitButton = tk.Button(self.cmdFrame, text = 'Quit', width = 15,
                                    command = self.close_windows)
        self.startButton.pack(side=tk.LEFT)
        self.stopButton.pack(side=tk.LEFT)
        self.quitButton.pack(side=tk.LEFT)
        
        self.frame.pack()
        self.buttonFrame.pack()
        self.cmdFrame.pack()

        # Draw matrices in canvas
        self.C = Matrix(self.M, self.N, self.MARGIN, self.MARGIN, self.PPE,
                        self.canvas, self.start_color)
        self.A = Matrix(self.M, self.P,
                        self.MARGIN + self.N * self.PPE + 2*self.SPACING,
                        self.MARGIN, self.PPE, self.canvas, self.default_color)
        self.B = Matrix(self.P, self.N, self.MARGIN
                        + (self.N + self.P) * self.PPE + 3*self.SPACING,
                        self.MARGIN, self.PPE, self.canvas, self.default_color)

    def increment(self, i, j, k):
        """Increment i, j, and k according to the selected mode"""
        if self.mode == 'ijk':
            k = k + 1
            if k >= self.P:
                k = 0
                j = j + 1
            if j >= self.N:
                j = 0
                i = i + 1
            if i >= self.M:
                i = j = k = 0
        elif self.mode == 'jik':
            k = k + 1
            if k >= self.P:
                k = 0
                i = i + 1
            if i >= self.M:
                i = 0
                j = j + 1
            if j >= self.N:
                i = j = k = 0
        elif self.mode == 'ikj':
            j = j + 1
            if j >= self.N:
                j = 0
                k = k + 1
            if k >= self.P:
                k = 0
                i = i + 1
            if i >= self.M:
                i = j = k = 0
        elif self.mode == 'kij':
            j = j + 1
            if j >= self.N:
                j = 0
                i = i + 1
            if i >= self.M:
                i = 0
                k = k + 1
            if k >= self.P:
                i = j = k = 0
        elif self.mode == 'jki':
            i = i + 1
            if i >= self.M:
                i = 0
                k = k + 1
            if k >= self.P:
                k = 0
                j = j + 1
            if j >= self.N:
                i = j = k = 0
        elif self.mode == 'kji':
            i = i + 1
            if i >= self.M:
                i = 0
                j = j + 1
            if j >= self.N:
                j = 0
                k = k + 1
            if k >= self.P:
                i = j = k = 0
        else:
            print('Invalid Mode')
        return (i, j, k)

    def tick(self):
        """Update matrices on canvas by doing one increment of product"""

        # Get local copies of i, j, k
        i, j, k = self.i, self.j, self.k

        self.A.setEntry(i, k, self.default_color)
        self.B.setEntry(k, j, self.default_color)
        
        i, j, k = self.increment(i, j, k)
    
        if (i == 0 and j == 0 and k == 0):
            self.animate = False
            self.C.setColor(self.start_color)
            return

        # Get current color and add delta color to it
        rgb = hex2rgb(self.C.getEntry(i, j)) + self.delta_rgb

        # Set updated color in C and make entries in A and B as accessed
        self.C.setEntry(i, j, rgb2hex(rgb))
        self.A.setEntry(i, k, self.accessed_color)
        self.B.setEntry(k, j, self.accessed_color)

        # save updated values of i, j, k
        self.i, self.j, self.k = i, j, k

        if self.animate and (i != 0 or j != 0 or k != 0):
            self.canvas.after(self.dt, self.tick)

    def reset_matrices(self):
        self.C.setColor(self.start_color)
        self.A.setColor(self.default_color)
        self.B.setColor(self.default_color)
        self.i = self.j = self.k = 0

    def set_mode_ijk(self):
        self.modeLabel.configure(text = 'ijk')
        self.mode = 'ijk'
        self.reset_matrices()
        
    def set_mode_ikj(self):
        self.modeLabel.configure(text = 'ikj')
        self.mode = 'ikj'
        self.reset_matrices()

    def set_mode_jik(self):
        self.modeLabel.configure(text = 'jik')
        self.mode = 'jik'
        self.reset_matrices()

    def set_mode_jki(self):
        self.modeLabel.configure(text = 'jki')
        self.mode = 'jki'
        self.reset_matrices()

    def set_mode_kij(self):
        self.modeLabel.configure(text = 'kij')
        self.mode = 'kij'
        self.reset_matrices()

    def set_mode_kji(self):
        self.modeLabel.configure(text = 'kji')
        self.mode = 'kji'
        self.reset_matrices()

    def start_multiply(self):
        if self.animate:
            return
        self.animate = True
        self.stopButton.configure(text = 'Stop')
        i, j, k = self.i, self.j, self.k
        if i*j*k == 0:
            rgb = hex2rgb(self.C.getEntry(i, j)) + self.delta_rgb
            self.canvas.itemconfig(self.C.X[i][j], fill = rgb2hex(rgb))
            self.canvas.itemconfig(self.A.X[i][k], fill = self.accessed_color)
            self.canvas.itemconfig(self.B.X[k][j], fill = self.accessed_color)
        self.dt = 100
        sleep(self.dt * 0.001)
        self.tick()

    def stop_multiply(self):
        if self.animate:
            self.animate = False
            self.stopButton.configure(text = 'Reset')
        else:
            self.reset_matrices()

    def close_windows(self):
        self.window.destroy()

if __name__ == '__main__':
    root = tk.Tk()
    app = MatrixProductHeatMap(root)
    root.mainloop()
