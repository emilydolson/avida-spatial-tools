import Tkinter as tk
import patch_analysis

class cell_picker():
    def __init__(self, rows=60, cols=60, selected=[]):
        self.begin_drag = None
        self.rows = rows
        self.cols = cols

        # Create a grid of None to store the references to the tiles
        self.tiles = [[None for _ in range(self.cols)] for _ in range(self.rows)]
        self.cells = []

        self.root = tk.Tk()
        self.c = tk.Canvas(self.root, width=self.rows*10, height=self.cols*10, borderwidth=0, background='white')
        self.c.pack()

        # Create the window, a canvas and the mouse click event binding
        self.c.bind("<Button-1>", self.callback)
        #self.c.bind("<B1-Motion>", self.dragstart)
        self.c.bind("<ButtonRelease-1>", self.dragend)
        self.initial = selected


    def callback(self, event):
        self.init_width()
        
        if len(self.initial) > 0:
            for cell in self.initial:
                self.color_square(cell[0], cell[1], True)
            self.initial = []
        self.begin_drag = event
        self.color_square(event.x, event.y)

    def init_width(self):
        # Get rectangle diameters
        self.col_width = self.c.winfo_width()/self.cols
        self.row_height = self.c.winfo_height()/self.rows

    def color_square(self, x, y, unit_coords=False):
        # Calculate column and row number

        if unit_coords:
            col = x
            row = y
        else:
            col = x//self.col_width
            row = y//self.row_height

        # If the tile is not filled, create a rectangle
        if not self.tiles[row][col]:
            self.tiles[row][col] = self.c.create_rectangle(col*self.col_width, row*self.row_height, (col+1)*self.col_width, (row+1)*self.row_height, fill="black")
            self.cells.append(row*self.cols + col)
        # If the tile is filled, delete the rectangle and clear the reference
        else:
            self.c.delete(self.tiles[row][col])
            self.tiles[row][col] = None
            self.cells.remove(row*self.cols + col)

    def dragstart(self, event):
        begin_drag = event

    def dragend(self, event):
        #self.init_width()
        x_range = (self.begin_drag.x//self.col_width, event.x//self.col_width)
        y_range = (self.begin_drag.y//self.row_height, event.y//self.row_height)
        for x in range(min(x_range), max(x_range)+1):
            for y in range(min(y_range), max(y_range)+1):
                if x == self.begin_drag.x//self.col_width and y == self.begin_drag.y//self.row_height:
                    continue
                self.color_square(x*self.col_width, y*self.row_height)
        self.begin_drag = None

        print len(self.cells)
        print patch_analysis.shape_index([[c%60, c/60] for c in list(set(self.cells))])

    def run(self):
        self.root.mainloop()
        return self.cells

#print cell_picker().run()
