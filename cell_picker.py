#The cell_picker class provides a graphical interface for selecting cells in
#an environment.

#Based on code from http://stackoverflow.com/questions/26988204/using-2d-array-to-create-clickable-tkinter-canvas

import Tkinter as tk

class cell_picker():
    """
    The cell_picker class produces a graphical window containing a grid of
    squares. Clicking a cell will turn it black and add it to the list of
    selected cells. Clicking again will turn it white and remove it. When
    the window is closed, the full list of currently selected cells will be
    returned.

    The cell_picker class also supports selecting and de-selecting cells by
    clicking and dragging across the grid. After every click, the number of
    cells currently selected is printed to the terminal.

    By default, the grid is 60x60 (the default in Avida), but the number of
    rows and columns can be specified as arguments to __init__. A third
    argument, a list of currently selected cells, can also be passed.

    Generally, you want to run the cell_picker as soon as it is created:
    >>> selected_cells = cell_picker().run()
    """
    def __init__(self, rows=60, cols=60, selected=[]):
        """
        Create the cell_picker - arguments explained in class docstring.
        """
        self.begin_drag = None
        self.rows = rows
        self.cols = cols

        # Create a grid of None to store the references to the tiles
        self.tiles = [[None for _ in range(self.cols)] \
                      for _ in range(self.rows)]
        self.cells = []

        self.root = tk.Tk()
        self.c = tk.Canvas(self.root, width=self.rows*10, \
                        height=self.cols*10, borderwidth=0, background='white')
        self.c.pack()

        # Create the window, a canvas and the mouse click event binding
        self.c.bind("<Button-1>", self.callback)
        self.c.bind("<ButtonRelease-1>", self.dragend)
        self.initial = selected


    def callback(self, event):
        """
        Selects cells on click.
        """
        self.init_width()
        
        if len(self.initial) > 0:
            for cell in self.initial:
                self.color_square(cell[0], cell[1], True)
            self.initial = []
        self.begin_drag = event
        self.color_square(event.x, event.y)

    def init_width(self):
        """
        Get rectangle diameters
        """
        self.col_width = self.c.winfo_width()/self.cols
        self.row_height = self.c.winfo_height()/self.rows

    def color_square(self, x, y, unit_coords=False):
        """
        Handles actually coloring the squares
        """
        # Calculate column and row number
        if unit_coords:
            col = x
            row = y
        else:
            col = x//self.col_width
            row = y//self.row_height

        # If the tile is not filled, create a rectangle
        if not self.tiles[row][col]:
            self.tiles[row][col] = \
            self.c.create_rectangle(col*self.col_width, row*self.row_height,\
                (col+1)*self.col_width, (row+1)*self.row_height, fill="black")
            self.cells.append(row*self.cols + col)
        
        # If the tile is filled, delete the rectangle and clear the reference
        else:
            self.c.delete(self.tiles[row][col])
            self.tiles[row][col] = None
            self.cells.remove(row*self.cols + col)

    def dragstart(self, event):
        """
        Handles the beggining of a drag action.
        """
        begin_drag = event

    def dragend(self, event):
        """
        Handles the end of a drag action.
        """
        x_range = [self.begin_drag.x//self.col_width, event.x//self.col_width]
        y_range = [self.begin_drag.y//self.row_height, event.y//self.row_height]

        #Check bounds
        for i in range(2):
            for ls in [x_range, y_range]:
                if ls[i] < 0:
                    ls[i] = 0
                if ls[i] >= self.rows:
                    ls[i] = self.rows-1
            
        for x in range(min(x_range), max(x_range)+1):
            for y in range(min(y_range), max(y_range)+1):
                if x == self.begin_drag.x//self.col_width and \
                   y == self.begin_drag.y//self.row_height:
                    continue
                self.color_square(x*self.col_width, y*self.row_height)
        self.begin_drag = None

        print len(self.cells), "cells selected"

    def run(self):
        """
        Makes the cell_picker actually do stuff and returns the result.
        """
        self.root.mainloop()
        return self.cells

if __name__ == "__main__":
    print cell_picker(60, 60).run()
