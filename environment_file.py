class EnvironmentFile:
    """
    A class to hold data related to a single environment file
    """

    def __init__(self, grid, resources, size, name):
        self.grid = grid
        self.resources = resources
        self.size = size
        self.name = name.split("/")[-1] #Extract filename from path
        self.name = self.name[:-4] #Remove file extension

    def __getitem__(self, index):
        return self.grid[index]

    def __len__(self):
        return len(self.grid)
    
