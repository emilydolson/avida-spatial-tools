import seaborn as sns

class EnvironmentFile:
    """
    A class to hold data related to a single environment file
    """

    def __init__(self, grid, resources, size, name, tasks):
        self.grid = grid
        self.resources = resources
        self.size = size
        self.name = name.split("/")[-1] #Extract filename from path
        self.name = self.name[:-4] #Remove file extension
        self.tasks = tasks
        
        if len(self.resources) == 1:
            #Yay, we can make the background nice and simple!
            self.resource_palette = sns.dark_palette("black", 1)
            #If we're only using two colors, it's better for them to be red
            #and yellow than red and green
            self.task_palette = sns.hls_palette(max(len(tasks), 4), s=1)
        
        elif set(self.tasks) == set(self.resources):
            #Yay, the tasks and resources correspond to each other!
            #If we're only using two colors, it's better for them to be red
            #and yellow than red and green
            self.resource_palette = sns.hls_palette(max(len(tasks), 4), s=1)
            self.task_palette = sns.hls_palette(max(len(tasks), 4), s=1)
        
        elif len(self.tasks) == 1:
            #Yay, we can display the tasks easily!
            self.task_palette = sns.dark_palette("black", 1)
            #If we're only using two colors, it's better for them to be red
            #and yellow than red and green
            self.resource_palette = sns.hls_palette(max(len(tasks), 4), s=1)

        else:
            #Since we don't know anything in particular about the tasks
            #and resources, let's be color-blind friendly and keep the
            #background simple. The user probably wants to assign new
            #palettes based on their data.
            self.task_palette = sns.color_palette("colorblind", len(tasks))
            self.resource_palette = sns.color_palette("bone", len(resources)) 

    def __getitem__(self, index):
        return self.grid[index]

    def __len__(self):
        return len(self.grid)
    
