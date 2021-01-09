class DataExcist(Exception):

    def __init__(self, path, message="File excist"):
        self.path = path
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "{}: {}".format(self.message, self.path)