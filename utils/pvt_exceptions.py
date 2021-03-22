class DataExcist(Exception):
    def __init__(self, path, message="File excist"):
        self.path = path
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "{}: {}".format(self.message, self.path)

class StationsTooClose(Exception):
    def __init__(self,n1,s1,n2,s2,message="Station too close:"):
        self.n1 = n1
        self.s1 = s1
        self.n2 = n2
        self.s2 = s2
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "{}: {}.{}-{}.{}".format(
            self.message,
            self.n1,
            self.s1,
            self.n2,
            self.s2
        )
    
class StationsTooFar(Exception):
    def __init__(self,n1,s1,n2,s2,message="Station too far:"):
        self.n1 = n1
        self.s1 = s1
        self.n2 = n2
        self.s2 = s2
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "{}: {}.{}-{}.{}".format(
            self.message,
            self.n1,
            self.s1,
            self.n2,
            self.s2
        )