import numpy as np

class Mesh:
    '''
    Contains informations about the discretization
    Elements: number of elements
    Nodes: number of nodes
    NodesElement: number of elements for each element
    dofsNode: degree of freedom per node
    d: dimension
    dofs: degree of freedom
    elements: array with caracteristcs of elements
    coordinates: coordinates of nodes

    '''
    def __init__(self):
        self.Elements = None
        self.Nodes = None
        self.NodesElement = None
        self.dofsNode = None
        self.d = None
        self.dofs = None
        self.elements = None
        self.coordinates = None
        self.points = None

    def nodesInElement(self, elementNumber):
        '''
        Evaluation of the number of nodes in each element
        elementNumber: current number of the element
        elements: contain the elements of the structure
        '''
        return self.elements[elementNumber][2:]

    def elementType(self, elementNumber):
        '''
        Definition of the element type
        elementNumber: current number of the element
        elements: contain the elements of the structure
        '''
        if self.elements[elementNumber][0] == 1:
            return "spring"

        elif self.elements[elementNumber][0] == 2:
            return "bar"

        elif self.elements[elementNumber,0] == 3:
            return "2D_bar"

        elif self.elements[elementNumber,0] == 4:
            return "triangle"

        else:
            print("elementType(): I don't know elementType =", self.elements[elementNumber][0])
            sys.exit()

    def element_coordinates(self, elementNumber):
        """Returns an array with the element nodal coordinates.
        elementNumber: current number of element
        nodesInElement: return just the nodes of each elements
        coordinates: an array with as rows the componets
        """
        return self.coordinates[self.nodesInElement(elementNumber)]