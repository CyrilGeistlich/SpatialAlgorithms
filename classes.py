## LIBRARIES ##

from numpy import sqrt, radians, arcsin, sin, cos
import json

## DEFINE CLASSES HERE ##

## POINT CLASS ##

class Point():
    _id_counter = 0  # Class-level attribute to track the ID

    # initialise
    def __init__(self, x=None, y=None,name = None):
        self.name = name
        self.x = x
        self.y = y
        # SET ID
        type(self)._id_counter += 1
        self.id = self._id_counter
    
    # representation
    def __repr__(self):
        return f'Point(x={self.x}, y={self.y})'

        # Test for equality between Points
    def __eq__(self, other): 
        if not isinstance(other, Point):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.x == other.x and self.y == other.y
    # We need this method so that the class will behave sensibly in sets and dictionaries
    def __hash__(self):
        return hash((self.x, self.y))
    
    # calculate Euclidean distance between two points
    def distEuclidean(self, other):
        return sqrt((self.x-other.x)**2 + (self.y-other.y)**2)
    
    # Calculate determinant with respect to three points. Note the order matters here - we use it to work out left/ right in the next method
    def __det(self, p1, p2):
        det = (self.x-p1.x)*(p2.y-p1.y)-(p2.x-p1.x)*(self.y-p1.y)       
        return det
    
    def leftRight(self, p1, p2):
    # based on GIS Algorithms, Ch2, p11-12, by Ningchuan Xiao, publ. 2016
    # -ve: this point is on the left side of a line connecting p1 and p2
    #   0: this point is collinear
    # +ve: this point is on the right side of the line
        side = int(self.__det(p1, p2))
        if side != 0:
            side = side/abs(side)  # will return 0 if collinear, -1 for left, 1 for right
        return side
    
    def set_polygon_id(self, polygon):
        self.polygon_id = polygon.id


## POINTGROUPS CLASS ##

# data provided should be as an array of points with x, y coordinates.
# class for a group of Points, assumes initial data is unsorted, spatially

class PointGroup(): 
    # initialise
    def __init__(self, data=None, xcol=None, ycol=None):
        self.points = []
        self.size = len(data)
        for d in data:
            self.points.append(Point(d[xcol], d[ycol]))
    
    # representation
    def __repr__(self):
        return f'PointGroup containing {self.size} points' 
 
    # create index of points in group for referencing
    def __getitem__(self, key):
        return self.points[key]
    
    # Calculate determinant with respect to three points. Note the order matters here - we use it to work out left/ right in the next method
    def __det(self, p1, p2):
        det = (self.x-p1.x)*(p2.y-p1.y)-(p2.x-p1.x)*(self.y-p1.y)       
        return det

    def leftRight(self, p1, p2):
        # based on GIS Algorithms, Ch2, p11-12, by Ningchuan Xiao, publ. 2016
        # -ve: this point is on the left side of a line connecting p1 and p2
        #   0: this point is collinear
        # +ve: this point is on the right side of the line
        side = int(self.__det(p1, p2))
        if side != 0:
            side = side/abs(side)  # will return 0 if collinear, -1 for left, 1 for right
        return side
      
## SEGMENT (LINE) CLASS ## 

class Segment():    
    # initialise
    def __init__(self, p0, p1):
        # A segement is made up of two points
        self.start = p0
        self.end = p1
        self.length = 0
    
    # representation
    def __repr__(self):
        return f'Segment with start {self.start} and end {self.end}.' 

    # Test for equality between Segments - we treat segments going in opposite directions as equal here (so direction doesn't matter)
    def __eq__(self, other): 
        if (self.start == other.start or self.start == other.end) and (self.end == other.end or self.end == other.start):
            return True
        else:
            return False
            # We need this method so that the class will behave sensibly in sets and dictionaries
    
    def __hash__(self):
        return hash((self.start.x, self.start.y,self.end.x, self.end.y))    
        
    # determine if intersects with another segment
    # - should we incorporate testing for identical segments and non-zero lengths?

    def intersects(self, other):        
        
        # first test if projections overlap 
        
        mnx = min(other.start.x, other.end.x)
        mny = min(other.start.y, other.end.y)
        mxx = max(other.start.x, other.end.x)
        mxy = max(other.start.y, other.end.y)
        
        if max(self.start.x,self.end.x) < mnx:
            return False
        if max(self.start.y,self.end.y) < mny: 
            return False
        if min(self.start.x,self.end.x) > mxx:
            return False
        if min(self.start.y,self.end.y) > mxy: 
            return False
        
        # If we get here, projections do overlap

        # find side for each point of each line in relation to the points of the other line
        # Positive det means both left or both right -> no intersection
        abp = self.start.leftRight(self.end, other.start)
        abq = self.start.leftRight(self.end, other.end)
        if abp * abq > 0:
            return False

        # We only calculate these determinants if we need to
        pqa = other.start.leftRight(other.end, self.start)
        pqb = other.start.leftRight(other.end, self.end)           
        if pqa * pqb > 0:
            return False
        
        return True

    def intersection(self, other):
        """Calculate intersection point between two line segments."""
        a = self.start
        b = self.end
        c = other.start
        d = other.end
        div = (a.x - b.x) * (c.y - d.y) - (a.y - b.y) * (c.x - d.x)
        if div == 0:
            print("you fucked up in function -> intersection")
            print(self)
            print(other)
            return f"doesn't work: {a}, {b},{c},{d} \n"  # Lines are parallel
        else:
            x = ((a.x * b.y - a.y * b.x) * (c.x - d.x) - (a.x - b.x) * (c.x * d.y - c.y * d.x)) / div
            y = ((a.x * b.y - a.y * b.x) * (c.y - d.y) - (a.y - b.y) * (c.x * d.y - c.y * d.x)) / div
            return Point(x, y)

## POLYGON CLASS ## 

class Vertex(Point):
    def __init__(self, x,y, name = None, intersect = False, alpha = 0.0):
        super().__init__(x,y,name)
        self.intersect = intersect
        self.entry_exit = None
        self.alpha = alpha
        self.next = None
        self.prev = None
        self.link = None ## Pointer to the Vertex in the OTHER Polygon, but same Intersection
        self.processed = False

    def __repr__(self):
        return f'Vertex(x={self.x}, y={self.y}, intersect = {self.intersect}, processed = {self.processed}, {self.entry_exit})'

    def next_vertex(self, next_original = True, unprocessed = None):
        """Return the next non intersecting vertex after the one specified."""
        start = self
        if next_original:
            c = self.next
            while c.intersect:
                c = c.next
                if c == start:
                    return c
            return c
        elif unprocessed: ## find next intersection
            c = self.next
            while not (c.intersect and not c.processed):
                c = c.next
                if c == start:
                    return c
            return c
        else:
            c = self.next
            while not c.intersect:
                c = c.next
                if c == start:
                    return c
            return c

    def mark_processed(self):
        self.processed = True
        if self.link and not self.link.processed:
            self.link.mark_processed()

    def perturb(self):
        self.x += 0.1
        self.y += 0.1

class Polygon(PointGroup):  
    _id_counter = 0  # Class-level attribute to track the ID

    # initialise
    def __init__(self, data=None, xcol=None, ycol=None,name=None):
        self.name = name
        self.points = []
        self.first = None
        # SET ID COUNTER
        type(self)._id_counter += 1
        self.id = self._id_counter

        if data:
            for i, d in enumerate(data):
                self.points.append(Point(name, d[xcol], d[ycol]))
                self.add(Vertex(d[xcol], d[ycol]))
            self.removeDuplicates()
            self.bbox = Bbox(self)


    # representation
    def __repr__(self):
        return f'Polygon PointGroup containing {self.size} points' 

    @property
    def size(self):
        """ Achtung, .first wird nur einmal gezählt! """
        return sum(1 for _ in self)

    def __iter__(self):
        """ Make the polygon iterable, traversing from head to the end. """
        current = self.first
        while current:
            yield current
            current = current.next
            if current == self.first: break

    def __getitem__(self, index):
        """ Retrieve the element at the specified index. """
        if index < 0 or index >= self.size:
            raise IndexError("Index out of range")
        current = self.first
        for _ in range(index):
            current = current.next
        return current

    def add(self, point):
        if not self.first:
            self.first = point
            self.first.next = point
            self.first.prev = point
        else:
            next = self.first
            prev = next.prev
            next.prev = point
            point.next = next
            point.prev = prev
            prev.next = point
    
    def replace(self, old, new):
        new.next = old.next
        new.prev = old.prev
        old.prev.next = new
        old.next.prev = new
        if old is self.first:
            self.first = new

    def insert(self, vertex, edge):
        """Insert and sort a vertex between a specified pair of vertices.

        This function inserts a vertex (most likely an intersection point)
        between two other vertices (start and end). These other vertices
        cannot be intersections (that is, they must be actual vertices of
        the original polygon). If there are multiple intersection points
        between the two vertices, then the new vertex is inserted based on
        its alpha value.
        """
        curr = edge.start
        while curr != edge.end and curr.alpha < vertex.alpha:
            curr = curr.next

        vertex.next = curr
        prev = curr.prev
        vertex.prev = prev
        prev.next = vertex
        curr.prev = vertex


    def remove(self, vertex):
        next = vertex.next
        previous = vertex.prev

        vertex.prev.next = next
        vertex.next.previous = previous

    def pop_vertices(self): ## Only used to display
        edge_ends = [i.end for i in self.edges()]
        edge_ends.insert(0, self.first)
        return edge_ends

    # test if polygon is closed: first and last point should be identical
    def isClosed(self):
        start = self.points[0]
        end = self.points[-1]
        return start == end

    def removeDuplicates(self):
        oldn = self.size
        points = list(dict.fromkeys(self)) # Get rid of the duplicates
        self.first = None
        for p in points:
            self.add(Vertex(p.x,p.y,p.name,p.intersect,p.alpha))
        n = self.size
        
        # find area and centre of the polygon
    # - based on GIS Algorithms, Ch.2 p9-10, by Ningchuan Xiao, publ. 2016 
    
    def __signedArea(self):  # used for both area and centre calculations - this is a private method (only used within the class)
        a = 0
        xmean = 0
        ymean = 0
        for i in range(0, self.size-1):
            ai = self[i].x * self[i+1].y - self[i+1].x * self[i].y
            a += ai
            xmean += (self[i+1].x + self[i].x) * ai
            ymean += (self[i+1].y + self[i].y) * ai

        a = a/2.0   # signed area of polygon (can be a negative)
    
        return a, xmean, ymean #Notice this method returns three values, which we use in area and centre
    
    def area(self):
        a, xmean, ymean = self.__signedArea()
        area = abs(a)   # absolute area of polygon
        
        return area

    def centre(self):
        a, xmean, ymean = self.__signedArea() # note we use the signed area here
        centre = Point(xmean/(6*a), ymean/(6*a)) # centre of polygon 
        return centre

    def containsPoint(self, p):
        if (self.bbox.containsPoint(p) == False):
            return False
        
        # Solution, as discussed in lecture, added here
        ray = Segment(p, Point(self.bbox.ur.x+1, p.y))
        count = 0
        
        for i in range(0, self.size-1):
            start = self[i]
            end = self[i+1]
            s = Segment(start, end)
            if s.intersects(ray):
                if (p.y != min(start.y, end.y)):
                    count = count + 1

        #print(f'count: {count}')
        #print(p)
        if (count%2 == 0):
            return False           

        return True

    def perturb(self, redo = False):
        """Perturb the coordinates of the polygon by adding 0.001 to both X and Y coordinates.""" ## Needed because Grainer Horman algorithm doesn't work when two points are the same
        if redo:
            for point in self:
                point.x -= 0.1
                #point.y -= 0.1
            self.perturbed = False
        else:
            for point in self:
                point.x += 0.1
                #point.y += 0.01
            self.perturbed = True

    
    def clip(self, clip_polygon):
        """Clip this polygon with another polygon using Greiner-Hormann algorithm."""

        if not self.bbox.intersects(clip_polygon.bbox):
            print("Bbox test failed, no polygon intersection")
            return 0

        #self.perturb()

        for i, s in enumerate(self):
            if not s.intersect:
                for j, c in enumerate(clip_polygon):
                    if not c.intersect:
                        if Segment(s,s.next_vertex()).intersects(Segment(c,c.next_vertex())):

                            if s.leftRight(c,c.next_vertex()) == 0: ## If Point in Polygon A is colinear with Edge in B
                                s.perturb()
                            if c.leftRight(s,s.next_vertex()) == 0:
                                c.perturb()

                            intersection_point = Segment(s,s.next_vertex()).intersection(Segment(c,c.next_vertex()))

                            c_vertex = Vertex(intersection_point.x, intersection_point.y, intersect = True,
                                alpha = intersection_point.distEuclidean(c))
                            
                            s_vertex = Vertex(intersection_point.x, intersection_point.y, intersect = True,
                                alpha = intersection_point.distEuclidean(s))

                            c_vertex.link = s_vertex
                            s_vertex.link = c_vertex

                            clip_polygon.insert(c_vertex, Segment(c,c.next_vertex()))
                            self.insert(s_vertex, Segment(s,s.next_vertex()))
        

        
        ## Phase 2:
        self.update_entry_exit(clip_polygon)
        clip_polygon.update_entry_exit(self)


    def unclip(self):
        current = self.first
        initial_first = self.first

        while True:
            next_node = current.next
            if current.intersect:
                current.prev.next = current.next
                current.next.prev = current.prev
                if current == self.first:
                    self.first = current.next

            if next_node == initial_first:
                break
            current = next_node

    def intersection(self, other):

        self.clip(other)
        current = self.first.next_vertex(next_original=False)
        result = []

        while current.intersect and not current.processed:
            #current = self.first.next_vertex(next_original=False, unprocessed=True)
            current.mark_processed()
            clipped = Polygon()
            clipped.add(Vertex(current.x,current.y))

            while True:
                current.mark_processed()
                if current.entry_exit == "entry":
                    while True:
                        current = current.next
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break
                #if current.entry_exit == "entry":
                else:
                    while True:
                        current = current.prev
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break

                current = current.link

                if current.processed:
                    result.append(clipped)
                    current = self.first.next_vertex(next_original=False, unprocessed=True)
                    break
        
        self.unclip()
       # self.perturb(redo=True)
        other.unclip()
        return result


    def union(self, other):
        self.clip(other)
        current = self.first.next_vertex(next_original=False)
        result = []

        while current.intersect and not current.processed:
            #current = self.first.next_vertex(next_original=False, unprocessed=True)
            current.mark_processed()
            clipped = Polygon()
            clipped.add(Vertex(current.x,current.y))

            while True:
                current.mark_processed()
                if current.entry_exit == "entry":
                    while True:
                        current = current.prev
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break
                #if current.entry_exit == "entry":
                else:
                    while True:
                        current = current.next
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break

                current = current.link

                if current.processed:
                    result.append(clipped)
                    current = self.first.next_vertex(next_original=False, unprocessed=True)
                    break
        
        self.unclip()
       # self.perturb(redo=True)
        other.unclip()
        return result


    def difference(self, other):
        self.clip(other)
        current = self.first.next_vertex(next_original=False)
        result = []

        while current.intersect and not current.processed:
            #current = self.first.next_vertex(next_original=False, unprocessed=True)
            current.mark_processed()
            clipped = Polygon()
            clipped.add(Vertex(current.x,current.y))

            while True:
                current.mark_processed()
                if current.entry_exit == "entry" and not (current.prev in self): ##yes
                    current = current.link
                    while True:
                        current = current.next
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break

                elif current.entry_exit == "entry" and current.prev in self: ##yes
                    current = current.link
                    while True:
                        current = current.prev
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break

                elif current.entry_exit == "exit" and current.prev in self:
                    while True:
                        current = current.next
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break

                else: # current.entry_exit == "exit" and not (current.prev in self):
                    while True: ## I think I never tested this part?
                        current = current.next
                        clipped.add(Vertex(current.x,current.y))
                        if current.intersect:
                            break

                if current.processed:
                    result.append(clipped)
                    current = self.first.next_vertex(next_original=False, unprocessed=True)
                    break
        
        self.unclip()
       # self.perturb(redo=True)
        other.unclip()
        return result
            

    def update_entry_exit(self, other):
        
        if other.containsPoint(self.first):
            status = "exit"
        else: status = "entry"

        for vertex in self:
            if vertex.intersect:
                vertex.entry_exit = status
                if status == "entry": status = "exit"
                else: status = "entry"
    def edges(self):
        """Generate edges of the polygon."""
        start = self.first
        segments = []
        while True:
            segment = Segment(start, start.next)
            segments.append(segment)
            start = start.next
            if start == self.first: break

        return segments

    def viz(self):
        x = [i.x for i in self]
        x.append(self.first.x)
        y = [i.y for i in self]
        y.append(self.first.y)
        plt.plot(x, y, linestyle='dashed')
        plt.scatter(x, y)


## BBOX CLASS ##

class Bbox():
    
    # initialise
    def __init__(self, data):
    # using built-in `isinstance` to test what class has been used to initialise the object   

        # for Segment objects 
        if isinstance(data, Segment) == True:
            x = [data.start.x, data.end.x]
            y = [data.start.y, data.end.y]
        
        # for PointGroup objects (including Polygon)
        else:      
            x = [i.x for i in data]   # extract all x coords as a list
            y = [i.y for i in data]   # extract all y coords as a list

        # determine corners, calculate centre and area
        self.ll = Point(min(x), min(y))    # lower-left corner (min x, min y)
        self.ur = Point(max(x), max(y))    # upper-right corner (max x, max y)
        self.ctr = Point((max(x)-min(x))/2, (max(y)-min(y))/2)   # centre of box
        self.area = (abs(max(x)-min(x)))*abs((max(y)-min(y)))    # area of box
           
    # representation
    def __repr__(self):
        return f'Bounding box with lower-left {self.ll} and upper-right {self.ur}' 
    
    def __eq__(self, other): 
        if (self.ll == other.ll and self.ur == other.ur):
            return True
        else:
            return False
            # We need this method so that the class will behave sensibly in sets and dictionaries
    
    def __hash__(self):
        return hash(self.ll, other.ur)  
        
    # test for overlap between two bounding boxes
    def intersects(self, other):       
        # test if bboxes overlap (touching is not enough to be compatible with the approach to segments)
        if (self.ur.x > other.ll.x and other.ur.x > self.ll.x and
            self.ur.y > other.ll.y and other.ur.y > self.ll.y):
            return True
        else:
            return False
        
    def containsPoint(self, p):
        if (self.ur.x > p.x and p.x > self.ll.x and
            self.ur.y > p.y and p.y > self.ll.y):
            return True
        return False

def process_json_file(json_files_data):

    all_data = {}
    for json_file, return_name in json_files_data:
        if return_name not in all_data:
            all_data[return_name] = []

        with open(json_file) as f:
            data = json.load(f)

            for feature in data['features']:
                if feature['geometry']['type'] == 'Point':
                    name = feature['properties']['NAME']
                    coordinates = feature['geometry']['coordinates']
                    x, y = coordinates[:2] 
                    point = Point(x, y, name)
                    all_data[return_name].append(point)
                elif feature['geometry']['type'] == 'Polygon':
                    name = feature['properties']['name']
                    points = feature['geometry']['coordinates'][0]
                    x_coords, y_coords = zip(*points)
                    polygon = Polygon(list(zip(x_coords, y_coords)), 0, 1,name)
                    all_data[return_name].append(polygon)

    return all_data

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    sample1 = [[10,0], 
             [13,10], [6, 53], [10,0]]

    sample2 = [[-1,5], 
             [26, 20],[25, 42], [0,40], [-1, 5]]

    sample4 = [[2,10] ,[20, 10], [20, 50],
             [2,50],  [2, 10]]

    sample3 = [[0,10], [5,0], [10,10], [15,0], [20,10], [25, 0],
             [30, 20], [35, 15], [45, 0], [50, 50], [45, 40], 
             [40, 50], [30, 45], [25, 40], [20, 30], [15, 50],
             [10,35],[5,50], [0, 10]]

    sample5 = [[30, 45], [25,40], [20,30], [30,45]]


    poly1 = Polygon(sample1, xcol=0, ycol=1)
    poly2 = Polygon(sample2, xcol=0, ycol=1)
    poly3 = Polygon(sample3, 0,1)
    poly4 = Polygon(sample4, 0,1)
    poly5 = Polygon(sample5, 0,1)

    #print(f"Polygon 1 is closed: {samplePolygon1.isClosed()}")
   # print(f"Polygon 2 is closed: {poly2.isClosed()}")

    #print(poly3.containsPoint(poly2.first.next))
    plt.figure()
    poly2.viz()
    poly3.viz()

    diff1 = poly3.difference(poly2)
    diff2 = poly2.difference(poly3)
    intersection = poly3.intersection(poly2)
    union = poly3.union(poly2)

    plt.figure()
    for p in diff1: p.viz()
    plt.figure()
    for p in diff2: p.viz()
    plt.figure()
    for p in intersection: p.viz()
    plt.figure()
    for p in union: p.viz()
    #for p in diff: p.viz()

    path = "/Users/sebastiangmur/Projekte/Uzh/GEO877_Spatial_Algorithms/SpatialAlgorithms/data/swissnames_points_json_export.geojson"
    json_files_data = [(path, "name")]

    test_polys = process_json_file(json_files_data)