import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

# Define the classes
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Segment:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        
    def intersects(self, other):
        def ccw(A, B, C):
            return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x)
        
        A, B = self.start, self.end
        C, D = other.start, other.end
        return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)

class BoundingBox:
    def __init__(self, ll, ur):
        self.ll = ll  # lower-left point
        self.ur = ur  # upper-right point
    
    def containsPoint(self, p):
        return self.ll.x <= p.x <= self.ur.x and self.ll.y <= p.y <= self.ur.y

class Polygon:
    def __init__(self, vertices):
        self.vertices = vertices
        self.size = len(vertices)
        self.bbox = self.compute_bounding_box()
        
    def compute_bounding_box(self):
        min_x = min(v.x for v in self.vertices)
        max_x = max(v.x for v in self.vertices)
        min_y = min(v.y for v in self.vertices)
        max_y = max(v.y for v in self.vertices)
        return BoundingBox(Point(min_x, min_y), Point(max_x, max_y))
        
    def containsPoint(self, p):
        if not self.bbox.containsPoint(p):
            return False
        
        ray = Segment(p, Point(self.bbox.ur.x + 1, p.y))
        count = 0
        intersections = []
        
        for i in range(self.size):
            start = self.vertices[i]
            end = self.vertices[(i + 1) % self.size]
            segment = Segment(start, end)
            if segment.intersects(ray):
                if p.y != min(start.y, end.y):
                    count += 1
                    intersections.append((start, end))
        
        return count % 2 == 1, intersections

# Define a polygon and a point
polygon_vertices = [Point(2, 1), Point(5, 1), Point(6, 4), Point(3, 5), Point(1, 3)]
polygon = Polygon(polygon_vertices)
test_point = Point(4, 3)

# Set up the plot
fig, ax = plt.subplots()
ax.set_xlim(0, 7)
ax.set_ylim(0, 6)

# Plot the polygon
polygon_patch = patches.Polygon([(v.x, v.y) for v in polygon_vertices], closed=True, fill=True, color='lightblue')
ax.add_patch(polygon_patch)

# Plot the point
ax.plot(test_point.x, test_point.y, 'ro')  # red point

# Function to update each frame
ray_end = Point(polygon.bbox.ur.x + 1, test_point.y)
ray = Segment(test_point, ray_end)
_, intersections = polygon.containsPoint(test_point)

def update(num):
    ax.clear()
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 6)
    ax.add_patch(polygon_patch)
    ax.plot(test_point.x, test_point.y, 'ro')
    ax.plot([test_point.x, ray_end.x], [test_point.y, ray_end.y], 'r--')
    
    for i in range(num):
        start, end = intersections[i]
        ax.plot([start.x, end.x], [start.y, end.y], 'g', linewidth=2.5)

# Create animation
ani = FuncAnimation(fig, update, frames=np.arange(1, len(intersections) + 1), repeat=False)

# Save animation as GIF
writer = PillowWriter(fps=1)
ani.save("ray_casting.gif", writer=writer)

plt.close()  # Close the plot as we have saved it to file
