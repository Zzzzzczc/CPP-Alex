# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 17:00:55 2018

@author: owner
"""

#from scipy.spatial import ConvexHull
from scipy.optimize import differential_evolution
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import LinearRing

from matplotlib import pyplot
#from descartes import PolygonPatch
from shapely.ops import triangulate
import matplotlib.pyplot as plt
import pyvisgraph as vg
import random

import math


################## Functions #############################
# change (#,#) into Point(#,#)
def P2list(P):
    vc = []
    for p in P:
        vc.append(vg.Point(p[0],p[1]))
    return vc


## plot
def plot_poly(polylist, FileName):
    i = 0
    for t in polylist:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], 'k', linewidth= 1.0)
        #plt.text(t.centroid.x, t.centroid.y, i)
        
    plt.axis('scaled')
    plt.savefig(FileName, dpi = 600)

def Plot_Poly(polylist, polylist2, FileName):
    i = 0
    for t in polylist:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], 'k', linewidth= 1.0)
    i = 0
    for t in polylist2:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], 'r', linewidth= 1.0)
        #plt.text(t.centroid.x, t.centroid.y, i)
        
    plt.axis('scaled')
    plt.savefig(FileName, dpi = 600)
        
def Plot_Poly_(polylist, ls, FileName):
    i = 0
    for t in polylist:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], 'k', linewidth= 1.0)
    
    coords = np.array(ls)
    plt.plot(coords[:,0], coords[:,1], 'r', linewidth= 1.0)
    
        #plt.text(t.centroid.x, t.centroid.y, i) 
    plt.axis('scaled')
    plt.savefig(FileName, dpi = 600)
# change Point(#,#) into (#,#)     
def sp2list(sp):
    spn = []
    for p in sp:
        spn.append((p.x, p.y))
    return spn


def get_linesc(x, y, r):
    
    return 

def Plot_Polyv3(polylist, polylist2, listp, FileName):
    i = 0
    for t in polylist:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], 'k', linewidth= 1.0)
    i = 0
    for t in polylist2:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], color = (random.random(),random.random(),random.random()), linewidth= 1.0, alpha = 0.75)
        #plt.text(t.centroid.x, t.centroid.y, i)
       
    coords = np.array(listp)
    plt.plot(coords[:,0], coords[:,1], 'ro', linewidth= 1.0, markersize=1.5)
    
    plt.axis('scaled')
    plt.savefig(FileName, dpi = 600)
    
    
def Plot_Polyv4(polylist, polylist2, listp, ls, FileName):
    i = 0
    for t in polylist:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], 'k', linewidth= 1.0)
    i = 0
    for t in polylist2:
        i = i + 1
        coords = np.array(t.exterior.coords)
        plt.plot(coords[:,0], coords[:,1], color = (random.random(),random.random(),random.random()), linewidth= 1.0, alpha = 0.75)
        #plt.text(t.centroid.x, t.centroid.y, i)
       
    coords = np.array(listp)
    plt.plot(coords[:,0], coords[:,1], 'ro', linewidth= 1.0, markersize=1.5)

    coords = np.array(ls)
    plt.plot(coords[:,0], coords[:,1], 'r', linewidth= 1.0)
    
    plt.axis('scaled')
    plt.savefig(FileName, dpi = 600)


################# End of Functions ######################


INIT = 100
P = [None]*INIT
Upoly = [None]*INIT

#P[0] = [(98.273808,116.32738),(92.982141,144.29762),(145.14286,144.29762),(145.14286,113.30357)]
#P[1] = [(120.95238,182.85119),(122.46428,245.59523),(128.5119,239.54762),(128.5119,148.07738),(119.44047,148.83333)]
#P[2] = [(167.82143,161.68452),(145.14286,220.64881),(185.20833,256.17857)]
#P[3] = [(89.202381,163.95238),(80.130951,244.83928),(104.32143,249.375)]
#P[4] = [(20.410714,94.404761),(20.410714,266.00594),(27.214285,266.00594),(27.214285,100.45238),(185.20833,100.45238),(185.20833,256.17857),(67.27976,256.17857),(67.27976,142.78571),(57.452381,142.78571),(57.452381,266.00594),(190.5,266.00594),(190.5,90.625)]
P[0] = [(18.898809,96.672618),(179.91666,96.672618),(179.91666,207.79762),(161.77381,207.79762),(161.77381,121.61905),(34.017857,121.61905),(34.017857,207.04167),(17.386905,207.04167)]
P[1] = [(69.547617,145.05357),(104.32143,145.05357),(59.720238,204.77381)]
P[2] = [(120.95238,146.56548),(86.934523,201.75),(147.41071,202.50595)]
# get num of polys
P = [x for x in P if x is not None]
NUMPOLY = len(P)

# obstacles
po_ = [Polygon(p) for p in P]
plot_poly(po_, 'Mypolymap')

po = Polygon()
for p in range(NUMPOLY):
    Upoly[p] = Polygon(P[p])
    po = po.union(Upoly[p])

#polys = [P2list(p.exterior.coords) for p in P]
polys = [P2list(p) for p in P]

g = vg.VisGraph()
g.build(polys)
#plot(po, [], '')

def DifferPo(Rx, Ry, Radius):
    xpoly_Temp = Point(Rx,Ry).buffer(Radius)
    xpoly= xpoly_Temp.difference(po)
    DifferPoly = xpoly.difference(xpoly.difference(po.convex_hull))
    #DifferPoly = xpoly.symmetric_difference(po)
    return DifferPoly
    
def Intersect(Rx, Ry, Rbuffer):
    xpoly = Point(Rx,Ry).buffer(Rbuffer)
    IntersectPoly = xpoly.intersection(po)
    return IntersectPoly

def LRF_(Rx, Ry, Radius):
    d = DifferPo(Rx, Ry, Radius)
    db = d.boundary
    
    step = 2*math.pi/500
    Po = Point(Rx, Ry)
    an = 0
    #L = []
    Lp = []
    rm = 5
    while an < 2*math.pi:
        Pe = Point(Po.x + (Radius + rm)*(math.cos(an)), Po.y + (Radius + rm)*(math.sin(an))) 
        l = LineString([Po, Pe])
        #L.append(l)
        p = db.intersection(l)
        if p.type is 'MultiPoint':
            md = Po.distance(p[0])
            i = 0
            ix = 0
            while i < (len(p) - 1):
                i = i + 1
                di = Po.distance(p[i]) 
                if di < md:
                    md = di
                    ix = i
            p = p[ix]
        Lp.append(p)
        an = an + step
    
    return Polygon(LineString(Lp))

#Create points randomly inside the map (Outside any polygon)
def CreatePs():
    all_points = []
    Counterforall_points = 0
    for i in range(5):
        while (Counterforall_points < 5):    
            #Set up the boundary of the random points
            min_x = po.bounds[0]
            max_x = po.bounds[2]
            min_y = po.bounds[1]
            max_y = po.bounds[3]
            generateddata=[min_x + (max_x - min_x)*random.random(), min_y + (max_y - min_y)*random.random()]
            if po.convex_hull.contains(Point(generateddata)) == True:
                if g.point_in_polygon(vg.Point(generateddata[0], generateddata[1])) == -1:
                   if not generateddata in all_points:
                           all_points.append(generateddata)
                           Counterforall_points = Counterforall_points + 1
    return all_points     

#Computing free space (fs)
fs = po.convex_hull.difference(po)

#A: The intersect part of path coverage (include polygons)
def getIntersection(co):
    NP = len(co)
    Counter_for_randomPoly = 0
    Counter_for_other_randomPoly = 1
    #Inter_Scan = ScanofUnion_po[Counter_for_randomPoly].intersection(ScanofUnion_po[Counter_for_other_randomPoly])
    Inter_Scan = Polygon()
    SInter_Scan = Inter_Scan
    
    while(Counter_for_randomPoly < NP-1):
        while (Counter_for_other_randomPoly < NP):
            Inter_Scan = co[Counter_for_randomPoly].intersection(co[Counter_for_other_randomPoly])
            SInter_Scan = SInter_Scan.union(Inter_Scan)
            Counter_for_other_randomPoly = Counter_for_other_randomPoly + 1
        Counter_for_randomPoly = Counter_for_randomPoly + 1
        Counter_for_other_randomPoly = Counter_for_randomPoly + 1
    
    return SInter_Scan

#B: The rest of the map's area except path coverage
def getUnscanned(co):
    #S_unscanned = po.convex_hull.difference(DUnion_po.union(po))
    NP = len(co)
    # union of all coverage
    DUnion_po = Polygon()
    for i in range(NP):
        DUnion_po = DUnion_po.union(co[i])
    
    S_unscanned = fs.difference(DUnion_po)
    return S_unscanned
#C: The distance between two adjacent points

#Function for optimization - Fun(A + B + C)
#Minf = A + B + C
Radius = 400  
min_x = po.bounds[0]
max_x = po.bounds[2]
min_y = po.bounds[1]
max_y = po.bounds[3]
bounds = [(min_x,max_x),(min_y,max_y),(min_x,max_x),(min_y,max_y),(min_x,max_x),(min_y,max_y),(min_x,max_x),(min_y,max_y),(min_x,max_x),(min_y,max_y)]
def Op_ABC(C_all_points):
    #len(C_all_points) = 10 (2 * Num_of_Points)
#    for i in range(len(C_all_points)/2):
#        i = i * 2
#        print C_all_points[i],',',C_all_points[i+1]
#    print ' '
    
    #Point in Polygon?
    for i in range(len(C_all_points)/2):
        i = i * 2
        if g.point_in_polygon(vg.Point(C_all_points[i],C_all_points[i+1])) == -1:
            pass
        else:
            return float('inf')
    
    #change[a,b,c,d] to [[a,b],[c,d]]
    CABC_all_points = []
    for i in range(len(C_all_points)/2):
        i = i * 2
        Temp_CABC = [C_all_points[i],C_all_points[i+1]]
        CABC_all_points.append(Temp_CABC)
        
    ap = [Point(p) for p in CABC_all_points]
    co = [LRF_(p.x, p.y, Radius) for p in ap]
    #A: The intersect part of path coverage (Enclude polygons)
    Sa = getIntersection(co).area
    #B: The rest of the map's area except path coverage
    Sb = getUnscanned(co).area   
#    #C: The shortest distance
#    def Dis_PosC(x):
#        Dist_ = 0
#        Store_Comb = []
#        list_Num = 5
#        Permutation_list = []
#        list_ = [0, 1, 2, 3, 4]
#        x_list = []    
#        
#        for i in range(len(C_all_points)/2):
#            x[i] = math.floor(x[i])
#            x_list.append(int(round(x[i])))
#    
#        for i in range(list_Num):
#            Permutation_list.append(list_[x_list[i]])
#            del list_[x_list[i]]
#        
##        for i in range(len(C_all_points)/2):
##            print Permutation_list[i]
##        print ' '
#        for i in range(len(C_all_points)/2):
#            Store_Comb.append(CABC_all_points[Permutation_list[i]])
#            
#        for i in range(len(C_all_points)/2 - 1):
#            Dis_0 = vg.Point(Store_Comb[i][0], Store_Comb[i][1])
#            Dis_1 = vg.Point(Store_Comb[i+1][0], Store_Comb[i+1][1])
#            Dist = g.shortest_path(Dis_0 , Dis_1)
#            _Dist_ = LineString(sp2list(Dist)).length
#            Dist_ = Dist_ + _Dist_
#        return Dist_
#
#    bounds_C = [(0.1,4.9),(0.1,3.9),(0.1,2.9),(0.1,1.9),(0.1,0.9)]
#    Nearest_Dis = differential_evolution(Dis_PosC, bounds_C, maxiter=5, polish = False)
#    print CABC_all_points
    #A + B + C
    Cal_Op_ABC = Sa + Sb #+ Nearest_Dis.fun
    return Cal_Op_ABC
        
COp_ABC = differential_evolution(Op_ABC, bounds, maxiter=5, polish = False)
#COp_ABC.x = C_all_points (The best one)
COp_ABC.x
COp_ABC.fun     

#Plot
COp_all_points = []
Store_Dist = []
for i in range(len(COp_ABC.x)/2):
    i = i * 2
    Temp_CABC = [COp_ABC.x[i],COp_ABC.x[i+1]]
    COp_all_points.append(Temp_CABC)

for i in range(len(COp_ABC.x)/2 - 1):
    Dis_0 = vg.Point(COp_all_points[i][0], COp_all_points[i][1])
    Dis_1 = vg.Point(COp_all_points[i+1][0], COp_all_points[i+1][1])
    Dist = g.shortest_path(Dis_0, Dis_1)
    del Dist[len(Dist)-1]
    for p in Dist:
        Store_Dist.append((p.x, p.y))
Store_Dist.append(COp_all_points[len(COp_ABC.x)/2 - 1])
   
Plot_Poly_(po, Store_Dist, 'ttt')
Plot_Poly_(po, COp_all_points, 'xxx')
ap = [Point(p_Op) for p_Op in COp_all_points]
co = [LRF_(p_Op.x, p_Op.y, Radius) for p_Op in ap]
Plot_Polyv3(po, co, COp_all_points, 'yyy')











