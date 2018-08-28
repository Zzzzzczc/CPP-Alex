# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 17:00:55 2018

@author: owner
"""

#from scipy.spatial import ConvexHull
from scipy.optimize import differential_evolution
from platypus import NSGAII, Problem, Real
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
import time
from mpl_toolkits.mplot3d import Axes3D

#print time
Time_Store = []
Time_Store.append('Starting Time')
Time_Store.append(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
Time_Store.append('               ')
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
#P[0] = [(18.898809,96.672618),(179.91666,96.672618),(179.91666,207.79762),(161.77381,207.79762),(161.77381,121.61905),(34.017857,121.61905),(34.017857,207.04167),(17.386905,207.04167)]
#P[1] = [(69.547617,145.05357),(104.32143,145.05357),(59.720238,204.77381)]
#P[2] = [(120.95238,146.56548),(86.934523,201.75),(147.41071,202.50595)]
P[0] = [(891.40431,852.57523),(891.40431,1013.8064),(919.20278,1013.8064),(919.20278,852.57523)]
P[1] = [(736.42946,650.573),(736.42946,700.61025),(788.31994,700.61025),(788.31994,650.573)]
P[2] = [(997.0385,1080.5227),(997.0385,1136.1196),(1067.4613,1136.1196),(1067.4613,1074.963)]
P[3] = [(174.20375,767.32658),(174.20375,715.4361),(372.49951,715.4361),(372.49951,769.17981)]
P[4] = [(646.77776,565.32436),(646.77776,641.30685),(817.27505,641.30685),(817.27505,569.03082)]
P[5] = [(172.35052,1228.7812),(172.35052,1267.6991),(213.12161,1267.6991),(213.12161,1228.7812)]
P[6] = [(724.61348,1193.5698),(748.70549,1230.6344),(787.62335,1204.6892),(767.2378,1173.1843)]
P[7] = [(444.77554,348.49629),(392.88506,348.49629),(392.88506,150.20053),(542.9968,150.20053),(542.9968,179.85223),(574.50174,179.85223),(574.50174,231.74271),(546.70326,231.74271),(546.70326,198.38454),(441.06907,198.38454)]
P[8] = [(974.79972,720.9958),(974.79972,774.73951),(1156.4164,774.73951),(1156.4164,717.28934)]
P[9] = [(620.83252,863.69462),(620.83252,984.15466),(704.22793,984.15466),(704.22793,865.54785)]
P[10] = [(1093.4065,502.31449),(1093.4065,669.10532),(1160.1229,669.10532),(1160.1229,500.46126)]
P[11] = [(750.55872,819.21706),(652.33745,819.21706),(652.33745,728.40872),(839.51383,728.40872),(839.51383,1008.2467),(748.70549,1008.2467)]
P[12] = [(631.95191,1023.0725),(648.14333,1083.4997),(713.49409,1063.8436),(689.40208,1000.8337)]
P[13] = [(980.35941,550.49851),(1082.2871,550.49851),(1082.2871,504.16772),(982.21264,504.16772)]
P[14] = [(736.42946,510.573),(736.42946,560.61025),(788.31994,560.61025),(788.31994,510.573)]
P[15] = [(1156.4164,1054.5775),(1130.4712,1054.5775),(1130.4712,1015.6596),(974.79972,1015.6596),(974.79972,847.01553),(1126.7647,847.01553),(1126.7647,811.80414),(1161.9761,811.80414)]
P[16] = [(290.95733,518.99357),(290.95733,620.9213),(259.4524,620.9213),(259.4524,515.28711)]
P[17] = [(448.482,772.88628),(418.8303,772.88628),(418.8303,663.54562),(261.30563,663.54562),(261.30563,628.33423),(450.33523,628.33423)]
P[18] = [(248.33301,1054.5775),(248.33301,1102.7615),(298.37026,1102.7615),(298.37026,1052.7242)]
P[19] = [(222.38777,835.89614),(222.38777,1010.0999),(322.46226,1010.0999),(322.46226,830.33645)]
P[20] = [(1098.9662,1091.6421),(1119.3518,1091.6421),(1119.3518,1123.147),(1098.9662,1123.147)]
P[21] = [(1130.4712,1074.963),(1130.4712,1134.2664),(1156.4164,1134.2664),(1156.4164,1076.8162)]
P[22] = [(967.38679,402.24),(967.38679,452.27724),(1154.5632,452.27724),(1154.5632,405.94646)]
P[23] = [(802.4492,153.90699),(802.4492,239.15563),(759.82488,239.15563),(759.82488,283.63319),(837.66059,283.63319),(837.66059,155.76022)]
P[24] = [(1156.4164,387.41415),(1126.7647,387.41415),(1126.7647,346.64305),(971.09326,346.64305),(971.09326,196.53131),(1128.6179,196.53131),(1128.6179,161.31991),(1160.1229,161.31991)]
P[25] = [(922.90924,322.55105),(922.90924,215.06362),(883.99138,215.06362),(883.99138,328.11074)]
P[26] = [(676.42946,650.573),(676.42946,700.61025),(728.31994,700.61025),(728.31994,650.573)]
P[27] = [(696.81501,152.05376),(696.81501,250.27502),(644.92453,250.27502),(644.92453,152.05376)]
P[28] = [(676.42946,510.573),(676.42946,560.61025),(728.31994,560.61025),(728.31994,510.573)]
P[29] = [(391.03183,1199.1295),(359.52689,1199.1295),(359.52689,1175.0375),(335.43488,1175.0375),(335.43488,1141.6793),(378.05921,1141.6793),(378.05921,843.30907),(415.12383,843.30907),(415.12383,1139.8261),(472.57401,1139.8261),(472.57401,1176.8907),(389.17859,1176.8907)]
P[30] = [(1006.3047,607.94868),(1006.3047,669.10532),(1061.9016,669.10532),(1061.9016,609.80191)]
P[31] = [(461.45462,850.722),(461.45462,1126.8535),(424.38999,1126.8535),(424.38999,845.1623)]
P[32] = [(793.18304,1230.6344),(837.66059,1236.1941),(837.66059,1263.9926),(765.38457,1262.1394)]
P[33] = [(209.41515,63.09865),(209.41515,33.446947),(137.13912,33.446947),(137.13912,487.48864),(224.241,487.48864),(224.241,676.51824),(142.69882,676.51824),(142.69882,1299.204),(878.43169,1299.204),(878.43169,1176.8907),(1197.1875,1176.8907),(1197.1875,116.84236),(385.47213,116.84236),(385.47213,37.15341),(311.34288,37.15341),(311.34288,61.245418),(355.82043,61.245418),(355.82043,493.04834),(524.46449,493.04834),(524.46449,459.69017),(379.91244,459.69017),(379.91244,142.7876),(548.5565,142.7876),(548.5565,168.73284),(587.47435,168.73284),(587.47435,237.3024),(550.40973,237.3024),(550.40973,250.27502),(598.59374,250.27502),(598.59374,365.17537),(633.80514,365.17537),(633.80514,179.85223),(613.41959,179.85223),(613.41959,142.7876),(848.77998,142.7876),(848.77998,318.84458),(720.90702,318.84458),(728.10492,345.70753),(702.3747,345.70753),(702.3747,361.46891),(848.77998,361.46891),(848.77998,381.85445),(876.57845,381.85445),(876.57845,148.34729),(1169.389,148.34729),(1169.389,467.1031),(878.43169,467.1031),(878.43169,446.71755),(833.95413,446.71755),(833.95413,457.83694),(596.74051,457.83694),(596.74051,491.1951),(845.07352,491.1951),(845.07352,1058.2839),(883.99138,1058.2839),(883.99138,676.51824),(861.7526,676.51824),(861.7526,493.04834),(876.57845,493.04834),(876.57845,478.22248),(909.93662,478.22248),(909.93662,500.46126),(948.85448,500.46126),(948.85448,593.12283),(965.53356,593.12283),(965.53356,491.1951),(1169.389,491.1951),(1169.389,676.51824),(943.29478,676.51824),(943.29478,743.23457),(934.02863,743.23457),(934.02863,767.32658),(961.8271,767.32658),(961.8271,706.16995),(1167.5358,706.16995),(1167.5358,1149.0923),(869.16553,1149.0923),(869.16553,1123.147),(843.22029,1123.147),(843.22029,1145.3858),(785.77012,1145.3858),(785.77012,1180.5972),(856.19291,1180.5972),(856.19291,1275.112),(163.08436,1275.112),(163.08436,1176.8907),(205.70869,1176.8907),(205.70869,1147.239),(166.79083,1147.239),(166.79083,708.02318),(379.91244,708.02318),(379.91244,776.59274),(413.2706,776.59274),(413.2706,676.51824),(252.03947,676.51824),(252.03947,493.04834),(277.98471,493.04834),(277.98471,465.24986),(168.64406,465.24986),(168.64406,63.09865)]

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
    
#    xpoly = Point(Rx,Ry).buffer(Radius)
#    DifferPoly = xpoly.difference(po)
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
#fs = po.convex_hull.difference(po)
a = po.convex_hull.difference(po)
po33 = po[33]
d = po33.convex_hull.difference(po33)
fs = d[3].difference(po)

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

#all_points = CreatePs()
#ap = [Point(p) for p in all_points]
#co = [LRF_(p.x, p.y, 400) for p in ap]
#Sb = getUnscanned(co).area  
#C: The distance between two adjacent points

#Function for optimization - Fun(A + B + C)
#Minf = A + B + C
Radius = 200
Adjust_x = []

def N_Rp_Op():
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
                return float('inf'), float('inf'), float('inf')
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
        #B:  The rest of the map's area except path coverage
        Sb = getUnscanned(co).area  
        #C: The shortest distance
        def Dis_PosC(x):
            Dist_ = 0
            Store_Comb = []
            Permutation_list = []
            list_ = []
            for i in range(len(C_all_points)/2):
                list_.append(i)
            list_Num = len(list_)
            x_list = []    
    
            for i in range(len(C_all_points)/2):
                x[i] = math.floor(x[i])
                x_list.append(int(round(x[i])))
    
            for i in range(list_Num):
                Permutation_list.append(list_[x_list[i]])
                del list_[x_list[i]]
    
    #        for i in range(len(C_all_points)/2):
    #            print Permutation_list[i]
    #        print ' '
            for i in range(len(C_all_points)/2):
                Store_Comb.append(CABC_all_points[Permutation_list[i]])                
    
            for i in range(len(C_all_points)/2 - 1):
                Dis_0 = vg.Point(Store_Comb[i][0], Store_Comb[i][1])
                Dis_1 = vg.Point(Store_Comb[i+1][0], Store_Comb[i+1][1])
                Dist = g.shortest_path(Dis_0 , Dis_1)
                _Dist_ = LineString(sp2list(Dist)).length
                Dist_ = Dist_ + _Dist_
            return Dist_

        bounds_C = []
        for i in range(len(C_all_points)/2):
            bounds_C.append((0.1,len(C_all_points)/2 - i - 0.1))
        Nearest_Dis = differential_evolution(Dis_PosC, bounds_C, maxiter=5, polish = False)
        Dc = Nearest_Dis.fun
#        for i in range(len(C_all_points)/2):
#            New_Nearest_Dis.append(Nearest_Dis.x[i])
        
        #A - Sa; B - Sb; C - Dc
        return Sa, Sb, Dc
       
    min_x = po.bounds[0] + 90
    max_x = po.bounds[2] - 10
    min_y = po.bounds[1] + 80
    max_y = po.bounds[3] - 123
    Store_Adjust_fun = []
    Store_Adjust_x = []
    
    for Num_of_Random_ in range(4,20):
        Num_of_Random = Num_of_Random_
        print Num_of_Random
        bounds_ABC = []
        for i in range(Num_of_Random * 2):
            if (i % 2 == 0):
                bounds_ABC.append((min_x,max_x))
            else:
                bounds_ABC.append((min_y,max_y))
        
        problem = Problem(Num_of_Random_ * 2, 3)
        problem.types = []
        for i in range(Num_of_Random_ * 2):
            if (i % 2 == 0):
                problem.types.append(Real(min_x,max_x))
            else:
                problem.types.append(Real(min_y,max_y))
            
        problem.function = Op_ABC
        #Execute MO-NSGAII 
        algorithm = NSGAII(problem)
        algorithm.run(1000)
        print len(algorithm.result)
        # display the results
        fig = plt.figure()    
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.scatter([s.objectives[0] for s in algorithm.result],
                   [s.objectives[1] for s in algorithm.result],
                   [s.objectives[2] for s in algorithm.result])
        ax.set_title(algorithm)
        ax.set_xlim([0, 10000])
        ax.set_ylim([0, 10000])
        ax.set_zlim([0, 2000])
        ax.view_init(elev=30.0, azim=45.0)
        ax.locator_params(nbins=12)
        plt.show()
        # display the results
        Result_Random_fun = []
        Result_Random_x = []
        for solution in algorithm.result:
            Result_Random_fun.append(solution.objectives)
            Result_Random_x.append(solution.variables)

        Show_Result_Random_fun = []
        Add_Result_Random_fun = []
        for i in range(len(algorithm.result)):
            Show_Result_Random_fun.append([Result_Random_fun[i][0],Result_Random_fun[i][1],Result_Random_fun[i][2]])
            Add_Result_Random_fun.append(Result_Random_fun[i][0]+Result_Random_fun[i][1]+Result_Random_fun[i][2])
        #Find the minimum
        Min_Result_Random_x = Add_Result_Random_fun.index(min(Add_Result_Random_fun))
        Store_Adjust_x.append(Result_Random_x[Min_Result_Random_x])
        Store_Adjust_fun.append(Add_Result_Random_fun[Min_Result_Random_x])
    #Find the Global minimum
    GlobalMin_Result_Random_x = Store_Adjust_fun.index(min(Store_Adjust_fun))
    Adjust_x.append(Store_Adjust_x[GlobalMin_Result_Random_x])

    print Adjust_x
    
    return Adjust_x

#Execute
N_Rp_Op()

CABC_all_points = []
for i in range(len(Adjust_x[0])/2):
    i = i * 2
    Temp_CABC = [Adjust_x[0][i],Adjust_x[0][i+1]]
    CABC_all_points.append(Temp_CABC)
def Dis_PosC(x):
    Dist_ = 0
    Store_Comb = []
    Permutation_list = []
    list_ = []
    for i in range(len(Adjust_x[0])/2):
        list_.append(i)
    list_Num = len(list_)
    x_list = []    

    for i in range(len(Adjust_x[0])/2):
        x[i] = math.floor(x[i])
        x_list.append(int(round(x[i])))

    for i in range(list_Num):
        Permutation_list.append(list_[x_list[i]])
        del list_[x_list[i]]

#        for i in range(len(C_all_points)/2):
#            print Permutation_list[i]
#        print ' '
    for i in range(len(Adjust_x[0])/2):
        Store_Comb.append(CABC_all_points[Permutation_list[i]])                

    for i in range(len(Adjust_x[0])/2 - 1):
        Dis_0 = vg.Point(Store_Comb[i][0], Store_Comb[i][1])
        Dis_1 = vg.Point(Store_Comb[i+1][0], Store_Comb[i+1][1])
        Dist = g.shortest_path(Dis_0 , Dis_1)
        _Dist_ = LineString(sp2list(Dist)).length
        Dist_ = Dist_ + _Dist_
    return Dist_

bounds_Ax = []
for i in range(len(Adjust_x[0])/2):
    bounds_Ax.append((0.1,len(Adjust_x[0])/2 - i - 0.1))

Nearest_Dis_Ax = differential_evolution(Dis_PosC, bounds_Ax, maxiter=5, polish = False)

#Plot
Dist_ = 0
Store_Comb = []
Permutation_list = []
list_ = []
for i in range(len(Adjust_x[0])/2):
    list_.append(i)
list_Num = len(list_)
x_list = []    

for i in range(len(Adjust_x[0])/2):
    Temp_x = math.floor(Nearest_Dis_Ax.x[i])
    x_list.append(int(round(Temp_x)))

for i in range(list_Num):
    Permutation_list.append(list_[x_list[i]])
    del list_[x_list[i]]

for i in range(len(Adjust_x[0])/2):
    Store_Comb.append(CABC_all_points[Permutation_list[i]])
    
Plot_Store_Dist = []
for i in range(len(Store_Comb) - 1):
    Dis_0 = vg.Point(Store_Comb[i][0], Store_Comb[i][1])
    Dis_1 = vg.Point(Store_Comb[i+1][0], Store_Comb[i+1][1])
    Dist = g.shortest_path(Dis_0, Dis_1)
    del Dist[len(Dist)-1]
    for p in Dist:
        Plot_Store_Dist.append((p.x, p.y))
Plot_Store_Dist.append(Store_Comb[len(Store_Comb) - 1])
   
Plot_Poly_(po, Plot_Store_Dist, 'ttt')
ap = [Point(p_Op) for p_Op in Store_Comb]
co = [LRF_(p_Op.x, p_Op.y, Radius) for p_Op in ap]
Plot_Polyv3(po, co, Plot_Store_Dist, 'yyy')

#print time
Time_Store.append('Finishing Time')
Time_Store.append(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
print Time_Store


