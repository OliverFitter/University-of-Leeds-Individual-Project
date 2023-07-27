from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import itertools
import numpy as np
from collections import defaultdict
import csv
# 1: Change work directory (lines 18-26)
os.chdir(r"C:\Users\olive\OneDrive\Documents\3rd Year\Individual Project\Chapter 1") # this needs changing if you want to run script
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
BaseDir = os.getcwd()
Foldername = 1
file_location = "/home/men-comp2/mn19of/Chapter 2/" + str(Foldername) + "/"# this needs changing if you want to run script
while os.path.exists(BaseDir+"/"+str(Foldername))==True:
    Foldername = Foldername + 1
os.mkdir(BaseDir+"/"+str(Foldername))
os.chdir(BaseDir+"/"+str(Foldername))
Load_Factor=4.9 # 1.0 to 4.9 in 0.5 increments (g)
global_mesh_size=0.008 #(m)
# 2: Creates the initital model (lines 31-9699)
def CreateStandardModel():
    Mdb()

    def CreateGeometry():
        def CreateRib(scaling_factor, name,Spar_position):
            points =[(1.0,0.0),
            (0.97004,0.00571),
            (0.94818,0.01036),
            (0.92127,0.01615),
            (0.88967,0.02289),
            (0.85381,0.03049),
            (0.81427,0.03885),
            (0.77166,0.04778),
            (0.63154,0.07517),
            (0.58254,0.08347),
            (0.53322,0.09089),
            (0.48402,0.09723),
            (0.4354,0.10229),
            (0.370592592592593,0.10592),
            (0.278,0.10592),
            (0.185407407407407,0.10592),
            (0.14432,0.095), #0.0934
            (0.11363,0.08562),
            (0.08637,0.07647),
            (0.06269,0.06616),
            (0.04272,0.05491),
            (0.0265,0.043),
            (0.0141,0.03081),
            (0.00562,0.01871),
            (0.00103,0.00711),
            (0.0,0.0),
            (0.00425,-0.01115),
            (0.01414,-0.01701),
            (0.02991,-0.02205),
            (0.05099,-0.02623),
            (0.0771,-0.02948),
            (0.10798,-0.03179),
            (0.14328,-0.03319),
            (0.185407407407407,-0.03327),
            (0.278,-0.03327),
            (0.370592592592593,-0.03327),
            (0.42547,-0.0265), #-0.02426
            (0.47979,-0.02064),
            (0.53477,-0.0168),
            (0.58975,-0.01295),
            (0.64403,-0.00927),
            (0.79557,-0.00079),
            (0.8399,0.00087),
            (0.87998,0.00189),
            (0.91518,0.0023),
            (0.94491,0.00216),
            (0.96864,0.00163),
            (1.0,0.0)]

            mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=5.0)
            for i in range(len(points)-1):
                mdb.models['Model-1'].sketches['__profile__'].Line(point1=(points[i][0]*scaling_factor, points[i][1]*scaling_factor), 
                point2=(points[i+1][0]*scaling_factor, points[i+1][1]*scaling_factor))

          
            mdb.models['Model-1'].Part(dimensionality=THREE_D, name=name, type=
                DEFORMABLE_BODY)
            mdb.models['Model-1'].parts[name].BaseShell(sketch=
                mdb.models['Model-1'].sketches['__profile__'])
            del mdb.models['Model-1'].sketches['__profile__']
            
            partition_location = 0.278*scaling_factor
            
            mdb.models['Model-1'].parts[name].DatumPlaneByPrincipalPlane(offset=partition_location+0.025, 
                principalPlane=YZPLANE)
            mdb.models['Model-1'].parts[name].DatumPlaneByPrincipalPlane(offset=partition_location, 
                principalPlane=YZPLANE)
            mdb.models['Model-1'].parts[name].DatumPlaneByPrincipalPlane(offset=partition_location-0.025, 
                principalPlane=YZPLANE)
            #second spar
            mdb.models['Model-1'].parts[name].DatumPlaneByPrincipalPlane(offset=Spar_position, 
                principalPlane=YZPLANE)
            #extra partitions
            
            mdb.models['Model-1'].parts[name].DatumPlaneByPrincipalPlane(offset=0.185407407407407*scaling_factor, 
                principalPlane=YZPLANE)
            mdb.models['Model-1'].parts[name].DatumPlaneByPrincipalPlane(offset=0.370592592592593*scaling_factor, 
                principalPlane=YZPLANE)
                
            mdb.models['Model-1'].parts[name].PartitionFaceByDatumPlane(datumPlane=
                mdb.models['Model-1'].parts[name].datums[2], faces=
                mdb.models['Model-1'].parts[name].faces)
            mdb.models['Model-1'].parts[name].PartitionFaceByDatumPlane(datumPlane=
                mdb.models['Model-1'].parts[name].datums[3], faces=
                mdb.models['Model-1'].parts[name].faces)
            mdb.models['Model-1'].parts[name].PartitionFaceByDatumPlane(datumPlane=
                mdb.models['Model-1'].parts[name].datums[4], faces=
                mdb.models['Model-1'].parts[name].faces)
            mdb.models['Model-1'].parts[name].PartitionFaceByDatumPlane(datumPlane=
                mdb.models['Model-1'].parts[name].datums[5], faces=
                mdb.models['Model-1'].parts[name].faces)
            mdb.models['Model-1'].parts[name].PartitionFaceByDatumPlane(datumPlane=
                mdb.models['Model-1'].parts[name].datums[6], faces=
                mdb.models['Model-1'].parts[name].faces)
            mdb.models['Model-1'].parts[name].PartitionFaceByDatumPlane(datumPlane=
                mdb.models['Model-1'].parts[name].datums[7], faces=
                mdb.models['Model-1'].parts[name].faces)
                #####
           

        Scaling_factor_list = [(0.59,'Rib 5_5',0.435),
            (0.43,'Rib 21_5',0.326613441487559),
            (0.59,'Rib 1',0.435),
            (0.59,'Rib 2',0.435),
            (0.59,'Rib 3',0.435),
            (0.59,'Rib 4',0.435),
            (0.59,'Rib 5',0.435),
            (0.584838709677419,'Rib 6',0.431503659402824),
            (0.578645161290323,'Rib 7',0.427308050686214),
            (0.568322580645161,'Rib 8',0.420315369491863),
            (0.558,'Rib 9',0.413322688297512),
            (0.547677419354839,'Rib 10',0.406330007103161),
            (0.537354838709677,'Rib 11',0.39933732590881),
            (0.527032258064516,'Rib 12',0.392344644714459),
            (0.516709677419355,'Rib 13',0.385351963520108),
            (0.506387096774194,'Rib 14',0.378359282325757),
            (0.496064516129032,'Rib 15',0.371366601131406),
            (0.485741935483871,'Rib 16',0.364373919937054),
            (0.47541935483871,'Rib 17',0.357381238742703),
            (0.465096774193548,'Rib 18',0.350388557548352),
            (0.454774193548387,'Rib 19',0.343395876354001),
            (0.444451612903226,'Rib 20',0.33640319515965),
            (0.434129032258065,'Rib 21',0.329410513965299),
            (0.404153846153846,'Rib 22',0.302371867937449),
            (0.361076923076923,'Rib 23',0.261969245353934),
            (0.318,'Rib 24',0.221566622770418),
            (0.29,'Rib 25',0.195304918091133),
            (0.59,'Base Rib',0.435)]


        for i in Scaling_factor_list:
            CreateRib(i[0], i[1], i[2])

        #start of cut
        def CreateCutout(scaling_factor, name, Spar_position):
            fillet=0.005
            #front cutout
            mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.006, name='__profile__', 
                sheetSize=0.273, transform=
                mdb.models['Model-1'].parts[name].MakeSketchTransform(
                sketchPlane=mdb.models['Model-1'].parts[name].faces[1], 
                sketchPlaneSide=SIDE1, 
                sketchUpEdge=mdb.models['Model-1'].parts[name].edges[4], 
                sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
            mdb.models['Model-1'].sketches['__profile__'].sketchOptions.setValues(
                decimalPlaces=3)
            mdb.models['Model-1'].parts[name].projectReferencesOntoSketch(filter=
                COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
            mdb.models['Model-1'].parts[name].projectEdgesOntoSketch(
                constrainToBackground=False, edges=(
                mdb.models['Model-1'].parts[name].edges[4], 
                mdb.models['Model-1'].parts[name].edges[5], 
                mdb.models['Model-1'].parts[name].edges[6], 
                mdb.models['Model-1'].parts[name].edges[7], 
                mdb.models['Model-1'].parts[name].edges[8], 
                mdb.models['Model-1'].parts[name].edges[9], 
                mdb.models['Model-1'].parts[name].edges[10], 
                mdb.models['Model-1'].parts[name].edges[11], 
                mdb.models['Model-1'].parts[name].edges[12], 
                mdb.models['Model-1'].parts[name].edges[13], 
                mdb.models['Model-1'].parts[name].edges[14], 
                mdb.models['Model-1'].parts[name].edges[15], 
                mdb.models['Model-1'].parts[name].edges[16], 
                mdb.models['Model-1'].parts[name].edges[17], 
                mdb.models['Model-1'].parts[name].edges[18], 
                mdb.models['Model-1'].parts[name].edges[19], 
                mdb.models['Model-1'].parts[name].edges[20], 
                mdb.models['Model-1'].parts[name].edges[21], 
                mdb.models['Model-1'].parts[name].edges[22]), sketch=
                mdb.models['Model-1'].sketches['__profile__'])
            mdb.models['Model-1'].sketches['__profile__'].scale(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[6], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[7], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[8], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[9], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[10], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[11], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[12], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[13], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[14], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[15], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[16], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[17], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[18], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[19], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[20], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[21], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[22], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[23], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[24], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[39], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[67], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[68], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[69], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[70], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[71], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[72], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[73], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[74], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[75], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[76], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[77], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[78], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[79], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[80], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[81], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[82], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[83], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[84], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[85]), scaleCenter=(
                0.0, 0.0), scaleValue=0.6)
            mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[17], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[16]))
            mdb.models['Model-1'].sketches['__profile__'].move(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[67], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[68], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[69], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[70], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[71], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[72], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[73], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[82], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[83], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[84], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[85], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[80], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[81], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[74], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[75], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[76], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[78], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[79], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[77]), vector=(
                (0.0233657777777779/0.59)*scaling_factor, (0.00814095/0.59)*scaling_factor))
                  
            mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
                mdb.models['Model-1'].sketches['__profile__'].geometry[68], curve2=
                mdb.models['Model-1'].sketches['__profile__'].geometry[67], nearPoint1=(
                (0.0797769501805305/0.59)*scaling_factor, (0.0418906882405281/0.59)*scaling_factor), nearPoint2=((0.0886033922433853/0.59)*scaling_factor, 
                (0.0392181575298309/0.59)*scaling_factor), radius=fillet)
            mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
                mdb.models['Model-1'].sketches['__profile__'].geometry[67], curve2=
                mdb.models['Model-1'].sketches['__profile__'].geometry[85], nearPoint1=(
                (0.088445782661438/0.59)*scaling_factor, (0.00683338567614555/0.59)*scaling_factor), nearPoint2=((0.0845054015517235/0.59)*scaling_factor, 
                (-0.00338511168956757/0.59)*scaling_factor), radius=fillet)
                
            mdb.models['Model-1'].parts[name].CutExtrude(flipExtrudeDirection=OFF, 
                sketch=mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=
                RIGHT, sketchPlane=mdb.models['Model-1'].parts[name].faces[1], 
                sketchPlaneSide=SIDE1, sketchUpEdge=
                mdb.models['Model-1'].parts[name].edges[4])
            del mdb.models['Model-1'].sketches['__profile__']
            
            
            
            #rear cut out
            mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.02, name='__profile__', 
                sheetSize=0.88, transform=
                mdb.models['Model-1'].parts[name].MakeSketchTransform(
                sketchPlane=mdb.models['Model-1'].parts[name].faces[2], 
                sketchPlaneSide=SIDE1, 
                sketchUpEdge=mdb.models['Model-1'].parts[name].edges[50], 
                sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
            mdb.models['Model-1'].sketches['__profile__'].sketchOptions.setValues(
                decimalPlaces=3)
            mdb.models['Model-1'].parts[name].projectReferencesOntoSketch(filter=
                COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
            mdb.models['Model-1'].parts[name].projectEdgesOntoSketch(
                constrainToBackground=False, edges=(
                mdb.models['Model-1'].parts[name].edges[0], 
                mdb.models['Model-1'].parts[name].edges[44], 
                mdb.models['Model-1'].parts[name].edges[45], 
                mdb.models['Model-1'].parts[name].edges[46], 
                mdb.models['Model-1'].parts[name].edges[47], 
                mdb.models['Model-1'].parts[name].edges[48], 
                mdb.models['Model-1'].parts[name].edges[49], 
                mdb.models['Model-1'].parts[name].edges[50], 
                mdb.models['Model-1'].parts[name].edges[51], 
                mdb.models['Model-1'].parts[name].edges[52], 
                mdb.models['Model-1'].parts[name].edges[53], 
                mdb.models['Model-1'].parts[name].edges[54], 
                mdb.models['Model-1'].parts[name].edges[55], 
                mdb.models['Model-1'].parts[name].edges[56]), sketch=
                mdb.models['Model-1'].sketches['__profile__'])
            mdb.models['Model-1'].sketches['__profile__'].scale(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[2], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[46], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[47], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[48], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[49], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[50], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[51], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[52], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[53], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[54], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[55], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[56], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[57], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[58], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[59], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[72], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[88], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[89], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[90], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[91], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[92], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[93], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[94], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[95], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[96], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[97], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[98], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[99], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[100], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[101]), scaleCenter=(
                (0.21864962962963/0.59)*scaling_factor, (0.02143175/0.59)*scaling_factor), scaleValue=0.6)
            mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[2], ))
            mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[46], ))
            mdb.models['Model-1'].sketches['__profile__'].move(objectList=(
                mdb.models['Model-1'].sketches['__profile__'].geometry[88], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[89], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[90], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[91], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[92], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[97], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[98], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[99], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[100], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[101], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[93], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[94], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[96], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[97], 
                mdb.models['Model-1'].sketches['__profile__'].geometry[95]), vector=(
                (0.0253503703703701/0.59)*scaling_factor, (-0.000431749999999998/0.59)*scaling_factor))
                
            
            
        # top left
            mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
                mdb.models['Model-1'].sketches['__profile__'].geometry[101], curve2=
                mdb.models['Model-1'].sketches['__profile__'].geometry[88], nearPoint1=(
                (0.253383994102478/0.59)*scaling_factor, (0.0456404611468315/0.59)*scaling_factor), nearPoint2=((0.245227456092834/0.59)*scaling_factor, 
                (0.0359553918242455/0.59)*scaling_factor), radius=fillet)
                
                
        # bottom left
            mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
                mdb.models['Model-1'].sketches['__profile__'].geometry[88], curve2=
                mdb.models['Model-1'].sketches['__profile__'].geometry[89], nearPoint1=(
                (0.243285417556763/0.59)*scaling_factor, (0.00883720815181732/0.59)*scaling_factor), nearPoint2=((0.248723119497299/0.59)*scaling_factor, 
                (-0.00162266939878464/0.59)*scaling_factor), radius=fillet)
                
                
            mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
                mdb.models['Model-1'].sketches['__profile__'].geometry[96], curve2=
                mdb.models['Model-1'].sketches['__profile__'].geometry[95], nearPoint1=(
                (0.363302975893021/0.59)*scaling_factor, (0.0305317640304565/0.59)*scaling_factor), nearPoint2=((0.371459513902664/0.59)*scaling_factor, 
                (0.023171104490757/0.59)*scaling_factor), radius=fillet)
                
                
                
            mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
                mdb.models['Model-1'].sketches['__profile__'].geometry[95], curve2=
                mdb.models['Model-1'].sketches['__profile__'].geometry[94], nearPoint1=(
                (0.372624725103378/0.59)*scaling_factor, (0.0154230520129204/0.59)*scaling_factor), nearPoint2=((0.364856600761414/0.59)*scaling_factor, 
                (0.00457577407360077/0.59)*scaling_factor), radius=fillet)
                
                
            mdb.models['Model-1'].parts[name].CutExtrude(flipExtrudeDirection=OFF, 
                sketch=mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=
                RIGHT, sketchPlane=mdb.models['Model-1'].parts[name].faces[2], 
                sketchPlaneSide=SIDE1, sketchUpEdge=
                mdb.models['Model-1'].parts[name].edges[50])
            del mdb.models['Model-1'].sketches['__profile__']


        #end of cut

        for i in Scaling_factor_list[:-2]:
            CreateCutout(i[0], i[1], i[2])


        Translation_coordinates_list=[(0,0,0.1,'Rib 1','Rib 1_1',0),
            (0,0,0.2,'Rib 2','Rib 2_1',0),
            (0,0,0.3,'Rib 3','Rib 3_1',0),
            (0,0,0.4,'Rib 4','Rib 4_1',0),
            (0,0,0.46,'Rib 5','Rib 5_1',0),
            (0.00349634059717552,0,0.56,'Rib 6','Rib 6_1',0),
            (0.00769194931378615,0,0.62,'Rib 7','Rib 7_1',0),
            (0.0146846305081372,0,0.72,'Rib 8','Rib 8_1',0),
            (0.0216773117024882,0,0.82,'Rib 9','Rib 9_1',0),
            (0.0286699928968393,0,0.92,'Rib 10','Rib 10_1',0),
            (0.0356626740911903,0,1.02,'Rib 11','Rib 11_1',0),
            (0.0426553552855414,0,1.12,'Rib 12','Rib 12_1',0),
            (0.0496480364798924,0,1.22,'Rib 13','Rib 13_1',0),
            (0.0566407176742434,0,1.32,'Rib 14','Rib 14_1',0),
            (0.0636333988685945,0,1.42,'Rib 15','Rib 15_1',0),
            (0.0706260800629455,0,1.52,'Rib 16','Rib 16_1',0),
            (0.0776187612572966,0,1.62,'Rib 17','Rib 17_1',0),
            (0.0846114424516476,0,1.72,'Rib 18','Rib 18_1',0),
            (0.0916041236459986,0,1.82,'Rib 19','Rib 19_1',0),
            (0.0985968048403497,0,1.92,'Rib 20','Rib 20_1',0),
            (0.105589486034701,0,2.02,'Rib 21','Rib 21_1',0),
            (0.132628132062551,-0.0116628185482631,2.12,'Rib 22','Rib 22_1',11),
            (0.173030754646066,-0.031100849462035,2.22,'Rib 23','Rib 23_1',11),
            (0.213433377229582,-0.0505388803758068,2.32,'Rib 24','Rib 24_1',11),
            (0.239695081908867,-0.0631736004697585,2.385,'Rib 25','Rib 25_1',11),
            (0,0,0.51,'Rib 5_5','Rib 5_5_1',0),
            (0.108386558512441,0,2.06,'Rib 21_5','Rib 21_5_1',0)]

        #cutout
        mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=
            mdb.models['Model-1'].parts['Rib 24'].features['Cut extrude-2'].sketch)
        mdb.models['Model-1'].parts['Rib 24'].projectReferencesOntoSketch(filter=
            COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__edit__'], 
            upToFeature=
            mdb.models['Model-1'].parts['Rib 24'].features['Cut extrude-2'])
        mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
            mdb.models['Model-1'].sketches['__edit__'].geometry[104], 
            mdb.models['Model-1'].sketches['__edit__'].geometry[105]))
        mdb.models['Model-1'].sketches['__edit__'].delete(objectList=(
            mdb.models['Model-1'].sketches['__edit__'].geometry[95], ))
        mdb.models['Model-1'].sketches['__edit__'].Line(point1=(0.196126460245368, 
            0.0158321863000336), point2=(0.195854514229442, 0.00330018128747136))
        mdb.models['Model-1'].sketches['__edit__'].FilletByRadius(curve1=
            mdb.models['Model-1'].sketches['__edit__'].geometry[96], curve2=
            mdb.models['Model-1'].sketches['__edit__'].geometry[192], nearPoint1=(
            0.194400727748871, 0.016193425282836), nearPoint2=(0.195895239710808, 
            0.0147812236100435), radius=0.005)
        mdb.models['Model-1'].sketches['__edit__'].FilletByRadius(curve1=
            mdb.models['Model-1'].sketches['__edit__'].geometry[192], curve2=
            mdb.models['Model-1'].sketches['__edit__'].geometry[94], nearPoint1=(
            0.195816576480865, 0.00779868103563786), nearPoint2=(0.1931421905756, 
            0.00324825849384069), radius=0.005)
        mdb.models['Model-1'].parts['Rib 24'].features['Cut extrude-2'].setValues(
            sketch=mdb.models['Model-1'].sketches['__edit__'])
        del mdb.models['Model-1'].sketches['__edit__']
        mdb.models['Model-1'].parts['Rib 24'].regenerate()
        mdb.models['Model-1'].parts['Rib 24'].regenerate()








        mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
        # def AssembleRib(names,name_1):
            # mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=name_1, part=
                # mdb.models['Model-1'].parts[names])
        mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Base Rib_1', part=
                mdb.models['Model-1'].parts['Base Rib'])#change back 
        def RotateRib(name_1,Angle):
            mdb.models['Model-1'].rootAssembly.rotate(angle=Angle, axisDirection=(10.0, 0.0, 
                0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=(name_1, ))


        def TranslateRibs(x_coord, y_coord, z_coord, name_1):
            mdb.models['Model-1'].rootAssembly.translate(instanceList=(name_1, ), 
                vector=(x_coord, y_coord, z_coord))
        for i in Translation_coordinates_list:
            # AssembleRib(i[3], i[4])
            mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=i[4], part=
                mdb.models['Model-1'].parts[i[3]])
        for i in Translation_coordinates_list:
            RotateRib(i[4], i[5])
        for i in Translation_coordinates_list:
            TranslateRibs(i[0], i[1], i[2], i[4])
        #creates ribs part from assembly
        mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
            instances=(mdb.models['Model-1'].rootAssembly.instances['Rib 1_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 2_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 3_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 4_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 5_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 6_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 7_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 8_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 9_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 10_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 11_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 12_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 13_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 14_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 15_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 16_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 17_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 18_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 19_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 20_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 21_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 22_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 23_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 24_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 25_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 5_5_1'], 
            mdb.models['Model-1'].rootAssembly.instances['Rib 21_5_1'],
            mdb.models['Model-1'].rootAssembly.instances['Base Rib_1']), name='Wing', 
            originalInstances=SUPPRESS)

        #spar 1
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2644], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1311], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1237], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1470], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1237], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1213], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1139], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1572], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1139], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1115], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1041], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2654], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1041], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1774], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1041], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1017], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1647], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1705], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1551], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1454], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1357], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1260], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1163], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1066], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[969], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[968], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1895], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[1044], ), (
            mdb.models['Model-1'].parts['Wing'].edges[943], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[943], ), (
            mdb.models['Model-1'].parts['Wing'].edges[922], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[845], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1999], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[845], ), (
            mdb.models['Model-1'].parts['Wing'].edges[824], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[747], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2101], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[747], ), (
            mdb.models['Model-1'].parts['Wing'].edges[726], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[649], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2203], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[649], ), (
            mdb.models['Model-1'].parts['Wing'].edges[628], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[551], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2305], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[551], ), (
            mdb.models['Model-1'].parts['Wing'].edges[530], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[453], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2407], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[453], ), (
            mdb.models['Model-1'].parts['Wing'].edges[432], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[355], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2509], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[355], ), (
            mdb.models['Model-1'].parts['Wing'].edges[334], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[257], ), (
            mdb.models['Model-1'].parts['Wing'].edges[79], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[159], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2613], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[159], ), (
            mdb.models['Model-1'].parts['Wing'].edges[138], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[62], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2716], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[61], ), (
            mdb.models['Model-1'].parts['Wing'].edges[40], )), startCondition=NONE)
        #end of spar 1
        #spar 2
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2420], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2516], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2364], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2268], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2172], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2076], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1980], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1884], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1788], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1668], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2631], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1568], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1471], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1374], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1277], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1180], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[1083], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[986], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[889], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[792], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[695], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[598], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[501], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[404], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[307], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[2758], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[209], )), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], ), (
            mdb.models['Model-1'].parts['Wing'].edges[112], )), startCondition=NONE)
        #end of spar




        #skin
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2451], 
            mdb.models['Model-1'].parts['Wing'].edges[2452], 
            mdb.models['Model-1'].parts['Wing'].edges[2454], 
            mdb.models['Model-1'].parts['Wing'].edges[2455], 
            mdb.models['Model-1'].parts['Wing'].edges[2456], 
            mdb.models['Model-1'].parts['Wing'].edges[2457], 
            mdb.models['Model-1'].parts['Wing'].edges[2458], 
            mdb.models['Model-1'].parts['Wing'].edges[2459], 
            mdb.models['Model-1'].parts['Wing'].edges[2460], 
            mdb.models['Model-1'].parts['Wing'].edges[2461], 
            mdb.models['Model-1'].parts['Wing'].edges[2462], 
            mdb.models['Model-1'].parts['Wing'].edges[2463], 
            mdb.models['Model-1'].parts['Wing'].edges[2464], 
            mdb.models['Model-1'].parts['Wing'].edges[2465], 
            mdb.models['Model-1'].parts['Wing'].edges[2466], 
            mdb.models['Model-1'].parts['Wing'].edges[2467], 
            mdb.models['Model-1'].parts['Wing'].edges[2468], 
            mdb.models['Model-1'].parts['Wing'].edges[2469], 
            mdb.models['Model-1'].parts['Wing'].edges[2470], 
            mdb.models['Model-1'].parts['Wing'].edges[2471], 
            mdb.models['Model-1'].parts['Wing'].edges[2472], 
            mdb.models['Model-1'].parts['Wing'].edges[2473], 
            mdb.models['Model-1'].parts['Wing'].edges[2474], 
            mdb.models['Model-1'].parts['Wing'].edges[2475], 
            mdb.models['Model-1'].parts['Wing'].edges[2476], 
            mdb.models['Model-1'].parts['Wing'].edges[2477], 
            mdb.models['Model-1'].parts['Wing'].edges[2478], 
            mdb.models['Model-1'].parts['Wing'].edges[2479], 
            mdb.models['Model-1'].parts['Wing'].edges[2480], 
            mdb.models['Model-1'].parts['Wing'].edges[2481], 
            mdb.models['Model-1'].parts['Wing'].edges[2482], 
            mdb.models['Model-1'].parts['Wing'].edges[2483], 
            mdb.models['Model-1'].parts['Wing'].edges[2484], 
            mdb.models['Model-1'].parts['Wing'].edges[2485], 
            mdb.models['Model-1'].parts['Wing'].edges[2486], 
            mdb.models['Model-1'].parts['Wing'].edges[2488], 
            mdb.models['Model-1'].parts['Wing'].edges[2489], 
            mdb.models['Model-1'].parts['Wing'].edges[2490], 
            mdb.models['Model-1'].parts['Wing'].edges[2491], 
            mdb.models['Model-1'].parts['Wing'].edges[2492], 
            mdb.models['Model-1'].parts['Wing'].edges[2493], 
            mdb.models['Model-1'].parts['Wing'].edges[2494], 
            mdb.models['Model-1'].parts['Wing'].edges[2495], 
            mdb.models['Model-1'].parts['Wing'].edges[2496], 
            mdb.models['Model-1'].parts['Wing'].edges[2497], 
            mdb.models['Model-1'].parts['Wing'].edges[2498], 
            mdb.models['Model-1'].parts['Wing'].edges[2499], 
            mdb.models['Model-1'].parts['Wing'].edges[2500], 
            mdb.models['Model-1'].parts['Wing'].edges[2501], 
            mdb.models['Model-1'].parts['Wing'].edges[2502], 
            mdb.models['Model-1'].parts['Wing'].edges[2503], 
            mdb.models['Model-1'].parts['Wing'].edges[2504], 
            mdb.models['Model-1'].parts['Wing'].edges[2505]), (
            mdb.models['Model-1'].parts['Wing'].edges[2507], 
            mdb.models['Model-1'].parts['Wing'].edges[2508], 
            mdb.models['Model-1'].parts['Wing'].edges[2531], 
            mdb.models['Model-1'].parts['Wing'].edges[2532], 
            mdb.models['Model-1'].parts['Wing'].edges[2533], 
            mdb.models['Model-1'].parts['Wing'].edges[2534], 
            mdb.models['Model-1'].parts['Wing'].edges[2535], 
            mdb.models['Model-1'].parts['Wing'].edges[2536], 
            mdb.models['Model-1'].parts['Wing'].edges[2537], 
            mdb.models['Model-1'].parts['Wing'].edges[2538], 
            mdb.models['Model-1'].parts['Wing'].edges[2539], 
            mdb.models['Model-1'].parts['Wing'].edges[2540], 
            mdb.models['Model-1'].parts['Wing'].edges[2541], 
            mdb.models['Model-1'].parts['Wing'].edges[2542], 
            mdb.models['Model-1'].parts['Wing'].edges[2543], 
            mdb.models['Model-1'].parts['Wing'].edges[2544], 
            mdb.models['Model-1'].parts['Wing'].edges[2545], 
            mdb.models['Model-1'].parts['Wing'].edges[2546], 
            mdb.models['Model-1'].parts['Wing'].edges[2547], 
            mdb.models['Model-1'].parts['Wing'].edges[2548], 
            mdb.models['Model-1'].parts['Wing'].edges[2567], 
            mdb.models['Model-1'].parts['Wing'].edges[2568], 
            mdb.models['Model-1'].parts['Wing'].edges[2569], 
            mdb.models['Model-1'].parts['Wing'].edges[2570], 
            mdb.models['Model-1'].parts['Wing'].edges[2571], 
            mdb.models['Model-1'].parts['Wing'].edges[2572], 
            mdb.models['Model-1'].parts['Wing'].edges[2573], 
            mdb.models['Model-1'].parts['Wing'].edges[2574], 
            mdb.models['Model-1'].parts['Wing'].edges[2575], 
            mdb.models['Model-1'].parts['Wing'].edges[2576], 
            mdb.models['Model-1'].parts['Wing'].edges[2577], 
            mdb.models['Model-1'].parts['Wing'].edges[2578], 
            mdb.models['Model-1'].parts['Wing'].edges[2579], 
            mdb.models['Model-1'].parts['Wing'].edges[2580], 
            mdb.models['Model-1'].parts['Wing'].edges[2581], 
            mdb.models['Model-1'].parts['Wing'].edges[2583], 
            mdb.models['Model-1'].parts['Wing'].edges[2584], 
            mdb.models['Model-1'].parts['Wing'].edges[2585], 
            mdb.models['Model-1'].parts['Wing'].edges[2586], 
            mdb.models['Model-1'].parts['Wing'].edges[2587], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2508], 
            mdb.models['Model-1'].parts['Wing'].edges[2509], 
            mdb.models['Model-1'].parts['Wing'].edges[2532], 
            mdb.models['Model-1'].parts['Wing'].edges[2533], 
            mdb.models['Model-1'].parts['Wing'].edges[2534], 
            mdb.models['Model-1'].parts['Wing'].edges[2535], 
            mdb.models['Model-1'].parts['Wing'].edges[2536], 
            mdb.models['Model-1'].parts['Wing'].edges[2537], 
            mdb.models['Model-1'].parts['Wing'].edges[2538], 
            mdb.models['Model-1'].parts['Wing'].edges[2539], 
            mdb.models['Model-1'].parts['Wing'].edges[2540], 
            mdb.models['Model-1'].parts['Wing'].edges[2541], 
            mdb.models['Model-1'].parts['Wing'].edges[2542], 
            mdb.models['Model-1'].parts['Wing'].edges[2543], 
            mdb.models['Model-1'].parts['Wing'].edges[2544], 
            mdb.models['Model-1'].parts['Wing'].edges[2545], 
            mdb.models['Model-1'].parts['Wing'].edges[2546], 
            mdb.models['Model-1'].parts['Wing'].edges[2547], 
            mdb.models['Model-1'].parts['Wing'].edges[2548], 
            mdb.models['Model-1'].parts['Wing'].edges[2549], 
            mdb.models['Model-1'].parts['Wing'].edges[2568], 
            mdb.models['Model-1'].parts['Wing'].edges[2569], 
            mdb.models['Model-1'].parts['Wing'].edges[2570], 
            mdb.models['Model-1'].parts['Wing'].edges[2571], 
            mdb.models['Model-1'].parts['Wing'].edges[2572], 
            mdb.models['Model-1'].parts['Wing'].edges[2573], 
            mdb.models['Model-1'].parts['Wing'].edges[2574], 
            mdb.models['Model-1'].parts['Wing'].edges[2575], 
            mdb.models['Model-1'].parts['Wing'].edges[2576], 
            mdb.models['Model-1'].parts['Wing'].edges[2577], 
            mdb.models['Model-1'].parts['Wing'].edges[2578], 
            mdb.models['Model-1'].parts['Wing'].edges[2579], 
            mdb.models['Model-1'].parts['Wing'].edges[2580], 
            mdb.models['Model-1'].parts['Wing'].edges[2581], 
            mdb.models['Model-1'].parts['Wing'].edges[2582], 
            mdb.models['Model-1'].parts['Wing'].edges[2584], 
            mdb.models['Model-1'].parts['Wing'].edges[2585], 
            mdb.models['Model-1'].parts['Wing'].edges[2586], 
            mdb.models['Model-1'].parts['Wing'].edges[2587], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2512], 
            mdb.models['Model-1'].parts['Wing'].edges[2513], 
            mdb.models['Model-1'].parts['Wing'].edges[2536], 
            mdb.models['Model-1'].parts['Wing'].edges[2537], 
            mdb.models['Model-1'].parts['Wing'].edges[2538], 
            mdb.models['Model-1'].parts['Wing'].edges[2539], 
            mdb.models['Model-1'].parts['Wing'].edges[2540], 
            mdb.models['Model-1'].parts['Wing'].edges[2541], 
            mdb.models['Model-1'].parts['Wing'].edges[2542], 
            mdb.models['Model-1'].parts['Wing'].edges[2543], 
            mdb.models['Model-1'].parts['Wing'].edges[2544], 
            mdb.models['Model-1'].parts['Wing'].edges[2545], 
            mdb.models['Model-1'].parts['Wing'].edges[2546], 
            mdb.models['Model-1'].parts['Wing'].edges[2547], 
            mdb.models['Model-1'].parts['Wing'].edges[2548], 
            mdb.models['Model-1'].parts['Wing'].edges[2549], 
            mdb.models['Model-1'].parts['Wing'].edges[2550], 
            mdb.models['Model-1'].parts['Wing'].edges[2551], 
            mdb.models['Model-1'].parts['Wing'].edges[2552], 
            mdb.models['Model-1'].parts['Wing'].edges[2553], 
            mdb.models['Model-1'].parts['Wing'].edges[2572], 
            mdb.models['Model-1'].parts['Wing'].edges[2573], 
            mdb.models['Model-1'].parts['Wing'].edges[2574], 
            mdb.models['Model-1'].parts['Wing'].edges[2575], 
            mdb.models['Model-1'].parts['Wing'].edges[2576], 
            mdb.models['Model-1'].parts['Wing'].edges[2577], 
            mdb.models['Model-1'].parts['Wing'].edges[2578], 
            mdb.models['Model-1'].parts['Wing'].edges[2579], 
            mdb.models['Model-1'].parts['Wing'].edges[2580], 
            mdb.models['Model-1'].parts['Wing'].edges[2581], 
            mdb.models['Model-1'].parts['Wing'].edges[2582], 
            mdb.models['Model-1'].parts['Wing'].edges[2583], 
            mdb.models['Model-1'].parts['Wing'].edges[2584], 
            mdb.models['Model-1'].parts['Wing'].edges[2585], 
            mdb.models['Model-1'].parts['Wing'].edges[2586], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2603], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2516], 
            mdb.models['Model-1'].parts['Wing'].edges[2517], 
            mdb.models['Model-1'].parts['Wing'].edges[2540], 
            mdb.models['Model-1'].parts['Wing'].edges[2541], 
            mdb.models['Model-1'].parts['Wing'].edges[2542], 
            mdb.models['Model-1'].parts['Wing'].edges[2543], 
            mdb.models['Model-1'].parts['Wing'].edges[2544], 
            mdb.models['Model-1'].parts['Wing'].edges[2545], 
            mdb.models['Model-1'].parts['Wing'].edges[2546], 
            mdb.models['Model-1'].parts['Wing'].edges[2547], 
            mdb.models['Model-1'].parts['Wing'].edges[2548], 
            mdb.models['Model-1'].parts['Wing'].edges[2549], 
            mdb.models['Model-1'].parts['Wing'].edges[2550], 
            mdb.models['Model-1'].parts['Wing'].edges[2551], 
            mdb.models['Model-1'].parts['Wing'].edges[2552], 
            mdb.models['Model-1'].parts['Wing'].edges[2553], 
            mdb.models['Model-1'].parts['Wing'].edges[2554], 
            mdb.models['Model-1'].parts['Wing'].edges[2555], 
            mdb.models['Model-1'].parts['Wing'].edges[2556], 
            mdb.models['Model-1'].parts['Wing'].edges[2557], 
            mdb.models['Model-1'].parts['Wing'].edges[2576], 
            mdb.models['Model-1'].parts['Wing'].edges[2577], 
            mdb.models['Model-1'].parts['Wing'].edges[2578], 
            mdb.models['Model-1'].parts['Wing'].edges[2579], 
            mdb.models['Model-1'].parts['Wing'].edges[2580], 
            mdb.models['Model-1'].parts['Wing'].edges[2581], 
            mdb.models['Model-1'].parts['Wing'].edges[2582], 
            mdb.models['Model-1'].parts['Wing'].edges[2583], 
            mdb.models['Model-1'].parts['Wing'].edges[2584], 
            mdb.models['Model-1'].parts['Wing'].edges[2585], 
            mdb.models['Model-1'].parts['Wing'].edges[2586], 
            mdb.models['Model-1'].parts['Wing'].edges[2587], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2603], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605], 
            mdb.models['Model-1'].parts['Wing'].edges[2606], 
            mdb.models['Model-1'].parts['Wing'].edges[2607], 
            mdb.models['Model-1'].parts['Wing'].edges[2608], 
            mdb.models['Model-1'].parts['Wing'].edges[2609])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2520], 
            mdb.models['Model-1'].parts['Wing'].edges[2521], 
            mdb.models['Model-1'].parts['Wing'].edges[2544], 
            mdb.models['Model-1'].parts['Wing'].edges[2545], 
            mdb.models['Model-1'].parts['Wing'].edges[2546], 
            mdb.models['Model-1'].parts['Wing'].edges[2547], 
            mdb.models['Model-1'].parts['Wing'].edges[2548], 
            mdb.models['Model-1'].parts['Wing'].edges[2549], 
            mdb.models['Model-1'].parts['Wing'].edges[2550], 
            mdb.models['Model-1'].parts['Wing'].edges[2551], 
            mdb.models['Model-1'].parts['Wing'].edges[2552], 
            mdb.models['Model-1'].parts['Wing'].edges[2553], 
            mdb.models['Model-1'].parts['Wing'].edges[2554], 
            mdb.models['Model-1'].parts['Wing'].edges[2555], 
            mdb.models['Model-1'].parts['Wing'].edges[2556], 
            mdb.models['Model-1'].parts['Wing'].edges[2557], 
            mdb.models['Model-1'].parts['Wing'].edges[2558], 
            mdb.models['Model-1'].parts['Wing'].edges[2559], 
            mdb.models['Model-1'].parts['Wing'].edges[2560], 
            mdb.models['Model-1'].parts['Wing'].edges[2561], 
            mdb.models['Model-1'].parts['Wing'].edges[2580], 
            mdb.models['Model-1'].parts['Wing'].edges[2581], 
            mdb.models['Model-1'].parts['Wing'].edges[2582], 
            mdb.models['Model-1'].parts['Wing'].edges[2583], 
            mdb.models['Model-1'].parts['Wing'].edges[2584], 
            mdb.models['Model-1'].parts['Wing'].edges[2585], 
            mdb.models['Model-1'].parts['Wing'].edges[2586], 
            mdb.models['Model-1'].parts['Wing'].edges[2587], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2603], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605], 
            mdb.models['Model-1'].parts['Wing'].edges[2606], 
            mdb.models['Model-1'].parts['Wing'].edges[2607], 
            mdb.models['Model-1'].parts['Wing'].edges[2608], 
            mdb.models['Model-1'].parts['Wing'].edges[2609], 
            mdb.models['Model-1'].parts['Wing'].edges[2610], 
            mdb.models['Model-1'].parts['Wing'].edges[2611], 
            mdb.models['Model-1'].parts['Wing'].edges[2612], 
            mdb.models['Model-1'].parts['Wing'].edges[2613])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2524], 
            mdb.models['Model-1'].parts['Wing'].edges[2525], 
            mdb.models['Model-1'].parts['Wing'].edges[2548], 
            mdb.models['Model-1'].parts['Wing'].edges[2549], 
            mdb.models['Model-1'].parts['Wing'].edges[2550], 
            mdb.models['Model-1'].parts['Wing'].edges[2551], 
            mdb.models['Model-1'].parts['Wing'].edges[2552], 
            mdb.models['Model-1'].parts['Wing'].edges[2553], 
            mdb.models['Model-1'].parts['Wing'].edges[2554], 
            mdb.models['Model-1'].parts['Wing'].edges[2555], 
            mdb.models['Model-1'].parts['Wing'].edges[2556], 
            mdb.models['Model-1'].parts['Wing'].edges[2557], 
            mdb.models['Model-1'].parts['Wing'].edges[2558], 
            mdb.models['Model-1'].parts['Wing'].edges[2559], 
            mdb.models['Model-1'].parts['Wing'].edges[2560], 
            mdb.models['Model-1'].parts['Wing'].edges[2561], 
            mdb.models['Model-1'].parts['Wing'].edges[2562], 
            mdb.models['Model-1'].parts['Wing'].edges[2563], 
            mdb.models['Model-1'].parts['Wing'].edges[2564], 
            mdb.models['Model-1'].parts['Wing'].edges[2565], 
            mdb.models['Model-1'].parts['Wing'].edges[2584], 
            mdb.models['Model-1'].parts['Wing'].edges[2585], 
            mdb.models['Model-1'].parts['Wing'].edges[2586], 
            mdb.models['Model-1'].parts['Wing'].edges[2587], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2603], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605], 
            mdb.models['Model-1'].parts['Wing'].edges[2606], 
            mdb.models['Model-1'].parts['Wing'].edges[2607], 
            mdb.models['Model-1'].parts['Wing'].edges[2608], 
            mdb.models['Model-1'].parts['Wing'].edges[2609], 
            mdb.models['Model-1'].parts['Wing'].edges[2610], 
            mdb.models['Model-1'].parts['Wing'].edges[2611], 
            mdb.models['Model-1'].parts['Wing'].edges[2612], 
            mdb.models['Model-1'].parts['Wing'].edges[2613], 
            mdb.models['Model-1'].parts['Wing'].edges[2614], 
            mdb.models['Model-1'].parts['Wing'].edges[2615], 
            mdb.models['Model-1'].parts['Wing'].edges[2616], 
            mdb.models['Model-1'].parts['Wing'].edges[2617])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2528], 
            mdb.models['Model-1'].parts['Wing'].edges[2529], 
            mdb.models['Model-1'].parts['Wing'].edges[2552], 
            mdb.models['Model-1'].parts['Wing'].edges[2553], 
            mdb.models['Model-1'].parts['Wing'].edges[2554], 
            mdb.models['Model-1'].parts['Wing'].edges[2555], 
            mdb.models['Model-1'].parts['Wing'].edges[2556], 
            mdb.models['Model-1'].parts['Wing'].edges[2557], 
            mdb.models['Model-1'].parts['Wing'].edges[2558], 
            mdb.models['Model-1'].parts['Wing'].edges[2559], 
            mdb.models['Model-1'].parts['Wing'].edges[2560], 
            mdb.models['Model-1'].parts['Wing'].edges[2561], 
            mdb.models['Model-1'].parts['Wing'].edges[2562], 
            mdb.models['Model-1'].parts['Wing'].edges[2563], 
            mdb.models['Model-1'].parts['Wing'].edges[2564], 
            mdb.models['Model-1'].parts['Wing'].edges[2565], 
            mdb.models['Model-1'].parts['Wing'].edges[2566], 
            mdb.models['Model-1'].parts['Wing'].edges[2567], 
            mdb.models['Model-1'].parts['Wing'].edges[2568], 
            mdb.models['Model-1'].parts['Wing'].edges[2569], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605], 
            mdb.models['Model-1'].parts['Wing'].edges[2606], 
            mdb.models['Model-1'].parts['Wing'].edges[2607], 
            mdb.models['Model-1'].parts['Wing'].edges[2608], 
            mdb.models['Model-1'].parts['Wing'].edges[2609], 
            mdb.models['Model-1'].parts['Wing'].edges[2610], 
            mdb.models['Model-1'].parts['Wing'].edges[2611], 
            mdb.models['Model-1'].parts['Wing'].edges[2612], 
            mdb.models['Model-1'].parts['Wing'].edges[2613], 
            mdb.models['Model-1'].parts['Wing'].edges[2614], 
            mdb.models['Model-1'].parts['Wing'].edges[2615], 
            mdb.models['Model-1'].parts['Wing'].edges[2616], 
            mdb.models['Model-1'].parts['Wing'].edges[2617], 
            mdb.models['Model-1'].parts['Wing'].edges[2618], 
            mdb.models['Model-1'].parts['Wing'].edges[2619], 
            mdb.models['Model-1'].parts['Wing'].edges[2620], 
            mdb.models['Model-1'].parts['Wing'].edges[2621])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2532], 
            mdb.models['Model-1'].parts['Wing'].edges[2533], 
            mdb.models['Model-1'].parts['Wing'].edges[2556], 
            mdb.models['Model-1'].parts['Wing'].edges[2557], 
            mdb.models['Model-1'].parts['Wing'].edges[2558], 
            mdb.models['Model-1'].parts['Wing'].edges[2559], 
            mdb.models['Model-1'].parts['Wing'].edges[2560], 
            mdb.models['Model-1'].parts['Wing'].edges[2561], 
            mdb.models['Model-1'].parts['Wing'].edges[2562], 
            mdb.models['Model-1'].parts['Wing'].edges[2563], 
            mdb.models['Model-1'].parts['Wing'].edges[2564], 
            mdb.models['Model-1'].parts['Wing'].edges[2565], 
            mdb.models['Model-1'].parts['Wing'].edges[2566], 
            mdb.models['Model-1'].parts['Wing'].edges[2567], 
            mdb.models['Model-1'].parts['Wing'].edges[2568], 
            mdb.models['Model-1'].parts['Wing'].edges[2569], 
            mdb.models['Model-1'].parts['Wing'].edges[2570], 
            mdb.models['Model-1'].parts['Wing'].edges[2571], 
            mdb.models['Model-1'].parts['Wing'].edges[2572], 
            mdb.models['Model-1'].parts['Wing'].edges[2573], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2601], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2603], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605], 
            mdb.models['Model-1'].parts['Wing'].edges[2606], 
            mdb.models['Model-1'].parts['Wing'].edges[2608], 
            mdb.models['Model-1'].parts['Wing'].edges[2609], 
            mdb.models['Model-1'].parts['Wing'].edges[2610], 
            mdb.models['Model-1'].parts['Wing'].edges[2611], 
            mdb.models['Model-1'].parts['Wing'].edges[2612], 
            mdb.models['Model-1'].parts['Wing'].edges[2613], 
            mdb.models['Model-1'].parts['Wing'].edges[2614], 
            mdb.models['Model-1'].parts['Wing'].edges[2615], 
            mdb.models['Model-1'].parts['Wing'].edges[2616], 
            mdb.models['Model-1'].parts['Wing'].edges[2617], 
            mdb.models['Model-1'].parts['Wing'].edges[2618], 
            mdb.models['Model-1'].parts['Wing'].edges[2619], 
            mdb.models['Model-1'].parts['Wing'].edges[2620], 
            mdb.models['Model-1'].parts['Wing'].edges[2621], 
            mdb.models['Model-1'].parts['Wing'].edges[2622], 
            mdb.models['Model-1'].parts['Wing'].edges[2623], 
            mdb.models['Model-1'].parts['Wing'].edges[2624], 
            mdb.models['Model-1'].parts['Wing'].edges[2625])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2527], 
            mdb.models['Model-1'].parts['Wing'].edges[2529], 
            mdb.models['Model-1'].parts['Wing'].edges[2552], 
            mdb.models['Model-1'].parts['Wing'].edges[2553], 
            mdb.models['Model-1'].parts['Wing'].edges[2554], 
            mdb.models['Model-1'].parts['Wing'].edges[2555], 
            mdb.models['Model-1'].parts['Wing'].edges[2556], 
            mdb.models['Model-1'].parts['Wing'].edges[2557], 
            mdb.models['Model-1'].parts['Wing'].edges[2558], 
            mdb.models['Model-1'].parts['Wing'].edges[2559], 
            mdb.models['Model-1'].parts['Wing'].edges[2560], 
            mdb.models['Model-1'].parts['Wing'].edges[2561], 
            mdb.models['Model-1'].parts['Wing'].edges[2562], 
            mdb.models['Model-1'].parts['Wing'].edges[2563], 
            mdb.models['Model-1'].parts['Wing'].edges[2564], 
            mdb.models['Model-1'].parts['Wing'].edges[2565], 
            mdb.models['Model-1'].parts['Wing'].edges[2566], 
            mdb.models['Model-1'].parts['Wing'].edges[2567], 
            mdb.models['Model-1'].parts['Wing'].edges[2568], 
            mdb.models['Model-1'].parts['Wing'].edges[2569], 
            mdb.models['Model-1'].parts['Wing'].edges[2588], 
            mdb.models['Model-1'].parts['Wing'].edges[2589], 
            mdb.models['Model-1'].parts['Wing'].edges[2590], 
            mdb.models['Model-1'].parts['Wing'].edges[2591], 
            mdb.models['Model-1'].parts['Wing'].edges[2592], 
            mdb.models['Model-1'].parts['Wing'].edges[2593], 
            mdb.models['Model-1'].parts['Wing'].edges[2594], 
            mdb.models['Model-1'].parts['Wing'].edges[2595], 
            mdb.models['Model-1'].parts['Wing'].edges[2596], 
            mdb.models['Model-1'].parts['Wing'].edges[2597], 
            mdb.models['Model-1'].parts['Wing'].edges[2598], 
            mdb.models['Model-1'].parts['Wing'].edges[2599], 
            mdb.models['Model-1'].parts['Wing'].edges[2600], 
            mdb.models['Model-1'].parts['Wing'].edges[2602], 
            mdb.models['Model-1'].parts['Wing'].edges[2603], 
            mdb.models['Model-1'].parts['Wing'].edges[2604], 
            mdb.models['Model-1'].parts['Wing'].edges[2605], 
            mdb.models['Model-1'].parts['Wing'].edges[2606], 
            mdb.models['Model-1'].parts['Wing'].edges[2607], 
            mdb.models['Model-1'].parts['Wing'].edges[2608], 
            mdb.models['Model-1'].parts['Wing'].edges[2609], 
            mdb.models['Model-1'].parts['Wing'].edges[2610], 
            mdb.models['Model-1'].parts['Wing'].edges[2611], 
            mdb.models['Model-1'].parts['Wing'].edges[2612], 
            mdb.models['Model-1'].parts['Wing'].edges[2613], 
            mdb.models['Model-1'].parts['Wing'].edges[2614], 
            mdb.models['Model-1'].parts['Wing'].edges[2615], 
            mdb.models['Model-1'].parts['Wing'].edges[2616], 
            mdb.models['Model-1'].parts['Wing'].edges[2617], 
            mdb.models['Model-1'].parts['Wing'].edges[2618], 
            mdb.models['Model-1'].parts['Wing'].edges[2619], 
            mdb.models['Model-1'].parts['Wing'].edges[2620], 
            mdb.models['Model-1'].parts['Wing'].edges[2621])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[4], 
            mdb.models['Model-1'].parts['Wing'].edges[7], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[55], 
            mdb.models['Model-1'].parts['Wing'].edges[58], 
            mdb.models['Model-1'].parts['Wing'].edges[61], 
            mdb.models['Model-1'].parts['Wing'].edges[64], 
            mdb.models['Model-1'].parts['Wing'].edges[67], 
            mdb.models['Model-1'].parts['Wing'].edges[70], 
            mdb.models['Model-1'].parts['Wing'].edges[73], 
            mdb.models['Model-1'].parts['Wing'].edges[76], 
            mdb.models['Model-1'].parts['Wing'].edges[79], 
            mdb.models['Model-1'].parts['Wing'].edges[82], 
            mdb.models['Model-1'].parts['Wing'].edges[85], 
            mdb.models['Model-1'].parts['Wing'].edges[88], 
            mdb.models['Model-1'].parts['Wing'].edges[91], 
            mdb.models['Model-1'].parts['Wing'].edges[94], 
            mdb.models['Model-1'].parts['Wing'].edges[97], 
            mdb.models['Model-1'].parts['Wing'].edges[100], 
            mdb.models['Model-1'].parts['Wing'].edges[103], 
            mdb.models['Model-1'].parts['Wing'].edges[106], 
            mdb.models['Model-1'].parts['Wing'].edges[109], 
            mdb.models['Model-1'].parts['Wing'].edges[112], 
            mdb.models['Model-1'].parts['Wing'].edges[115], 
            mdb.models['Model-1'].parts['Wing'].edges[117], 
            mdb.models['Model-1'].parts['Wing'].edges[120], 
            mdb.models['Model-1'].parts['Wing'].edges[123], 
            mdb.models['Model-1'].parts['Wing'].edges[126], 
            mdb.models['Model-1'].parts['Wing'].edges[129], 
            mdb.models['Model-1'].parts['Wing'].edges[132], 
            mdb.models['Model-1'].parts['Wing'].edges[135], 
            mdb.models['Model-1'].parts['Wing'].edges[138], 
            mdb.models['Model-1'].parts['Wing'].edges[141], 
            mdb.models['Model-1'].parts['Wing'].edges[144], 
            mdb.models['Model-1'].parts['Wing'].edges[147], 
            mdb.models['Model-1'].parts['Wing'].edges[150], 
            mdb.models['Model-1'].parts['Wing'].edges[153], 
            mdb.models['Model-1'].parts['Wing'].edges[156], 
            mdb.models['Model-1'].parts['Wing'].edges[159], 
            mdb.models['Model-1'].parts['Wing'].edges[162]), (
            mdb.models['Model-1'].parts['Wing'].edges[3029], 
            mdb.models['Model-1'].parts['Wing'].edges[3031], 
            mdb.models['Model-1'].parts['Wing'].edges[3054], 
            mdb.models['Model-1'].parts['Wing'].edges[3055], 
            mdb.models['Model-1'].parts['Wing'].edges[3056], 
            mdb.models['Model-1'].parts['Wing'].edges[3057], 
            mdb.models['Model-1'].parts['Wing'].edges[3058], 
            mdb.models['Model-1'].parts['Wing'].edges[3059], 
            mdb.models['Model-1'].parts['Wing'].edges[3060], 
            mdb.models['Model-1'].parts['Wing'].edges[3061], 
            mdb.models['Model-1'].parts['Wing'].edges[3062], 
            mdb.models['Model-1'].parts['Wing'].edges[3063], 
            mdb.models['Model-1'].parts['Wing'].edges[3064], 
            mdb.models['Model-1'].parts['Wing'].edges[3065], 
            mdb.models['Model-1'].parts['Wing'].edges[3066], 
            mdb.models['Model-1'].parts['Wing'].edges[3067], 
            mdb.models['Model-1'].parts['Wing'].edges[3068], 
            mdb.models['Model-1'].parts['Wing'].edges[3069], 
            mdb.models['Model-1'].parts['Wing'].edges[3070], 
            mdb.models['Model-1'].parts['Wing'].edges[3071], 
            mdb.models['Model-1'].parts['Wing'].edges[3090], 
            mdb.models['Model-1'].parts['Wing'].edges[3091], 
            mdb.models['Model-1'].parts['Wing'].edges[3092], 
            mdb.models['Model-1'].parts['Wing'].edges[3093], 
            mdb.models['Model-1'].parts['Wing'].edges[3094], 
            mdb.models['Model-1'].parts['Wing'].edges[3095], 
            mdb.models['Model-1'].parts['Wing'].edges[3096], 
            mdb.models['Model-1'].parts['Wing'].edges[3097], 
            mdb.models['Model-1'].parts['Wing'].edges[3098], 
            mdb.models['Model-1'].parts['Wing'].edges[3099], 
            mdb.models['Model-1'].parts['Wing'].edges[3100], 
            mdb.models['Model-1'].parts['Wing'].edges[3101], 
            mdb.models['Model-1'].parts['Wing'].edges[3102], 
            mdb.models['Model-1'].parts['Wing'].edges[3104], 
            mdb.models['Model-1'].parts['Wing'].edges[3105], 
            mdb.models['Model-1'].parts['Wing'].edges[3106], 
            mdb.models['Model-1'].parts['Wing'].edges[3107], 
            mdb.models['Model-1'].parts['Wing'].edges[3108], 
            mdb.models['Model-1'].parts['Wing'].edges[3109], 
            mdb.models['Model-1'].parts['Wing'].edges[3110], 
            mdb.models['Model-1'].parts['Wing'].edges[3111], 
            mdb.models['Model-1'].parts['Wing'].edges[3112], 
            mdb.models['Model-1'].parts['Wing'].edges[3113], 
            mdb.models['Model-1'].parts['Wing'].edges[3114], 
            mdb.models['Model-1'].parts['Wing'].edges[3115], 
            mdb.models['Model-1'].parts['Wing'].edges[3116], 
            mdb.models['Model-1'].parts['Wing'].edges[3117], 
            mdb.models['Model-1'].parts['Wing'].edges[3118], 
            mdb.models['Model-1'].parts['Wing'].edges[3119], 
            mdb.models['Model-1'].parts['Wing'].edges[3120], 
            mdb.models['Model-1'].parts['Wing'].edges[3121], 
            mdb.models['Model-1'].parts['Wing'].edges[3122], 
            mdb.models['Model-1'].parts['Wing'].edges[3123])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2634], 
            mdb.models['Model-1'].parts['Wing'].edges[2636], 
            mdb.models['Model-1'].parts['Wing'].edges[2659], 
            mdb.models['Model-1'].parts['Wing'].edges[2660], 
            mdb.models['Model-1'].parts['Wing'].edges[2661], 
            mdb.models['Model-1'].parts['Wing'].edges[2662], 
            mdb.models['Model-1'].parts['Wing'].edges[2663], 
            mdb.models['Model-1'].parts['Wing'].edges[2664], 
            mdb.models['Model-1'].parts['Wing'].edges[2665], 
            mdb.models['Model-1'].parts['Wing'].edges[2666], 
            mdb.models['Model-1'].parts['Wing'].edges[2667], 
            mdb.models['Model-1'].parts['Wing'].edges[2668], 
            mdb.models['Model-1'].parts['Wing'].edges[2669], 
            mdb.models['Model-1'].parts['Wing'].edges[2670], 
            mdb.models['Model-1'].parts['Wing'].edges[2671], 
            mdb.models['Model-1'].parts['Wing'].edges[2672], 
            mdb.models['Model-1'].parts['Wing'].edges[2673], 
            mdb.models['Model-1'].parts['Wing'].edges[2674], 
            mdb.models['Model-1'].parts['Wing'].edges[2675], 
            mdb.models['Model-1'].parts['Wing'].edges[2676], 
            mdb.models['Model-1'].parts['Wing'].edges[2695], 
            mdb.models['Model-1'].parts['Wing'].edges[2696], 
            mdb.models['Model-1'].parts['Wing'].edges[2697], 
            mdb.models['Model-1'].parts['Wing'].edges[2698], 
            mdb.models['Model-1'].parts['Wing'].edges[2699], 
            mdb.models['Model-1'].parts['Wing'].edges[2700], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2702], 
            mdb.models['Model-1'].parts['Wing'].edges[2703], 
            mdb.models['Model-1'].parts['Wing'].edges[2704], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2707], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2713], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2639], 
            mdb.models['Model-1'].parts['Wing'].edges[2641], 
            mdb.models['Model-1'].parts['Wing'].edges[2664], 
            mdb.models['Model-1'].parts['Wing'].edges[2665], 
            mdb.models['Model-1'].parts['Wing'].edges[2666], 
            mdb.models['Model-1'].parts['Wing'].edges[2667], 
            mdb.models['Model-1'].parts['Wing'].edges[2668], 
            mdb.models['Model-1'].parts['Wing'].edges[2669], 
            mdb.models['Model-1'].parts['Wing'].edges[2670], 
            mdb.models['Model-1'].parts['Wing'].edges[2671], 
            mdb.models['Model-1'].parts['Wing'].edges[2672], 
            mdb.models['Model-1'].parts['Wing'].edges[2673], 
            mdb.models['Model-1'].parts['Wing'].edges[2674], 
            mdb.models['Model-1'].parts['Wing'].edges[2675], 
            mdb.models['Model-1'].parts['Wing'].edges[2676], 
            mdb.models['Model-1'].parts['Wing'].edges[2677], 
            mdb.models['Model-1'].parts['Wing'].edges[2678], 
            mdb.models['Model-1'].parts['Wing'].edges[2679], 
            mdb.models['Model-1'].parts['Wing'].edges[2680], 
            mdb.models['Model-1'].parts['Wing'].edges[2681], 
            mdb.models['Model-1'].parts['Wing'].edges[2700], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2702], 
            mdb.models['Model-1'].parts['Wing'].edges[2703], 
            mdb.models['Model-1'].parts['Wing'].edges[2704], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2707], 
            mdb.models['Model-1'].parts['Wing'].edges[2708], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2644], 
            mdb.models['Model-1'].parts['Wing'].edges[2646], 
            mdb.models['Model-1'].parts['Wing'].edges[2669], 
            mdb.models['Model-1'].parts['Wing'].edges[2670], 
            mdb.models['Model-1'].parts['Wing'].edges[2671], 
            mdb.models['Model-1'].parts['Wing'].edges[2672], 
            mdb.models['Model-1'].parts['Wing'].edges[2673], 
            mdb.models['Model-1'].parts['Wing'].edges[2674], 
            mdb.models['Model-1'].parts['Wing'].edges[2675], 
            mdb.models['Model-1'].parts['Wing'].edges[2676], 
            mdb.models['Model-1'].parts['Wing'].edges[2677], 
            mdb.models['Model-1'].parts['Wing'].edges[2678], 
            mdb.models['Model-1'].parts['Wing'].edges[2679], 
            mdb.models['Model-1'].parts['Wing'].edges[2680], 
            mdb.models['Model-1'].parts['Wing'].edges[2681], 
            mdb.models['Model-1'].parts['Wing'].edges[2682], 
            mdb.models['Model-1'].parts['Wing'].edges[2683], 
            mdb.models['Model-1'].parts['Wing'].edges[2684], 
            mdb.models['Model-1'].parts['Wing'].edges[2685], 
            mdb.models['Model-1'].parts['Wing'].edges[2686], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2707], 
            mdb.models['Model-1'].parts['Wing'].edges[2708], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2713], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2649], 
            mdb.models['Model-1'].parts['Wing'].edges[2651], 
            mdb.models['Model-1'].parts['Wing'].edges[2674], 
            mdb.models['Model-1'].parts['Wing'].edges[2675], 
            mdb.models['Model-1'].parts['Wing'].edges[2676], 
            mdb.models['Model-1'].parts['Wing'].edges[2677], 
            mdb.models['Model-1'].parts['Wing'].edges[2678], 
            mdb.models['Model-1'].parts['Wing'].edges[2679], 
            mdb.models['Model-1'].parts['Wing'].edges[2680], 
            mdb.models['Model-1'].parts['Wing'].edges[2681], 
            mdb.models['Model-1'].parts['Wing'].edges[2682], 
            mdb.models['Model-1'].parts['Wing'].edges[2683], 
            mdb.models['Model-1'].parts['Wing'].edges[2684], 
            mdb.models['Model-1'].parts['Wing'].edges[2685], 
            mdb.models['Model-1'].parts['Wing'].edges[2686], 
            mdb.models['Model-1'].parts['Wing'].edges[2687], 
            mdb.models['Model-1'].parts['Wing'].edges[2688], 
            mdb.models['Model-1'].parts['Wing'].edges[2689], 
            mdb.models['Model-1'].parts['Wing'].edges[2690], 
            mdb.models['Model-1'].parts['Wing'].edges[2691], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2713], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2743])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2654], 
            mdb.models['Model-1'].parts['Wing'].edges[2656], 
            mdb.models['Model-1'].parts['Wing'].edges[2679], 
            mdb.models['Model-1'].parts['Wing'].edges[2680], 
            mdb.models['Model-1'].parts['Wing'].edges[2681], 
            mdb.models['Model-1'].parts['Wing'].edges[2682], 
            mdb.models['Model-1'].parts['Wing'].edges[2683], 
            mdb.models['Model-1'].parts['Wing'].edges[2684], 
            mdb.models['Model-1'].parts['Wing'].edges[2685], 
            mdb.models['Model-1'].parts['Wing'].edges[2686], 
            mdb.models['Model-1'].parts['Wing'].edges[2687], 
            mdb.models['Model-1'].parts['Wing'].edges[2688], 
            mdb.models['Model-1'].parts['Wing'].edges[2689], 
            mdb.models['Model-1'].parts['Wing'].edges[2690], 
            mdb.models['Model-1'].parts['Wing'].edges[2691], 
            mdb.models['Model-1'].parts['Wing'].edges[2692], 
            mdb.models['Model-1'].parts['Wing'].edges[2693], 
            mdb.models['Model-1'].parts['Wing'].edges[2694], 
            mdb.models['Model-1'].parts['Wing'].edges[2695], 
            mdb.models['Model-1'].parts['Wing'].edges[2696], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2743], 
            mdb.models['Model-1'].parts['Wing'].edges[2744], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2748])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2659], 
            mdb.models['Model-1'].parts['Wing'].edges[2661], 
            mdb.models['Model-1'].parts['Wing'].edges[2684], 
            mdb.models['Model-1'].parts['Wing'].edges[2685], 
            mdb.models['Model-1'].parts['Wing'].edges[2686], 
            mdb.models['Model-1'].parts['Wing'].edges[2687], 
            mdb.models['Model-1'].parts['Wing'].edges[2688], 
            mdb.models['Model-1'].parts['Wing'].edges[2689], 
            mdb.models['Model-1'].parts['Wing'].edges[2690], 
            mdb.models['Model-1'].parts['Wing'].edges[2691], 
            mdb.models['Model-1'].parts['Wing'].edges[2692], 
            mdb.models['Model-1'].parts['Wing'].edges[2693], 
            mdb.models['Model-1'].parts['Wing'].edges[2694], 
            mdb.models['Model-1'].parts['Wing'].edges[2695], 
            mdb.models['Model-1'].parts['Wing'].edges[2696], 
            mdb.models['Model-1'].parts['Wing'].edges[2697], 
            mdb.models['Model-1'].parts['Wing'].edges[2698], 
            mdb.models['Model-1'].parts['Wing'].edges[2699], 
            mdb.models['Model-1'].parts['Wing'].edges[2700], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2743], 
            mdb.models['Model-1'].parts['Wing'].edges[2744], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2748], 
            mdb.models['Model-1'].parts['Wing'].edges[2749], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2753])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2664], 
            mdb.models['Model-1'].parts['Wing'].edges[2666], 
            mdb.models['Model-1'].parts['Wing'].edges[2689], 
            mdb.models['Model-1'].parts['Wing'].edges[2690], 
            mdb.models['Model-1'].parts['Wing'].edges[2691], 
            mdb.models['Model-1'].parts['Wing'].edges[2692], 
            mdb.models['Model-1'].parts['Wing'].edges[2693], 
            mdb.models['Model-1'].parts['Wing'].edges[2694], 
            mdb.models['Model-1'].parts['Wing'].edges[2695], 
            mdb.models['Model-1'].parts['Wing'].edges[2696], 
            mdb.models['Model-1'].parts['Wing'].edges[2697], 
            mdb.models['Model-1'].parts['Wing'].edges[2698], 
            mdb.models['Model-1'].parts['Wing'].edges[2699], 
            mdb.models['Model-1'].parts['Wing'].edges[2700], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2702], 
            mdb.models['Model-1'].parts['Wing'].edges[2703], 
            mdb.models['Model-1'].parts['Wing'].edges[2704], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2743], 
            mdb.models['Model-1'].parts['Wing'].edges[2744], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2748], 
            mdb.models['Model-1'].parts['Wing'].edges[2749], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2753], 
            mdb.models['Model-1'].parts['Wing'].edges[2754], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2758])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2669], 
            mdb.models['Model-1'].parts['Wing'].edges[2671], 
            mdb.models['Model-1'].parts['Wing'].edges[2694], 
            mdb.models['Model-1'].parts['Wing'].edges[2695], 
            mdb.models['Model-1'].parts['Wing'].edges[2696], 
            mdb.models['Model-1'].parts['Wing'].edges[2697], 
            mdb.models['Model-1'].parts['Wing'].edges[2698], 
            mdb.models['Model-1'].parts['Wing'].edges[2699], 
            mdb.models['Model-1'].parts['Wing'].edges[2700], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2702], 
            mdb.models['Model-1'].parts['Wing'].edges[2703], 
            mdb.models['Model-1'].parts['Wing'].edges[2704], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2707], 
            mdb.models['Model-1'].parts['Wing'].edges[2708], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2744], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2748], 
            mdb.models['Model-1'].parts['Wing'].edges[2749], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2753], 
            mdb.models['Model-1'].parts['Wing'].edges[2754], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2758], 
            mdb.models['Model-1'].parts['Wing'].edges[2759], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2763])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2674], 
            mdb.models['Model-1'].parts['Wing'].edges[2676], 
            mdb.models['Model-1'].parts['Wing'].edges[2699], 
            mdb.models['Model-1'].parts['Wing'].edges[2700], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2702], 
            mdb.models['Model-1'].parts['Wing'].edges[2703], 
            mdb.models['Model-1'].parts['Wing'].edges[2704], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2707], 
            mdb.models['Model-1'].parts['Wing'].edges[2708], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2713], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2743], 
            mdb.models['Model-1'].parts['Wing'].edges[2744], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2749], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2753], 
            mdb.models['Model-1'].parts['Wing'].edges[2754], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2758], 
            mdb.models['Model-1'].parts['Wing'].edges[2759], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2763], 
            mdb.models['Model-1'].parts['Wing'].edges[2764], 
            mdb.models['Model-1'].parts['Wing'].edges[2765], 
            mdb.models['Model-1'].parts['Wing'].edges[2766], 
            mdb.models['Model-1'].parts['Wing'].edges[2767], 
            mdb.models['Model-1'].parts['Wing'].edges[2768])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2679], 
            mdb.models['Model-1'].parts['Wing'].edges[2681], 
            mdb.models['Model-1'].parts['Wing'].edges[2704], 
            mdb.models['Model-1'].parts['Wing'].edges[2705], 
            mdb.models['Model-1'].parts['Wing'].edges[2706], 
            mdb.models['Model-1'].parts['Wing'].edges[2707], 
            mdb.models['Model-1'].parts['Wing'].edges[2708], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2713], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2742], 
            mdb.models['Model-1'].parts['Wing'].edges[2743], 
            mdb.models['Model-1'].parts['Wing'].edges[2744], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2748], 
            mdb.models['Model-1'].parts['Wing'].edges[2749], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2754], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2758], 
            mdb.models['Model-1'].parts['Wing'].edges[2759], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2763], 
            mdb.models['Model-1'].parts['Wing'].edges[2764], 
            mdb.models['Model-1'].parts['Wing'].edges[2765], 
            mdb.models['Model-1'].parts['Wing'].edges[2766], 
            mdb.models['Model-1'].parts['Wing'].edges[2767], 
            mdb.models['Model-1'].parts['Wing'].edges[2768], 
            mdb.models['Model-1'].parts['Wing'].edges[2769], 
            mdb.models['Model-1'].parts['Wing'].edges[2770], 
            mdb.models['Model-1'].parts['Wing'].edges[2771], 
            mdb.models['Model-1'].parts['Wing'].edges[2772], 
            mdb.models['Model-1'].parts['Wing'].edges[2773])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2684], 
            mdb.models['Model-1'].parts['Wing'].edges[2686], 
            mdb.models['Model-1'].parts['Wing'].edges[2709], 
            mdb.models['Model-1'].parts['Wing'].edges[2710], 
            mdb.models['Model-1'].parts['Wing'].edges[2711], 
            mdb.models['Model-1'].parts['Wing'].edges[2712], 
            mdb.models['Model-1'].parts['Wing'].edges[2713], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2745], 
            mdb.models['Model-1'].parts['Wing'].edges[2746], 
            mdb.models['Model-1'].parts['Wing'].edges[2747], 
            mdb.models['Model-1'].parts['Wing'].edges[2748], 
            mdb.models['Model-1'].parts['Wing'].edges[2749], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2753], 
            mdb.models['Model-1'].parts['Wing'].edges[2754], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2759], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2763], 
            mdb.models['Model-1'].parts['Wing'].edges[2764], 
            mdb.models['Model-1'].parts['Wing'].edges[2765], 
            mdb.models['Model-1'].parts['Wing'].edges[2766], 
            mdb.models['Model-1'].parts['Wing'].edges[2767], 
            mdb.models['Model-1'].parts['Wing'].edges[2768], 
            mdb.models['Model-1'].parts['Wing'].edges[2769], 
            mdb.models['Model-1'].parts['Wing'].edges[2770], 
            mdb.models['Model-1'].parts['Wing'].edges[2771], 
            mdb.models['Model-1'].parts['Wing'].edges[2772], 
            mdb.models['Model-1'].parts['Wing'].edges[2773], 
            mdb.models['Model-1'].parts['Wing'].edges[2774], 
            mdb.models['Model-1'].parts['Wing'].edges[2775], 
            mdb.models['Model-1'].parts['Wing'].edges[2776], 
            mdb.models['Model-1'].parts['Wing'].edges[2777], 
            mdb.models['Model-1'].parts['Wing'].edges[2778])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2689], 
            mdb.models['Model-1'].parts['Wing'].edges[2691], 
            mdb.models['Model-1'].parts['Wing'].edges[2714], 
            mdb.models['Model-1'].parts['Wing'].edges[2715], 
            mdb.models['Model-1'].parts['Wing'].edges[2716], 
            mdb.models['Model-1'].parts['Wing'].edges[2717], 
            mdb.models['Model-1'].parts['Wing'].edges[2718], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2750], 
            mdb.models['Model-1'].parts['Wing'].edges[2751], 
            mdb.models['Model-1'].parts['Wing'].edges[2752], 
            mdb.models['Model-1'].parts['Wing'].edges[2753], 
            mdb.models['Model-1'].parts['Wing'].edges[2754], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2758], 
            mdb.models['Model-1'].parts['Wing'].edges[2759], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2764], 
            mdb.models['Model-1'].parts['Wing'].edges[2765], 
            mdb.models['Model-1'].parts['Wing'].edges[2766], 
            mdb.models['Model-1'].parts['Wing'].edges[2767], 
            mdb.models['Model-1'].parts['Wing'].edges[2768], 
            mdb.models['Model-1'].parts['Wing'].edges[2769], 
            mdb.models['Model-1'].parts['Wing'].edges[2770], 
            mdb.models['Model-1'].parts['Wing'].edges[2771], 
            mdb.models['Model-1'].parts['Wing'].edges[2772], 
            mdb.models['Model-1'].parts['Wing'].edges[2773], 
            mdb.models['Model-1'].parts['Wing'].edges[2774], 
            mdb.models['Model-1'].parts['Wing'].edges[2775], 
            mdb.models['Model-1'].parts['Wing'].edges[2776], 
            mdb.models['Model-1'].parts['Wing'].edges[2777], 
            mdb.models['Model-1'].parts['Wing'].edges[2778], 
            mdb.models['Model-1'].parts['Wing'].edges[2779], 
            mdb.models['Model-1'].parts['Wing'].edges[2780], 
            mdb.models['Model-1'].parts['Wing'].edges[2781], 
            mdb.models['Model-1'].parts['Wing'].edges[2782], 
            mdb.models['Model-1'].parts['Wing'].edges[2783])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2694], 
            mdb.models['Model-1'].parts['Wing'].edges[2696], 
            mdb.models['Model-1'].parts['Wing'].edges[2719], 
            mdb.models['Model-1'].parts['Wing'].edges[2720], 
            mdb.models['Model-1'].parts['Wing'].edges[2721], 
            mdb.models['Model-1'].parts['Wing'].edges[2722], 
            mdb.models['Model-1'].parts['Wing'].edges[2723], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2755], 
            mdb.models['Model-1'].parts['Wing'].edges[2756], 
            mdb.models['Model-1'].parts['Wing'].edges[2757], 
            mdb.models['Model-1'].parts['Wing'].edges[2758], 
            mdb.models['Model-1'].parts['Wing'].edges[2759], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2763], 
            mdb.models['Model-1'].parts['Wing'].edges[2764], 
            mdb.models['Model-1'].parts['Wing'].edges[2765], 
            mdb.models['Model-1'].parts['Wing'].edges[2766], 
            mdb.models['Model-1'].parts['Wing'].edges[2767], 
            mdb.models['Model-1'].parts['Wing'].edges[2769], 
            mdb.models['Model-1'].parts['Wing'].edges[2770], 
            mdb.models['Model-1'].parts['Wing'].edges[2771], 
            mdb.models['Model-1'].parts['Wing'].edges[2772], 
            mdb.models['Model-1'].parts['Wing'].edges[2773], 
            mdb.models['Model-1'].parts['Wing'].edges[2774], 
            mdb.models['Model-1'].parts['Wing'].edges[2775], 
            mdb.models['Model-1'].parts['Wing'].edges[2776], 
            mdb.models['Model-1'].parts['Wing'].edges[2777], 
            mdb.models['Model-1'].parts['Wing'].edges[2778], 
            mdb.models['Model-1'].parts['Wing'].edges[2779], 
            mdb.models['Model-1'].parts['Wing'].edges[2780], 
            mdb.models['Model-1'].parts['Wing'].edges[2781], 
            mdb.models['Model-1'].parts['Wing'].edges[2782], 
            mdb.models['Model-1'].parts['Wing'].edges[2783], 
            mdb.models['Model-1'].parts['Wing'].edges[2784], 
            mdb.models['Model-1'].parts['Wing'].edges[2785], 
            mdb.models['Model-1'].parts['Wing'].edges[2786], 
            mdb.models['Model-1'].parts['Wing'].edges[2787], 
            mdb.models['Model-1'].parts['Wing'].edges[2788])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2699], 
            mdb.models['Model-1'].parts['Wing'].edges[2701], 
            mdb.models['Model-1'].parts['Wing'].edges[2724], 
            mdb.models['Model-1'].parts['Wing'].edges[2725], 
            mdb.models['Model-1'].parts['Wing'].edges[2726], 
            mdb.models['Model-1'].parts['Wing'].edges[2727], 
            mdb.models['Model-1'].parts['Wing'].edges[2728], 
            mdb.models['Model-1'].parts['Wing'].edges[2729], 
            mdb.models['Model-1'].parts['Wing'].edges[2730], 
            mdb.models['Model-1'].parts['Wing'].edges[2731], 
            mdb.models['Model-1'].parts['Wing'].edges[2732], 
            mdb.models['Model-1'].parts['Wing'].edges[2733], 
            mdb.models['Model-1'].parts['Wing'].edges[2734], 
            mdb.models['Model-1'].parts['Wing'].edges[2735], 
            mdb.models['Model-1'].parts['Wing'].edges[2736], 
            mdb.models['Model-1'].parts['Wing'].edges[2737], 
            mdb.models['Model-1'].parts['Wing'].edges[2738], 
            mdb.models['Model-1'].parts['Wing'].edges[2739], 
            mdb.models['Model-1'].parts['Wing'].edges[2740], 
            mdb.models['Model-1'].parts['Wing'].edges[2741], 
            mdb.models['Model-1'].parts['Wing'].edges[2760], 
            mdb.models['Model-1'].parts['Wing'].edges[2761], 
            mdb.models['Model-1'].parts['Wing'].edges[2762], 
            mdb.models['Model-1'].parts['Wing'].edges[2763], 
            mdb.models['Model-1'].parts['Wing'].edges[2764], 
            mdb.models['Model-1'].parts['Wing'].edges[2765], 
            mdb.models['Model-1'].parts['Wing'].edges[2766], 
            mdb.models['Model-1'].parts['Wing'].edges[2767], 
            mdb.models['Model-1'].parts['Wing'].edges[2768], 
            mdb.models['Model-1'].parts['Wing'].edges[2769], 
            mdb.models['Model-1'].parts['Wing'].edges[2770], 
            mdb.models['Model-1'].parts['Wing'].edges[2771], 
            mdb.models['Model-1'].parts['Wing'].edges[2772], 
            mdb.models['Model-1'].parts['Wing'].edges[2774], 
            mdb.models['Model-1'].parts['Wing'].edges[2775], 
            mdb.models['Model-1'].parts['Wing'].edges[2776], 
            mdb.models['Model-1'].parts['Wing'].edges[2777], 
            mdb.models['Model-1'].parts['Wing'].edges[2778], 
            mdb.models['Model-1'].parts['Wing'].edges[2779], 
            mdb.models['Model-1'].parts['Wing'].edges[2780], 
            mdb.models['Model-1'].parts['Wing'].edges[2781], 
            mdb.models['Model-1'].parts['Wing'].edges[2782], 
            mdb.models['Model-1'].parts['Wing'].edges[2783], 
            mdb.models['Model-1'].parts['Wing'].edges[2784], 
            mdb.models['Model-1'].parts['Wing'].edges[2785], 
            mdb.models['Model-1'].parts['Wing'].edges[2786], 
            mdb.models['Model-1'].parts['Wing'].edges[2787], 
            mdb.models['Model-1'].parts['Wing'].edges[2788], 
            mdb.models['Model-1'].parts['Wing'].edges[2789], 
            mdb.models['Model-1'].parts['Wing'].edges[2790], 
            mdb.models['Model-1'].parts['Wing'].edges[2791], 
            mdb.models['Model-1'].parts['Wing'].edges[2792], 
            mdb.models['Model-1'].parts['Wing'].edges[2793])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[3860], 
            mdb.models['Model-1'].parts['Wing'].edges[3862], 
            mdb.models['Model-1'].parts['Wing'].edges[3885], 
            mdb.models['Model-1'].parts['Wing'].edges[3886], 
            mdb.models['Model-1'].parts['Wing'].edges[3887], 
            mdb.models['Model-1'].parts['Wing'].edges[3888], 
            mdb.models['Model-1'].parts['Wing'].edges[3889], 
            mdb.models['Model-1'].parts['Wing'].edges[3890], 
            mdb.models['Model-1'].parts['Wing'].edges[3891], 
            mdb.models['Model-1'].parts['Wing'].edges[3892], 
            mdb.models['Model-1'].parts['Wing'].edges[3893], 
            mdb.models['Model-1'].parts['Wing'].edges[3894], 
            mdb.models['Model-1'].parts['Wing'].edges[3895], 
            mdb.models['Model-1'].parts['Wing'].edges[3896], 
            mdb.models['Model-1'].parts['Wing'].edges[3897], 
            mdb.models['Model-1'].parts['Wing'].edges[3898], 
            mdb.models['Model-1'].parts['Wing'].edges[3899], 
            mdb.models['Model-1'].parts['Wing'].edges[3900], 
            mdb.models['Model-1'].parts['Wing'].edges[3901], 
            mdb.models['Model-1'].parts['Wing'].edges[3902], 
            mdb.models['Model-1'].parts['Wing'].edges[3921], 
            mdb.models['Model-1'].parts['Wing'].edges[3922], 
            mdb.models['Model-1'].parts['Wing'].edges[3923], 
            mdb.models['Model-1'].parts['Wing'].edges[3924], 
            mdb.models['Model-1'].parts['Wing'].edges[3925], 
            mdb.models['Model-1'].parts['Wing'].edges[3926], 
            mdb.models['Model-1'].parts['Wing'].edges[3927], 
            mdb.models['Model-1'].parts['Wing'].edges[3928], 
            mdb.models['Model-1'].parts['Wing'].edges[3929], 
            mdb.models['Model-1'].parts['Wing'].edges[3930], 
            mdb.models['Model-1'].parts['Wing'].edges[3931], 
            mdb.models['Model-1'].parts['Wing'].edges[3932], 
            mdb.models['Model-1'].parts['Wing'].edges[3933], 
            mdb.models['Model-1'].parts['Wing'].edges[3935], 
            mdb.models['Model-1'].parts['Wing'].edges[3936], 
            mdb.models['Model-1'].parts['Wing'].edges[3937], 
            mdb.models['Model-1'].parts['Wing'].edges[3938], 
            mdb.models['Model-1'].parts['Wing'].edges[3939], 
            mdb.models['Model-1'].parts['Wing'].edges[3940], 
            mdb.models['Model-1'].parts['Wing'].edges[3941], 
            mdb.models['Model-1'].parts['Wing'].edges[3942], 
            mdb.models['Model-1'].parts['Wing'].edges[3943], 
            mdb.models['Model-1'].parts['Wing'].edges[3944], 
            mdb.models['Model-1'].parts['Wing'].edges[3945], 
            mdb.models['Model-1'].parts['Wing'].edges[3946], 
            mdb.models['Model-1'].parts['Wing'].edges[3947], 
            mdb.models['Model-1'].parts['Wing'].edges[3948], 
            mdb.models['Model-1'].parts['Wing'].edges[3949], 
            mdb.models['Model-1'].parts['Wing'].edges[3950], 
            mdb.models['Model-1'].parts['Wing'].edges[3951], 
            mdb.models['Model-1'].parts['Wing'].edges[3952], 
            mdb.models['Model-1'].parts['Wing'].edges[3953], 
            mdb.models['Model-1'].parts['Wing'].edges[3954])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[105], 
            mdb.models['Model-1'].parts['Wing'].edges[108], 
            mdb.models['Model-1'].parts['Wing'].edges[111], 
            mdb.models['Model-1'].parts['Wing'].edges[114], 
            mdb.models['Model-1'].parts['Wing'].edges[117], 
            mdb.models['Model-1'].parts['Wing'].edges[120], 
            mdb.models['Model-1'].parts['Wing'].edges[123], 
            mdb.models['Model-1'].parts['Wing'].edges[126], 
            mdb.models['Model-1'].parts['Wing'].edges[129], 
            mdb.models['Model-1'].parts['Wing'].edges[132], 
            mdb.models['Model-1'].parts['Wing'].edges[135], 
            mdb.models['Model-1'].parts['Wing'].edges[138], 
            mdb.models['Model-1'].parts['Wing'].edges[141], 
            mdb.models['Model-1'].parts['Wing'].edges[144], 
            mdb.models['Model-1'].parts['Wing'].edges[147], 
            mdb.models['Model-1'].parts['Wing'].edges[150], 
            mdb.models['Model-1'].parts['Wing'].edges[153], 
            mdb.models['Model-1'].parts['Wing'].edges[156], 
            mdb.models['Model-1'].parts['Wing'].edges[159]), (
            mdb.models['Model-1'].parts['Wing'].edges[2805], 
            mdb.models['Model-1'].parts['Wing'].edges[2807], 
            mdb.models['Model-1'].parts['Wing'].edges[2830], 
            mdb.models['Model-1'].parts['Wing'].edges[2831], 
            mdb.models['Model-1'].parts['Wing'].edges[2832], 
            mdb.models['Model-1'].parts['Wing'].edges[2833], 
            mdb.models['Model-1'].parts['Wing'].edges[2834], 
            mdb.models['Model-1'].parts['Wing'].edges[2835], 
            mdb.models['Model-1'].parts['Wing'].edges[2836], 
            mdb.models['Model-1'].parts['Wing'].edges[2837], 
            mdb.models['Model-1'].parts['Wing'].edges[2838], 
            mdb.models['Model-1'].parts['Wing'].edges[2839], 
            mdb.models['Model-1'].parts['Wing'].edges[2840], 
            mdb.models['Model-1'].parts['Wing'].edges[2841], 
            mdb.models['Model-1'].parts['Wing'].edges[2842], 
            mdb.models['Model-1'].parts['Wing'].edges[2843], 
            mdb.models['Model-1'].parts['Wing'].edges[2844], 
            mdb.models['Model-1'].parts['Wing'].edges[2845], 
            mdb.models['Model-1'].parts['Wing'].edges[2846], 
            mdb.models['Model-1'].parts['Wing'].edges[2847], 
            mdb.models['Model-1'].parts['Wing'].edges[2866], 
            mdb.models['Model-1'].parts['Wing'].edges[2867], 
            mdb.models['Model-1'].parts['Wing'].edges[2868], 
            mdb.models['Model-1'].parts['Wing'].edges[2869], 
            mdb.models['Model-1'].parts['Wing'].edges[2870], 
            mdb.models['Model-1'].parts['Wing'].edges[2871], 
            mdb.models['Model-1'].parts['Wing'].edges[2872], 
            mdb.models['Model-1'].parts['Wing'].edges[2873], 
            mdb.models['Model-1'].parts['Wing'].edges[2874], 
            mdb.models['Model-1'].parts['Wing'].edges[2875], 
            mdb.models['Model-1'].parts['Wing'].edges[2876], 
            mdb.models['Model-1'].parts['Wing'].edges[2877], 
            mdb.models['Model-1'].parts['Wing'].edges[2878], 
            mdb.models['Model-1'].parts['Wing'].edges[2880], 
            mdb.models['Model-1'].parts['Wing'].edges[2881], 
            mdb.models['Model-1'].parts['Wing'].edges[2882], 
            mdb.models['Model-1'].parts['Wing'].edges[2883], 
            mdb.models['Model-1'].parts['Wing'].edges[2884], 
            mdb.models['Model-1'].parts['Wing'].edges[2885], 
            mdb.models['Model-1'].parts['Wing'].edges[2886], 
            mdb.models['Model-1'].parts['Wing'].edges[2887], 
            mdb.models['Model-1'].parts['Wing'].edges[2888], 
            mdb.models['Model-1'].parts['Wing'].edges[2889], 
            mdb.models['Model-1'].parts['Wing'].edges[2890], 
            mdb.models['Model-1'].parts['Wing'].edges[2891], 
            mdb.models['Model-1'].parts['Wing'].edges[2892], 
            mdb.models['Model-1'].parts['Wing'].edges[2893], 
            mdb.models['Model-1'].parts['Wing'].edges[2894], 
            mdb.models['Model-1'].parts['Wing'].edges[2895], 
            mdb.models['Model-1'].parts['Wing'].edges[2896], 
            mdb.models['Model-1'].parts['Wing'].edges[2897], 
            mdb.models['Model-1'].parts['Wing'].edges[2898], 
            mdb.models['Model-1'].parts['Wing'].edges[2899])), startCondition=NONE)
        mdb.models['Model-1'].parts['Wing'].ShellLoft(endCondition=NONE, loftsections=(
            (mdb.models['Model-1'].parts['Wing'].edges[2], 
            mdb.models['Model-1'].parts['Wing'].edges[5], 
            mdb.models['Model-1'].parts['Wing'].edges[8], 
            mdb.models['Model-1'].parts['Wing'].edges[11], 
            mdb.models['Model-1'].parts['Wing'].edges[14], 
            mdb.models['Model-1'].parts['Wing'].edges[17], 
            mdb.models['Model-1'].parts['Wing'].edges[20], 
            mdb.models['Model-1'].parts['Wing'].edges[23], 
            mdb.models['Model-1'].parts['Wing'].edges[26], 
            mdb.models['Model-1'].parts['Wing'].edges[29], 
            mdb.models['Model-1'].parts['Wing'].edges[32], 
            mdb.models['Model-1'].parts['Wing'].edges[35], 
            mdb.models['Model-1'].parts['Wing'].edges[38], 
            mdb.models['Model-1'].parts['Wing'].edges[41], 
            mdb.models['Model-1'].parts['Wing'].edges[44], 
            mdb.models['Model-1'].parts['Wing'].edges[47], 
            mdb.models['Model-1'].parts['Wing'].edges[50], 
            mdb.models['Model-1'].parts['Wing'].edges[53], 
            mdb.models['Model-1'].parts['Wing'].edges[56], 
            mdb.models['Model-1'].parts['Wing'].edges[59], 
            mdb.models['Model-1'].parts['Wing'].edges[62], 
            mdb.models['Model-1'].parts['Wing'].edges[65], 
            mdb.models['Model-1'].parts['Wing'].edges[68], 
            mdb.models['Model-1'].parts['Wing'].edges[71], 
            mdb.models['Model-1'].parts['Wing'].edges[74], 
            mdb.models['Model-1'].parts['Wing'].edges[77], 
            mdb.models['Model-1'].parts['Wing'].edges[80], 
            mdb.models['Model-1'].parts['Wing'].edges[83], 
            mdb.models['Model-1'].parts['Wing'].edges[86], 
            mdb.models['Model-1'].parts['Wing'].edges[89], 
            mdb.models['Model-1'].parts['Wing'].edges[92], 
            mdb.models['Model-1'].parts['Wing'].edges[95], 
            mdb.models['Model-1'].parts['Wing'].edges[98], 
            mdb.models['Model-1'].parts['Wing'].edges[101], 
            mdb.models['Model-1'].parts['Wing'].edges[104], 
            mdb.models['Model-1'].parts['Wing'].edges[107], 
            mdb.models['Model-1'].parts['Wing'].edges[110], 
            mdb.models['Model-1'].parts['Wing'].edges[113], 
            mdb.models['Model-1'].parts['Wing'].edges[116], 
            mdb.models['Model-1'].parts['Wing'].edges[119], 
            mdb.models['Model-1'].parts['Wing'].edges[122], 
            mdb.models['Model-1'].parts['Wing'].edges[125], 
            mdb.models['Model-1'].parts['Wing'].edges[128], 
            mdb.models['Model-1'].parts['Wing'].edges[131], 
            mdb.models['Model-1'].parts['Wing'].edges[134], 
            mdb.models['Model-1'].parts['Wing'].edges[137], 
            mdb.models['Model-1'].parts['Wing'].edges[140], 
            mdb.models['Model-1'].parts['Wing'].edges[143], 
            mdb.models['Model-1'].parts['Wing'].edges[146], 
            mdb.models['Model-1'].parts['Wing'].edges[149], 
            mdb.models['Model-1'].parts['Wing'].edges[152], 
            mdb.models['Model-1'].parts['Wing'].edges[155], 
            mdb.models['Model-1'].parts['Wing'].edges[158]), (
            mdb.models['Model-1'].parts['Wing'].edges[2849], 
            mdb.models['Model-1'].parts['Wing'].edges[2851], 
            mdb.models['Model-1'].parts['Wing'].edges[2853], 
            mdb.models['Model-1'].parts['Wing'].edges[2854], 
            mdb.models['Model-1'].parts['Wing'].edges[2855], 
            mdb.models['Model-1'].parts['Wing'].edges[2856], 
            mdb.models['Model-1'].parts['Wing'].edges[2857], 
            mdb.models['Model-1'].parts['Wing'].edges[2858], 
            mdb.models['Model-1'].parts['Wing'].edges[2859], 
            mdb.models['Model-1'].parts['Wing'].edges[2860], 
            mdb.models['Model-1'].parts['Wing'].edges[2861], 
            mdb.models['Model-1'].parts['Wing'].edges[2862], 
            mdb.models['Model-1'].parts['Wing'].edges[2863], 
            mdb.models['Model-1'].parts['Wing'].edges[2864], 
            mdb.models['Model-1'].parts['Wing'].edges[2865], 
            mdb.models['Model-1'].parts['Wing'].edges[2866], 
            mdb.models['Model-1'].parts['Wing'].edges[2867], 
            mdb.models['Model-1'].parts['Wing'].edges[2868], 
            mdb.models['Model-1'].parts['Wing'].edges[2869], 
            mdb.models['Model-1'].parts['Wing'].edges[2870], 
            mdb.models['Model-1'].parts['Wing'].edges[2871], 
            mdb.models['Model-1'].parts['Wing'].edges[2872], 
            mdb.models['Model-1'].parts['Wing'].edges[2873], 
            mdb.models['Model-1'].parts['Wing'].edges[2874], 
            mdb.models['Model-1'].parts['Wing'].edges[2875], 
            mdb.models['Model-1'].parts['Wing'].edges[2876], 
            mdb.models['Model-1'].parts['Wing'].edges[2877], 
            mdb.models['Model-1'].parts['Wing'].edges[2878], 
            mdb.models['Model-1'].parts['Wing'].edges[2879], 
            mdb.models['Model-1'].parts['Wing'].edges[2880], 
            mdb.models['Model-1'].parts['Wing'].edges[2881], 
            mdb.models['Model-1'].parts['Wing'].edges[2882], 
            mdb.models['Model-1'].parts['Wing'].edges[2883], 
            mdb.models['Model-1'].parts['Wing'].edges[2885], 
            mdb.models['Model-1'].parts['Wing'].edges[2886], 
            mdb.models['Model-1'].parts['Wing'].edges[2887], 
            mdb.models['Model-1'].parts['Wing'].edges[2888], 
            mdb.models['Model-1'].parts['Wing'].edges[2889], 
            mdb.models['Model-1'].parts['Wing'].edges[2890], 
            mdb.models['Model-1'].parts['Wing'].edges[2891], 
            mdb.models['Model-1'].parts['Wing'].edges[2892], 
            mdb.models['Model-1'].parts['Wing'].edges[2893], 
            mdb.models['Model-1'].parts['Wing'].edges[2894], 
            mdb.models['Model-1'].parts['Wing'].edges[2895], 
            mdb.models['Model-1'].parts['Wing'].edges[2896], 
            mdb.models['Model-1'].parts['Wing'].edges[2897], 
            mdb.models['Model-1'].parts['Wing'].edges[2898], 
            mdb.models['Model-1'].parts['Wing'].edges[2899], 
            mdb.models['Model-1'].parts['Wing'].edges[2900], 
            mdb.models['Model-1'].parts['Wing'].edges[2901], 
            mdb.models['Model-1'].parts['Wing'].edges[2902], 
            mdb.models['Model-1'].parts['Wing'].edges[2903], 
            mdb.models['Model-1'].parts['Wing'].edges[2904])), startCondition=NONE)
            
            
        mdb.models['Model-1'].parts['Wing'].RemoveFaces(deleteCells=False, faceList=(
            mdb.models['Model-1'].parts['Wing'].faces.findAt((0.262803, 0.025595, 
            2.06), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((0.178205, 
            0.004937, 2.06), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((
            0.294514, -0.004752, 2.06), ), 
            mdb.models['Model-1'].parts['Wing'].faces.findAt((0.19305, 0.005644, 2.06), 
            ), mdb.models['Model-1'].parts['Wing'].faces.findAt((0.21126, 0.005644, 
            2.06), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((0.23626, 
            0.005644, 2.06), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((
            0.436733, 0.013634, 2.06), )))
        mdb.models['Model-1'].parts['Wing'].RemoveFaces(deleteCells=False, faceList=(
            mdb.models['Model-1'].parts['Wing'].faces.findAt((0.208773, 0.035119, 
            0.51), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((0.095797, 
            0.006959, 0.51), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((
            0.254689, -0.006606, 0.51), ), 
            mdb.models['Model-1'].parts['Wing'].faces.findAt((0.119267, 0.007745, 
            0.51), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((0.147353, 
            0.007745, 0.51), ), mdb.models['Model-1'].parts['Wing'].faces.findAt((
            0.172353, 0.007745, 0.51), ), 
            mdb.models['Model-1'].parts['Wing'].faces.findAt((0.44176, 0.019318, 0.51), 
            )))


    CreateGeometry()

    #sets of components
    def createcomponentsets():
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.34655, -0.046229, 
            2.388294), ), ((0.289491, -0.047265, 2.388092), ), ((0.352472, -0.058794, 
            2.385851), ), ((0.294081, -0.059437, 2.385726), ), ((0.303648, -0.059437, 
            2.385726), ), ((0.328648, -0.059437, 2.385726), ), ((0.444492, -0.053008, 
            2.386976), ), ((0.3298, -0.031958, 2.323612), ), ((0.265066, -0.047082, 
            2.320672), ), ((0.351569, -0.05393, 2.319341), ), ((0.273874, -0.046441, 
            2.320796), ), ((0.285171, -0.046441, 2.320796), ), ((0.310171, -0.046441, 
            2.320796), ), ((0.44294, -0.03973, 2.322101), ), ((0.278264, 0.011952, 
            2.12459), ), ((0.19825, -0.007137, 2.12088), ), ((0.307681, -0.016034, 
            2.11915), ), ((0.211702, -0.006455, 2.121012), ), ((0.228316, -0.006455, 
            2.121012), ), ((0.253316, -0.006455, 2.121012), ), ((0.438166, 0.001126, 
            2.122486), ), ((0.261409, 0.025841, 2.02), ), ((0.176078, 0.004989, 2.02), 
            ), ((0.293487, -0.0048, 2.02), ), ((0.191146, 0.005699, 2.02), ), ((
            0.209611, 0.005699, 2.02), ), ((0.234611, 0.005699, 2.02), ), ((0.436863, 
            0.01378, 2.02), ), ((0.257923, 0.026455, 1.92), ), ((0.170761, 0.00512, 
            1.92), ), ((0.290917, -0.00492, 1.92), ), ((0.186386, 0.005834, 1.92), ), (
            (0.205488, 0.005834, 1.92), ), ((0.230488, 0.005834, 1.92), ), ((0.437187, 
            0.014147, 1.92), ), ((0.254437, 0.02707, 1.82), ), ((0.165445, 0.00525, 
            1.82), ), ((0.288348, -0.005039, 1.82), ), ((0.181626, 0.00597, 1.82), ), (
            (0.201365, 0.00597, 1.82), ), ((0.226365, 0.00597, 1.82), ), ((0.437512, 
            0.014514, 1.82), ), ((0.250951, 0.027684, 1.72), ), ((0.160128, 0.005381, 
            1.72), ), ((0.285778, -0.005159, 1.72), ), ((0.176865, 0.006105, 1.72), ), 
            ((0.197242, 0.006105, 1.72), ), ((0.222242, 0.006105, 1.72), ), ((0.437836, 
            0.014881, 1.72), ), ((0.247466, 0.028299, 1.62), ), ((0.154811, 0.005511, 
            1.62), ), ((0.283209, -0.005279, 1.62), ), ((0.172105, 0.006241, 1.62), ), 
            ((0.193119, 0.006241, 1.62), ), ((0.218119, 0.006241, 1.62), ), ((0.43816, 
            0.015247, 1.62), ), ((0.24398, 0.028913, 1.52), ), ((0.149495, 0.005642, 
            1.52), ), ((0.28064, -0.005398, 1.52), ), ((0.167345, 0.006376, 1.52), ), (
            (0.188996, 0.006376, 1.52), ), ((0.213996, 0.006376, 1.52), ), ((0.438485, 
            0.015614, 1.52), ), ((0.240494, 0.029527, 1.42), ), ((0.144178, 0.005772, 
            1.42), ), ((0.27807, -0.005518, 1.42), ), ((0.162585, 0.006512, 1.42), ), (
            (0.184873, 0.006512, 1.42), ), ((0.209873, 0.006512, 1.42), ), ((0.438809, 
            0.015981, 1.42), ), ((0.237008, 0.030142, 1.32), ), ((0.138861, 0.005902, 
            1.32), ), ((0.275501, -0.005637, 1.32), ), ((0.157825, 0.006647, 1.32), ), 
            ((0.18075, 0.006647, 1.32), ), ((0.20575, 0.006647, 1.32), ), ((0.439133, 
            0.016348, 1.32), ), ((0.233522, 0.030756, 1.22), ), ((0.133545, 0.006033, 
            1.22), ), ((0.272931, -0.005757, 1.22), ), ((0.153064, 0.006783, 1.22), ), 
            ((0.176627, 0.006783, 1.22), ), ((0.201627, 0.006783, 1.22), ), ((0.439457, 
            0.016714, 1.22), ), ((0.230037, 0.031371, 1.12), ), ((0.128228, 0.006163, 
            1.12), ), ((0.270362, -0.005877, 1.12), ), ((0.148304, 0.006918, 1.12), ), 
            ((0.172504, 0.006918, 1.12), ), ((0.197504, 0.006918, 1.12), ), ((0.439782, 
            0.017081, 1.12), ), ((0.226551, 0.031985, 1.02), ), ((0.122912, 0.006294, 
            1.02), ), ((0.267793, -0.005996, 1.02), ), ((0.143544, 0.007054, 1.02), ), 
            ((0.168381, 0.007054, 1.02), ), ((0.193381, 0.007054, 1.02), ), ((0.440106, 
            0.017448, 1.02), ), ((0.223065, 0.0326, 0.92), ), ((0.117595, 0.006424, 
            0.92), ), ((0.265223, -0.006116, 0.92), ), ((0.138784, 0.007189, 0.92), ), 
            ((0.164258, 0.007189, 0.92), ), ((0.189258, 0.007189, 0.92), ), ((0.44043, 
            0.017814, 0.92), ), ((0.216093, 0.033828, 0.72), ), ((0.106962, 0.006685, 
            0.72), ), ((0.260084, -0.006355, 0.72), ), ((0.129263, 0.00746, 0.72), ), (
            (0.156012, 0.00746, 0.72), ), ((0.181012, 0.00746, 0.72), ), ((0.441079, 
            0.018548, 0.72), ), ((0.212607, 0.034443, 0.62), ), ((0.101645, 0.006816, 
            0.62), ), ((0.257515, -0.006475, 0.62), ), ((0.124503, 0.007596, 0.62), ), 
            ((0.151889, 0.007596, 0.62), ), ((0.176889, 0.007596, 0.62), ), ((0.441403, 
            0.018915, 0.62), ), ((0.210516, 0.034812, 0.56), ), ((0.098455, 0.006894, 
            0.56), ), ((0.255973, -0.006547, 0.56), ), ((0.121647, 0.007677, 0.56), ), 
            ((0.149415, 0.007677, 0.56), ), ((0.174415, 0.007677, 0.56), ), ((0.441598, 
            0.019135, 0.56), ), ((0.208773, 0.035119, 0.46), ), ((0.095797, 0.006959, 
            0.46), ), ((0.254689, -0.006606, 0.46), ), ((0.119267, 0.007745, 0.46), ), 
            ((0.147353, 0.007745, 0.46), ), ((0.172353, 0.007745, 0.46), ), ((0.44176, 
            0.019318, 0.46), ), ((0.208773, 0.035119, 0.4), ), ((0.095797, 0.006959, 
            0.4), ), ((0.254689, -0.006606, 0.4), ), ((0.119267, 0.007745, 0.4), ), ((
            0.147353, 0.007745, 0.4), ), ((0.172353, 0.007745, 0.4), ), ((0.44176, 
            0.019318, 0.4), ), ((0.208773, 0.035119, 0.3), ), ((0.095797, 0.006959, 
            0.3), ), ((0.254689, -0.006606, 0.3), ), ((0.119267, 0.007745, 0.3), ), ((
            0.147353, 0.007745, 0.3), ), ((0.172353, 0.007745, 0.3), ), ((0.44176, 
            0.019318, 0.3), ), ((0.208773, 0.035119, 0.2), ), ((0.095797, 0.006959, 
            0.2), ), ((0.254689, -0.006606, 0.2), ), ((0.119267, 0.007745, 0.2), ), ((
            0.147353, 0.007745, 0.2), ), ((0.172353, 0.007745, 0.2), ), ((0.44176, 
            0.019318, 0.2), ), ((0.208773, 0.035119, 0.1), ), ((0.095797, 0.006959, 
            0.1), ), ((0.254689, -0.006606, 0.1), ), ((0.119267, 0.007745, 0.1), ), ((
            0.147353, 0.007745, 0.1), ), ((0.172353, 0.007745, 0.1), ), ((0.44176, 
            0.019318, 0.1), ), ((0.219579, 0.033214, 0.82), ), ((0.112278, 0.006555, 
            0.82), ), ((0.262654, -0.006236, 0.82), ), ((0.134024, 0.007325, 0.82), ), 
            ((0.160135, 0.007325, 0.82), ), ((0.185135, 0.007325, 0.82), ), ((0.440755, 
            0.018181, 0.82), ), ((0.304032, -0.010003, 2.224101), ), ((0.231658, 
            -0.02711, 2.220776), ), ((0.329625, -0.034982, 2.219246), ), ((0.242788, 
            -0.026448, 2.220904), ), ((0.256743, -0.026448, 2.220904), ), ((0.281743, 
            -0.026448, 2.220904), ), ((0.440553, -0.019302, 2.222294), ), ), name=
            'Ribs')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.208773, 0.035119, 
            0.0), ), ((0.10131, 0.032971, 0.0), ), ((0.229442, 0.009076, 0.0), ), ((
            0.119267, 0.007745, 0.0), ), ((0.147353, 0.007745, 0.0), ), ((0.172353, 
            0.007745, 0.0), ), ((0.44176, 0.019318, 0.0), ), ), name='Base Rib')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.186596, 0.034078, 
            0.653333), ), ((0.153263, 0.034078, 0.653333), ), ((0.191906, 0.034593, 
            0.58), ), ((0.189707, 0.034937, 0.526667), ), ((0.18902, 0.007745, 
            0.476667), ), ((0.18902, 0.007745, 0.42), ), ((0.18902, 0.007745, 
            0.333333), ), ((0.18902, 0.007745, 0.233333), ), ((0.18902, 0.007745, 
            0.133333), ), ((0.18902, 0.007745, 0.033333), ), ((0.141906, 0.034593, 
            0.58), ), ((0.139707, 0.034937, 0.526667), ), ((0.13902, 0.007745, 
            0.476667), ), ((0.13902, 0.007745, 0.42), ), ((0.13902, 0.007745, 
            0.333333), ), ((0.13902, 0.007745, 0.233333), ), ((0.13902, 0.007745, 
            0.133333), ), ((0.13902, 0.007745, 0.033333), ), ), name=
            'Spar web-before boom')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.307997, -0.03714, 
            2.34509), ), ((0.282886, -0.017976, 2.257144), ), ((0.254459, 0.008379, 
            2.158489), ), ((0.233911, 0.020538, 2.082645), ), ((0.226827, 0.025695, 
            2.033333), ), ((0.223529, 0.026091, 1.953333), ), ((0.219406, 0.026705, 
            1.853333), ), ((0.215283, 0.02732, 1.753333), ), ((0.21116, 0.027934, 
            1.653333), ), ((0.207037, 0.028549, 1.553333), ), ((0.202914, 0.029163, 
            1.453333), ), ((0.198791, 0.029777, 1.353333), ), ((0.194668, 0.030392, 
            1.253333), ), ((0.190545, 0.031006, 1.153333), ), ((0.186422, 0.031621, 
            1.053333), ), ((0.182299, 0.032235, 0.953333), ), ((0.178176, 0.03285, 
            0.853333), ), ((0.174053, 0.033464, 0.753333), ), ), name=
            'Spar Web-after boom')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.435, 0.019939, 
            0.653333), ), ((0.435, 0.020288, 0.58), ), ((0.435, 0.020526, 0.526667), ), 
            ((0.435, 0.009125, 0.476667), ), ((0.435, 0.009125, 0.42), ), ((0.435, 
            0.009125, 0.333333), ), ((0.435, 0.009125, 0.233333), ), ((0.435, 0.009125, 
            0.133333), ), ((0.435, 0.009125, 0.033333), ), ), name=
            'Rear Spar-before boom')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.435, -0.042558, 
            2.344037), ), ((0.435, -0.024931, 2.255792), ), ((0.435, -0.004961, 
            2.155896), ), ((0.435, 0.009824, 2.081346), ), ((0.435, 0.014054, 
            2.033333), ), ((0.435, 0.014347, 1.953333), ), ((0.435, 0.014777, 
            1.853333), ), ((0.435, 0.015207, 1.753333), ), ((0.435, 0.015637, 
            1.653333), ), ((0.435, 0.016067, 1.553333), ), ((0.435, 0.016498, 
            1.453333), ), ((0.435, 0.016928, 1.353333), ), ((0.435, 0.017358, 
            1.253333), ), ((0.435, 0.017788, 1.153333), ), ((0.435, 0.018218, 
            1.053333), ), ((0.435, 0.018648, 0.953333), ), ((0.435, 0.019078, 
            0.853333), ), ((0.435, 0.019508, 0.753333), ), ), name=
            'Rear Spar-after boom')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.005096, 0.001386, 
            0.58), (-0.987308, 0.143028, 0.069039)), ((0.00639, 0.006405, 0.58), (
            -0.927883, 0.367154, 0.065055)), ((0.009823, 0.013263, 0.58), (-0.817545, 
            0.572957, 0.057801)), ((0.015529, 0.020332, 0.58), (-0.700158, 0.71222, 
            0.050206)), ((0.0235, 0.027381, 0.58), (-0.591304, 0.805285, 0.043305)), ((
            0.033684, 0.034193, 0.58), (-0.490483, 0.870662, 0.03707)), ((0.046045, 
            0.040566, 0.58), (-0.398994, 0.91641, 0.031577)), ((0.060543, 0.046349, 
            0.58), (-0.318094, 0.947678, 0.026888)), ((0.077098, 0.051726, 0.58), (
            -0.292195, 0.95602, 0.025454)), ((0.097011, 0.057492, 0.58), (-0.256786, 
            0.96618, 0.023606)), ((0.122663, 0.061727, 0.58), (0.0, 0.99994, 
            0.010933)), ((0.201623, 0.061727, 0.58), (0.0, 0.99994, 0.010933)), ((
            0.233501, 0.06102, 0.58), (0.055922, 0.998393, 0.009145)), ((0.268113, 
            0.058626, 0.58), (0.10351, 0.994597, 0.007916)), ((0.296561, 0.055427, 
            0.58), (0.127802, 0.991772, 0.007403)), ((0.325256, 0.051522, 0.58), (
            0.148768, 0.988847, 0.007063)), ((0.353937, 0.047026, 0.58), (0.167005, 
            0.985932, 0.006859)), ((0.393658, 0.039757, 0.58), (0.19184, 0.981403, 
            0.006707)), ((0.448001, 0.029135, 0.58), (0.19184, 0.981403, 0.006707)), ((
            0.471124, 0.024382, 0.58), (0.205114, 0.978714, 0.006823)), ((0.494765, 
            0.019399, 0.58), (0.206854, 0.978348, 0.006846)), ((0.516381, 0.014821, 
            0.58), (0.207325, 0.978248, 0.006854)), ((0.535627, 0.010726, 0.58), (
            0.208594, 0.977978, 0.006881)), ((0.552224, 0.007166, 0.58), (0.210343, 
            0.977603, 0.006925)), ((0.565948, 0.004234, 0.58), (0.208057, 0.978093, 
            0.006861)), ((0.581828, 0.001113, 0.58), (0.187214, 0.982299, 0.006234)), (
            (0.581556, 0.000318, 0.58), (-0.051907, -0.99865, -0.001728)), ((0.564767, 
            0.001053, 0.58), (-0.022329, -0.99975, -0.000839)), ((0.549768, 0.001286, 
            0.58), (-0.004709, -0.999989, -0.000353)), ((0.531376, 0.00126, 0.58), (
            0.011647, -0.999932, 4.8e-05)), ((0.509911, 0.000903, 0.58), (0.025441, 
            -0.999676, 0.000337)), ((0.485725, 0.000183, 0.58), (0.03742, -0.999299, 
            0.000538)), ((0.457274, -0.00109, 0.58), (0.055871, -0.998438, 0.000763)), 
            ((0.398502, -0.004379, 0.58), (0.055871, -0.998438, 0.000763)), ((0.359168, 
            -0.00683, 0.58), (0.067641, -0.997709, 0.000722)), ((0.327263, -0.00904, 
            0.58), (0.069854, -0.997557, 0.000701)), ((0.295222, -0.01128, 0.58), (
            0.069674, -0.99757, 0.000704)), ((0.263437, -0.014301, 0.58), (0.107257, 
            -0.994231, -7e-05)), ((0.231565, -0.018069, 0.58), (0.122438, -0.992476, 
            -0.000469)), ((0.201623, -0.019389, 0.58), (0.0, -0.999994, 0.003434)), ((
            0.122663, -0.019389, 0.58), (0.0, -0.999994, 0.003434)), ((0.096607, 
            -0.019358, 0.58), (-0.001899, -0.999992, 0.003531)), ((0.074704, -0.018799, 
            0.58), (-0.039628, -0.999199, 0.005608)), ((0.055847, -0.017631, 0.58), (
            -0.074595, -0.997185, 0.007657)), ((0.039701, -0.01592, 0.58), (-0.123513, 
            -0.992286, 0.010673)), ((0.026435, -0.013665, 0.58), (-0.194483, -0.980788, 
            0.015231)), ((0.01621, -0.010896, 0.58), (-0.304348, -0.952295, 0.02251)), 
            ((0.0093, -0.00764, 0.58), (-0.509417, -0.85975, 0.036388)), ((0.005723, 
            -0.002174, 0.58), (-0.932433, -0.355412, 0.065202)), ((0.583278, 0.000321, 
            0.526667), (-0.051907, -0.99865, -0.001728)), ((0.56633, 0.001063, 
            0.526667), (-0.022329, -0.99975, -0.000839)), ((0.55119, 0.001298, 
            0.526667), (-0.004709, -0.999989, -0.000353)), ((0.532624, 0.001272, 
            0.526667), (0.011647, -0.999932, 4.8e-05)), ((0.510957, 0.000911, 
            0.526667), (0.025441, -0.999676, 0.000337)), ((0.486543, 0.000185, 
            0.526667), (0.03742, -0.999299, 0.000538)), ((0.457721, -0.001106, 
            0.526667), (0.055871, -0.998438, 0.000763)), ((0.398376, -0.004427, 
            0.526667), (0.055871, -0.998438, 0.000763)), ((0.358778, -0.006894, 
            0.526667), (0.067641, -0.997709, 0.000722)), ((0.326572, -0.009126, 
            0.526667), (0.069854, -0.997557, 0.000701)), ((0.294229, -0.011387, 
            0.526667), (0.069674, -0.99757, 0.000704)), ((0.262144, -0.014437, 
            0.526667), (0.107257, -0.994231, -7e-05)), ((0.22997, -0.018241, 0.526667), 
            (0.122438, -0.992476, -0.000469)), ((0.199584, -0.019572, 0.526667), (0.0, 
            -0.999994, 0.003434)), ((0.120113, -0.019572, 0.526667), (0.0, -0.999994, 
            0.003434)), ((0.093739, -0.019541, 0.526667), (-0.001899, -0.999992, 
            0.003531)), ((0.07163, -0.018977, 0.526667), (-0.039628, -0.999199, 
            0.005608)), ((0.052595, -0.017797, 0.526667), (-0.074595, -0.997185, 
            0.007657)), ((0.036297, -0.01607, 0.526667), (-0.123513, -0.992286, 
            0.010673)), ((0.022907, -0.013794, 0.526667), (-0.194483, -0.980788, 
            0.015231)), ((0.012585, -0.010998, 0.526667), (-0.304348, -0.952295, 
            0.02251)), ((0.005611, -0.007712, 0.526667), (-0.509417, -0.85975, 
            0.036388)), ((0.002001, -0.002193, 0.526667), (-0.932433, -0.355412, 
            0.065202)), ((0.001368, 0.001398, 0.526667), (-0.987308, 0.143028, 
            0.069039)), ((0.002674, 0.006464, 0.526667), (-0.927883, 0.367154, 
            0.065055)), ((0.006139, 0.013386, 0.526667), (-0.817545, 0.572957, 
            0.057801)), ((0.011899, 0.020522, 0.526667), (-0.700158, 0.71222, 
            0.050206)), ((0.019945, 0.027638, 0.526667), (-0.591304, 0.805285, 
            0.043305)), ((0.030224, 0.034515, 0.526667), (-0.490483, 0.870662, 
            0.03707)), ((0.042702, 0.040948, 0.526667), (-0.398994, 0.91641, 
            0.031577)), ((0.057336, 0.046785, 0.526667), (-0.318094, 0.947678, 
            0.026888)), ((0.074047, 0.052213, 0.526667), (-0.292195, 0.95602, 
            0.025454)), ((0.094146, 0.058034, 0.526667), (-0.256786, 0.96618, 
            0.023606)), ((0.120113, 0.062311, 0.526667), (0.0, 0.99994, 0.010933)), ((
            0.199584, 0.062311, 0.526667), (0.0, 0.99994, 0.010933)), ((0.231923, 
            0.061597, 0.526667), (0.055922, 0.998393, 0.009145)), ((0.266864, 0.05918, 
            0.526667), (0.10351, 0.994597, 0.007916)), ((0.295581, 0.055952, 0.526667), 
            (0.127802, 0.991772, 0.007403)), ((0.324547, 0.052009, 0.526667), (
            0.148768, 0.988847, 0.007063)), ((0.353498, 0.047471, 0.526667), (0.167005, 
            0.985932, 0.006859)), ((0.393485, 0.040156, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.448357, 0.029429, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.471804, 0.024611, 0.526667), (0.205114, 0.978714, 
            0.006823)), ((0.495668, 0.019581, 0.526667), (0.206854, 0.978348, 
            0.006846)), ((0.517488, 0.01496, 0.526667), (0.207325, 0.978248, 
            0.006854)), ((0.536915, 0.010826, 0.526667), (0.208594, 0.977978, 
            0.006881)), ((0.553668, 0.007233, 0.526667), (0.210343, 0.977603, 
            0.006925)), ((0.567521, 0.004274, 0.526667), (0.208057, 0.978093, 
            0.006861)), ((0.583553, 0.001123, 0.526667), (0.187214, 0.982299, 
            0.006234)), ((0.000203, 0.001398, 0.476667), (-0.989669, 0.14337, 0.0)), ((
            0.00151, 0.006476, 0.476667), (-0.929852, 0.367933, 0.0)), ((0.004984, 
            0.013419, 0.476667), (-0.818914, 0.573916, 0.0)), ((0.010758, 0.020575, 
            0.476667), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.476667), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.476667), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.476667), (-0.399193, 0.916867, 
            0.0)), ((0.056319, 0.046917, 0.476667), (-0.318209, 0.94802, 0.0)), ((
            0.073077, 0.052361, 0.476667), (-0.29229, 0.95633, 0.0)), ((0.093229, 
            0.058198, 0.476667), (-0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 
            0.476667), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.476667), (0.0, 1.0, 
            0.0)), ((0.231395, 0.061779, 0.476667), (0.055924, 0.998435, 0.0)), ((
            0.266448, 0.059356, 0.476667), (0.103513, 0.994628, 0.0)), ((0.295248, 
            0.056119, 0.476667), (0.127805, 0.991799, 0.0)), ((0.324299, 0.052166, 
            0.476667), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 0.476667), (
            0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.476667), (0.191844, 
            0.981425, 0.0)), ((0.44176, 0.030833, 0.476667), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.476667), (0.205119, 0.978737, 0.0)), ((
            0.488195, 0.021277, 0.476667), (0.206858, 0.978371, 0.0)), ((0.5108, 
            0.016494, 0.476667), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.476667), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.476667), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.476667), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.476667), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.476667), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.476667), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.476667), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.476667), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.476667), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.476667), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.476667), (0.055871, -0.998438, 
            0.0)), ((0.416659, -0.003417, 0.476667), (0.055871, -0.998438, 0.0)), ((
            0.369303, -0.006193, 0.476667), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.476667), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.476667), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.476667), (
            0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.476667), (0.122438, 
            -0.992476, 0.0)), ((0.208773, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((
            0.129143, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 
            0.476667), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.476667), 
            (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.476667), (-0.074597, 
            -0.997214, 0.0)), ((0.040354, -0.016754, 0.476667), (-0.12352, -0.992342, 
            0.0)), ((0.025938, -0.014654, 0.476667), (-0.194505, -0.980901, 0.0)), ((
            0.014545, -0.012018, 0.476667), (-0.304425, -0.952536, 0.0)), ((0.006398, 
            -0.008883, 0.476667), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 
            0.476667), (-0.934421, -0.35617, 0.0)), ((0.583833, 0.000321, 0.42), (
            -0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 0.42), (-0.022329, 
            -0.999751, 0.0)), ((0.55165, 0.001302, 0.42), (-0.004709, -0.999989, 0.0)), 
            ((0.533034, 0.001276, 0.42), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.42), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.42), 
            (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.42), (0.055871, 
            -0.998438, 0.0)), ((0.416659, -0.003417, 0.42), (0.055871, -0.998438, 
            0.0)), ((0.369303, -0.006193, 0.42), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.42), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.42), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.42), (0.107257, 
            -0.994231, 0.0)), ((0.240235, -0.016966, 0.42), (0.122438, -0.992476, 
            0.0)), ((0.208773, -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.129143, 
            -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.42), (
            -0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.42), (-0.039629, 
            -0.999214, 0.0)), ((0.057635, -0.018302, 0.42), (-0.074597, -0.997214, 
            0.0)), ((0.040354, -0.016754, 0.42), (-0.12352, -0.992342, 0.0)), ((
            0.025938, -0.014654, 0.42), (-0.194505, -0.980901, 0.0)), ((0.014545, 
            -0.012018, 0.42), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 
            0.42), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.42), (
            -0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 0.42), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.42), (-0.929852, 0.367933, 0.0)), ((
            0.004984, 0.013419, 0.42), (-0.818914, 0.573916, 0.0)), ((0.010758, 
            0.020575, 0.42), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.42), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.42), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.42), (-0.399193, 0.916867, 0.0)), 
            ((0.056319, 0.046917, 0.42), (-0.318209, 0.94802, 0.0)), ((0.073077, 
            0.052361, 0.42), (-0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.42), (
            -0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 0.42), (0.0, 1.0, 0.0)), 
            ((0.198897, 0.062493, 0.42), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.42), 
            (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.42), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.42), (0.127805, 0.991799, 0.0)), (
            (0.324299, 0.052166, 0.42), (0.148772, 0.988872, 0.0)), ((0.353335, 
            0.047615, 0.42), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.42), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.42), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.42), (0.205119, 0.978737, 0.0)), ((0.488195, 
            0.021277, 0.42), (0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.42), (
            0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 0.42), (0.208599, 0.978001, 
            0.0)), ((0.548842, 0.00839, 0.42), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.42), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.42), (
            0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.333333), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.333333), (-0.929852, 0.367933, 
            0.0)), ((0.004984, 0.013419, 0.333333), (-0.818914, 0.573916, 0.0)), ((
            0.010758, 0.020575, 0.333333), (-0.701042, 0.71312, 0.0)), ((0.018825, 
            0.027712, 0.333333), (-0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 
            0.333333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 0.333333), (
            -0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.333333), (-0.318209, 
            0.94802, 0.0)), ((0.073077, 0.052361, 0.333333), (-0.29229, 0.95633, 0.0)), 
            ((0.093229, 0.058198, 0.333333), (-0.256858, 0.966449, 0.0)), ((0.119267, 
            0.062493, 0.333333), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.333333), (
            0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.333333), (0.055924, 0.998435, 
            0.0)), ((0.266448, 0.059356, 0.333333), (0.103513, 0.994628, 0.0)), ((
            0.295248, 0.056119, 0.333333), (0.127805, 0.991799, 0.0)), ((0.324299, 
            0.052166, 0.333333), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 
            0.333333), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.333333), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.333333), (0.191844, 
            0.981425, 0.0)), ((0.463659, 0.026434, 0.333333), (0.205119, 0.978737, 
            0.0)), ((0.488195, 0.021277, 0.333333), (0.206858, 0.978371, 0.0)), ((
            0.5108, 0.016494, 0.333333), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.333333), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.333333), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.333333), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.333333), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.333333), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.333333), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.333333), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.333333), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.333333), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.333333), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.333333), (0.055871, -0.998438, 
            0.0)), ((0.416659, -0.003417, 0.333333), (0.055871, -0.998438, 0.0)), ((
            0.369303, -0.006193, 0.333333), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.333333), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.333333), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.333333), (
            0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.333333), (0.122438, 
            -0.992476, 0.0)), ((0.208773, -0.019629, 0.333333), (0.0, -1.0, 0.0)), ((
            0.129143, -0.019629, 0.333333), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 
            0.333333), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.333333), 
            (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.333333), (-0.074597, 
            -0.997214, 0.0)), ((0.040354, -0.016754, 0.333333), (-0.12352, -0.992342, 
            0.0)), ((0.025938, -0.014654, 0.333333), (-0.194505, -0.980901, 0.0)), ((
            0.014545, -0.012018, 0.333333), (-0.304425, -0.952536, 0.0)), ((0.006398, 
            -0.008883, 0.333333), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 
            0.333333), (-0.934421, -0.35617, 0.0)), ((0.583833, 0.000321, 0.233333), (
            -0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 0.233333), (-0.022329, 
            -0.999751, 0.0)), ((0.55165, 0.001302, 0.233333), (-0.004709, -0.999989, 
            0.0)), ((0.533034, 0.001276, 0.233333), (0.011647, -0.999932, 0.0)), ((
            0.511306, 0.000915, 0.233333), (0.025441, -0.999676, 0.0)), ((0.486823, 
            0.000187, 0.233333), (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 
            0.233333), (0.055871, -0.998438, 0.0)), ((0.416659, -0.003417, 0.233333), (
            0.055871, -0.998438, 0.0)), ((0.369303, -0.006193, 0.233333), (0.067641, 
            -0.99771, 0.0)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((
            0.33714, -0.008398, 0.233333), (0.069854, -0.997557, 0.0)), ((0.304702, 
            -0.010667, 0.233333), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 
            0.233333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.233333), (
            0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.233333), (0.0, -1.0, 
            0.0)), ((0.129143, -0.019629, 0.233333), (0.0, -1.0, 0.0)), ((0.101105, 
            -0.019614, 0.233333), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 
            0.233333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.233333), 
            (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.233333), (-0.12352, 
            -0.992342, 0.0)), ((0.025938, -0.014654, 0.233333), (-0.194505, -0.980901, 
            0.0)), ((0.014545, -0.012018, 0.233333), (-0.304425, -0.952536, 0.0)), ((
            0.006398, -0.008883, 0.233333), (-0.509755, -0.86032, 0.0)), ((0.001672, 
            -0.004386, 0.233333), (-0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 
            0.233333), (-0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.233333), (
            -0.929852, 0.367933, 0.0)), ((0.004984, 0.013419, 0.233333), (-0.818914, 
            0.573916, 0.0)), ((0.010758, 0.020575, 0.233333), (-0.701042, 0.71312, 
            0.0)), ((0.018825, 0.027712, 0.233333), (-0.591859, 0.806041, 0.0)), ((
            0.029132, 0.034609, 0.233333), (-0.49082, 0.871261, 0.0)), ((0.041644, 
            0.041062, 0.233333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 
            0.233333), (-0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.233333), (
            -0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.233333), (-0.256858, 
            0.966449, 0.0)), ((0.119267, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((
            0.198897, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 
            0.233333), (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.233333), (
            0.103513, 0.994628, 0.0)), ((0.295248, 0.056119, 0.233333), (0.127805, 
            0.991799, 0.0)), ((0.324299, 0.052166, 0.233333), (0.148772, 0.988872, 
            0.0)), ((0.353335, 0.047615, 0.233333), (0.167009, 0.985955, 0.0)), ((
            0.393406, 0.040285, 0.233333), (0.191844, 0.981425, 0.0)), ((0.44176, 
            0.030833, 0.233333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 
            0.233333), (0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.233333), (
            0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.233333), (0.20733, 
            0.978271, 0.0)), ((0.53112, 0.01218, 0.233333), (0.208599, 0.978001, 0.0)), 
            ((0.548842, 0.00839, 0.233333), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.233333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 
            0.233333), (0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.133333), (
            -0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.133333), (-0.929852, 
            0.367933, 0.0)), ((0.004984, 0.013419, 0.133333), (-0.818914, 0.573916, 
            0.0)), ((0.010758, 0.020575, 0.133333), (-0.701042, 0.71312, 0.0)), ((
            0.018825, 0.027712, 0.133333), (-0.591859, 0.806041, 0.0)), ((0.029132, 
            0.034609, 0.133333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 
            0.133333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.133333), (
            -0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.133333), (-0.29229, 
            0.95633, 0.0)), ((0.093229, 0.058198, 0.133333), (-0.256858, 0.966449, 
            0.0)), ((0.119267, 0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.198897, 
            0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.133333), (
            0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.133333), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.133333), (0.127805, 0.991799, 
            0.0)), ((0.324299, 0.052166, 0.133333), (0.148772, 0.988872, 0.0)), ((
            0.353335, 0.047615, 0.133333), (0.167009, 0.985955, 0.0)), ((0.393406, 
            0.040285, 0.133333), (0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 
            0.133333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 0.133333), (
            0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.133333), (0.206858, 
            0.978371, 0.0)), ((0.5108, 0.016494, 0.133333), (0.20733, 0.978271, 0.0)), 
            ((0.53112, 0.01218, 0.133333), (0.208599, 0.978001, 0.0)), ((0.548842, 
            0.00839, 0.133333), (0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 
            0.133333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.133333), (
            0.187218, 0.982318, 0.0)), ((0.583833, 0.000321, 0.133333), (-0.051907, 
            -0.998652, 0.0)), ((0.566831, 0.001066, 0.133333), (-0.022329, -0.999751, 
            0.0)), ((0.55165, 0.001302, 0.133333), (-0.004709, -0.999989, 0.0)), ((
            0.533034, 0.001276, 0.133333), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.133333), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 
            0.133333), (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.133333), (
            0.055871, -0.998438, 0.0)), ((0.416659, -0.003417, 0.133333), (0.055871, 
            -0.998438, 0.0)), ((0.369303, -0.006193, 0.133333), (0.067641, -0.99771, 
            0.0)), ((0.33714, -0.008398, 0.133333), (0.069854, -0.997557, 0.0)), ((
            0.304702, -0.010667, 0.133333), (0.069674, -0.99757, 0.0)), ((0.272393, 
            -0.01333, 0.133333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 
            0.133333), (0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.133333), (
            0.0, -1.0, 0.0)), ((0.129143, -0.019629, 0.133333), (0.0, -1.0, 0.0)), ((
            0.101105, -0.019614, 0.133333), (-0.001899, -0.999998, 0.0)), ((0.077593, 
            -0.019307, 0.133333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 
            0.133333), (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.133333), 
            (-0.12352, -0.992342, 0.0)), ((0.025938, -0.014654, 0.133333), (-0.194505, 
            -0.980901, 0.0)), ((0.014545, -0.012018, 0.133333), (-0.304425, -0.952536, 
            0.0)), ((0.006398, -0.008883, 0.133333), (-0.509755, -0.86032, 0.0)), ((
            0.001672, -0.004386, 0.133333), (-0.934421, -0.35617, 0.0)), ((0.583833, 
            0.000321, 0.033333), (-0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 
            0.033333), (-0.022329, -0.999751, 0.0)), ((0.55165, 0.001302, 0.033333), (
            -0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 0.033333), (0.011647, 
            -0.999932, 0.0)), ((0.511306, 0.000915, 0.033333), (0.025441, -0.999676, 
            0.0)), ((0.486823, 0.000187, 0.033333), (0.03742, -0.9993, 0.0)), ((
            0.457924, -0.001108, 0.033333), (0.055871, -0.998438, 0.0)), ((0.416659, 
            -0.003417, 0.033333), (0.055871, -0.998438, 0.0)), ((0.369303, -0.006193, 
            0.033333), (0.067641, -0.99771, 0.0)), ((0.33714, -0.008398, 0.033333), (
            0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 0.033333), (0.069674, 
            -0.99757, 0.0)), ((0.272393, -0.01333, 0.033333), (0.107257, -0.994231, 
            0.0)), ((0.240235, -0.016966, 0.033333), (0.122438, -0.992476, 0.0)), ((
            0.208773, -0.019629, 0.033333), (0.0, -1.0, 0.0)), ((0.129143, -0.019629, 
            0.033333), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.033333), (-0.001899, 
            -0.999998, 0.0)), ((0.077593, -0.019307, 0.033333), (-0.039629, -0.999214, 
            0.0)), ((0.057635, -0.018302, 0.033333), (-0.074597, -0.997214, 0.0)), ((
            0.040354, -0.016754, 0.033333), (-0.12352, -0.992342, 0.0)), ((0.025938, 
            -0.014654, 0.033333), (-0.194505, -0.980901, 0.0)), ((0.014545, -0.012018, 
            0.033333), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 0.033333), 
            (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.033333), (-0.934421, 
            -0.35617, 0.0)), ((0.000203, 0.001398, 0.033333), (-0.989669, 0.14337, 
            0.0)), ((0.00151, 0.006476, 0.033333), (-0.929852, 0.367933, 0.0)), ((
            0.004984, 0.013419, 0.033333), (-0.818914, 0.573916, 0.0)), ((0.010758, 
            0.020575, 0.033333), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 
            0.033333), (-0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.033333), (
            -0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 0.033333), (-0.399193, 
            0.916867, 0.0)), ((0.056319, 0.046917, 0.033333), (-0.318209, 0.94802, 
            0.0)), ((0.073077, 0.052361, 0.033333), (-0.29229, 0.95633, 0.0)), ((
            0.093229, 0.058198, 0.033333), (-0.256858, 0.966449, 0.0)), ((0.119267, 
            0.062493, 0.033333), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.033333), (
            0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.033333), (0.055924, 0.998435, 
            0.0)), ((0.266448, 0.059356, 0.033333), (0.103513, 0.994628, 0.0)), ((
            0.295248, 0.056119, 0.033333), (0.127805, 0.991799, 0.0)), ((0.324299, 
            0.052166, 0.033333), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 
            0.033333), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.033333), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.033333), (0.191844, 
            0.981425, 0.0)), ((0.463659, 0.026434, 0.033333), (0.205119, 0.978737, 
            0.0)), ((0.488195, 0.021277, 0.033333), (0.206858, 0.978371, 0.0)), ((
            0.5108, 0.016494, 0.033333), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.033333), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.033333), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.033333), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.033333), (0.187218, 0.982318, 
            0.0)), ), name='Skin-root')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.556624, 0.000275, 
            1.353333), (-0.051907, -0.99865, -0.001728)), ((0.54214, 0.000909, 
            1.353333), (-0.022329, -0.99975, -0.000839)), ((0.529192, 0.00111, 
            1.353333), (-0.004709, -0.999989, -0.000353)), ((0.513316, 0.001088, 
            1.353333), (0.011647, -0.999932, 4.8e-05)), ((0.494789, 0.000778, 
            1.353333), (0.025441, -0.999676, 0.000337)), ((0.473913, 0.000157, 
            1.353333), (0.03742, -0.999299, 0.000538)), ((0.450932, -0.000854, 
            1.353333), (0.055871, -0.998438, 0.000763)), ((0.400294, -0.003688, 
            1.353333), (0.055871, -0.998438, 0.000763)), ((0.364746, -0.005892, 
            1.353333), (0.067641, -0.997709, 0.000722)), ((0.337213, -0.0078, 
            1.353333), (0.069854, -0.997557, 0.000701)), ((0.309561, -0.009733, 
            1.353333), (0.069674, -0.99757, 0.000704)), ((0.282129, -0.012339, 
            1.353333), (0.107257, -0.994231, -7e-05)), ((0.254623, -0.01559, 1.353333), 
            (0.122438, -0.992476, -0.000469)), ((0.231087, -0.016733, 1.353333), (0.0, 
            -0.999994, 0.003434)), ((0.159517, -0.016733, 1.353333), (0.0, -0.999994, 
            0.003434)), ((0.138145, -0.016706, 1.353333), (-0.001899, -0.999992, 
            0.003531)), ((0.119238, -0.016225, 1.353333), (-0.039628, -0.999199, 
            0.005608)), ((0.102961, -0.015217, 1.353333), (-0.074595, -0.997185, 
            0.007657)), ((0.089024, -0.013741, 1.353333), (-0.123513, -0.992286, 
            0.010673)), ((0.077573, -0.011796, 1.353333), (-0.194483, -0.980788, 
            0.015231)), ((0.068745, -0.009406, 1.353333), (-0.304348, -0.952295, 
            0.02251)), ((0.062779, -0.006597, 1.353333), (-0.509417, -0.85975, 
            0.036388)), ((0.059689, -0.001882, 1.353333), (-0.932433, -0.355412, 
            0.065202)), ((0.059145, 0.0012, 1.353333), (-0.987308, 0.143028, 
            0.069039)), ((0.060264, 0.005534, 1.353333), (-0.927883, 0.367154, 
            0.065055)), ((0.06323, 0.011453, 1.353333), (-0.817545, 0.572957, 
            0.057801)), ((0.068156, 0.017553, 1.353333), (-0.700158, 0.71222, 
            0.050206)), ((0.075038, 0.023637, 1.353333), (-0.591304, 0.805285, 
            0.043305)), ((0.083828, 0.029516, 1.353333), (-0.490483, 0.870662, 
            0.03707)), ((0.094498, 0.035015, 1.353333), (-0.398994, 0.91641, 
            0.031577)), ((0.107012, 0.040005, 1.353333), (-0.318094, 0.947678, 
            0.026888)), ((0.121302, 0.044646, 1.353333), (-0.292195, 0.95602, 
            0.025454)), ((0.138492, 0.049623, 1.353333), (-0.256786, 0.96618, 
            0.023606)), ((0.159517, 0.053272, 1.353333), (0.0, 0.99994, 0.010933)), ((
            0.231087, 0.053272, 1.353333), (0.0, 0.99994, 0.010933)), ((0.256299, 
            0.052659, 1.353333), (0.055922, 0.998393, 0.009145)), ((0.286161, 0.050592, 
            1.353333), (0.10351, 0.994597, 0.007916)), ((0.310712, 0.047831, 1.353333), 
            (0.127802, 0.991772, 0.007403)), ((0.335478, 0.04446, 1.353333), (0.148768, 
            0.988847, 0.007063)), ((0.360229, 0.04058, 1.353333), (0.167005, 0.985932, 
            0.006859)), ((0.396121, 0.033991, 1.353333), (0.19184, 0.981403, 
            0.006707)), ((0.442942, 0.024839, 1.353333), (0.19184, 0.981403, 
            0.006707)), ((0.461313, 0.021047, 1.353333), (0.205114, 0.978714, 
            0.006823)), ((0.481718, 0.016746, 1.353333), (0.206854, 0.978348, 
            0.006846)), ((0.500375, 0.012795, 1.353333), (0.207325, 0.978248, 
            0.006854)), ((0.516987, 0.00926, 1.353333), (0.208594, 0.977978, 
            0.006881)), ((0.531313, 0.006188, 1.353333), (0.210343, 0.977603, 
            0.006925)), ((0.54316, 0.003657, 1.353333), (0.208057, 0.978093, 
            0.006861)), ((0.556861, 0.000964, 1.353333), (0.187214, 0.982299, 
            0.006234)), ((0.052156, 0.001225, 1.253333), (-0.987308, 0.143028, 
            0.069039)), ((0.053298, 0.005647, 1.253333), (-0.927883, 0.367154, 
            0.065055)), ((0.056324, 0.011687, 1.253333), (-0.817545, 0.572957, 
            0.057801)), ((0.061352, 0.017913, 1.253333), (-0.700158, 0.71222, 
            0.050206)), ((0.068374, 0.024122, 1.253333), (-0.591304, 0.805285, 
            0.043305)), ((0.077345, 0.030121, 1.253333), (-0.490483, 0.870662, 
            0.03707)), ((0.088234, 0.035734, 1.253333), (-0.398994, 0.91641, 
            0.031577)), ((0.101005, 0.040826, 1.253333), (-0.318094, 0.947678, 
            0.026888)), ((0.115588, 0.045562, 1.253333), (-0.292195, 0.95602, 
            0.025454)), ((0.133131, 0.050641, 1.253333), (-0.256786, 0.96618, 
            0.023606)), ((0.154757, 0.054365, 1.253333), (0.0, 0.99994, 0.010933)), ((
            0.227282, 0.054365, 1.253333), (0.0, 0.99994, 0.010933)), ((0.253355, 
            0.05374, 1.253333), (0.055922, 0.998393, 0.009145)), ((0.28383, 0.051631, 
            1.253333), (0.10351, 0.994597, 0.007916)), ((0.308885, 0.048813, 1.253333), 
            (0.127802, 0.991772, 0.007403)), ((0.334159, 0.045373, 1.253333), (
            0.148768, 0.988847, 0.007063)), ((0.359418, 0.041413, 1.253333), (0.167005, 
            0.985932, 0.006859)), ((0.395805, 0.034736, 1.253333), (0.19184, 0.981403, 
            0.006707)), ((0.443591, 0.025395, 1.253333), (0.19184, 0.981403, 
            0.006707)), ((0.462579, 0.021479, 1.253333), (0.205114, 0.978714, 
            0.006823)), ((0.483403, 0.017089, 1.253333), (0.206854, 0.978348, 
            0.006846)), ((0.502442, 0.013058, 1.253333), (0.207325, 0.978248, 
            0.006854)), ((0.519395, 0.00945, 1.253333), (0.208594, 0.977978, 
            0.006881)), ((0.534015, 0.006315, 1.253333), (0.210343, 0.977603, 
            0.006925)), ((0.546105, 0.003732, 1.253333), (0.208057, 0.978093, 
            0.006861)), ((0.560088, 0.000983, 1.253333), (0.187214, 0.982299, 
            0.006234)), ((0.559846, 0.000281, 1.253333), (-0.051907, -0.99865, 
            -0.001728)), ((0.545064, 0.000928, 1.253333), (-0.022329, -0.99975, 
            -0.000839)), ((0.531851, 0.001133, 1.253333), (-0.004709, -0.999989, 
            -0.000353)), ((0.51565, 0.00111, 1.253333), (0.011647, -0.999932, 
            4.8e-05)), ((0.496742, 0.000794, 1.253333), (0.025441, -0.999676, 
            0.000337)), ((0.475438, 0.000161, 1.253333), (0.03742, -0.999299, 
            0.000538)), ((0.451745, -0.000885, 1.253333), (0.055871, -0.998438, 
            0.000763)), ((0.400065, -0.003777, 1.253333), (0.055871, -0.998438, 
            0.000763)), ((0.364028, -0.006013, 1.253333), (0.067641, -0.997709, 
            0.000722)), ((0.335929, -0.00796, 1.253333), (0.069854, -0.997557, 
            0.000701)), ((0.30771, -0.009932, 1.253333), (0.069674, -0.99757, 
            0.000704)), ((0.279715, -0.012592, 1.253333), (0.107257, -0.994231, 
            -7e-05)), ((0.251644, -0.01591, 1.253333), (0.122438, -0.992476, 
            -0.000469)), ((0.227282, -0.017076, 1.253333), (0.0, -0.999994, 0.003434)), 
            ((0.154757, -0.017076, 1.253333), (0.0, -0.999994, 0.003434)), ((0.132776, 
            -0.017049, 1.253333), (-0.001899, -0.999992, 0.003531)), ((0.113482, 
            -0.016558, 1.253333), (-0.039628, -0.999199, 0.005608)), ((0.096871, 
            -0.015529, 1.253333), (-0.074595, -0.997185, 0.007657)), ((0.082648, 
            -0.014023, 1.253333), (-0.123513, -0.992286, 0.010673)), ((0.070962, 
            -0.012038, 1.253333), (-0.194483, -0.980788, 0.015231)), ((0.061953, 
            -0.009599, 1.253333), (-0.304348, -0.952295, 0.02251)), ((0.055864, 
            -0.006732, 1.253333), (-0.509417, -0.85975, 0.036388)), ((0.052711, 
            -0.00192, 1.253333), (-0.932433, -0.355412, 0.065202)), ((0.563068, 
            0.000286, 1.153333), (-0.051907, -0.99865, -0.001728)), ((0.547989, 
            0.000947, 1.153333), (-0.022329, -0.99975, -0.000839)), ((0.53451, 
            0.001156, 1.153333), (-0.004709, -0.999989, -0.000353)), ((0.517983, 
            0.001132, 1.153333), (0.011647, -0.999932, 4.8e-05)), ((0.498695, 0.00081, 
            1.153333), (0.025441, -0.999676, 0.000337)), ((0.476963, 0.000164, 
            1.153333), (0.03742, -0.999299, 0.000538)), ((0.452558, -0.000916, 
            1.153333), (0.055871, -0.998438, 0.000763)), ((0.399835, -0.003867, 
            1.153333), (0.055871, -0.998438, 0.000763)), ((0.36331, -0.006134, 
            1.153333), (0.067641, -0.997709, 0.000722)), ((0.334646, -0.00812, 
            1.153333), (0.069854, -0.997557, 0.000701)), ((0.305859, -0.010132, 
            1.153333), (0.069674, -0.99757, 0.000704)), ((0.277301, -0.012846, 
            1.153333), (0.107257, -0.994231, -7e-05)), ((0.248666, -0.016231, 
            1.153333), (0.122438, -0.992476, -0.000469)), ((0.223478, -0.01742, 
            1.153333), (0.0, -0.999994, 0.003434)), ((0.149997, -0.01742, 1.153333), (
            0.0, -0.999994, 0.003434)), ((0.127407, -0.017392, 1.153333), (-0.001899, 
            -0.999992, 0.003531)), ((0.107725, -0.016891, 1.153333), (-0.039628, 
            -0.999199, 0.005608)), ((0.09078, -0.015841, 1.153333), (-0.074595, 
            -0.997185, 0.007657)), ((0.076271, -0.014305, 1.153333), (-0.123513, 
            -0.992286, 0.010673)), ((0.06435, -0.01228, 1.153333), (-0.194483, 
            -0.980788, 0.015231)), ((0.05516, -0.009792, 1.153333), (-0.304348, 
            -0.952295, 0.02251)), ((0.048949, -0.006868, 1.153333), (-0.509417, 
            -0.85975, 0.036388)), ((0.045733, -0.001959, 1.153333), (-0.932433, 
            -0.355412, 0.065202)), ((0.045167, 0.001249, 1.153333), (-0.987308, 
            0.143028, 0.069039)), ((0.046332, 0.005761, 1.153333), (-0.927883, 
            0.367154, 0.065055)), ((0.049419, 0.011922, 1.153333), (-0.817545, 
            0.572957, 0.057801)), ((0.054547, 0.018273, 1.153333), (-0.700158, 0.71222, 
            0.050206)), ((0.061711, 0.024607, 1.153333), (-0.591304, 0.805285, 
            0.043305)), ((0.070862, 0.030727, 1.153333), (-0.490483, 0.870662, 
            0.03707)), ((0.08197, 0.036452, 1.153333), (-0.398994, 0.91641, 0.031577)), 
            ((0.094998, 0.041646, 1.153333), (-0.318094, 0.947678, 0.026888)), ((
            0.109873, 0.046478, 1.153333), (-0.292195, 0.95602, 0.025454)), ((0.127769, 
            0.05166, 1.153333), (-0.256786, 0.96618, 0.023606)), ((0.149997, 0.055459, 
            1.153333), (0.0, 0.99994, 0.010933)), ((0.223478, 0.055459, 1.153333), (
            0.0, 0.99994, 0.010933)), ((0.250411, 0.054821, 1.153333), (0.055922, 
            0.998393, 0.009145)), ((0.281499, 0.052669, 1.153333), (0.10351, 0.994597, 
            0.007916)), ((0.307058, 0.049795, 1.153333), (0.127802, 0.991772, 
            0.007403)), ((0.33284, 0.046286, 1.153333), (0.148768, 0.988847, 
            0.007063)), ((0.358607, 0.042246, 1.153333), (0.167005, 0.985932, 
            0.006859)), ((0.395489, 0.035481, 1.153333), (0.19184, 0.981403, 
            0.006707)), ((0.444239, 0.025952, 1.153333), (0.19184, 0.981403, 
            0.006707)), ((0.463845, 0.02191, 1.153333), (0.205114, 0.978714, 
            0.006823)), ((0.485088, 0.017433, 1.153333), (0.206854, 0.978348, 
            0.006846)), ((0.50451, 0.01332, 1.153333), (0.207325, 0.978248, 0.006854)), 
            ((0.521804, 0.00964, 1.153333), (0.208594, 0.977978, 0.006881)), ((
            0.536718, 0.006442, 1.153333), (0.210343, 0.977603, 0.006925)), ((0.549051, 
            0.003807, 1.153333), (0.208057, 0.978093, 0.006861)), ((0.563314, 0.001003, 
            1.153333), (0.187214, 0.982299, 0.006234)), ((0.038178, 0.001274, 
            1.053333), (-0.987308, 0.143028, 0.069039)), ((0.039366, 0.005874, 
            1.053333), (-0.927883, 0.367154, 0.065055)), ((0.042513, 0.012157, 
            1.053333), (-0.817545, 0.572957, 0.057801)), ((0.047743, 0.018633, 
            1.053333), (-0.700158, 0.71222, 0.050206)), ((0.055048, 0.025092, 
            1.053333), (-0.591304, 0.805285, 0.043305)), ((0.064379, 0.031332, 
            1.053333), (-0.490483, 0.870662, 0.03707)), ((0.075706, 0.03717, 1.053333), 
            (-0.398994, 0.91641, 0.031577)), ((0.08899, 0.042467, 1.053333), (
            -0.318094, 0.947678, 0.026888)), ((0.104159, 0.047394, 1.053333), (
            -0.292195, 0.95602, 0.025454)), ((0.122408, 0.052678, 1.053333), (
            -0.256786, 0.96618, 0.023606)), ((0.145237, 0.056552, 1.053333), (0.0, 
            0.99994, 0.010933)), ((0.219673, 0.056552, 1.053333), (0.0, 0.99994, 
            0.010933)), ((0.247466, 0.055902, 1.053333), (0.055922, 0.998393, 
            0.009145)), ((0.279168, 0.053708, 1.053333), (0.10351, 0.994597, 
            0.007916)), ((0.305231, 0.050777, 1.053333), (0.127802, 0.991772, 
            0.007403)), ((0.331521, 0.047198, 1.053333), (0.148768, 0.988847, 
            0.007063)), ((0.357797, 0.043079, 1.053333), (0.167005, 0.985932, 
            0.006859)), ((0.395174, 0.036226, 1.053333), (0.19184, 0.981403, 
            0.006707)), ((0.444888, 0.026508, 1.053333), (0.19184, 0.981403, 
            0.006707)), ((0.465111, 0.022342, 1.053333), (0.205114, 0.978714, 
            0.006823)), ((0.486772, 0.017776, 1.053333), (0.206854, 0.978348, 
            0.006846)), ((0.506578, 0.013583, 1.053333), (0.207325, 0.978248, 
            0.006854)), ((0.524212, 0.00983, 1.053333), (0.208594, 0.977978, 
            0.006881)), ((0.53942, 0.006568, 1.053333), (0.210343, 0.977603, 
            0.006925)), ((0.551996, 0.003882, 1.053333), (0.208057, 0.978093, 
            0.006861)), ((0.566541, 0.001023, 1.053333), (0.187214, 0.982299, 
            0.006234)), ((0.56629, 0.000292, 1.053333), (-0.051907, -0.99865, 
            -0.001728)), ((0.550914, 0.000965, 1.053333), (-0.022329, -0.99975, 
            -0.000839)), ((0.537169, 0.001178, 1.053333), (-0.004709, -0.999989, 
            -0.000353)), ((0.520316, 0.001155, 1.053333), (0.011647, -0.999932, 
            4.8e-05)), ((0.500648, 0.000826, 1.053333), (0.025441, -0.999676, 
            0.000337)), ((0.478488, 0.000167, 1.053333), (0.03742, -0.999299, 
            0.000538)), ((0.453371, -0.000947, 1.053333), (0.055871, -0.998438, 
            0.000763)), ((0.399605, -0.003956, 1.053333), (0.055871, -0.998438, 
            0.000763)), ((0.362592, -0.006255, 1.053333), (0.067641, -0.997709, 
            0.000722)), ((0.333363, -0.00828, 1.053333), (0.069854, -0.997557, 
            0.000701)), ((0.304008, -0.010332, 1.053333), (0.069674, -0.99757, 
            0.000704)), ((0.274888, -0.013099, 1.053333), (0.107257, -0.994231, 
            -7e-05)), ((0.245688, -0.016551, 1.053333), (0.122438, -0.992476, 
            -0.000469)), ((0.219673, -0.017763, 1.053333), (0.0, -0.999994, 0.003434)), 
            ((0.145237, -0.017763, 1.053333), (0.0, -0.999994, 0.003434)), ((0.122039, 
            -0.017735, 1.053333), (-0.001899, -0.999992, 0.003531)), ((0.101968, 
            -0.017224, 1.053333), (-0.039628, -0.999199, 0.005608)), ((0.08469, 
            -0.016154, 1.053333), (-0.074595, -0.997185, 0.007657)), ((0.069895, 
            -0.014587, 1.053333), (-0.123513, -0.992286, 0.010673)), ((0.057739, 
            -0.012522, 1.053333), (-0.194483, -0.980788, 0.015231)), ((0.048368, 
            -0.009985, 1.053333), (-0.304348, -0.952295, 0.02251)), ((0.042034, 
            -0.007003, 1.053333), (-0.509417, -0.85975, 0.036388)), ((0.038755, 
            -0.001997, 1.053333), (-0.932433, -0.355412, 0.065202)), ((0.569512, 
            0.000298, 0.953333), (-0.051907, -0.99865, -0.001728)), ((0.553838, 
            0.000984, 0.953333), (-0.022329, -0.99975, -0.000839)), ((0.539828, 
            0.001201, 0.953333), (-0.004709, -0.999989, -0.000353)), ((0.522649, 
            0.001177, 0.953333), (0.011647, -0.999932, 4.8e-05)), ((0.502601, 0.000842, 
            0.953333), (0.025441, -0.999676, 0.000337)), ((0.480012, 0.00017, 
            0.953333), (0.03742, -0.999299, 0.000538)), ((0.454184, -0.000978, 
            0.953333), (0.055871, -0.998438, 0.000763)), ((0.399375, -0.004045, 
            0.953333), (0.055871, -0.998438, 0.000763)), ((0.361874, -0.006376, 
            0.953333), (0.067641, -0.997709, 0.000722)), ((0.332079, -0.00844, 
            0.953333), (0.069854, -0.997557, 0.000701)), ((0.302157, -0.010532, 
            0.953333), (0.069674, -0.99757, 0.000704)), ((0.272474, -0.013352, 
            0.953333), (0.107257, -0.994231, -7e-05)), ((0.242709, -0.016871, 
            0.953333), (0.122438, -0.992476, -0.000469)), ((0.215869, -0.018107, 
            0.953333), (0.0, -0.999994, 0.003434)), ((0.140477, -0.018107, 0.953333), (
            0.0, -0.999994, 0.003434)), ((0.11667, -0.018078, 0.953333), (-0.001899, 
            -0.999992, 0.003531)), ((0.096212, -0.017557, 0.953333), (-0.039628, 
            -0.999199, 0.005608)), ((0.078599, -0.016466, 0.953333), (-0.074595, 
            -0.997185, 0.007657)), ((0.063518, -0.014869, 0.953333), (-0.123513, 
            -0.992286, 0.010673)), ((0.051127, -0.012764, 0.953333), (-0.194483, 
            -0.980788, 0.015231)), ((0.041575, -0.010178, 0.953333), (-0.304348, 
            -0.952295, 0.02251)), ((0.035119, -0.007138, 0.953333), (-0.509417, 
            -0.85975, 0.036388)), ((0.031777, -0.002036, 0.953333), (-0.932433, 
            -0.355412, 0.065202)), ((0.031189, 0.001298, 0.953333), (-0.987308, 
            0.143028, 0.069039)), ((0.032399, 0.005987, 0.953333), (-0.927883, 
            0.367154, 0.065055)), ((0.035608, 0.012392, 0.953333), (-0.817545, 
            0.572957, 0.057801)), ((0.040938, 0.018993, 0.953333), (-0.700158, 0.71222, 
            0.050206)), ((0.048384, 0.025576, 0.953333), (-0.591304, 0.805285, 
            0.043305)), ((0.057896, 0.031938, 0.953333), (-0.490483, 0.870662, 
            0.03707)), ((0.069442, 0.037889, 0.953333), (-0.398994, 0.91641, 
            0.031577)), ((0.082983, 0.043288, 0.953333), (-0.318094, 0.947678, 
            0.026888)), ((0.098445, 0.04831, 0.953333), (-0.292195, 0.95602, 
            0.025454)), ((0.117046, 0.053696, 0.953333), (-0.256786, 0.96618, 
            0.023606)), ((0.140477, 0.057646, 0.953333), (0.0, 0.99994, 0.010933)), ((
            0.215869, 0.057646, 0.953333), (0.0, 0.99994, 0.010933)), ((0.244522, 
            0.056983, 0.953333), (0.055922, 0.998393, 0.009145)), ((0.276838, 0.054746, 
            0.953333), (0.10351, 0.994597, 0.007916)), ((0.303404, 0.051759, 0.953333), 
            (0.127802, 0.991772, 0.007403)), ((0.330203, 0.048111, 0.953333), (
            0.148768, 0.988847, 0.007063)), ((0.356986, 0.043912, 0.953333), (0.167005, 
            0.985932, 0.006859)), ((0.394858, 0.036971, 0.953333), (0.19184, 0.981403, 
            0.006707)), ((0.445536, 0.027065, 0.953333), (0.19184, 0.981403, 
            0.006707)), ((0.466378, 0.022774, 0.953333), (0.205114, 0.978714, 
            0.006823)), ((0.488457, 0.01812, 0.953333), (0.206854, 0.978348, 
            0.006846)), ((0.508645, 0.013845, 0.953333), (0.207325, 0.978248, 
            0.006854)), ((0.526621, 0.01002, 0.953333), (0.208594, 0.977978, 
            0.006881)), ((0.542122, 0.006695, 0.953333), (0.210343, 0.977603, 
            0.006925)), ((0.554941, 0.003956, 0.953333), (0.208057, 0.978093, 
            0.006861)), ((0.569768, 0.001042, 0.953333), (0.187214, 0.982299, 
            0.006234)), ((0.0242, 0.001322, 0.853333), (-0.987308, 0.143028, 
            0.069039)), ((0.025433, 0.006101, 0.853333), (-0.927883, 0.367154, 
            0.065055)), ((0.028702, 0.012626, 0.853333), (-0.817545, 0.572957, 
            0.057801)), ((0.034134, 0.019353, 0.853333), (-0.700158, 0.71222, 
            0.050206)), ((0.041721, 0.026061, 0.853333), (-0.591304, 0.805285, 
            0.043305)), ((0.051413, 0.032543, 0.853333), (-0.490483, 0.870662, 
            0.03707)), ((0.063178, 0.038607, 0.853333), (-0.398994, 0.91641, 
            0.031577)), ((0.076976, 0.044109, 0.853333), (-0.318094, 0.947678, 
            0.026888)), ((0.092731, 0.049226, 0.853333), (-0.292195, 0.95602, 
            0.025454)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.111684, 
            0.054714, 0.853333), (-0.256786, 0.96618, 0.023606)), ((0.135716, 0.058739, 
            0.853333), (0.0, 0.99994, 0.010933)), ((0.212065, 0.058739, 0.853333), (
            0.0, 0.99994, 0.010933)), ((0.241578, 0.058064, 0.853333), (0.055922, 
            0.998393, 0.009145)), ((0.274507, 0.055785, 0.853333), (0.10351, 0.994597, 
            0.007916)), ((0.301577, 0.052741, 0.853333), (0.127802, 0.991772, 
            0.007403)), ((0.328884, 0.049024, 0.853333), (0.148768, 0.988847, 
            0.007063)), ((0.356175, 0.044745, 0.853333), (0.167005, 0.985932, 
            0.006859)), ((0.394542, 0.037717, 0.853333), (0.19184, 0.981403, 
            0.006707)), ((0.446185, 0.027622, 0.853333), (0.19184, 0.981403, 
            0.006707)), ((0.467644, 0.023206, 0.853333), (0.205114, 0.978714, 
            0.006823)), ((0.490142, 0.018463, 0.853333), (0.206854, 0.978348, 
            0.006846)), ((0.510713, 0.014107, 0.853333), (0.207325, 0.978248, 
            0.006854)), ((0.529029, 0.01021, 0.853333), (0.208594, 0.977978, 
            0.006881)), ((0.544825, 0.006822, 0.853333), (0.210343, 0.977603, 
            0.006925)), ((0.557887, 0.004031, 0.853333), (0.208057, 0.978093, 
            0.006861)), ((0.572995, 0.001062, 0.853333), (0.187214, 0.982299, 
            0.006234)), ((0.572734, 0.000303, 0.853333), (-0.051907, -0.99865, 
            -0.001728)), ((0.556763, 0.001003, 0.853333), (-0.022329, -0.99975, 
            -0.000839)), ((0.542487, 0.001224, 0.853333), (-0.004709, -0.999989, 
            -0.000353)), ((0.524982, 0.001199, 0.853333), (0.011647, -0.999932, 
            4.8e-05)), ((0.504554, 0.000858, 0.853333), (0.025441, -0.999676, 
            0.000337)), ((0.481537, 0.000174, 0.853333), (0.03742, -0.999299, 
            0.000538)), ((0.454997, -0.001009, 0.853333), (0.055871, -0.998438, 
            0.000763)), ((0.399146, -0.004134, 0.853333), (0.055871, -0.998438, 
            0.000763)), ((0.361156, -0.006497, 0.853333), (0.067641, -0.997709, 
            0.000722)), ((0.330796, -0.0086, 0.853333), (0.069854, -0.997557, 
            0.000701)), ((0.300306, -0.010732, 0.853333), (0.069674, -0.99757, 
            0.000704)), ((0.27006, -0.013606, 0.853333), (0.107257, -0.994231, 
            -7e-05)), ((0.239731, -0.017191, 0.853333), (0.122438, -0.992476, 
            -0.000469)), ((0.212065, -0.01845, 0.853333), (0.0, -0.999994, 0.003434)), 
            ((0.135716, -0.01845, 0.853333), (0.0, -0.999994, 0.003434)), ((0.111301, 
            -0.018421, 0.853333), (-0.001899, -0.999992, 0.003531)), ((0.090455, 
            -0.01789, 0.853333), (-0.039628, -0.999199, 0.005608)), ((0.072508, 
            -0.016778, 0.853333), (-0.074595, -0.997185, 0.007657)), ((0.057142, 
            -0.015151, 0.853333), (-0.123513, -0.992286, 0.010673)), ((0.044516, 
            -0.013006, 0.853333), (-0.194483, -0.980788, 0.015231)), ((0.034783, 
            -0.01037, 0.853333), (-0.304348, -0.952295, 0.02251)), ((0.028205, 
            -0.007273, 0.853333), (-0.509417, -0.85975, 0.036388)), ((0.024799, 
            -0.002074, 0.853333), (-0.932433, -0.355412, 0.065202)), ((0.575956, 
            0.000309, 0.753333), (-0.051907, -0.99865, -0.001728)), ((0.559687, 
            0.001021, 0.753333), (-0.022329, -0.99975, -0.000839)), ((0.545146, 
            0.001247, 0.753333), (-0.004709, -0.999989, -0.000353)), ((0.527316, 
            0.001222, 0.753333), (0.011647, -0.999932, 4.8e-05)), ((0.506507, 0.000874, 
            0.753333), (0.025441, -0.999676, 0.000337)), ((0.483062, 0.000177, 
            0.753333), (0.03742, -0.999299, 0.000538)), ((0.45581, -0.00104, 0.753333), 
            (0.055871, -0.998438, 0.000763)), ((0.398916, -0.004224, 0.753333), (
            0.055871, -0.998438, 0.000763)), ((0.360437, -0.006618, 0.753333), (
            0.067641, -0.997709, 0.000722)), ((0.329513, -0.008761, 0.753333), (
            0.069854, -0.997557, 0.000701)), ((0.298456, -0.010932, 0.753333), (
            0.069674, -0.99757, 0.000704)), ((0.267646, -0.013859, 0.753333), (
            0.107257, -0.994231, -7e-05)), ((0.236753, -0.017511, 0.753333), (0.122438, 
            -0.992476, -0.000469)), ((0.20826, -0.018794, 0.753333), (0.0, -0.999994, 
            0.003434)), ((0.130956, -0.018794, 0.753333), (0.0, -0.999994, 0.003434)), 
            ((0.105932, -0.018764, 0.753333), (-0.001899, -0.999992, 0.003531)), ((
            0.084699, -0.018223, 0.753333), (-0.039628, -0.999199, 0.005608)), ((
            0.066418, -0.01709, 0.753333), (-0.074595, -0.997185, 0.007657)), ((
            0.050765, -0.015433, 0.753333), (-0.123513, -0.992286, 0.010673)), ((
            0.037905, -0.013248, 0.753333), (-0.194483, -0.980788, 0.015231)), ((
            0.02799, -0.010563, 0.753333), (-0.304348, -0.952295, 0.02251)), ((0.02129, 
            -0.007409, 0.753333), (-0.509417, -0.85975, 0.036388)), ((0.017821, 
            -0.002112, 0.753333), (-0.932433, -0.355412, 0.065202)), ((0.017211, 
            0.001347, 0.753333), (-0.987308, 0.143028, 0.069039)), ((0.018467, 
            0.006214, 0.753333), (-0.927883, 0.367154, 0.065055)), ((0.021797, 
            0.012861, 0.753333), (-0.817545, 0.572957, 0.057801)), ((0.027329, 
            0.019713, 0.753333), (-0.700158, 0.71222, 0.050206)), ((0.035058, 0.026546, 
            0.753333), (-0.591304, 0.805285, 0.043305)), ((0.04493, 0.033149, 
            0.753333), (-0.490483, 0.870662, 0.03707)), ((0.056914, 0.039326, 
            0.753333), (-0.398994, 0.91641, 0.031577)), ((0.070969, 0.04493, 0.753333), 
            (-0.318094, 0.947678, 0.026888)), ((0.087017, 0.050142, 0.753333), (
            -0.292195, 0.95602, 0.025454)), ((0.106323, 0.055732, 0.753333), (
            -0.256786, 0.96618, 0.023606)), ((0.130956, 0.059832, 0.753333), (0.0, 
            0.99994, 0.010933)), ((0.20826, 0.059832, 0.753333), (0.0, 0.99994, 
            0.010933)), ((0.238634, 0.059145, 0.753333), (0.055922, 0.998393, 
            0.009145)), ((0.272176, 0.056823, 0.753333), (0.10351, 0.994597, 
            0.007916)), ((0.29975, 0.053722, 0.753333), (0.127802, 0.991772, 
            0.007403)), ((0.327565, 0.049936, 0.753333), (0.148768, 0.988847, 
            0.007063)), ((0.355364, 0.045578, 0.753333), (0.167005, 0.985932, 
            0.006859)), ((0.394227, 0.038462, 0.753333), (0.19184, 0.981403, 
            0.006707)), ((0.446833, 0.028178, 0.753333), (0.19184, 0.981403, 
            0.006707)), ((0.46891, 0.023637, 0.753333), (0.205114, 0.978714, 
            0.006823)), ((0.491827, 0.018807, 0.753333), (0.206854, 0.978348, 
            0.006846)), ((0.512781, 0.01437, 0.753333), (0.207325, 0.978248, 
            0.006854)), ((0.531438, 0.0104, 0.753333), (0.208594, 0.977978, 0.006881)), 
            ((0.547527, 0.006949, 0.753333), (0.210343, 0.977603, 0.006925)), ((
            0.560832, 0.004106, 0.753333), (0.208057, 0.978093, 0.006861)), ((0.576222, 
            0.001082, 0.753333), (0.187214, 0.982299, 0.006234)), ((0.579178, 0.000314, 
            0.653333), (-0.051907, -0.99865, -0.001728)), ((0.562612, 0.00104, 
            0.653333), (-0.022329, -0.99975, -0.000839)), ((0.547805, 0.001269, 
            0.653333), (-0.004709, -0.999989, -0.000353)), ((0.529649, 0.001244, 
            0.653333), (0.011647, -0.999932, 4.8e-05)), ((0.50846, 0.00089, 0.653333), 
            (0.025441, -0.999676, 0.000337)), ((0.484586, 0.00018, 0.653333), (0.03742, 
            -0.999299, 0.000538)), ((0.456623, -0.001071, 0.653333), (0.055871, 
            -0.998438, 0.000763)), ((0.398686, -0.004313, 0.653333), (0.055871, 
            -0.998438, 0.000763)), ((0.359719, -0.006739, 0.653333), (0.067641, 
            -0.997709, 0.000722)), ((0.328229, -0.008921, 0.653333), (0.069854, 
            -0.997557, 0.000701)), ((0.296605, -0.011132, 0.653333), (0.069674, 
            -0.99757, 0.000704)), ((0.265232, -0.014113, 0.653333), (0.107257, 
            -0.994231, -7e-05)), ((0.233774, -0.017831, 0.653333), (0.122438, 
            -0.992476, -0.000469)), ((0.204456, -0.019137, 0.653333), (0.0, -0.999994, 
            0.003434)), ((0.126196, -0.019137, 0.653333), (0.0, -0.999994, 0.003434)), 
            ((0.100564, -0.019106, 0.653333), (-0.001899, -0.999992, 0.003531)), ((
            0.078942, -0.018556, 0.653333), (-0.039628, -0.999199, 0.005608)), ((
            0.060327, -0.017403, 0.653333), (-0.074595, -0.997185, 0.007657)), ((
            0.044389, -0.015714, 0.653333), (-0.123513, -0.992286, 0.010673)), ((
            0.031293, -0.01349, 0.653333), (-0.194483, -0.980788, 0.015231)), ((
            0.021198, -0.010756, 0.653333), (-0.304348, -0.952295, 0.02251)), ((
            0.014375, -0.007544, 0.653333), (-0.509417, -0.85975, 0.036388)), ((
            0.010843, -0.002151, 0.653333), (-0.932433, -0.355412, 0.065202)), ((
            0.010222, 0.001371, 0.653333), (-0.987308, 0.143028, 0.069039)), ((
            0.011501, 0.006327, 0.653333), (-0.927883, 0.367154, 0.065055)), ((
            0.014891, 0.013096, 0.653333), (-0.817545, 0.572957, 0.057801)), ((
            0.020525, 0.020073, 0.653333), (-0.700158, 0.71222, 0.050206)), ((0.028394, 
            0.027031, 0.653333), (-0.591304, 0.805285, 0.043305)), ((0.038447, 
            0.033754, 0.653333), (-0.490483, 0.870662, 0.03707)), ((0.05065, 0.040044, 
            0.653333), (-0.398994, 0.91641, 0.031577)), ((0.064961, 0.045751, 
            0.653333), (-0.318094, 0.947678, 0.026888)), ((0.081303, 0.051058, 
            0.653333), (-0.292195, 0.95602, 0.025454)), ((0.100961, 0.056751, 
            0.653333), (-0.256786, 0.96618, 0.023606)), ((0.126196, 0.060926, 
            0.653333), (0.0, 0.99994, 0.010933)), ((0.204456, 0.060926, 0.653333), (
            0.0, 0.99994, 0.010933)), ((0.235689, 0.060225, 0.653333), (0.055922, 
            0.998393, 0.009145)), ((0.269845, 0.057862, 0.653333), (0.10351, 0.994597, 
            0.007916)), ((0.297923, 0.054704, 0.653333), (0.127802, 0.991772, 
            0.007403)), ((0.326246, 0.050849, 0.653333), (0.148768, 0.988847, 
            0.007063)), ((0.354554, 0.046411, 0.653333), (0.167005, 0.985932, 
            0.006859)), ((0.393911, 0.039207, 0.653333), (0.19184, 0.981403, 
            0.006707)), ((0.447482, 0.028735, 0.653333), (0.19184, 0.981403, 
            0.006707)), ((0.470176, 0.024069, 0.653333), (0.205114, 0.978714, 
            0.006823)), ((0.493511, 0.01915, 0.653333), (0.206854, 0.978348, 
            0.006846)), ((0.514848, 0.014632, 0.653333), (0.207325, 0.978248, 
            0.006854)), ((0.533846, 0.01059, 0.653333), (0.208594, 0.977978, 
            0.006881)), ((0.55023, 0.007076, 0.653333), (0.210343, 0.977603, 
            0.006925)), ((0.563778, 0.004181, 0.653333), (0.208057, 0.978093, 
            0.006861)), ((0.579448, 0.001101, 0.653333), (0.187214, 0.982299, 
            0.006234)), ), name='Skin-mid')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.222296, -0.054011, 
            2.34181), (-0.921246, 0.06129, 0.384121)), ((0.222992, -0.051389, 2.34232), 
            (-0.87222, 0.272656, 0.406069)), ((0.224821, -0.047822, 2.343013), (
            -0.778269, 0.476044, 0.409487)), ((0.227854, -0.044147, 2.343728), (
            -0.674525, 0.621555, 0.398353)), ((0.232086, -0.040482, 2.34444), (
            -0.575154, 0.723888, 0.38103)), ((0.23749, -0.036942, 2.345128), (
            -0.480724, 0.799158, 0.360904)), ((0.244048, -0.033631, 2.345772), (
            -0.393325, 0.854192, 0.340076)), ((0.251736, -0.030628, 2.346356), (
            -0.314905, 0.893516, 0.3201)), ((0.260514, -0.027832, 2.346899), (
            -0.289604, 0.904327, 0.313563)), ((0.271089, -0.02483, 2.347483), (
            -0.254874, 0.917694, 0.304757)), ((0.280898, -0.022657, 2.347905), (0.0, 
            0.972107, 0.23454)), ((0.334478, -0.022657, 2.347905), (0.0, 0.972107, 
            0.23454)), ((0.343446, -0.023035, 2.347832), (0.055897, 0.973641, 
            0.22113)), ((0.361734, -0.024284, 2.347589), (0.10349, 0.972127, 
            0.210377)), ((0.376803, -0.02595, 2.347265), (0.127788, 0.970349, 
            0.205166)), ((0.392002, -0.027983, 2.34687), (0.14876, 0.968233, 
            0.200986)), ((0.407192, -0.030323, 2.346415), (0.167001, 0.965944, 
            0.197645)), ((0.424035, -0.033301, 2.345836), (0.191841, 0.962164, 
            0.193486)), ((0.44294, -0.036928, 2.345131), (0.191841, 0.962164, 
            0.193486)), ((0.46489, -0.041203, 2.3443), (0.205116, 0.95972, 0.192005)), 
            ((0.477717, -0.043849, 2.343786), (0.206856, 0.959379, 0.19184)), ((
            0.489531, -0.046303, 2.343309), (0.207327, 0.959285, 0.191803)), ((
            0.500148, -0.048516, 2.342878), (0.208596, 0.959025, 0.191724)), ((
            0.509405, -0.05046, 2.342501), (0.210345, 0.958661, 0.191637)), ((0.519289, 
            -0.052536, 2.342097), (0.208059, 0.959142, 0.191723)), ((0.527678, 
            -0.054156, 2.341782), (0.187215, 0.963321, 0.192257)), ((0.52753, 
            -0.054581, 2.3417), (-0.051907, -0.980043, -0.191889)), ((0.518659, 
            -0.054201, 2.341773), (-0.022329, -0.981195, -0.191722)), ((0.507824, 
            -0.054068, 2.341799), (-0.004709, -0.981439, -0.191717)), ((0.497539, 
            -0.054135, 2.341786), (0.011647, -0.981353, -0.19186)), ((0.485685, 
            -0.054381, 2.341739), (0.025441, -0.981037, -0.192143)), ((0.472452, 
            -0.054817, 2.341654), (0.03742, -0.980572, -0.192558)), ((0.445475, 
            -0.056214, 2.341382), (0.055871, -0.979511, -0.193484)), ((0.426566, 
            -0.057252, 2.34118), (0.055871, -0.979511, -0.193484)), ((0.409977, 
            -0.058291, 2.340978), (0.067641, -0.978511, -0.194782)), ((0.393081, 
            -0.05944, 2.340755), (0.069854, -0.978298, -0.195074)), ((0.37611, 
            -0.060605, 2.340529), (0.069673, -0.978316, -0.195046)), ((0.359274, 
            -0.06217, 2.340224), (0.107249, -0.973584, -0.201572)), ((0.342394, 
            -0.064127, 2.339844), (0.122423, -0.971186, -0.204476)), ((0.334478, 
            -0.064831, 2.339707), (0.0, -0.984214, -0.176982)), ((0.280898, -0.064831, 
            2.339707), (0.0, -0.984214, -0.176982)), ((0.270879, -0.064815, 2.33971), (
            -0.001899, -0.984319, -0.176387)), ((0.259259, -0.064528, 2.339766), (
            -0.039614, -0.985704, -0.163761)), ((0.249259, -0.063923, 2.339884), (
            -0.074539, -0.985673, -0.151302)), ((0.240694, -0.063036, 2.340056), (
            -0.123318, -0.983434, -0.132851)), ((0.233654, -0.061866, 2.340283), (
            -0.193825, -0.975453, -0.104511)), ((0.228223, -0.060429, 2.340563), (
            -0.302046, -0.951537, -0.05784)), ((0.224547, -0.058739, 2.340891), (
            -0.499339, -0.865632, 0.036618)), ((0.222638, -0.055911, 2.341441), (
            -0.876188, -0.394143, 0.277389)), ((0.529442, -0.037388, 2.253371), (
            -0.051907, -0.980043, -0.191889)), ((0.519487, -0.036963, 2.253453), (
            -0.022329, -0.981195, -0.191722)), ((0.507386, -0.036814, 2.253482), (
            -0.004709, -0.981439, -0.191717)), ((0.49584, -0.036888, 2.253468), (
            0.011647, -0.981353, -0.19186)), ((0.482531, -0.037164, 2.253414), (
            0.025441, -0.981037, -0.192143)), ((0.467672, -0.037653, 2.253319), (
            0.03742, -0.980572, -0.192558)), ((0.443431, -0.038888, 2.253079), (
            0.055871, -0.979511, -0.193484)), ((0.419603, -0.040196, 2.252825), (
            0.055871, -0.979511, -0.193484)), ((0.397508, -0.041553, 2.252561), (
            0.067641, -0.978511, -0.194782)), ((0.37853, -0.042843, 2.25231), (
            0.069854, -0.978298, -0.195074)), ((0.359467, -0.044151, 2.252056), (
            0.069673, -0.978316, -0.195046)), ((0.340554, -0.045907, 2.251715), (
            0.107249, -0.973584, -0.201572)), ((0.321594, -0.048104, 2.251288), (
            0.122423, -0.971186, -0.204476)), ((0.310697, -0.048904, 2.251132), (0.0, 
            -0.984214, -0.176982)), ((0.253593, -0.048904, 2.251132), (0.0, -0.984214, 
            -0.176982)), ((0.241246, -0.048886, 2.251136), (-0.001899, -0.984319, 
            -0.176387)), ((0.228186, -0.048565, 2.251198), (-0.039614, -0.985704, 
            -0.163761)), ((0.216947, -0.047887, 2.25133), (-0.074539, -0.985673, 
            -0.151302)), ((0.20732, -0.046891, 2.251523), (-0.123318, -0.983434, 
            -0.132851)), ((0.199406, -0.045579, 2.251779), (-0.193825, -0.975453, 
            -0.104511)), ((0.193299, -0.043965, 2.252092), (-0.302046, -0.951537, 
            -0.05784)), ((0.189162, -0.042067, 2.252461), (-0.499339, -0.865632, 
            0.036618)), ((0.18701, -0.038898, 2.253077), (-0.876188, -0.394143, 
            0.277389)), ((0.186622, -0.03674, 2.253497), (-0.921246, 0.06129, 
            0.384121)), ((0.187408, -0.03379, 2.25407), (-0.87222, 0.272656, 
            0.406069)), ((0.189468, -0.029783, 2.254849), (-0.778269, 0.476044, 
            0.409487)), ((0.192879, -0.025654, 2.255652), (-0.674525, 0.621555, 
            0.398353)), ((0.197639, -0.021538, 2.256452), (-0.575154, 0.723888, 
            0.38103)), ((0.203714, -0.017563, 2.257224), (-0.480724, 0.799158, 
            0.360904)), ((0.211084, -0.013845, 2.257947), (-0.393325, 0.854192, 
            0.340076)), ((0.219725, -0.010473, 2.258603), (-0.314905, 0.893516, 
            0.3201)), ((0.22959, -0.007331, 2.259213), (-0.289604, 0.904327, 
            0.313563)), ((0.241482, -0.003957, 2.259869), (-0.254874, 0.917694, 
            0.304757)), ((0.253593, -0.001531, 2.260341), (0.0, 0.972107, 0.23454)), ((
            0.310697, -0.001531, 2.260341), (0.0, 0.972107, 0.23454)), ((0.32279, 
            -0.001959, 2.260257), (0.055897, 0.973641, 0.22113)), ((0.343311, 
            -0.003364, 2.259984), (0.10349, 0.972127, 0.210377)), ((0.360238, 
            -0.005237, 2.25962), (0.127788, 0.970349, 0.205166)), ((0.377311, 
            -0.007523, 2.259176), (0.14876, 0.968233, 0.200986)), ((0.394373, 
            -0.010152, 2.258665), (0.167001, 0.965944, 0.197645)), ((0.416776, 
            -0.014167, 2.257884), (0.191841, 0.962164, 0.193486)), ((0.440553, 
            -0.018729, 2.256998), (0.191841, 0.962164, 0.193486)), ((0.459175, 
            -0.022373, 2.256289), (0.205116, 0.95972, 0.192005)), ((0.473579, 
            -0.025345, 2.255712), (0.206856, 0.959379, 0.19184)), ((0.486846, 
            -0.028101, 2.255176), (0.207327, 0.959285, 0.191803)), ((0.498766, 
            -0.030586, 2.254693), (0.208596, 0.959025, 0.191724)), ((0.509158, 
            -0.032768, 2.254269), (0.210345, 0.958661, 0.191637)), ((0.520198, 
            -0.035087, 2.253818), (0.208059, 0.959142, 0.191723)), ((0.52961, 
            -0.036906, 2.253464), (0.187215, 0.963321, 0.192257)), ((0.146234, 
            -0.017202, 2.153516), (-0.921246, 0.06129, 0.384121)), ((0.147115, 
            -0.013888, 2.15416), (-0.87222, 0.272656, 0.406069)), ((0.149429, 
            -0.009383, 2.155036), (-0.778269, 0.476044, 0.409487)), ((0.153262, 
            -0.004741, 2.155938), (-0.674525, 0.621555, 0.398353)), ((0.15861, 
            -0.000114, 2.156838), (-0.575154, 0.723888, 0.38103)), ((0.165438, 
            0.004356, 2.157707), (-0.480724, 0.799158, 0.360904)), ((0.173722, 
            0.008536, 2.158519), (-0.393325, 0.854192, 0.340076)), ((0.183435, 
            0.012328, 2.159256), (-0.314905, 0.893516, 0.3201)), ((0.194523, 0.015859, 
            2.159943), (-0.289604, 0.904327, 0.313563)), ((0.207886, 0.019652, 
            2.16068), (-0.254874, 0.917694, 0.304757)), ((0.222507, 0.022386, 
            2.161211), (0.0, 0.972107, 0.23454)), ((0.283599, 0.022386, 2.161211), (
            0.0, 0.972107, 0.23454)), ((0.299281, 0.021906, 2.161118), (0.055897, 
            0.973641, 0.22113)), ((0.322362, 0.020328, 2.160811), (0.10349, 0.972127, 
            0.210377)), ((0.341392, 0.018223, 2.160402), (0.127788, 0.970349, 
            0.205166)), ((0.360586, 0.015654, 2.159903), (0.14876, 0.968233, 
            0.200986)), ((0.379768, 0.012699, 2.159328), (0.167001, 0.965944, 
            0.197645)), ((0.407978, 0.007606, 2.158338), (0.191841, 0.962164, 
            0.193486)), ((0.438166, 0.001813, 2.157212), (0.191841, 0.962164, 
            0.193486)), ((0.452625, -0.001041, 2.156657), (0.205116, 0.95972, 
            0.192005)), ((0.468821, -0.004382, 2.156008), (0.206856, 0.959379, 
            0.19184)), ((0.483737, -0.007481, 2.155406), (0.207327, 0.959285, 
            0.191803)), ((0.497142, -0.010275, 2.154862), (0.208596, 0.959025, 
            0.191724)), ((0.508827, -0.012728, 2.154386), (0.210345, 0.958661, 
            0.191637)), ((0.521267, -0.015342, 2.153877), (0.208059, 0.959142, 
            0.191723)), ((0.531854, -0.017387, 2.15348), (0.187215, 0.963321, 
            0.192257)), ((0.531666, -0.017927, 2.153375), (-0.051907, -0.980043, 
            -0.191889)), ((0.52047, -0.017448, 2.153468), (-0.022329, -0.981195, 
            -0.191722)), ((0.506833, -0.017281, 2.153501), (-0.004709, -0.981439, 
            -0.191717)), ((0.493849, -0.017365, 2.153484), (0.011647, -0.981353, 
            -0.19186)), ((0.478884, -0.017674, 2.153424), (0.025441, -0.981037, 
            -0.192143)), ((0.459047, -0.01834, 2.153295), (0.03742, -0.980572, 
            -0.192558)), ((0.441387, -0.019258, 2.153116), (0.055871, -0.979511, 
            -0.193484)), ((0.411163, -0.020919, 2.152794), (0.055871, -0.979511, 
            -0.193484)), ((0.38329, -0.022611, 2.152465), (0.067641, -0.978511, 
            -0.194782)), ((0.361953, -0.024061, 2.152183), (0.069854, -0.978298, 
            -0.195074)), ((0.340522, -0.025532, 2.151897), (0.069673, -0.978316, 
            -0.195046)), ((0.31926, -0.027507, 2.151513), (0.107249, -0.973584, 
            -0.201572)), ((0.297944, -0.029977, 2.151033), (0.122423, -0.971186, 
            -0.204476)), ((0.283599, -0.030872, 2.150859), (0.0, -0.984214, 
            -0.176982)), ((0.222507, -0.030872, 2.150859), (0.0, -0.984214, 
            -0.176982)), ((0.207621, -0.030852, 2.150863), (-0.001899, -0.984319, 
            -0.176387)), ((0.192941, -0.030491, 2.150933), (-0.039614, -0.985704, 
            -0.163761)), ((0.180309, -0.029728, 2.151081), (-0.074539, -0.985673, 
            -0.151302)), ((0.169489, -0.028608, 2.151299), (-0.123318, -0.983434, 
            -0.132851)), ((0.160594, -0.027132, 2.151586), (-0.193825, -0.975453, 
            -0.104511)), ((0.153732, -0.025317, 2.151939), (-0.302046, -0.951537, 
            -0.05784)), ((0.149085, -0.023183, 2.152353), (-0.499339, -0.865632, 
            0.036618)), ((0.146668, -0.019617, 2.153047), (-0.876188, -0.394143, 
            0.277389)), ((0.533357, -0.003653, 2.08), (-0.051251, -0.980093, 
            -0.191811)), ((0.521236, -0.003128, 2.080042), (-0.022048, -0.981247, 
            -0.19149)), ((0.506371, -0.002944, 2.080059), (-0.00465, -0.981495, 
            -0.19143)), ((0.492322, -0.003037, 2.080049), (0.0115, -0.981406, 
            -0.191598)), ((0.476133, -0.003376, 2.080023), (0.025119, -0.981078, 
            -0.191977)), ((0.458062, -0.003979, 2.07998), (0.036944, -0.980588, 
            -0.192566)), ((0.44016, -0.004856, 2.079912), (0.055216, -0.979493, 
            -0.193767)), ((0.404411, -0.006839, 2.079763), (0.055109, -0.979381, 
            -0.194358)), ((0.372759, -0.008782, 2.079668), (0.066745, -0.978281, 
            -0.196243)), ((0.349691, -0.01037, 2.079569), (0.068922, -0.977971, 
            -0.197034)), ((0.326523, -0.01198, 2.07947), (0.068739, -0.977886, 
            -0.197517)), ((0.303539, -0.014143, 2.07932), (0.105703, -0.973066, 
            -0.204866)), ((0.280494, -0.016848, 2.079146), (0.120592, -0.970497, 
            -0.208791)), ((0.26355, -0.017825, 2.079145), (0.0, -0.983527, -0.180763)), 
            ((0.199533, -0.017825, 2.079145), (0.0, -0.983527, -0.180763)), ((0.182881, 
            -0.017803, 2.079147), (-0.001868, -0.983637, -0.18015)), ((0.167028, 
            -0.017406, 2.079183), (-0.039247, -0.985181, -0.166968)), ((0.153382, 
            -0.016569, 2.079242), (-0.073977, -0.985317, -0.153878)), ((0.141696, 
            -0.015343, 2.079326), (-0.122696, -0.983282, -0.134545)), ((0.132092, 
            -0.013724, 2.079433), (-0.193559, -0.975457, -0.104969)), ((0.124686, 
            -0.011734, 2.079563), (-0.303316, -0.95122, -0.05638)), ((0.119679, 
            -0.009393, 2.079713), (-0.505936, -0.861534, 0.042276)), ((0.117104, 
            -0.005474, 2.079991), (-0.887331, -0.358133, 0.290489)), ((0.116775, 
            -0.001931, 2.080174), (-0.914895, 0.109414, 0.388581)), ((0.117527, 
            0.000766, 2.080198), (-0.858275, 0.316949, 0.403619)), ((0.120023, 
            0.005703, 2.080496), (-0.760791, 0.511954, 0.398873)), ((0.124164, 
            0.010788, 2.080805), (-0.657481, 0.649404, 0.382091)), ((0.129942, 
            0.015855, 2.081116), (-0.560265, 0.745443, 0.361133)), ((0.137319, 
            0.020748, 2.08142), (-0.468525, 0.815872, 0.338875)), ((0.14627, 0.025323, 
            2.081707), (-0.383785, 0.867286, 0.317056)), ((0.156764, 0.029473, 
            2.081971), (-0.307708, 0.903978, 0.296882)), ((0.168743, 0.033339, 
            2.082206), (-0.283224, 0.91445, 0.289078)), ((0.183166, 0.037491, 
            2.082447), (-0.249498, 0.927212, 0.279334)), ((0.199533, 0.040483, 
            2.082723), (0.0, 0.97607, 0.217455)), ((0.26355, 0.040483, 2.082723), (0.0, 
            0.97607, 0.217455)), ((0.281918, 0.03996, 2.082722), (0.055107, 0.977006, 
            0.205969)), ((0.306906, 0.038233, 2.082628), (0.102133, 0.975058, 
            0.197057)), ((0.327477, 0.035929, 2.082498), (0.12617, 0.973036, 
            0.193084)), ((0.348226, 0.033117, 2.082335), (0.146928, 0.970706, 
            0.190108)), ((0.368963, 0.029882, 2.082144), (0.164987, 0.968222, 
            0.187948)), ((0.400935, 0.023997, 2.081929), (0.189443, 0.964183, 
            0.18564)), ((0.436733, 0.017072, 2.081346), (0.190042, 0.963722, 
            0.187414)), ((0.447739, 0.014841, 2.081227), (0.202669, 0.961281, 
            0.186717)), ((0.465254, 0.011182, 2.080997), (0.204365, 0.960759, 
            0.187555)), ((0.481388, 0.007789, 2.080782), (0.204806, 0.960488, 
            0.188459)), ((0.495888, 0.00473, 2.080587), (0.206038, 0.96007, 0.189246)), 
            ((0.508532, 0.002043, 2.080414), (0.207745, 0.959566, 0.189934)), ((
            0.519148, -0.000218, 2.080265), (0.205467, 0.95991, 0.190675)), ((0.533558, 
            -0.003064, 2.080001), (0.184843, 0.963823, 0.192039)), ((0.106671, 
            0.001029, 2.033333), (-0.987308, 0.143028, 0.069039)), ((0.107632, 
            0.004756, 2.033333), (-0.927883, 0.367154, 0.065055)), ((0.110181, 
            0.009848, 2.033333), (-0.817545, 0.572957, 0.057801)), ((0.114418, 
            0.015097, 2.033333), (-0.700158, 0.71222, 0.050206)), ((0.120337, 0.020332, 
            2.033333), (-0.591304, 0.805285, 0.043305)), ((0.127899, 0.02539, 
            2.033333), (-0.490483, 0.870662, 0.03707)), ((0.137078, 0.030123, 
            2.033333), (-0.398994, 0.91641, 0.031577)), ((0.147843, 0.034417, 
            2.033333), (-0.318094, 0.947678, 0.026888)), ((0.160137, 0.03841, 
            2.033333), (-0.292195, 0.95602, 0.025454)), ((0.174922, 0.042692, 
            2.033333), (-0.256786, 0.96618, 0.023606)), ((0.191823, 0.045837, 
            2.033333), (0.0, 0.99994, 0.010933)), ((0.256893, 0.045837, 2.033333), (
            0.0, 0.99994, 0.010933)), ((0.276275, 0.045312, 2.033333), (0.055922, 
            0.998393, 0.009145)), ((0.301978, 0.043534, 2.033333), (0.10351, 0.994597, 
            0.007916)), ((0.323103, 0.041159, 2.033333), (0.127802, 0.991772, 
            0.007403)), ((0.344411, 0.038259, 2.033333), (0.148768, 0.988847, 
            0.007063)), ((0.365708, 0.034921, 2.033333), (0.167005, 0.985932, 
            0.006859)), ((0.398236, 0.028931, 2.033333), (0.19184, 0.981403, 
            0.006707)), ((0.438597, 0.021041, 2.033333), (0.19184, 0.981403, 
            0.006707)), ((0.452733, 0.018105, 2.033333), (0.205114, 0.978714, 
            0.006823)), ((0.470289, 0.014404, 2.033333), (0.206854, 0.978348, 
            0.006846)), ((0.48634, 0.011006, 2.033333), (0.207325, 0.978248, 
            0.006854)), ((0.500631, 0.007964, 2.033333), (0.208594, 0.977978, 
            0.006881)), ((0.512955, 0.005321, 2.033333), (0.210343, 0.977603, 
            0.006925)), ((0.523146, 0.003144, 2.033333), (0.208057, 0.978093, 
            0.006861)), ((0.534939, 0.000826, 2.033333), (0.187214, 0.982299, 
            0.006234)), ((0.534736, 0.000236, 2.033333), (-0.051907, -0.99865, 
            -0.001728)), ((0.522269, 0.000782, 2.033333), (-0.022329, -0.99975, 
            -0.000839)), ((0.511132, 0.000955, 2.033333), (-0.004709, -0.999989, 
            -0.000353)), ((0.497475, 0.000936, 2.033333), (0.011647, -0.999932, 
            4.8e-05)), ((0.481536, 0.00067, 2.033333), (0.025441, -0.999676, 
            0.000337)), ((0.463576, 0.000136, 2.033333), (0.03742, -0.999299, 
            0.000538)), ((0.445484, -0.00064, 2.033333), (0.055871, -0.998438, 
            0.000763)), ((0.401834, -0.003082, 2.033333), (0.055871, -0.998438, 
            0.000763)), ((0.369593, -0.005072, 2.033333), (0.067641, -0.997709, 
            0.000722)), ((0.345901, -0.006713, 2.033333), (0.069854, -0.997557, 
            0.000701)), ((0.322108, -0.008376, 2.033333), (0.069674, -0.99757, 
            0.000704)), ((0.298506, -0.01062, 2.033333), (0.107257, -0.994231, 
            -7e-05)), ((0.274838, -0.013418, 2.033333), (0.122438, -0.992476, 
            -0.000469)), ((0.256893, -0.014398, 2.033333), (0.0, -0.999994, 0.003434)), 
            ((0.191823, -0.014398, 2.033333), (0.0, -0.999994, 0.003434)), ((0.174623, 
            -0.014375, 2.033333), (-0.001899, -0.999992, 0.003531)), ((0.158359, 
            -0.01396, 2.033333), (-0.039628, -0.999199, 0.005608)), ((0.144356, 
            -0.013092, 2.033333), (-0.074595, -0.997185, 0.007657)), ((0.132366, 
            -0.011821, 2.033333), (-0.123513, -0.992286, 0.010673)), ((0.122516, 
            -0.010147, 2.033333), (-0.194483, -0.980788, 0.015231)), ((0.114923, 
            -0.00809, 2.033333), (-0.304348, -0.952295, 0.02251)), ((0.109792, 
            -0.005673, 2.033333), (-0.509417, -0.85975, 0.036388)), ((0.107137, 
            -0.001614, 2.033333), (-0.932433, -0.355412, 0.065202)), ((0.537292, 
            0.000241, 1.953333), (-0.051907, -0.99865, -0.001728)), ((0.524593, 
            0.000797, 1.953333), (-0.022329, -0.99975, -0.000839)), ((0.513239, 
            0.000973, 1.953333), (-0.004709, -0.999989, -0.000353)), ((0.499317, 
            0.000954, 1.953333), (0.011647, -0.999932, 4.8e-05)), ((0.48307, 0.000682, 
            1.953333), (0.025441, -0.999676, 0.000337)), ((0.464765, 0.000138, 
            1.953333), (0.03742, -0.999299, 0.000538)), ((0.446053, -0.000669, 
            1.953333), (0.055871, -0.998438, 0.000763)), ((0.401673, -0.003152, 
            1.953333), (0.055871, -0.998438, 0.000763)), ((0.369055, -0.005166, 
            1.953333), (0.067641, -0.997709, 0.000722)), )+\
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.344912, -0.006839, 
            1.953333), (0.069854, -0.997557, 0.000701)), ((0.320666, -0.008534, 
            1.953333), (0.069674, -0.99757, 0.000704)), ((0.296612, -0.010819, 
            1.953333), (0.107257, -0.994231, -7e-05)), ((0.272493, -0.013669, 
            1.953333), (0.122438, -0.992476, -0.000469)), ((0.253913, -0.014672, 
            1.953333), (0.0, -0.999994, 0.003434)), ((0.188079, -0.014672, 1.953333), (
            0.0, -0.999994, 0.003434)), ((0.170357, -0.014649, 1.953333), (-0.001899, 
            -0.999992, 0.003531)), ((0.153778, -0.014227, 1.953333), (-0.039628, 
            -0.999199, 0.005608)), ((0.139505, -0.013343, 1.953333), (-0.074595, 
            -0.997185, 0.007657)), ((0.127283, -0.012049, 1.953333), (-0.123513, 
            -0.992286, 0.010673)), ((0.117241, -0.010344, 1.953333), (-0.194483, 
            -0.980788, 0.015231)), ((0.1095, -0.008248, 1.953333), (-0.304348, 
            -0.952295, 0.02251)), ((0.104267, -0.005785, 1.953333), (-0.509417, 
            -0.85975, 0.036388)), ((0.101557, -0.001652, 1.953333), (-0.932433, 
            -0.355412, 0.065202)), ((0.10108, 0.001053, 1.953333), (-0.987308, 
            0.143028, 0.069039)), ((0.102062, 0.004854, 1.953333), (-0.927883, 
            0.367154, 0.065055)), ((0.104662, 0.010044, 1.953333), (-0.817545, 
            0.572957, 0.057801)), ((0.108983, 0.015393, 1.953333), (-0.700158, 0.71222, 
            0.050206)), ((0.115017, 0.020728, 1.953333), (-0.591304, 0.805285, 
            0.043305)), ((0.122726, 0.025883, 1.953333), (-0.490483, 0.870662, 
            0.03707)), ((0.132083, 0.030705, 1.953333), (-0.398994, 0.91641, 
            0.031577)), ((0.143056, 0.03508, 1.953333), (-0.318094, 0.947678, 
            0.026888)), ((0.155586, 0.039149, 1.953333), (-0.292195, 0.95602, 
            0.025454)), ((0.170661, 0.043514, 1.953333), (-0.256786, 0.96618, 
            0.023606)), ((0.188079, 0.046712, 1.953333), (0.0, 0.99994, 0.010933)), ((
            0.253913, 0.046712, 1.953333), (0.0, 0.99994, 0.010933)), ((0.273964, 
            0.046174, 1.953333), (0.055922, 0.998393, 0.009145)), ((0.300147, 0.044361, 
            1.953333), (0.10351, 0.994597, 0.007916)), ((0.321675, 0.04194, 1.953333), 
            (0.127802, 0.991772, 0.007403)), ((0.34339, 0.038984, 1.953333), (0.148768, 
            0.988847, 0.007063)), ((0.365093, 0.035582, 1.953333), (0.167005, 0.985932, 
            0.006859)), ((0.398015, 0.029521, 1.953333), (0.19184, 0.981403, 
            0.006707)), ((0.439051, 0.021499, 1.953333), (0.19184, 0.981403, 
            0.006707)), ((0.453717, 0.018456, 1.953333), (0.205114, 0.978714, 
            0.006823)), ((0.471609, 0.014685, 1.953333), (0.206854, 0.978348, 
            0.006846)), ((0.487969, 0.011221, 1.953333), (0.207325, 0.978248, 
            0.006854)), ((0.502536, 0.008121, 1.953333), (0.208594, 0.977978, 
            0.006881)), ((0.515099, 0.005427, 1.953333), (0.210343, 0.977603, 
            0.006925)), ((0.525487, 0.003207, 1.953333), (0.208057, 0.978093, 
            0.006861)), ((0.5375, 0.000846, 1.953333), (0.187214, 0.982299, 0.006234)), 
            ((0.094091, 0.001078, 1.853333), (-0.987308, 0.143028, 0.069039)), ((
            0.095096, 0.004967, 1.853333), (-0.927883, 0.367154, 0.065055)), ((
            0.097757, 0.010279, 1.853333), (-0.817545, 0.572957, 0.057801)), ((
            0.102179, 0.015753, 1.853333), (-0.700158, 0.71222, 0.050206)), ((0.108354, 
            0.021213, 1.853333), (-0.591304, 0.805285, 0.043305)), ((0.116243, 
            0.026488, 1.853333), (-0.490483, 0.870662, 0.03707)), ((0.125819, 0.031423, 
            1.853333), (-0.398994, 0.91641, 0.031577)), ((0.137049, 0.035901, 
            1.853333), (-0.318094, 0.947678, 0.026888)), ((0.149872, 0.040065, 
            1.853333), (-0.292195, 0.95602, 0.025454)), ((0.1653, 0.044532, 1.853333), 
            (-0.256786, 0.96618, 0.023606)), ((0.183318, 0.047805, 1.853333), (0.0, 
            0.99994, 0.010933)), ((0.250109, 0.047805, 1.853333), (0.0, 0.99994, 
            0.010933)), ((0.27102, 0.047255, 1.853333), (0.055922, 0.998393, 
            0.009145)), ((0.297816, 0.0454, 1.853333), (0.10351, 0.994597, 0.007916)), 
            ((0.319848, 0.042922, 1.853333), (0.127802, 0.991772, 0.007403)), ((
            0.342071, 0.039897, 1.853333), (0.148768, 0.988847, 0.007063)), ((0.364283, 
            0.036415, 1.853333), (0.167005, 0.985932, 0.006859)), ((0.397699, 0.030266, 
            1.853333), (0.19184, 0.981403, 0.006707)), ((0.439699, 0.022056, 1.853333), 
            (0.19184, 0.981403, 0.006707)), ((0.454983, 0.018888, 1.853333), (0.205114, 
            0.978714, 0.006823)), ((0.473294, 0.015028, 1.853333), (0.206854, 0.978348, 
            0.006846)), ((0.490037, 0.011483, 1.853333), (0.207325, 0.978248, 
            0.006854)), ((0.504945, 0.008311, 1.853333), (0.208594, 0.977978, 
            0.006881)), ((0.517801, 0.005554, 1.853333), (0.210343, 0.977603, 
            0.006925)), ((0.528433, 0.003282, 1.853333), (0.208057, 0.978093, 
            0.006861)), ((0.540727, 0.000866, 1.853333), (0.187214, 0.982299, 
            0.006234)), ((0.540514, 0.000247, 1.853333), (-0.051907, -0.99865, 
            -0.001728)), ((0.527517, 0.000816, 1.853333), (-0.022329, -0.99975, 
            -0.000839)), ((0.515898, 0.000996, 1.853333), (-0.004709, -0.999989, 
            -0.000353)), ((0.50165, 0.000976, 1.853333), (0.011647, -0.999932, 
            4.8e-05)), ((0.485024, 0.000698, 1.853333), (0.025441, -0.999676, 
            0.000337)), ((0.46629, 0.000141, 1.853333), (0.03742, -0.999299, 
            0.000538)), ((0.446866, -0.0007, 1.853333), (0.055871, -0.998438, 
            0.000763)), ((0.401443, -0.003242, 1.853333), (0.055871, -0.998438, 
            0.000763)), ((0.368337, -0.005287, 1.853333), (0.067641, -0.997709, 
            0.000722)), ((0.343629, -0.006999, 1.853333), (0.069854, -0.997557, 
            0.000701)), ((0.318815, -0.008733, 1.853333), (0.069674, -0.99757, 
            0.000704)), ((0.294198, -0.011072, 1.853333), (0.107257, -0.994231, 
            -7e-05)), ((0.269515, -0.01399, 1.853333), (0.122438, -0.992476, 
            -0.000469)), ((0.250109, -0.015016, 1.853333), (0.0, -0.999994, 0.003434)), 
            ((0.183318, -0.015016, 1.853333), (0.0, -0.999994, 0.003434)), ((0.164988, 
            -0.014992, 1.853333), (-0.001899, -0.999992, 0.003531)), ((0.148021, 
            -0.01456, 1.853333), (-0.039628, -0.999199, 0.005608)), ((0.133414, 
            -0.013655, 1.853333), (-0.074595, -0.997185, 0.007657)), ((0.120907, 
            -0.012331, 1.853333), (-0.123513, -0.992286, 0.010673)), ((0.11063, 
            -0.010586, 1.853333), (-0.194483, -0.980788, 0.015231)), ((0.102707, 
            -0.008441, 1.853333), (-0.304348, -0.952295, 0.02251)), ((0.097352, 
            -0.005921, 1.853333), (-0.509417, -0.85975, 0.036388)), ((0.094579, 
            -0.00169, 1.853333), (-0.932433, -0.355412, 0.065202)), ((0.543736, 
            0.000253, 1.753333), (-0.051907, -0.99865, -0.001728)), ((0.530442, 
            0.000835, 1.753333), (-0.022329, -0.99975, -0.000839)), ((0.518557, 
            0.001019, 1.753333), (-0.004709, -0.999989, -0.000353)), ((0.503983, 
            0.000998, 1.753333), (0.011647, -0.999932, 4.8e-05)), ((0.486977, 0.000714, 
            1.753333), (0.025441, -0.999676, 0.000337)), ((0.467815, 0.000144, 
            1.753333), (0.03742, -0.999299, 0.000538)), ((0.447679, -0.000731, 
            1.753333), (0.055871, -0.998438, 0.000763)), ((0.401213, -0.003331, 
            1.753333), (0.055871, -0.998438, 0.000763)), ((0.367619, -0.005408, 
            1.753333), (0.067641, -0.997709, 0.000722)), ((0.342346, -0.007159, 
            1.753333), (0.069854, -0.997557, 0.000701)), ((0.316964, -0.008933, 
            1.753333), (0.069674, -0.99757, 0.000704)), ((0.291784, -0.011325, 
            1.753333), (0.107257, -0.994231, -7e-05)), ((0.266536, -0.01431, 1.753333), 
            (0.122438, -0.992476, -0.000469)), ((0.246304, -0.015359, 1.753333), (0.0, 
            -0.999994, 0.003434)), ((0.178558, -0.015359, 1.753333), (0.0, -0.999994, 
            0.003434)), ((0.15962, -0.015335, 1.753333), (-0.001899, -0.999992, 
            0.003531)), ((0.142265, -0.014893, 1.753333), (-0.039628, -0.999199, 
            0.005608)), ((0.127323, -0.013968, 1.753333), (-0.074595, -0.997185, 
            0.007657)), ((0.11453, -0.012613, 1.753333), (-0.123513, -0.992286, 
            0.010673)), ((0.104019, -0.010828, 1.753333), (-0.194483, -0.980788, 
            0.015231)), ((0.095915, -0.008634, 1.753333), (-0.304348, -0.952295, 
            0.02251)), ((0.090438, -0.006056, 1.753333), (-0.509417, -0.85975, 
            0.036388)), ((0.087601, -0.001729, 1.753333), (-0.932433, -0.355412, 
            0.065202)), ((0.087102, 0.001102, 1.753333), (-0.987308, 0.143028, 
            0.069039)), ((0.088129, 0.005081, 1.753333), (-0.927883, 0.367154, 
            0.065055)), ((0.090852, 0.010513, 1.753333), (-0.817545, 0.572957, 
            0.057801)), ((0.095374, 0.016113, 1.753333), (-0.700158, 0.71222, 
            0.050206)), ((0.101691, 0.021698, 1.753333), (-0.591304, 0.805285, 
            0.043305)), ((0.10976, 0.027094, 1.753333), (-0.490483, 0.870662, 
            0.03707)), ((0.119555, 0.032142, 1.753333), (-0.398994, 0.91641, 
            0.031577)), ((0.131042, 0.036721, 1.753333), (-0.318094, 0.947678, 
            0.026888)), ((0.144158, 0.040981, 1.753333), (-0.292195, 0.95602, 
            0.025454)), ((0.159938, 0.04555, 1.753333), (-0.256786, 0.96618, 
            0.023606)), ((0.178558, 0.048899, 1.753333), (0.0, 0.99994, 0.010933)), ((
            0.246304, 0.048899, 1.753333), (0.0, 0.99994, 0.010933)), ((0.268076, 
            0.048336, 1.753333), (0.055922, 0.998393, 0.009145)), ((0.295485, 0.046438, 
            1.753333), (0.10351, 0.994597, 0.007916)), ((0.318021, 0.043904, 1.753333), 
            (0.127802, 0.991772, 0.007403)), ((0.340753, 0.04081, 1.753333), (0.148768, 
            0.988847, 0.007063)), ((0.363472, 0.037248, 1.753333), (0.167005, 0.985932, 
            0.006859)), ((0.397384, 0.031011, 1.753333), (0.19184, 0.981403, 
            0.006707)), ((0.440348, 0.022612, 1.753333), (0.19184, 0.981403, 
            0.006707)), ((0.456249, 0.01932, 1.753333), (0.205114, 0.978714, 
            0.006823)), ((0.474979, 0.015372, 1.753333), (0.206854, 0.978348, 
            0.006846)), ((0.492104, 0.011746, 1.753333), (0.207325, 0.978248, 
            0.006854)), ((0.507353, 0.008501, 1.753333), (0.208594, 0.977978, 
            0.006881)), ((0.520503, 0.00568, 1.753333), (0.210343, 0.977603, 
            0.006925)), ((0.531378, 0.003357, 1.753333), (0.208057, 0.978093, 
            0.006861)), ((0.543953, 0.000885, 1.753333), (0.187214, 0.982299, 
            0.006234)), ((0.080113, 0.001127, 1.653333), (-0.987308, 0.143028, 
            0.069039)), ((0.081163, 0.005194, 1.653333), (-0.927883, 0.367154, 
            0.065055)), ((0.083946, 0.010748, 1.653333), (-0.817545, 0.572957, 
            0.057801)), ((0.08857, 0.016473, 1.653333), (-0.700158, 0.71222, 
            0.050206)), ((0.095028, 0.022182, 1.653333), (-0.591304, 0.805285, 
            0.043305)), ((0.103277, 0.027699, 1.653333), (-0.490483, 0.870662, 
            0.03707)), ((0.113291, 0.03286, 1.653333), (-0.398994, 0.91641, 0.031577)), 
            ((0.125034, 0.037542, 1.653333), (-0.318094, 0.947678, 0.026888)), ((
            0.138444, 0.041897, 1.653333), (-0.292195, 0.95602, 0.025454)), ((0.154577, 
            0.046568, 1.653333), (-0.256786, 0.96618, 0.023606)), ((0.173798, 0.049992, 
            1.653333), (0.0, 0.99994, 0.010933)), ((0.2425, 0.049992, 1.653333), (0.0, 
            0.99994, 0.010933)), ((0.265132, 0.049417, 1.653333), (0.055922, 0.998393, 
            0.009145)), ((0.293154, 0.047477, 1.653333), (0.10351, 0.994597, 
            0.007916)), ((0.316194, 0.044886, 1.653333), (0.127802, 0.991772, 
            0.007403)), ((0.339434, 0.041722, 1.653333), (0.148768, 0.988847, 
            0.007063)), ((0.362661, 0.038081, 1.653333), (0.167005, 0.985932, 
            0.006859)), ((0.397068, 0.031756, 1.653333), (0.19184, 0.981403, 
            0.006707)), ((0.440996, 0.023169, 1.653333), (0.19184, 0.981403, 
            0.006707)), ((0.457515, 0.019752, 1.653333), (0.205114, 0.978714, 
            0.006823)), ((0.476664, 0.015715, 1.653333), (0.206854, 0.978348, 
            0.006846)), ((0.494172, 0.012008, 1.653333), (0.207325, 0.978248, 
            0.006854)), ((0.509762, 0.008691, 1.653333), (0.208594, 0.977978, 
            0.006881)), ((0.523206, 0.005807, 1.653333), (0.210343, 0.977603, 
            0.006925)), ((0.534323, 0.003432, 1.653333), (0.208057, 0.978093, 
            0.006861)), ((0.54718, 0.000905, 1.653333), (0.187214, 0.982299, 
            0.006234)), ((0.546958, 0.000258, 1.653333), (-0.051907, -0.99865, 
            -0.001728)), ((0.533366, 0.000853, 1.653333), (-0.022329, -0.99975, 
            -0.000839)), ((0.521215, 0.001042, 1.653333), (-0.004709, -0.999989, 
            -0.000353)), ((0.506317, 0.001021, 1.653333), (0.011647, -0.999932, 
            4.8e-05)), ((0.48893, 0.00073, 1.653333), (0.025441, -0.999676, 0.000337)), 
            ((0.469339, 0.000148, 1.653333), (0.03742, -0.999299, 0.000538)), ((
            0.448492, -0.000762, 1.653333), (0.055871, -0.998438, 0.000763)), ((
            0.400984, -0.00342, 1.653333), (0.055871, -0.998438, 0.000763)), ((
            0.366901, -0.005529, 1.653333), (0.067641, -0.997709, 0.000722)), ((
            0.341062, -0.007319, 1.653333), (0.069854, -0.997557, 0.000701)), ((
            0.315113, -0.009133, 1.653333), (0.069674, -0.99757, 0.000704)), ((
            0.289371, -0.011579, 1.653333), (0.107257, -0.994231, -7e-05)), ((0.263558, 
            -0.01463, 1.653333), (0.122438, -0.992476, -0.000469)), ((0.2425, 
            -0.015703, 1.653333), (0.0, -0.999994, 0.003434)), ((0.173798, -0.015703, 
            1.653333), (0.0, -0.999994, 0.003434)), ((0.154251, -0.015678, 1.653333), (
            -0.001899, -0.999992, 0.003531)), ((0.136508, -0.015226, 1.653333), (
            -0.039628, -0.999199, 0.005608)), ((0.121233, -0.01428, 1.653333), (
            -0.074595, -0.997185, 0.007657)), ((0.108154, -0.012895, 1.653333), (
            -0.123513, -0.992286, 0.010673)), ((0.097407, -0.01107, 1.653333), (
            -0.194483, -0.980788, 0.015231)), ((0.089123, -0.008827, 1.653333), (
            -0.304348, -0.952295, 0.02251)), ((0.083523, -0.006191, 1.653333), (
            -0.509417, -0.85975, 0.036388)), ((0.080623, -0.001767, 1.653333), (
            -0.932433, -0.355412, 0.065202)), ((0.55018, 0.000264, 1.553333), (
            -0.051907, -0.99865, -0.001728)), ((0.536291, 0.000872, 1.553333), (
            -0.022329, -0.99975, -0.000839)), ((0.523874, 0.001064, 1.553333), (
            -0.004709, -0.999989, -0.000353)), ((0.50865, 0.001043, 1.553333), (
            0.011647, -0.999932, 4.8e-05)), ((0.490883, 0.000746, 1.553333), (0.025441, 
            -0.999676, 0.000337)), ((0.470864, 0.000151, 1.553333), (0.03742, 
            -0.999299, 0.000538)), ((0.449305, -0.000793, 1.553333), (0.055871, 
            -0.998438, 0.000763)), ((0.400754, -0.00351, 1.553333), (0.055871, 
            -0.998438, 0.000763)), ((0.366183, -0.00565, 1.553333), (0.067641, 
            -0.997709, 0.000722)), ((0.339779, -0.007479, 1.553333), (0.069854, 
            -0.997557, 0.000701)), ((0.313262, -0.009333, 1.553333), (0.069674, 
            -0.99757, 0.000704)), ((0.286957, -0.011832, 1.553333), (0.107257, 
            -0.994231, -7e-05)), ((0.26058, -0.01495, 1.553333), (0.122438, -0.992476, 
            -0.000469)), ((0.238695, -0.016046, 1.553333), (0.0, -0.999994, 0.003434)), 
            ((0.169038, -0.016046, 1.553333), (0.0, -0.999994, 0.003434)), ((0.148882, 
            -0.016021, 1.553333), (-0.001899, -0.999992, 0.003531)), ((0.130751, 
            -0.015559, 1.553333), (-0.039628, -0.999199, 0.005608)), ((0.115142, 
            -0.014592, 1.553333), (-0.074595, -0.997185, 0.007657)), ((0.101777, 
            -0.013177, 1.553333), (-0.123513, -0.992286, 0.010673)), ((0.090796, 
            -0.011312, 1.553333), (-0.194483, -0.980788, 0.015231)), ((0.08233, 
            -0.00902, 1.553333), (-0.304348, -0.952295, 0.02251)), ((0.076608, 
            -0.006326, 1.553333), (-0.509417, -0.85975, 0.036388)), ((0.073645, 
            -0.001805, 1.553333), (-0.932433, -0.355412, 0.065202)), ((0.073124, 
            0.001151, 1.553333), (-0.987308, 0.143028, 0.069039)), ((0.074197, 
            0.005307, 1.553333), (-0.927883, 0.367154, 0.065055)), ((0.077041, 
            0.010983, 1.553333), (-0.817545, 0.572957, 0.057801)), ((0.081765, 
            0.016833, 1.553333), (-0.700158, 0.71222, 0.050206)), ((0.088364, 0.022667, 
            1.553333), (-0.591304, 0.805285, 0.043305)), ((0.096794, 0.028305, 
            1.553333), (-0.490483, 0.870662, 0.03707)), ((0.107027, 0.033578, 
            1.553333), (-0.398994, 0.91641, 0.031577)), ((0.119027, 0.038363, 
            1.553333), (-0.318094, 0.947678, 0.026888)), ((0.13273, 0.042813, 
            1.553333), (-0.292195, 0.95602, 0.025454)), ((0.149215, 0.047587, 
            1.553333), (-0.256786, 0.96618, 0.023606)), ((0.169038, 0.051085, 
            1.553333), (0.0, 0.99994, 0.010933)), ((0.238695, 0.051085, 1.553333), (
            0.0, 0.99994, 0.010933)), ((0.262187, 0.050498, 1.553333), (0.055922, 
            0.998393, 0.009145)), ((0.290823, 0.048515, 1.553333), (0.10351, 0.994597, 
            0.007916)), ((0.314367, 0.045868, 1.553333), (0.127802, 0.991772, 
            0.007403)), ((0.338115, 0.042635, 1.553333), (0.148768, 0.988847, 
            0.007063)), ((0.36185, 0.038914, 1.553333), (0.167005, 0.985932, 
            0.006859)), ((0.396752, 0.032501, 1.553333), (0.19184, 0.981403, 
            0.006707)), ((0.441645, 0.023725, 1.553333), (0.19184, 0.981403, 
            0.006707)), ((0.458781, 0.020183, 1.553333), (0.205114, 0.978714, 
            0.006823)), ((0.478348, 0.016059, 1.553333), (0.206854, 0.978348, 
            0.006846)), ((0.49624, 0.01227, 1.553333), (0.207325, 0.978248, 0.006854)), 
            ((0.51217, 0.00888, 1.553333), (0.208594, 0.977978, 0.006881)), ((0.525908, 
            0.005934, 1.553333), (0.210343, 0.977603, 0.006925)), ((0.537269, 0.003507, 
            1.553333), (0.208057, 0.978093, 0.006861)), ((0.550407, 0.000925, 
            1.553333), (0.187214, 0.982299, 0.006234)), ((0.066135, 0.001176, 
            1.453333), (-0.987308, 0.143028, 0.069039)), ((0.067231, 0.005421, 
            1.453333), (-0.927883, 0.367154, 0.065055)), ((0.070135, 0.011218, 
            1.453333), (-0.817545, 0.572957, 0.057801)), ((0.074961, 0.017193, 
            1.453333), (-0.700158, 0.71222, 0.050206)), ((0.081701, 0.023152, 
            1.453333), (-0.591304, 0.805285, 0.043305)), ((0.090311, 0.02891, 
            1.453333), (-0.490483, 0.870662, 0.03707)), ((0.100762, 0.034297, 
            1.453333), (-0.398994, 0.91641, 0.031577)), ((0.11302, 0.039184, 1.453333), 
            (-0.318094, 0.947678, 0.026888)), ((0.127016, 0.043729, 1.453333), (
            -0.292195, 0.95602, 0.025454)), ((0.143854, 0.048605, 1.453333), (
            -0.256786, 0.96618, 0.023606)), ((0.164278, 0.052179, 1.453333), (0.0, 
            0.99994, 0.010933)), ((0.234891, 0.052179, 1.453333), (0.0, 0.99994, 
            0.010933)), ((0.259243, 0.051578, 1.453333), (0.055922, 0.998393, 
            0.009145)), ((0.288492, 0.049554, 1.453333), (0.10351, 0.994597, 
            0.007916)), ((0.312539, 0.046849, 1.453333), (0.127802, 0.991772, 
            0.007403)), ((0.336796, 0.043548, 1.453333), (0.148768, 0.988847, 
            0.007063)), ((0.36104, 0.039747, 1.453333), (0.167005, 0.985932, 
            0.006859)), ((0.396437, 0.033246, 1.453333), (0.19184, 0.981403, 
            0.006707)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.442293, 
            0.024282, 1.453333), (0.19184, 0.981403, 0.006707)), ((0.460047, 0.020615, 
            1.453333), (0.205114, 0.978714, 0.006823)), ((0.480033, 0.016402, 
            1.453333), (0.206854, 0.978348, 0.006846)), ((0.498307, 0.012533, 
            1.453333), (0.207325, 0.978248, 0.006854)), ((0.514578, 0.00907, 1.453333), 
            (0.208594, 0.977978, 0.006881)), ((0.528611, 0.006061, 1.453333), (
            0.210343, 0.977603, 0.006925)), ((0.540214, 0.003582, 1.453333), (0.208057, 
            0.978093, 0.006861)), ((0.553634, 0.000944, 1.453333), (0.187214, 0.982299, 
            0.006234)), ((0.553402, 0.00027, 1.453333), (-0.051907, -0.99865, 
            -0.001728)), ((0.539215, 0.000891, 1.453333), (-0.022329, -0.99975, 
            -0.000839)), ((0.526533, 0.001087, 1.453333), (-0.004709, -0.999989, 
            -0.000353)), ((0.510983, 0.001065, 1.453333), (0.011647, -0.999932, 
            4.8e-05)), ((0.492836, 0.000762, 1.453333), (0.025441, -0.999676, 
            0.000337)), ((0.472389, 0.000154, 1.453333), (0.03742, -0.999299, 
            0.000538)), ((0.450118, -0.000824, 1.453333), (0.055871, -0.998438, 
            0.000763)), ((0.400524, -0.003599, 1.453333), (0.055871, -0.998438, 
            0.000763)), ((0.365465, -0.005771, 1.453333), (0.067641, -0.997709, 
            0.000722)), ((0.338496, -0.007639, 1.453333), (0.069854, -0.997557, 
            0.000701)), ((0.311411, -0.009533, 1.453333), (0.069674, -0.99757, 
            0.000704)), ((0.284543, -0.012086, 1.453333), (0.107257, -0.994231, 
            -7e-05)), ((0.257601, -0.01527, 1.453333), (0.122438, -0.992476, 
            -0.000469)), ((0.234891, -0.01639, 1.453333), (0.0, -0.999994, 0.003434)), 
            ((0.164278, -0.01639, 1.453333), (0.0, -0.999994, 0.003434)), ((0.143513, 
            -0.016363, 1.453333), (-0.001899, -0.999992, 0.003531)), ((0.124995, 
            -0.015892, 1.453333), (-0.039628, -0.999199, 0.005608)), ((0.109052, 
            -0.014905, 1.453333), (-0.074595, -0.997185, 0.007657)), ((0.095401, 
            -0.013459, 1.453333), (-0.123513, -0.992286, 0.010673)), ((0.084184, 
            -0.011554, 1.453333), (-0.194483, -0.980788, 0.015231)), ((0.075538, 
            -0.009213, 1.453333), (-0.304348, -0.952295, 0.02251)), ((0.069693, 
            -0.006462, 1.453333), (-0.509417, -0.85975, 0.036388)), ((0.066667, 
            -0.001844, 1.453333), (-0.932433, -0.355412, 0.065202)), ), name=
            'Skin-tip')
            
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.29133, -0.022657, 
            2.347905), ), ((0.31633, -0.022657, 2.347905), ), ((0.266219, -0.001531, 
            2.260341), ), ((0.291219, -0.001531, 2.260341), ), ((0.237792, 0.022386, 
            2.161211), ), ((0.262792, 0.022386, 2.161211), ), ((0.216945, 0.040483, 
            2.082723), ), ((0.241945, 0.040483, 2.082723), ), ((0.21016, 0.045837, 
            2.033333), ), ((0.23516, 0.045837, 2.033333), ), ((0.206862, 0.046712, 
            1.953333), ), ((0.231862, 0.046712, 1.953333), ), ((0.202739, 0.047805, 
            1.853333), ), ((0.227739, 0.047805, 1.853333), ), ((0.198616, 0.048899, 
            1.753333), ), ((0.223616, 0.048899, 1.753333), ), ((0.194493, 0.049992, 
            1.653333), ), ((0.219493, 0.049992, 1.653333), ), ((0.19037, 0.051085, 
            1.553333), ), ((0.21537, 0.051085, 1.553333), ), ((0.186247, 0.052179, 
            1.453333), ), ((0.211247, 0.052179, 1.453333), ), ((0.182124, 0.053272, 
            1.353333), ), ((0.207124, 0.053272, 1.353333), ), ((0.178001, 0.054365, 
            1.253333), ), ((0.203001, 0.054365, 1.253333), ), ((0.173878, 0.055459, 
            1.153333), ), ((0.198878, 0.055459, 1.153333), ), ((0.169755, 0.056552, 
            1.053333), ), ((0.194755, 0.056552, 1.053333), ), ((0.165632, 0.057646, 
            0.953333), ), ((0.190632, 0.057646, 0.953333), ), ((0.161509, 0.058739, 
            0.853333), ), ((0.186509, 0.058739, 0.853333), ), ((0.157386, 0.059832, 
            0.753333), ), ((0.182386, 0.059832, 0.753333), ), ((0.161596, 0.060926, 
            0.653333), ), ((0.187971, 0.060561, 0.686667), ), ((0.154637, 0.060561, 
            0.686667), ), ((0.178263, 0.060926, 0.653333), ), ((0.150239, 0.061727, 
            0.58), ), ((0.175239, 0.061727, 0.58), ), ((0.14804, 0.062311, 0.526667), 
            ), ((0.17304, 0.062311, 0.526667), ), ((0.147353, 0.062493, 0.476667), ), (
            (0.172353, 0.062493, 0.476667), ), ((0.147353, 0.062493, 0.42), ), ((
            0.172353, 0.062493, 0.42), ), ((0.147353, 0.062493, 0.333333), ), ((
            0.172353, 0.062493, 0.333333), ), ((0.147353, 0.062493, 0.233333), ), ((
            0.172353, 0.062493, 0.233333), ), ((0.147353, 0.062493, 0.133333), ), ((
            0.172353, 0.062493, 0.133333), ), ((0.147353, 0.062493, 0.033333), ), ((
            0.172353, 0.062493, 0.033333), ), ), name='top flange')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.31633, -0.064831, 
            2.339707), ), ((0.29133, -0.064831, 2.339707), ), ((0.291219, -0.048904, 
            2.251132), ), ((0.266219, -0.048904, 2.251132), ), ((0.262792, -0.030872, 
            2.150859), ), ((0.237792, -0.030872, 2.150859), ), ((0.241945, -0.017825, 
            2.079145), ), ((0.216945, -0.017825, 2.079145), ), ((0.23516, -0.014398, 
            2.033333), ), ((0.21016, -0.014398, 2.033333), ), ((0.231862, -0.014672, 
            1.953333), ), ((0.206862, -0.014672, 1.953333), ), ((0.227739, -0.015016, 
            1.853333), ), ((0.202739, -0.015016, 1.853333), ), ((0.223616, -0.015359, 
            1.753333), ), ((0.198616, -0.015359, 1.753333), ), ((0.219493, -0.015703, 
            1.653333), ), ((0.194493, -0.015703, 1.653333), ), ((0.21537, -0.016046, 
            1.553333), ), ((0.19037, -0.016046, 1.553333), ), ((0.211247, -0.01639, 
            1.453333), ), ((0.186247, -0.01639, 1.453333), ), ((0.207124, -0.016733, 
            1.353333), ), ((0.182124, -0.016733, 1.353333), ), ((0.203001, -0.017076, 
            1.253333), ), ((0.178001, -0.017076, 1.253333), ), ((0.198878, -0.01742, 
            1.153333), ), ((0.173878, -0.01742, 1.153333), ), ((0.194755, -0.017763, 
            1.053333), ), ((0.169755, -0.017763, 1.053333), ), ((0.190632, -0.018107, 
            0.953333), ), ((0.165632, -0.018107, 0.953333), ), ((0.186509, -0.01845, 
            0.853333), ), ((0.161509, -0.01845, 0.853333), ), ((0.182386, -0.018794, 
            0.753333), ), ((0.157386, -0.018794, 0.753333), ), ((0.154637, -0.019023, 
            0.686667), ), ((0.178263, -0.019137, 0.653333), ), ((0.187971, -0.019023, 
            0.686667), ), ((0.161596, -0.019137, 0.653333), ), ((0.175239, -0.019389, 
            0.58), ), ((0.150239, -0.019389, 0.58), ), ((0.17304, -0.019572, 0.526667), 
            ), ((0.14804, -0.019572, 0.526667), ), ((0.180687, -0.019629, 0.476667), ), 
            ((0.155687, -0.019629, 0.476667), ), ((0.180687, -0.019629, 0.42), ), ((
            0.155687, -0.019629, 0.42), ), ((0.180687, -0.019629, 0.333333), ), ((
            0.155687, -0.019629, 0.333333), ), ((0.180687, -0.019629, 0.233333), ), ((
            0.155687, -0.019629, 0.233333), ), ((0.180687, -0.019629, 0.133333), ), ((
            0.155687, -0.019629, 0.133333), ), ((0.180687, -0.019629, 0.033333), ), ((
            0.155687, -0.019629, 0.033333), ), ), name='bottom flange')

     
    createcomponentsets()
     
    def createsurfaces():    
        #load surfaces
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-1', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.583833, 0.000321, 
            0.133333), ), ((0.566831, 0.001066, 0.133333), ), ((0.55165, 0.001302, 
            0.133333), ), ((0.533034, 0.001276, 0.133333), ), ((0.511306, 0.000915, 
            0.133333), ), ((0.486823, 0.000187, 0.133333), ), ((0.457924, -0.001108, 
            0.133333), ), ((0.416659, -0.003417, 0.133333), ), ((0.369303, -0.006193, 
            0.133333), ), ((0.33714, -0.008398, 0.133333), ), ((0.304702, -0.010667, 
            0.133333), ), ((0.272393, -0.01333, 0.133333), ), ((0.240235, -0.016966, 
            0.133333), ), ((0.208773, -0.019629, 0.133333), ), ((0.180687, -0.019629, 
            0.133333), ), ((0.155687, -0.019629, 0.133333), ), ((0.129143, -0.019629, 
            0.133333), ), ((0.101105, -0.019614, 0.133333), ), ((0.077593, -0.019307, 
            0.133333), ), ((0.057635, -0.018302, 0.133333), ), ((0.040354, -0.016754, 
            0.133333), ), ((0.025938, -0.014654, 0.133333), ), ((0.014545, -0.012018, 
            0.133333), ), ((0.006398, -0.008883, 0.133333), ), ((0.001672, -0.004386, 
            0.133333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-2', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.583833, 0.000321, 
            0.476667), ), ((0.566831, 0.001066, 0.476667), ), ((0.55165, 0.001302, 
            0.476667), ), ((0.533034, 0.001276, 0.476667), ), ((0.511306, 0.000915, 
            0.476667), ), ((0.486823, 0.000187, 0.476667), ), ((0.457924, -0.001108, 
            0.476667), ), ((0.416659, -0.003417, 0.476667), ), ((0.369303, -0.006193, 
            0.476667), ), ((0.33714, -0.008398, 0.476667), ), ((0.304702, -0.010667, 
            0.476667), ), ((0.272393, -0.01333, 0.476667), ), ((0.240235, -0.016966, 
            0.476667), ), ((0.208773, -0.019629, 0.476667), ), ((0.180687, -0.019629, 
            0.476667), ), ((0.155687, -0.019629, 0.476667), ), ((0.129143, -0.019629, 
            0.476667), ), ((0.101105, -0.019614, 0.476667), ), ((0.077593, -0.019307, 
            0.476667), ), ((0.057635, -0.018302, 0.476667), ), ((0.040354, -0.016754, 
            0.476667), ), ((0.025938, -0.014654, 0.476667), ), ((0.014545, -0.012018, 
            0.476667), ), ((0.006398, -0.008883, 0.476667), ), ((0.001672, -0.004386, 
            0.476667), ), ((0.583833, 0.000321, 0.42), ), ((0.566831, 0.001066, 0.42), 
            ), ((0.55165, 0.001302, 0.42), ), ((0.533034, 0.001276, 0.42), ), ((
            0.511306, 0.000915, 0.42), ), ((0.486823, 0.000187, 0.42), ), ((0.457924, 
            -0.001108, 0.42), ), ((0.416659, -0.003417, 0.42), ), ((0.369303, 
            -0.006193, 0.42), ), ((0.33714, -0.008398, 0.42), ), ((0.304702, -0.010667, 
            0.42), ), ((0.272393, -0.01333, 0.42), ), ((0.240235, -0.016966, 0.42), ), 
            ((0.208773, -0.019629, 0.42), ), ((0.180687, -0.019629, 0.42), ), ((
            0.155687, -0.019629, 0.42), ), ((0.129143, -0.019629, 0.42), ), ((0.101105, 
            -0.019614, 0.42), ), ((0.077593, -0.019307, 0.42), ), ((0.057635, 
            -0.018302, 0.42), ), ((0.040354, -0.016754, 0.42), ), ((0.025938, 
            -0.014654, 0.42), ), ((0.014545, -0.012018, 0.42), ), ((0.006398, 
            -0.008883, 0.42), ), ((0.001672, -0.004386, 0.42), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-3', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.575956, 0.000309, 
            0.753333), ), ((0.559687, 0.001021, 0.753333), ), ((0.545146, 0.001247, 
            0.753333), ), ((0.527316, 0.001222, 0.753333), ), ((0.506507, 0.000874, 
            0.753333), ), ((0.483062, 0.000177, 0.753333), ), ((0.45581, -0.00104, 
            0.753333), ), ((0.398916, -0.004224, 0.753333), ), ((0.360437, -0.006618, 
            0.753333), ), ((0.329513, -0.008761, 0.753333), ), ((0.298456, -0.010932, 
            0.753333), ), ((0.267646, -0.013859, 0.753333), ), ((0.236753, -0.017511, 
            0.753333), ), ((0.20826, -0.018794, 0.753333), ), ((0.182386, -0.018794, 
            0.753333), ), ((0.157386, -0.018794, 0.753333), ), ((0.130956, -0.018794, 
            0.753333), ), ((0.105932, -0.018764, 0.753333), ), ((0.084699, -0.018223, 
            0.753333), ), ((0.066418, -0.01709, 0.753333), ), ((0.050765, -0.015433, 
            0.753333), ), ((0.037905, -0.013248, 0.753333), ), ((0.02799, -0.010563, 
            0.753333), ), ((0.02129, -0.007409, 0.753333), ), ((0.017821, -0.002112, 
            0.753333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-4', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.56629, 0.000292, 
            1.053333), ), ((0.550914, 0.000965, 1.053333), ), ((0.537169, 0.001178, 
            1.053333), ), ((0.520316, 0.001155, 1.053333), ), ((0.500648, 0.000826, 
            1.053333), ), ((0.478488, 0.000167, 1.053333), ), ((0.453371, -0.000947, 
            1.053333), ), ((0.399605, -0.003956, 1.053333), ), ((0.362592, -0.006255, 
            1.053333), ), ((0.333363, -0.00828, 1.053333), ), ((0.304008, -0.010332, 
            1.053333), ), ((0.274888, -0.013099, 1.053333), ), ((0.245688, -0.016551, 
            1.053333), ), ((0.219673, -0.017763, 1.053333), ), ((0.194755, -0.017763, 
            1.053333), ), ((0.169755, -0.017763, 1.053333), ), ((0.145237, -0.017763, 
            1.053333), ), ((0.122039, -0.017735, 1.053333), ), ((0.101968, -0.017224, 
            1.053333), ), ((0.08469, -0.016154, 1.053333), ), ((0.069895, -0.014587, 
            1.053333), ), ((0.057739, -0.012522, 1.053333), ), ((0.048368, -0.009985, 
            1.053333), ), ((0.042034, -0.007003, 1.053333), ), ((0.038755, -0.001997, 
            1.053333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-5', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.556624, 0.000275, 
            1.353333), ), ((0.54214, 0.000909, 1.353333), ), ((0.529192, 0.00111, 
            1.353333), ), ((0.513316, 0.001088, 1.353333), ), ((0.494789, 0.000778, 
            1.353333), ), ((0.473913, 0.000157, 1.353333), ), ((0.450932, -0.000854, 
            1.353333), ), ((0.400294, -0.003688, 1.353333), ), ((0.364746, -0.005892, 
            1.353333), ), ((0.337213, -0.0078, 1.353333), ), ((0.309561, -0.009733, 
            1.353333), ), ((0.282129, -0.012339, 1.353333), ), ((0.254623, -0.01559, 
            1.353333), ), ((0.231087, -0.016733, 1.353333), ), ((0.207124, -0.016733, 
            1.353333), ), ((0.182124, -0.016733, 1.353333), ), ((0.159517, -0.016733, 
            1.353333), ), ((0.138145, -0.016706, 1.353333), ), ((0.119238, -0.016225, 
            1.353333), ), ((0.102961, -0.015217, 1.353333), ), ((0.089024, -0.013741, 
            1.353333), ), ((0.077573, -0.011796, 1.353333), ), ((0.068745, -0.009406, 
            1.353333), ), ((0.062779, -0.006597, 1.353333), ), ((0.059689, -0.001882, 
            1.353333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-6', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.546958, 0.000258, 
            1.653333), ), ((0.533366, 0.000853, 1.653333), ), ((0.521215, 0.001042, 
            1.653333), ), ((0.506317, 0.001021, 1.653333), ), ((0.48893, 0.00073, 
            1.653333), ), ((0.469339, 0.000148, 1.653333), ), ((0.448492, -0.000762, 
            1.653333), ), ((0.400984, -0.00342, 1.653333), ), ((0.366901, -0.005529, 
            1.653333), ), ((0.341062, -0.007319, 1.653333), ), ((0.315113, -0.009133, 
            1.653333), ), ((0.289371, -0.011579, 1.653333), ), ((0.263558, -0.01463, 
            1.653333), ), ((0.2425, -0.015703, 1.653333), ), ((0.219493, -0.015703, 
            1.653333), ), ((0.194493, -0.015703, 1.653333), ), ((0.173798, -0.015703, 
            1.653333), ), ((0.154251, -0.015678, 1.653333), ), ((0.136508, -0.015226, 
            1.653333), ), ((0.121233, -0.01428, 1.653333), ), ((0.108154, -0.012895, 
            1.653333), ), ((0.097407, -0.01107, 1.653333), ), ((0.089123, -0.008827, 
            1.653333), ), ((0.083523, -0.006191, 1.653333), ), ((0.080623, -0.001767, 
            1.653333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-7', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.537292, 0.000241, 
            1.953333), ), ((0.524593, 0.000797, 1.953333), ), ((0.513239, 0.000973, 
            1.953333), ), ((0.499317, 0.000954, 1.953333), ), ((0.48307, 0.000682, 
            1.953333), ), ((0.464765, 0.000138, 1.953333), ), ((0.446053, -0.000669, 
            1.953333), ), ((0.401673, -0.003152, 1.953333), ), ((0.369055, -0.005166, 
            1.953333), ), ((0.344912, -0.006839, 1.953333), ), ((0.320666, -0.008534, 
            1.953333), ), ((0.296612, -0.010819, 1.953333), ), ((0.272493, -0.013669, 
            1.953333), ), ((0.253913, -0.014672, 1.953333), ), ((0.231862, -0.014672, 
            1.953333), ), ((0.206862, -0.014672, 1.953333), ), ((0.188079, -0.014672, 
            1.953333), ), ((0.170357, -0.014649, 1.953333), ), ((0.153778, -0.014227, 
            1.953333), ), ((0.139505, -0.013343, 1.953333), ), ((0.127283, -0.012049, 
            1.953333), ), ((0.117241, -0.010344, 1.953333), ), ((0.1095, -0.008248, 
            1.953333), ), ((0.104267, -0.005785, 1.953333), ), ((0.101557, -0.001652, 
            1.953333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='load-8', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.529442, -0.037388, 
            2.253371), ), ((0.519487, -0.036963, 2.253453), ), ((0.507386, -0.036814, 
            2.253482), ), ((0.49584, -0.036888, 2.253468), ), ((0.482531, -0.037164, 
            2.253414), ), ((0.467672, -0.037653, 2.253319), ), ((0.443431, -0.038888, 
            2.253079), ), ((0.419603, -0.040196, 2.252825), ), ((0.397508, -0.041553, 
            2.252561), ), ((0.37853, -0.042843, 2.25231), ), ((0.359467, -0.044151, 
            2.252056), ), ((0.340554, -0.045907, 2.251715), ), ((0.321594, -0.048104, 
            2.251288), ), ((0.310697, -0.048904, 2.251132), ), ((0.291219, -0.048904, 
            2.251132), ), ((0.266219, -0.048904, 2.251132), ), ((0.253593, -0.048904, 
            2.251132), ), ((0.241246, -0.048886, 2.251136), ), ((0.228186, -0.048565, 
            2.251198), ), ((0.216947, -0.047887, 2.25133), ), ((0.20732, -0.046891, 
            2.251523), ), ((0.199406, -0.045579, 2.251779), ), ((0.193299, -0.043965, 
            2.252092), ), ((0.189162, -0.042067, 2.252461), ), ((0.18701, -0.038898, 
            2.253077), ), ))
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.150239, 0.061727, 
            0.58), ), ((0.175239, 0.061727, 0.58), ), ((0.14804, 0.062311, 0.526667), 
            ), ((0.17304, 0.062311, 0.526667), ), ((0.147353, 0.062493, 0.476667), ), (
            (0.172353, 0.062493, 0.476667), ), ((0.147353, 0.062493, 0.42), ), ((
            0.172353, 0.062493, 0.42), ), ((0.147353, 0.062493, 0.333333), ), ((
            0.172353, 0.062493, 0.333333), ), ((0.147353, 0.062493, 0.233333), ), ((
            0.172353, 0.062493, 0.233333), ), ((0.147353, 0.062493, 0.133333), ), ((
            0.172353, 0.062493, 0.133333), ), ((0.147353, 0.062493, 0.033333), ), ((
            0.172353, 0.062493, 0.033333), ), ), name='Spar cap-top root')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.182124, 0.053272, 
            1.353333), ), ((0.207124, 0.053272, 1.353333), ), ((0.178001, 0.054365, 
            1.253333), ), ((0.203001, 0.054365, 1.253333), ), ((0.173878, 0.055459, 
            1.153333), ), ((0.198878, 0.055459, 1.153333), ), ((0.169755, 0.056552, 
            1.053333), ), ((0.194755, 0.056552, 1.053333), ), ((0.165632, 0.057646, 
            0.953333), ), ((0.190632, 0.057646, 0.953333), ), ((0.161509, 0.058739, 
            0.853333), ), ((0.186509, 0.058739, 0.853333), ), ((0.157386, 0.059832, 
            0.753333), ), ((0.182386, 0.059832, 0.753333), ), ((0.161596, 0.060926, 
            0.653333), ), ((0.187971, 0.060561, 0.686667), ), ((0.154637, 0.060561, 
            0.686667), ), ((0.178263, 0.060926, 0.653333), ), ), name=
            'Spar cap-top mid')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.29133, -0.022657, 
            2.347905), ), ((0.31633, -0.022657, 2.347905), ), ((0.266219, -0.001531, 
            2.260341), ), ((0.291219, -0.001531, 2.260341), ), ((0.237792, 0.022386, 
            2.161211), ), ((0.262792, 0.022386, 2.161211), ), ((0.216945, 0.040483, 
            2.082723), ), ((0.241945, 0.040483, 2.082723), ), ((0.21016, 0.045837, 
            2.033333), ), ((0.23516, 0.045837, 2.033333), ), ((0.206862, 0.046712, 
            1.953333), ), ((0.231862, 0.046712, 1.953333), ), ((0.202739, 0.047805, 
            1.853333), ), ((0.227739, 0.047805, 1.853333), ), ((0.198616, 0.048899, 
            1.753333), ), ((0.223616, 0.048899, 1.753333), ), ((0.194493, 0.049992, 
            1.653333), ), ((0.219493, 0.049992, 1.653333), ), ((0.19037, 0.051085, 
            1.553333), ), ((0.21537, 0.051085, 1.553333), ), ((0.186247, 0.052179, 
            1.453333), ), ((0.211247, 0.052179, 1.453333), ), ), name=
            'Spar cap-top tip')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.175239, -0.019389, 
            0.58), ), ((0.150239, -0.019389, 0.58), ), ((0.17304, -0.019572, 0.526667), 
            ), ((0.14804, -0.019572, 0.526667), ), ((0.180687, -0.019629, 0.476667), ), 
            ((0.155687, -0.019629, 0.476667), ), ((0.180687, -0.019629, 0.42), ), ((
            0.155687, -0.019629, 0.42), ), ((0.180687, -0.019629, 0.333333), ), ((
            0.155687, -0.019629, 0.333333), ), ((0.180687, -0.019629, 0.233333), ), ((
            0.155687, -0.019629, 0.233333), ), ((0.180687, -0.019629, 0.133333), ), ((
            0.155687, -0.019629, 0.133333), ), ((0.180687, -0.019629, 0.033333), ), ((
            0.155687, -0.019629, 0.033333), ), ), name='Spar cap bottom root')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.207124, -0.016733, 
            1.353333), ), ((0.182124, -0.016733, 1.353333), ), ((0.203001, -0.017076, 
            1.253333), ), ((0.178001, -0.017076, 1.253333), ), ((0.198878, -0.01742, 
            1.153333), ), ((0.173878, -0.01742, 1.153333), ), ((0.194755, -0.017763, 
            1.053333), ), ((0.169755, -0.017763, 1.053333), ), ((0.190632, -0.018107, 
            0.953333), ), ((0.165632, -0.018107, 0.953333), ), ((0.186509, -0.01845, 
            0.853333), ), ((0.161509, -0.01845, 0.853333), ), ((0.182386, -0.018794, 
            0.753333), ), ((0.157386, -0.018794, 0.753333), ), ((0.154637, -0.019023, 
            0.686667), ), ((0.178263, -0.019137, 0.653333), ), ((0.187971, -0.019023, 
            0.686667), ), ((0.161596, -0.019137, 0.653333), ), ), name=
            'Spar cap bottom mid')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.31633, -0.064831, 
            2.339707), ), ((0.29133, -0.064831, 2.339707), ), ((0.291219, -0.048904, 
            2.251132), ), ((0.266219, -0.048904, 2.251132), ), ((0.262792, -0.030872, 
            2.150859), ), ((0.237792, -0.030872, 2.150859), ), ((0.241945, -0.017825, 
            2.079145), ), ((0.216945, -0.017825, 2.079145), ), ((0.23516, -0.014398, 
            2.033333), ), ((0.21016, -0.014398, 2.033333), ), ((0.231862, -0.014672, 
            1.953333), ), ((0.206862, -0.014672, 1.953333), ), ((0.227739, -0.015016, 
            1.853333), ), ((0.202739, -0.015016, 1.853333), ), ((0.223616, -0.015359, 
            1.753333), ), ((0.198616, -0.015359, 1.753333), ), ((0.219493, -0.015703, 
            1.653333), ), ((0.194493, -0.015703, 1.653333), ), ((0.21537, -0.016046, 
            1.553333), ), ((0.19037, -0.016046, 1.553333), ), ((0.211247, -0.01639, 
            1.453333), ), ((0.186247, -0.01639, 1.453333), ), ), name=
            'Spar cap-bottom tip')


        mdb.models['Model-1'].parts['Wing'].Surface(name='Rear Spar-after boom', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.435, 
            -0.042558, 2.344037), ), ((0.435, -0.024931, 2.255792), ), ((0.435, 
            -0.004961, 2.155896), ), ((0.435, 0.009824, 2.081346), ), ((0.435, 
            0.014054, 2.033333), ), ((0.435, 0.014347, 1.953333), ), ((0.435, 0.014777, 
            1.853333), ), ((0.435, 0.015207, 1.753333), ), ((0.435, 0.015637, 
            1.653333), ), ((0.435, 0.016067, 1.553333), ), ((0.435, 0.016498, 
            1.453333), ), ((0.435, 0.016928, 1.353333), ), ((0.435, 0.017358, 
            1.253333), ), ((0.435, 0.017788, 1.153333), ), ((0.435, 0.018218, 
            1.053333), ), ((0.435, 0.018648, 0.953333), ), ((0.435, 0.019078, 
            0.853333), ), ((0.435, 0.019508, 0.753333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='Rear Spar-before boom', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.435, 
            0.019939, 0.653333), ), ((0.435, 0.020288, 0.58), ), ((0.435, 0.020526, 
            0.526667), ), ((0.435, 0.009125, 0.476667), ), ((0.435, 0.009125, 0.42), ), 
            ((0.435, 0.009125, 0.333333), ), ((0.435, 0.009125, 0.233333), ), ((0.435, 
            0.009125, 0.133333), ), ((0.435, 0.009125, 0.033333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='Ribs', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.34655, -0.046229, 
            2.388294), ), ((0.289491, -0.047265, 2.388092), ), ((0.352472, -0.058794, 
            2.385851), ), ((0.294081, -0.059437, 2.385726), ), ((0.303648, -0.059437, 
            2.385726), ), ((0.328648, -0.059437, 2.385726), ), ((0.444492, -0.053008, 
            2.386976), ), ((0.3298, -0.031958, 2.323612), ), ((0.265066, -0.047082, 
            2.320672), ), ((0.351569, -0.05393, 2.319341), ), ((0.273874, -0.046441, 
            2.320796), ), ((0.285171, -0.046441, 2.320796), ), ((0.310171, -0.046441, 
            2.320796), ), ((0.44294, -0.03973, 2.322101), ), ((0.278264, 0.011952, 
            2.12459), ), ((0.19825, -0.007137, 2.12088), ), ((0.307681, -0.016034, 
            2.11915), ), ((0.211702, -0.006455, 2.121012), ), ((0.228316, -0.006455, 
            2.121012), ), ((0.253316, -0.006455, 2.121012), ), ((0.438166, 0.001126, 
            2.122486), ), ((0.261409, 0.025841, 2.02), ), ((0.176078, 0.004989, 2.02), 
            ), ((0.293487, -0.0048, 2.02), ), ((0.191146, 0.005699, 2.02), ), ((
            0.209611, 0.005699, 2.02), ), ((0.234611, 0.005699, 2.02), ), ((0.436863, 
            0.01378, 2.02), ), ((0.257923, 0.026455, 1.92), ), ((0.170761, 0.00512, 
            1.92), ), ((0.290917, -0.00492, 1.92), ), ((0.186386, 0.005834, 1.92), ), (
            (0.205488, 0.005834, 1.92), ), ((0.230488, 0.005834, 1.92), ), ((0.437187, 
            0.014147, 1.92), ), ((0.254437, 0.02707, 1.82), ), ((0.165445, 0.00525, 
            1.82), ), ((0.288348, -0.005039, 1.82), ), ((0.181626, 0.00597, 1.82), ), (
            (0.201365, 0.00597, 1.82), ), ((0.226365, 0.00597, 1.82), ), ((0.437512, 
            0.014514, 1.82), ), ((0.250951, 0.027684, 1.72), ), ((0.160128, 0.005381, 
            1.72), ), ((0.285778, -0.005159, 1.72), ), ((0.176865, 0.006105, 1.72), ), 
            ((0.197242, 0.006105, 1.72), ), ((0.222242, 0.006105, 1.72), ), ((0.437836, 
            0.014881, 1.72), ), ((0.247466, 0.028299, 1.62), ), ((0.154811, 0.005511, 
            1.62), ), ((0.283209, -0.005279, 1.62), ), ((0.172105, 0.006241, 1.62), ), 
            ((0.193119, 0.006241, 1.62), ), ((0.218119, 0.006241, 1.62), ), ((0.43816, 
            0.015247, 1.62), ), ((0.24398, 0.028913, 1.52), ), ((0.149495, 0.005642, 
            1.52), ), ((0.28064, -0.005398, 1.52), ), ((0.167345, 0.006376, 1.52), ), (
            (0.188996, 0.006376, 1.52), ), ((0.213996, 0.006376, 1.52), ), ((0.438485, 
            0.015614, 1.52), ), ((0.240494, 0.029527, 1.42), ), ((0.144178, 0.005772, 
            1.42), ), ((0.27807, -0.005518, 1.42), ), ((0.162585, 0.006512, 1.42), ), (
            (0.184873, 0.006512, 1.42), ), ((0.209873, 0.006512, 1.42), ), ((0.438809, 
            0.015981, 1.42), ), ((0.237008, 0.030142, 1.32), ), ((0.138861, 0.005902, 
            1.32), ), ((0.275501, -0.005637, 1.32), ), ((0.157825, 0.006647, 1.32), ), 
            ((0.18075, 0.006647, 1.32), ), ((0.20575, 0.006647, 1.32), ), ((0.439133, 
            0.016348, 1.32), ), ((0.233522, 0.030756, 1.22), ), ((0.133545, 0.006033, 
            1.22), ), ((0.272931, -0.005757, 1.22), ), ((0.153064, 0.006783, 1.22), ), 
            ((0.176627, 0.006783, 1.22), ), ((0.201627, 0.006783, 1.22), ), ((0.439457, 
            0.016714, 1.22), ), ((0.230037, 0.031371, 1.12), ), ((0.128228, 0.006163, 
            1.12), ), ((0.270362, -0.005877, 1.12), ), ((0.148304, 0.006918, 1.12), ), 
            ((0.172504, 0.006918, 1.12), ), ((0.197504, 0.006918, 1.12), ), ((0.439782, 
            0.017081, 1.12), ), ((0.226551, 0.031985, 1.02), ), ((0.122912, 0.006294, 
            1.02), ), ((0.267793, -0.005996, 1.02), ), ((0.143544, 0.007054, 1.02), ), 
            ((0.168381, 0.007054, 1.02), ), ((0.193381, 0.007054, 1.02), ), ((0.440106, 
            0.017448, 1.02), ), ((0.223065, 0.0326, 0.92), ), ((0.117595, 0.006424, 
            0.92), ), ((0.265223, -0.006116, 0.92), ), ((0.138784, 0.007189, 0.92), ), 
            ((0.164258, 0.007189, 0.92), ), ((0.189258, 0.007189, 0.92), ), ((0.44043, 
            0.017814, 0.92), ), ((0.216093, 0.033828, 0.72), ), ((0.106962, 0.006685, 
            0.72), ), ((0.260084, -0.006355, 0.72), ), ((0.129263, 0.00746, 0.72), ), (
            (0.156012, 0.00746, 0.72), ), ((0.181012, 0.00746, 0.72), ), ((0.441079, 
            0.018548, 0.72), ), ((0.212607, 0.034443, 0.62), ), ((0.101645, 0.006816, 
            0.62), ), ((0.257515, -0.006475, 0.62), ), ((0.124503, 0.007596, 0.62), ), 
            ((0.151889, 0.007596, 0.62), ), ((0.176889, 0.007596, 0.62), ), ((0.441403, 
            0.018915, 0.62), ), ((0.210516, 0.034812, 0.56), ), ((0.098455, 0.006894, 
            0.56), ), ((0.255973, -0.006547, 0.56), ), ((0.121647, 0.007677, 0.56), ), 
            ((0.149415, 0.007677, 0.56), ), ((0.174415, 0.007677, 0.56), ), ((0.441598, 
            0.019135, 0.56), ), ((0.208773, 0.035119, 0.46), ), ((0.095797, 0.006959, 
            0.46), ), ((0.254689, -0.006606, 0.46), ), ((0.119267, 0.007745, 0.46), ), 
            ((0.147353, 0.007745, 0.46), ), ((0.172353, 0.007745, 0.46), ), ((0.44176, 
            0.019318, 0.46), ), ((0.208773, 0.035119, 0.4), ), ((0.095797, 0.006959, 
            0.4), ), ((0.254689, -0.006606, 0.4), ), ((0.119267, 0.007745, 0.4), ), ((
            0.147353, 0.007745, 0.4), ), ((0.172353, 0.007745, 0.4), ), ((0.44176, 
            0.019318, 0.4), ), ((0.208773, 0.035119, 0.3), ), ((0.095797, 0.006959, 
            0.3), ), ((0.254689, -0.006606, 0.3), ), ((0.119267, 0.007745, 0.3), ), ((
            0.147353, 0.007745, 0.3), ), ((0.172353, 0.007745, 0.3), ), ((0.44176, 
            0.019318, 0.3), ), ((0.208773, 0.035119, 0.2), ), ((0.095797, 0.006959, 
            0.2), ), ((0.254689, -0.006606, 0.2), ), ((0.119267, 0.007745, 0.2), ), ((
            0.147353, 0.007745, 0.2), ), ((0.172353, 0.007745, 0.2), ), ((0.44176, 
            0.019318, 0.2), ), ((0.208773, 0.035119, 0.1), ), ((0.095797, 0.006959, 
            0.1), ), ((0.254689, -0.006606, 0.1), ), ((0.119267, 0.007745, 0.1), ), ((
            0.147353, 0.007745, 0.1), ), ((0.172353, 0.007745, 0.1), ), ((0.44176, 
            0.019318, 0.1), ), ((0.219579, 0.033214, 0.82), ), ((0.112278, 0.006555, 
            0.82), ), ((0.262654, -0.006236, 0.82), ), ((0.134024, 0.007325, 0.82), ), 
            ((0.160135, 0.007325, 0.82), ), ((0.185135, 0.007325, 0.82), ), ((0.440755, 
            0.018181, 0.82), ), ((0.304032, -0.010003, 2.224101), ), ((0.231658, 
            -0.02711, 2.220776), ), ((0.329625, -0.034982, 2.219246), ), ((0.242788, 
            -0.026448, 2.220904), ), ((0.256743, -0.026448, 2.220904), ), ((0.281743, 
            -0.026448, 2.220904), ), ((0.440553, -0.019302, 2.222294), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='skin-mid', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.556624, 0.000275, 
            1.353333), (-0.051907, -0.99865, -0.001728)), ((0.54214, 0.000909, 
            1.353333), (-0.022329, -0.99975, -0.000839)), ((0.529192, 0.00111, 
            1.353333), (-0.004709, -0.999989, -0.000353)), ((0.513316, 0.001088, 
            1.353333), (0.011647, -0.999932, 4.8e-05)), ((0.494789, 0.000778, 
            1.353333), (0.025441, -0.999676, 0.000337)), ((0.473913, 0.000157, 
            1.353333), (0.03742, -0.999299, 0.000538)), ((0.450932, -0.000854, 
            1.353333), (0.055871, -0.998438, 0.000763)), ((0.400294, -0.003688, 
            1.353333), (0.055871, -0.998438, 0.000763)), ((0.364746, -0.005892, 
            1.353333), (0.067641, -0.997709, 0.000722)), ((0.337213, -0.0078, 
            1.353333), (0.069854, -0.997557, 0.000701)), ((0.309561, -0.009733, 
            1.353333), (0.069674, -0.99757, 0.000704)), ((0.282129, -0.012339, 
            1.353333), (0.107257, -0.994231, -7e-05)), ((0.254623, -0.01559, 1.353333), 
            (0.122438, -0.992476, -0.000469)), ((0.231087, -0.016733, 1.353333), (0.0, 
            -0.999994, 0.003434)), ((0.159517, -0.016733, 1.353333), (0.0, -0.999994, 
            0.003434)), ((0.138145, -0.016706, 1.353333), (-0.001899, -0.999992, 
            0.003531)), ((0.119238, -0.016225, 1.353333), (-0.039628, -0.999199, 
            0.005608)), ((0.102961, -0.015217, 1.353333), (-0.074595, -0.997185, 
            0.007657)), ((0.089024, -0.013741, 1.353333), (-0.123513, -0.992286, 
            0.010673)), ((0.077573, -0.011796, 1.353333), (-0.194483, -0.980788, 
            0.015231)), ((0.068745, -0.009406, 1.353333), (-0.304348, -0.952295, 
            0.02251)), ((0.062779, -0.006597, 1.353333), (-0.509417, -0.85975, 
            0.036388)), ((0.059689, -0.001882, 1.353333), (-0.932433, -0.355412, 
            0.065202)), ((0.059145, 0.0012, 1.353333), (-0.987308, 0.143028, 
            0.069039)), ((0.060264, 0.005534, 1.353333), (-0.927883, 0.367154, 
            0.065055)), ((0.06323, 0.011453, 1.353333), (-0.817545, 0.572957, 
            0.057801)), ((0.068156, 0.017553, 1.353333), (-0.700158, 0.71222, 
            0.050206)), ((0.075038, 0.023637, 1.353333), (-0.591304, 0.805285, 
            0.043305)), ((0.083828, 0.029516, 1.353333), (-0.490483, 0.870662, 
            0.03707)), ((0.094498, 0.035015, 1.353333), (-0.398994, 0.91641, 
            0.031577)), ((0.107012, 0.040005, 1.353333), (-0.318094, 0.947678, 
            0.026888)), ((0.121302, 0.044646, 1.353333), (-0.292195, 0.95602, 
            0.025454)), ((0.138492, 0.049623, 1.353333), (-0.256786, 0.96618, 
            0.023606)), ((0.159517, 0.053272, 1.353333), (0.0, 0.99994, 0.010933)), ((
            0.231087, 0.053272, 1.353333), (0.0, 0.99994, 0.010933)), ((0.256299, 
            0.052659, 1.353333), (0.055922, 0.998393, 0.009145)), ((0.286161, 0.050592, 
            1.353333), (0.10351, 0.994597, 0.007916)), ((0.310712, 0.047831, 1.353333), 
            (0.127802, 0.991772, 0.007403)), ((0.335478, 0.04446, 1.353333), (0.148768, 
            0.988847, 0.007063)), ((0.360229, 0.04058, 1.353333), (0.167005, 0.985932, 
            0.006859)), ((0.396121, 0.033991, 1.353333), (0.19184, 0.981403, 
            0.006707)), ((0.442942, 0.024839, 1.353333), (0.19184, 0.981403, 
            0.006707)), ((0.461313, 0.021047, 1.353333), (0.205114, 0.978714, 
            0.006823)), ((0.481718, 0.016746, 1.353333), (0.206854, 0.978348, 
            0.006846)), ((0.500375, 0.012795, 1.353333), (0.207325, 0.978248, 
            0.006854)), ((0.516987, 0.00926, 1.353333), (0.208594, 0.977978, 
            0.006881)), ((0.531313, 0.006188, 1.353333), (0.210343, 0.977603, 
            0.006925)), ((0.54316, 0.003657, 1.353333), (0.208057, 0.978093, 
            0.006861)), ((0.556861, 0.000964, 1.353333), (0.187214, 0.982299, 
            0.006234)), ((0.052156, 0.001225, 1.253333), (-0.987308, 0.143028, 
            0.069039)), ((0.053298, 0.005647, 1.253333), (-0.927883, 0.367154, 
            0.065055)), ((0.056324, 0.011687, 1.253333), (-0.817545, 0.572957, 
            0.057801)), ((0.061352, 0.017913, 1.253333), (-0.700158, 0.71222, 
            0.050206)), ((0.068374, 0.024122, 1.253333), (-0.591304, 0.805285, 
            0.043305)), ((0.077345, 0.030121, 1.253333), (-0.490483, 0.870662, 
            0.03707)), ((0.088234, 0.035734, 1.253333), (-0.398994, 0.91641, 
            0.031577)), ((0.101005, 0.040826, 1.253333), (-0.318094, 0.947678, 
            0.026888)), ((0.115588, 0.045562, 1.253333), (-0.292195, 0.95602, 
            0.025454)), ((0.133131, 0.050641, 1.253333), (-0.256786, 0.96618, 
            0.023606)), ((0.154757, 0.054365, 1.253333), (0.0, 0.99994, 0.010933)), ((
            0.227282, 0.054365, 1.253333), (0.0, 0.99994, 0.010933)), ((0.253355, 
            0.05374, 1.253333), (0.055922, 0.998393, 0.009145)), ((0.28383, 0.051631, 
            1.253333), (0.10351, 0.994597, 0.007916)), ((0.308885, 0.048813, 1.253333), 
            (0.127802, 0.991772, 0.007403)), ((0.334159, 0.045373, 1.253333), (
            0.148768, 0.988847, 0.007063)), ((0.359418, 0.041413, 1.253333), (0.167005, 
            0.985932, 0.006859)), ((0.395805, 0.034736, 1.253333), (0.19184, 0.981403, 
            0.006707)), ((0.443591, 0.025395, 1.253333), (0.19184, 0.981403, 
            0.006707)), ((0.462579, 0.021479, 1.253333), (0.205114, 0.978714, 
            0.006823)), ((0.483403, 0.017089, 1.253333), (0.206854, 0.978348, 
            0.006846)), ((0.502442, 0.013058, 1.253333), (0.207325, 0.978248, 
            0.006854)), ((0.519395, 0.00945, 1.253333), (0.208594, 0.977978, 
            0.006881)), ((0.534015, 0.006315, 1.253333), (0.210343, 0.977603, 
            0.006925)), ((0.546105, 0.003732, 1.253333), (0.208057, 0.978093, 
            0.006861)), ((0.560088, 0.000983, 1.253333), (0.187214, 0.982299, 
            0.006234)), ((0.559846, 0.000281, 1.253333), (-0.051907, -0.99865, 
            -0.001728)), ((0.545064, 0.000928, 1.253333), (-0.022329, -0.99975, 
            -0.000839)), ((0.531851, 0.001133, 1.253333), (-0.004709, -0.999989, 
            -0.000353)), ((0.51565, 0.00111, 1.253333), (0.011647, -0.999932, 
            4.8e-05)), ((0.496742, 0.000794, 1.253333), (0.025441, -0.999676, 
            0.000337)), ((0.475438, 0.000161, 1.253333), (0.03742, -0.999299, 
            0.000538)), ((0.451745, -0.000885, 1.253333), (0.055871, -0.998438, 
            0.000763)), ((0.400065, -0.003777, 1.253333), (0.055871, -0.998438, 
            0.000763)), ((0.364028, -0.006013, 1.253333), (0.067641, -0.997709, 
            0.000722)), ((0.335929, -0.00796, 1.253333), (0.069854, -0.997557, 
            0.000701)), ((0.30771, -0.009932, 1.253333), (0.069674, -0.99757, 
            0.000704)), ((0.279715, -0.012592, 1.253333), (0.107257, -0.994231, 
            -7e-05)), ((0.251644, -0.01591, 1.253333), (0.122438, -0.992476, 
            -0.000469)), ((0.227282, -0.017076, 1.253333), (0.0, -0.999994, 0.003434)), 
            ((0.154757, -0.017076, 1.253333), (0.0, -0.999994, 0.003434)), ((0.132776, 
            -0.017049, 1.253333), (-0.001899, -0.999992, 0.003531)), ((0.113482, 
            -0.016558, 1.253333), (-0.039628, -0.999199, 0.005608)), ((0.096871, 
            -0.015529, 1.253333), (-0.074595, -0.997185, 0.007657)), ((0.082648, 
            -0.014023, 1.253333), (-0.123513, -0.992286, 0.010673)), ((0.070962, 
            -0.012038, 1.253333), (-0.194483, -0.980788, 0.015231)), ((0.061953, 
            -0.009599, 1.253333), (-0.304348, -0.952295, 0.02251)), ((0.055864, 
            -0.006732, 1.253333), (-0.509417, -0.85975, 0.036388)), ((0.052711, 
            -0.00192, 1.253333), (-0.932433, -0.355412, 0.065202)), ((0.563068, 
            0.000286, 1.153333), (-0.051907, -0.99865, -0.001728)), ((0.547989, 
            0.000947, 1.153333), (-0.022329, -0.99975, -0.000839)), ((0.53451, 
            0.001156, 1.153333), (-0.004709, -0.999989, -0.000353)), ((0.517983, 
            0.001132, 1.153333), (0.011647, -0.999932, 4.8e-05)), ((0.498695, 0.00081, 
            1.153333), (0.025441, -0.999676, 0.000337)), ((0.476963, 0.000164, 
            1.153333), (0.03742, -0.999299, 0.000538)), ((0.452558, -0.000916, 
            1.153333), (0.055871, -0.998438, 0.000763)), ((0.399835, -0.003867, 
            1.153333), (0.055871, -0.998438, 0.000763)), ((0.36331, -0.006134, 
            1.153333), (0.067641, -0.997709, 0.000722)), ((0.334646, -0.00812, 
            1.153333), (0.069854, -0.997557, 0.000701)), ((0.305859, -0.010132, 
            1.153333), (0.069674, -0.99757, 0.000704)), ((0.277301, -0.012846, 
            1.153333), (0.107257, -0.994231, -7e-05)), ((0.248666, -0.016231, 
            1.153333), (0.122438, -0.992476, -0.000469)), ((0.223478, -0.01742, 
            1.153333), (0.0, -0.999994, 0.003434)), ((0.149997, -0.01742, 1.153333), (
            0.0, -0.999994, 0.003434)), ((0.127407, -0.017392, 1.153333), (-0.001899, 
            -0.999992, 0.003531)), ((0.107725, -0.016891, 1.153333), (-0.039628, 
            -0.999199, 0.005608)), ((0.09078, -0.015841, 1.153333), (-0.074595, 
            -0.997185, 0.007657)), ((0.076271, -0.014305, 1.153333), (-0.123513, 
            -0.992286, 0.010673)), ((0.06435, -0.01228, 1.153333), (-0.194483, 
            -0.980788, 0.015231)), ((0.05516, -0.009792, 1.153333), (-0.304348, 
            -0.952295, 0.02251)), ((0.048949, -0.006868, 1.153333), (-0.509417, 
            -0.85975, 0.036388)), ((0.045733, -0.001959, 1.153333), (-0.932433, 
            -0.355412, 0.065202)), ((0.045167, 0.001249, 1.153333), (-0.987308, 
            0.143028, 0.069039)), ((0.046332, 0.005761, 1.153333), (-0.927883, 
            0.367154, 0.065055)), ((0.049419, 0.011922, 1.153333), (-0.817545, 
            0.572957, 0.057801)), ((0.054547, 0.018273, 1.153333), (-0.700158, 0.71222, 
            0.050206)), ((0.061711, 0.024607, 1.153333), (-0.591304, 0.805285, 
            0.043305)), ((0.070862, 0.030727, 1.153333), (-0.490483, 0.870662, 
            0.03707)), ((0.08197, 0.036452, 1.153333), (-0.398994, 0.91641, 0.031577)), 
            ((0.094998, 0.041646, 1.153333), (-0.318094, 0.947678, 0.026888)), ((
            0.109873, 0.046478, 1.153333), (-0.292195, 0.95602, 0.025454)), ((0.127769, 
            0.05166, 1.153333), (-0.256786, 0.96618, 0.023606)), ((0.149997, 0.055459, 
            1.153333), (0.0, 0.99994, 0.010933)), ((0.223478, 0.055459, 1.153333), (
            0.0, 0.99994, 0.010933)), ((0.250411, 0.054821, 1.153333), (0.055922, 
            0.998393, 0.009145)), ((0.281499, 0.052669, 1.153333), (0.10351, 0.994597, 
            0.007916)), ((0.307058, 0.049795, 1.153333), (0.127802, 0.991772, 
            0.007403)), ((0.33284, 0.046286, 1.153333), (0.148768, 0.988847, 
            0.007063)), ((0.358607, 0.042246, 1.153333), (0.167005, 0.985932, 
            0.006859)), ((0.395489, 0.035481, 1.153333), (0.19184, 0.981403, 
            0.006707)), ((0.444239, 0.025952, 1.153333), (0.19184, 0.981403, 
            0.006707)), ((0.463845, 0.02191, 1.153333), (0.205114, 0.978714, 
            0.006823)), ((0.485088, 0.017433, 1.153333), (0.206854, 0.978348, 
            0.006846)), ((0.50451, 0.01332, 1.153333), (0.207325, 0.978248, 0.006854)), 
            ((0.521804, 0.00964, 1.153333), (0.208594, 0.977978, 0.006881)), ((
            0.536718, 0.006442, 1.153333), (0.210343, 0.977603, 0.006925)), ((0.549051, 
            0.003807, 1.153333), (0.208057, 0.978093, 0.006861)), ((0.563314, 0.001003, 
            1.153333), (0.187214, 0.982299, 0.006234)), ((0.038178, 0.001274, 
            1.053333), (-0.987308, 0.143028, 0.069039)), ((0.039366, 0.005874, 
            1.053333), (-0.927883, 0.367154, 0.065055)), ((0.042513, 0.012157, 
            1.053333), (-0.817545, 0.572957, 0.057801)), ((0.047743, 0.018633, 
            1.053333), (-0.700158, 0.71222, 0.050206)), ((0.055048, 0.025092, 
            1.053333), (-0.591304, 0.805285, 0.043305)), ((0.064379, 0.031332, 
            1.053333), (-0.490483, 0.870662, 0.03707)), ((0.075706, 0.03717, 1.053333), 
            (-0.398994, 0.91641, 0.031577)), ((0.08899, 0.042467, 1.053333), (
            -0.318094, 0.947678, 0.026888)), ((0.104159, 0.047394, 1.053333), (
            -0.292195, 0.95602, 0.025454)), ((0.122408, 0.052678, 1.053333), (
            -0.256786, 0.96618, 0.023606)), ((0.145237, 0.056552, 1.053333), (0.0, 
            0.99994, 0.010933)), ((0.219673, 0.056552, 1.053333), (0.0, 0.99994, 
            0.010933)), ((0.247466, 0.055902, 1.053333), (0.055922, 0.998393, 
            0.009145)), ((0.279168, 0.053708, 1.053333), (0.10351, 0.994597, 
            0.007916)), ((0.305231, 0.050777, 1.053333), (0.127802, 0.991772, 
            0.007403)), ((0.331521, 0.047198, 1.053333), (0.148768, 0.988847, 
            0.007063)), ((0.357797, 0.043079, 1.053333), (0.167005, 0.985932, 
            0.006859)), ((0.395174, 0.036226, 1.053333), (0.19184, 0.981403, 
            0.006707)), ((0.444888, 0.026508, 1.053333), (0.19184, 0.981403, 
            0.006707)), ((0.465111, 0.022342, 1.053333), (0.205114, 0.978714, 
            0.006823)), ((0.486772, 0.017776, 1.053333), (0.206854, 0.978348, 
            0.006846)), ((0.506578, 0.013583, 1.053333), (0.207325, 0.978248, 
            0.006854)), ((0.524212, 0.00983, 1.053333), (0.208594, 0.977978, 
            0.006881)), ((0.53942, 0.006568, 1.053333), (0.210343, 0.977603, 
            0.006925)), ((0.551996, 0.003882, 1.053333), (0.208057, 0.978093, 
            0.006861)), ((0.566541, 0.001023, 1.053333), (0.187214, 0.982299, 
            0.006234)), ((0.56629, 0.000292, 1.053333), (-0.051907, -0.99865, 
            -0.001728)), ((0.550914, 0.000965, 1.053333), (-0.022329, -0.99975, 
            -0.000839)), ((0.537169, 0.001178, 1.053333), (-0.004709, -0.999989, 
            -0.000353)), ((0.520316, 0.001155, 1.053333), (0.011647, -0.999932, 
            4.8e-05)), ((0.500648, 0.000826, 1.053333), (0.025441, -0.999676, 
            0.000337)), ((0.478488, 0.000167, 1.053333), (0.03742, -0.999299, 
            0.000538)), ((0.453371, -0.000947, 1.053333), (0.055871, -0.998438, 
            0.000763)), ((0.399605, -0.003956, 1.053333), (0.055871, -0.998438, 
            0.000763)), ((0.362592, -0.006255, 1.053333), (0.067641, -0.997709, 
            0.000722)), ((0.333363, -0.00828, 1.053333), (0.069854, -0.997557, 
            0.000701)), ((0.304008, -0.010332, 1.053333), (0.069674, -0.99757, 
            0.000704)), ((0.274888, -0.013099, 1.053333), (0.107257, -0.994231, 
            -7e-05)), ((0.245688, -0.016551, 1.053333), (0.122438, -0.992476, 
            -0.000469)), ((0.219673, -0.017763, 1.053333), (0.0, -0.999994, 0.003434)), 
            ((0.145237, -0.017763, 1.053333), (0.0, -0.999994, 0.003434)), ((0.122039, 
            -0.017735, 1.053333), (-0.001899, -0.999992, 0.003531)), ((0.101968, 
            -0.017224, 1.053333), (-0.039628, -0.999199, 0.005608)), ((0.08469, 
            -0.016154, 1.053333), (-0.074595, -0.997185, 0.007657)), ((0.069895, 
            -0.014587, 1.053333), (-0.123513, -0.992286, 0.010673)), ((0.057739, 
            -0.012522, 1.053333), (-0.194483, -0.980788, 0.015231)), ((0.048368, 
            -0.009985, 1.053333), (-0.304348, -0.952295, 0.02251)), ((0.042034, 
            -0.007003, 1.053333), (-0.509417, -0.85975, 0.036388)), ((0.038755, 
            -0.001997, 1.053333), (-0.932433, -0.355412, 0.065202)), ((0.569512, 
            0.000298, 0.953333), (-0.051907, -0.99865, -0.001728)), ((0.553838, 
            0.000984, 0.953333), (-0.022329, -0.99975, -0.000839)), ((0.539828, 
            0.001201, 0.953333), (-0.004709, -0.999989, -0.000353)), ((0.522649, 
            0.001177, 0.953333), (0.011647, -0.999932, 4.8e-05)), ((0.502601, 0.000842, 
            0.953333), (0.025441, -0.999676, 0.000337)), ((0.480012, 0.00017, 
            0.953333), (0.03742, -0.999299, 0.000538)), ((0.454184, -0.000978, 
            0.953333), (0.055871, -0.998438, 0.000763)), ((0.399375, -0.004045, 
            0.953333), (0.055871, -0.998438, 0.000763)), ((0.361874, -0.006376, 
            0.953333), (0.067641, -0.997709, 0.000722)), ((0.332079, -0.00844, 
            0.953333), (0.069854, -0.997557, 0.000701)), ((0.302157, -0.010532, 
            0.953333), (0.069674, -0.99757, 0.000704)), ((0.272474, -0.013352, 
            0.953333), (0.107257, -0.994231, -7e-05)), ((0.242709, -0.016871, 
            0.953333), (0.122438, -0.992476, -0.000469)), ((0.215869, -0.018107, 
            0.953333), (0.0, -0.999994, 0.003434)), ((0.140477, -0.018107, 0.953333), (
            0.0, -0.999994, 0.003434)), ((0.11667, -0.018078, 0.953333), (-0.001899, 
            -0.999992, 0.003531)), ((0.096212, -0.017557, 0.953333), (-0.039628, 
            -0.999199, 0.005608)), ((0.078599, -0.016466, 0.953333), (-0.074595, 
            -0.997185, 0.007657)), ((0.063518, -0.014869, 0.953333), (-0.123513, 
            -0.992286, 0.010673)), ((0.051127, -0.012764, 0.953333), (-0.194483, 
            -0.980788, 0.015231)), ((0.041575, -0.010178, 0.953333), (-0.304348, 
            -0.952295, 0.02251)), ((0.035119, -0.007138, 0.953333), (-0.509417, 
            -0.85975, 0.036388)), ((0.031777, -0.002036, 0.953333), (-0.932433, 
            -0.355412, 0.065202)), ((0.031189, 0.001298, 0.953333), (-0.987308, 
            0.143028, 0.069039)), ((0.032399, 0.005987, 0.953333), (-0.927883, 
            0.367154, 0.065055)), ((0.035608, 0.012392, 0.953333), (-0.817545, 
            0.572957, 0.057801)), ((0.040938, 0.018993, 0.953333), (-0.700158, 0.71222, 
            0.050206)), ((0.048384, 0.025576, 0.953333), (-0.591304, 0.805285, 
            0.043305)), ((0.057896, 0.031938, 0.953333), (-0.490483, 0.870662, 
            0.03707)), ((0.069442, 0.037889, 0.953333), (-0.398994, 0.91641, 
            0.031577)), ((0.082983, 0.043288, 0.953333), (-0.318094, 0.947678, 
            0.026888)), ((0.098445, 0.04831, 0.953333), (-0.292195, 0.95602, 
            0.025454)), ((0.117046, 0.053696, 0.953333), (-0.256786, 0.96618, 
            0.023606)), ((0.140477, 0.057646, 0.953333), (0.0, 0.99994, 0.010933)), ((
            0.215869, 0.057646, 0.953333), (0.0, 0.99994, 0.010933)), ((0.244522, 
            0.056983, 0.953333), (0.055922, 0.998393, 0.009145)), ((0.276838, 0.054746, 
            0.953333), (0.10351, 0.994597, 0.007916)), ((0.303404, 0.051759, 0.953333), 
            (0.127802, 0.991772, 0.007403)), ((0.330203, 0.048111, 0.953333), (
            0.148768, 0.988847, 0.007063)), ((0.356986, 0.043912, 0.953333), (0.167005, 
            0.985932, 0.006859)), ((0.394858, 0.036971, 0.953333), (0.19184, 0.981403, 
            0.006707)), ((0.445536, 0.027065, 0.953333), (0.19184, 0.981403, 
            0.006707)), ((0.466378, 0.022774, 0.953333), (0.205114, 0.978714, 
            0.006823)), ((0.488457, 0.01812, 0.953333), (0.206854, 0.978348, 
            0.006846)), ((0.508645, 0.013845, 0.953333), (0.207325, 0.978248, 
            0.006854)), ((0.526621, 0.01002, 0.953333), (0.208594, 0.977978, 
            0.006881)), ((0.542122, 0.006695, 0.953333), (0.210343, 0.977603, 
            0.006925)), ((0.554941, 0.003956, 0.953333), (0.208057, 0.978093, 
            0.006861)), ((0.569768, 0.001042, 0.953333), (0.187214, 0.982299, 
            0.006234)), ((0.0242, 0.001322, 0.853333), (-0.987308, 0.143028, 
            0.069039)), ((0.025433, 0.006101, 0.853333), (-0.927883, 0.367154, 
            0.065055)), ((0.028702, 0.012626, 0.853333), (-0.817545, 0.572957, 
            0.057801)), ((0.034134, 0.019353, 0.853333), (-0.700158, 0.71222, 
            0.050206)), ((0.041721, 0.026061, 0.853333), (-0.591304, 0.805285, 
            0.043305)), ((0.051413, 0.032543, 0.853333), (-0.490483, 0.870662, 
            0.03707)), ((0.063178, 0.038607, 0.853333), (-0.398994, 0.91641, 
            0.031577)), ((0.076976, 0.044109, 0.853333), (-0.318094, 0.947678, 
            0.026888)), ((0.092731, 0.049226, 0.853333), (-0.292195, 0.95602, 
            0.025454)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.111684, 
            0.054714, 0.853333), (-0.256786, 0.96618, 0.023606)), ((0.135716, 0.058739, 
            0.853333), (0.0, 0.99994, 0.010933)), ((0.212065, 0.058739, 0.853333), (
            0.0, 0.99994, 0.010933)), ((0.241578, 0.058064, 0.853333), (0.055922, 
            0.998393, 0.009145)), ((0.274507, 0.055785, 0.853333), (0.10351, 0.994597, 
            0.007916)), ((0.301577, 0.052741, 0.853333), (0.127802, 0.991772, 
            0.007403)), ((0.328884, 0.049024, 0.853333), (0.148768, 0.988847, 
            0.007063)), ((0.356175, 0.044745, 0.853333), (0.167005, 0.985932, 
            0.006859)), ((0.394542, 0.037717, 0.853333), (0.19184, 0.981403, 
            0.006707)), ((0.446185, 0.027622, 0.853333), (0.19184, 0.981403, 
            0.006707)), ((0.467644, 0.023206, 0.853333), (0.205114, 0.978714, 
            0.006823)), ((0.490142, 0.018463, 0.853333), (0.206854, 0.978348, 
            0.006846)), ((0.510713, 0.014107, 0.853333), (0.207325, 0.978248, 
            0.006854)), ((0.529029, 0.01021, 0.853333), (0.208594, 0.977978, 
            0.006881)), ((0.544825, 0.006822, 0.853333), (0.210343, 0.977603, 
            0.006925)), ((0.557887, 0.004031, 0.853333), (0.208057, 0.978093, 
            0.006861)), ((0.572995, 0.001062, 0.853333), (0.187214, 0.982299, 
            0.006234)), ((0.572734, 0.000303, 0.853333), (-0.051907, -0.99865, 
            -0.001728)), ((0.556763, 0.001003, 0.853333), (-0.022329, -0.99975, 
            -0.000839)), ((0.542487, 0.001224, 0.853333), (-0.004709, -0.999989, 
            -0.000353)), ((0.524982, 0.001199, 0.853333), (0.011647, -0.999932, 
            4.8e-05)), ((0.504554, 0.000858, 0.853333), (0.025441, -0.999676, 
            0.000337)), ((0.481537, 0.000174, 0.853333), (0.03742, -0.999299, 
            0.000538)), ((0.454997, -0.001009, 0.853333), (0.055871, -0.998438, 
            0.000763)), ((0.399146, -0.004134, 0.853333), (0.055871, -0.998438, 
            0.000763)), ((0.361156, -0.006497, 0.853333), (0.067641, -0.997709, 
            0.000722)), ((0.330796, -0.0086, 0.853333), (0.069854, -0.997557, 
            0.000701)), ((0.300306, -0.010732, 0.853333), (0.069674, -0.99757, 
            0.000704)), ((0.27006, -0.013606, 0.853333), (0.107257, -0.994231, 
            -7e-05)), ((0.239731, -0.017191, 0.853333), (0.122438, -0.992476, 
            -0.000469)), ((0.212065, -0.01845, 0.853333), (0.0, -0.999994, 0.003434)), 
            ((0.135716, -0.01845, 0.853333), (0.0, -0.999994, 0.003434)), ((0.111301, 
            -0.018421, 0.853333), (-0.001899, -0.999992, 0.003531)), ((0.090455, 
            -0.01789, 0.853333), (-0.039628, -0.999199, 0.005608)), ((0.072508, 
            -0.016778, 0.853333), (-0.074595, -0.997185, 0.007657)), ((0.057142, 
            -0.015151, 0.853333), (-0.123513, -0.992286, 0.010673)), ((0.044516, 
            -0.013006, 0.853333), (-0.194483, -0.980788, 0.015231)), ((0.034783, 
            -0.01037, 0.853333), (-0.304348, -0.952295, 0.02251)), ((0.028205, 
            -0.007273, 0.853333), (-0.509417, -0.85975, 0.036388)), ((0.024799, 
            -0.002074, 0.853333), (-0.932433, -0.355412, 0.065202)), ((0.575956, 
            0.000309, 0.753333), (-0.051907, -0.99865, -0.001728)), ((0.559687, 
            0.001021, 0.753333), (-0.022329, -0.99975, -0.000839)), ((0.545146, 
            0.001247, 0.753333), (-0.004709, -0.999989, -0.000353)), ((0.527316, 
            0.001222, 0.753333), (0.011647, -0.999932, 4.8e-05)), ((0.506507, 0.000874, 
            0.753333), (0.025441, -0.999676, 0.000337)), ((0.483062, 0.000177, 
            0.753333), (0.03742, -0.999299, 0.000538)), ((0.45581, -0.00104, 0.753333), 
            (0.055871, -0.998438, 0.000763)), ((0.398916, -0.004224, 0.753333), (
            0.055871, -0.998438, 0.000763)), ((0.360437, -0.006618, 0.753333), (
            0.067641, -0.997709, 0.000722)), ((0.329513, -0.008761, 0.753333), (
            0.069854, -0.997557, 0.000701)), ((0.298456, -0.010932, 0.753333), (
            0.069674, -0.99757, 0.000704)), ((0.267646, -0.013859, 0.753333), (
            0.107257, -0.994231, -7e-05)), ((0.236753, -0.017511, 0.753333), (0.122438, 
            -0.992476, -0.000469)), ((0.20826, -0.018794, 0.753333), (0.0, -0.999994, 
            0.003434)), ((0.130956, -0.018794, 0.753333), (0.0, -0.999994, 0.003434)), 
            ((0.105932, -0.018764, 0.753333), (-0.001899, -0.999992, 0.003531)), ((
            0.084699, -0.018223, 0.753333), (-0.039628, -0.999199, 0.005608)), ((
            0.066418, -0.01709, 0.753333), (-0.074595, -0.997185, 0.007657)), ((
            0.050765, -0.015433, 0.753333), (-0.123513, -0.992286, 0.010673)), ((
            0.037905, -0.013248, 0.753333), (-0.194483, -0.980788, 0.015231)), ((
            0.02799, -0.010563, 0.753333), (-0.304348, -0.952295, 0.02251)), ((0.02129, 
            -0.007409, 0.753333), (-0.509417, -0.85975, 0.036388)), ((0.017821, 
            -0.002112, 0.753333), (-0.932433, -0.355412, 0.065202)), ((0.017211, 
            0.001347, 0.753333), (-0.987308, 0.143028, 0.069039)), ((0.018467, 
            0.006214, 0.753333), (-0.927883, 0.367154, 0.065055)), ((0.021797, 
            0.012861, 0.753333), (-0.817545, 0.572957, 0.057801)), ((0.027329, 
            0.019713, 0.753333), (-0.700158, 0.71222, 0.050206)), ((0.035058, 0.026546, 
            0.753333), (-0.591304, 0.805285, 0.043305)), ((0.04493, 0.033149, 
            0.753333), (-0.490483, 0.870662, 0.03707)), ((0.056914, 0.039326, 
            0.753333), (-0.398994, 0.91641, 0.031577)), ((0.070969, 0.04493, 0.753333), 
            (-0.318094, 0.947678, 0.026888)), ((0.087017, 0.050142, 0.753333), (
            -0.292195, 0.95602, 0.025454)), ((0.106323, 0.055732, 0.753333), (
            -0.256786, 0.96618, 0.023606)), ((0.130956, 0.059832, 0.753333), (0.0, 
            0.99994, 0.010933)), ((0.20826, 0.059832, 0.753333), (0.0, 0.99994, 
            0.010933)), ((0.238634, 0.059145, 0.753333), (0.055922, 0.998393, 
            0.009145)), ((0.272176, 0.056823, 0.753333), (0.10351, 0.994597, 
            0.007916)), ((0.29975, 0.053722, 0.753333), (0.127802, 0.991772, 
            0.007403)), ((0.327565, 0.049936, 0.753333), (0.148768, 0.988847, 
            0.007063)), ((0.355364, 0.045578, 0.753333), (0.167005, 0.985932, 
            0.006859)), ((0.394227, 0.038462, 0.753333), (0.19184, 0.981403, 
            0.006707)), ((0.446833, 0.028178, 0.753333), (0.19184, 0.981403, 
            0.006707)), ((0.46891, 0.023637, 0.753333), (0.205114, 0.978714, 
            0.006823)), ((0.491827, 0.018807, 0.753333), (0.206854, 0.978348, 
            0.006846)), ((0.512781, 0.01437, 0.753333), (0.207325, 0.978248, 
            0.006854)), ((0.531438, 0.0104, 0.753333), (0.208594, 0.977978, 0.006881)), 
            ((0.547527, 0.006949, 0.753333), (0.210343, 0.977603, 0.006925)), ((
            0.560832, 0.004106, 0.753333), (0.208057, 0.978093, 0.006861)), ((0.576222, 
            0.001082, 0.753333), (0.187214, 0.982299, 0.006234)), ((0.579178, 0.000314, 
            0.653333), (-0.051907, -0.99865, -0.001728)), ((0.562612, 0.00104, 
            0.653333), (-0.022329, -0.99975, -0.000839)), ((0.547805, 0.001269, 
            0.653333), (-0.004709, -0.999989, -0.000353)), ((0.529649, 0.001244, 
            0.653333), (0.011647, -0.999932, 4.8e-05)), ((0.50846, 0.00089, 0.653333), 
            (0.025441, -0.999676, 0.000337)), ((0.484586, 0.00018, 0.653333), (0.03742, 
            -0.999299, 0.000538)), ((0.456623, -0.001071, 0.653333), (0.055871, 
            -0.998438, 0.000763)), ((0.398686, -0.004313, 0.653333), (0.055871, 
            -0.998438, 0.000763)), ((0.359719, -0.006739, 0.653333), (0.067641, 
            -0.997709, 0.000722)), ((0.328229, -0.008921, 0.653333), (0.069854, 
            -0.997557, 0.000701)), ((0.296605, -0.011132, 0.653333), (0.069674, 
            -0.99757, 0.000704)), ((0.265232, -0.014113, 0.653333), (0.107257, 
            -0.994231, -7e-05)), ((0.233774, -0.017831, 0.653333), (0.122438, 
            -0.992476, -0.000469)), ((0.204456, -0.019137, 0.653333), (0.0, -0.999994, 
            0.003434)), ((0.126196, -0.019137, 0.653333), (0.0, -0.999994, 0.003434)), 
            ((0.100564, -0.019106, 0.653333), (-0.001899, -0.999992, 0.003531)), ((
            0.078942, -0.018556, 0.653333), (-0.039628, -0.999199, 0.005608)), ((
            0.060327, -0.017403, 0.653333), (-0.074595, -0.997185, 0.007657)), ((
            0.044389, -0.015714, 0.653333), (-0.123513, -0.992286, 0.010673)), ((
            0.031293, -0.01349, 0.653333), (-0.194483, -0.980788, 0.015231)), ((
            0.021198, -0.010756, 0.653333), (-0.304348, -0.952295, 0.02251)), ((
            0.014375, -0.007544, 0.653333), (-0.509417, -0.85975, 0.036388)), ((
            0.010843, -0.002151, 0.653333), (-0.932433, -0.355412, 0.065202)), ((
            0.010222, 0.001371, 0.653333), (-0.987308, 0.143028, 0.069039)), ((
            0.011501, 0.006327, 0.653333), (-0.927883, 0.367154, 0.065055)), ((
            0.014891, 0.013096, 0.653333), (-0.817545, 0.572957, 0.057801)), ((
            0.020525, 0.020073, 0.653333), (-0.700158, 0.71222, 0.050206)), ((0.028394, 
            0.027031, 0.653333), (-0.591304, 0.805285, 0.043305)), ((0.038447, 
            0.033754, 0.653333), (-0.490483, 0.870662, 0.03707)), ((0.05065, 0.040044, 
            0.653333), (-0.398994, 0.91641, 0.031577)), ((0.064961, 0.045751, 
            0.653333), (-0.318094, 0.947678, 0.026888)), ((0.081303, 0.051058, 
            0.653333), (-0.292195, 0.95602, 0.025454)), ((0.100961, 0.056751, 
            0.653333), (-0.256786, 0.96618, 0.023606)), ((0.126196, 0.060926, 
            0.653333), (0.0, 0.99994, 0.010933)), ((0.204456, 0.060926, 0.653333), (
            0.0, 0.99994, 0.010933)), ((0.235689, 0.060225, 0.653333), (0.055922, 
            0.998393, 0.009145)), ((0.269845, 0.057862, 0.653333), (0.10351, 0.994597, 
            0.007916)), ((0.297923, 0.054704, 0.653333), (0.127802, 0.991772, 
            0.007403)), ((0.326246, 0.050849, 0.653333), (0.148768, 0.988847, 
            0.007063)), ((0.354554, 0.046411, 0.653333), (0.167005, 0.985932, 
            0.006859)), ((0.393911, 0.039207, 0.653333), (0.19184, 0.981403, 
            0.006707)), ((0.447482, 0.028735, 0.653333), (0.19184, 0.981403, 
            0.006707)), ((0.470176, 0.024069, 0.653333), (0.205114, 0.978714, 
            0.006823)), ((0.493511, 0.01915, 0.653333), (0.206854, 0.978348, 
            0.006846)), ((0.514848, 0.014632, 0.653333), (0.207325, 0.978248, 
            0.006854)), ((0.533846, 0.01059, 0.653333), (0.208594, 0.977978, 
            0.006881)), ((0.55023, 0.007076, 0.653333), (0.210343, 0.977603, 
            0.006925)), ((0.563778, 0.004181, 0.653333), (0.208057, 0.978093, 
            0.006861)), ((0.579448, 0.001101, 0.653333), (0.187214, 0.982299, 
            0.006234)), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='skin-root', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.005096, 0.001386, 
            0.58), (-0.987308, 0.143028, 0.069039)), ((0.00639, 0.006405, 0.58), (
            -0.927883, 0.367154, 0.065055)), ((0.009823, 0.013263, 0.58), (-0.817545, 
            0.572957, 0.057801)), ((0.015529, 0.020332, 0.58), (-0.700158, 0.71222, 
            0.050206)), ((0.0235, 0.027381, 0.58), (-0.591304, 0.805285, 0.043305)), ((
            0.033684, 0.034193, 0.58), (-0.490483, 0.870662, 0.03707)), ((0.046045, 
            0.040566, 0.58), (-0.398994, 0.91641, 0.031577)), ((0.060543, 0.046349, 
            0.58), (-0.318094, 0.947678, 0.026888)), ((0.077098, 0.051726, 0.58), (
            -0.292195, 0.95602, 0.025454)), ((0.097011, 0.057492, 0.58), (-0.256786, 
            0.96618, 0.023606)), ((0.122663, 0.061727, 0.58), (0.0, 0.99994, 
            0.010933)), ((0.201623, 0.061727, 0.58), (0.0, 0.99994, 0.010933)), ((
            0.233501, 0.06102, 0.58), (0.055922, 0.998393, 0.009145)), ((0.268113, 
            0.058626, 0.58), (0.10351, 0.994597, 0.007916)), ((0.296561, 0.055427, 
            0.58), (0.127802, 0.991772, 0.007403)), ((0.325256, 0.051522, 0.58), (
            0.148768, 0.988847, 0.007063)), ((0.353937, 0.047026, 0.58), (0.167005, 
            0.985932, 0.006859)), ((0.393658, 0.039757, 0.58), (0.19184, 0.981403, 
            0.006707)), ((0.448001, 0.029135, 0.58), (0.19184, 0.981403, 0.006707)), ((
            0.471124, 0.024382, 0.58), (0.205114, 0.978714, 0.006823)), ((0.494765, 
            0.019399, 0.58), (0.206854, 0.978348, 0.006846)), ((0.516381, 0.014821, 
            0.58), (0.207325, 0.978248, 0.006854)), ((0.535627, 0.010726, 0.58), (
            0.208594, 0.977978, 0.006881)), ((0.552224, 0.007166, 0.58), (0.210343, 
            0.977603, 0.006925)), ((0.565948, 0.004234, 0.58), (0.208057, 0.978093, 
            0.006861)), ((0.581828, 0.001113, 0.58), (0.187214, 0.982299, 0.006234)), (
            (0.581556, 0.000318, 0.58), (-0.051907, -0.99865, -0.001728)), ((0.564767, 
            0.001053, 0.58), (-0.022329, -0.99975, -0.000839)), ((0.549768, 0.001286, 
            0.58), (-0.004709, -0.999989, -0.000353)), ((0.531376, 0.00126, 0.58), (
            0.011647, -0.999932, 4.8e-05)), ((0.509911, 0.000903, 0.58), (0.025441, 
            -0.999676, 0.000337)), ((0.485725, 0.000183, 0.58), (0.03742, -0.999299, 
            0.000538)), ((0.457274, -0.00109, 0.58), (0.055871, -0.998438, 0.000763)), 
            ((0.398502, -0.004379, 0.58), (0.055871, -0.998438, 0.000763)), ((0.359168, 
            -0.00683, 0.58), (0.067641, -0.997709, 0.000722)), ((0.327263, -0.00904, 
            0.58), (0.069854, -0.997557, 0.000701)), ((0.295222, -0.01128, 0.58), (
            0.069674, -0.99757, 0.000704)), ((0.263437, -0.014301, 0.58), (0.107257, 
            -0.994231, -7e-05)), ((0.231565, -0.018069, 0.58), (0.122438, -0.992476, 
            -0.000469)), ((0.201623, -0.019389, 0.58), (0.0, -0.999994, 0.003434)), ((
            0.122663, -0.019389, 0.58), (0.0, -0.999994, 0.003434)), ((0.096607, 
            -0.019358, 0.58), (-0.001899, -0.999992, 0.003531)), ((0.074704, -0.018799, 
            0.58), (-0.039628, -0.999199, 0.005608)), ((0.055847, -0.017631, 0.58), (
            -0.074595, -0.997185, 0.007657)), ((0.039701, -0.01592, 0.58), (-0.123513, 
            -0.992286, 0.010673)), ((0.026435, -0.013665, 0.58), (-0.194483, -0.980788, 
            0.015231)), ((0.01621, -0.010896, 0.58), (-0.304348, -0.952295, 0.02251)), 
            ((0.0093, -0.00764, 0.58), (-0.509417, -0.85975, 0.036388)), ((0.005723, 
            -0.002174, 0.58), (-0.932433, -0.355412, 0.065202)), ((0.583278, 0.000321, 
            0.526667), (-0.051907, -0.99865, -0.001728)), ((0.56633, 0.001063, 
            0.526667), (-0.022329, -0.99975, -0.000839)), ((0.55119, 0.001298, 
            0.526667), (-0.004709, -0.999989, -0.000353)), ((0.532624, 0.001272, 
            0.526667), (0.011647, -0.999932, 4.8e-05)), ((0.510957, 0.000911, 
            0.526667), (0.025441, -0.999676, 0.000337)), ((0.486543, 0.000185, 
            0.526667), (0.03742, -0.999299, 0.000538)), ((0.457721, -0.001106, 
            0.526667), (0.055871, -0.998438, 0.000763)), ((0.398376, -0.004427, 
            0.526667), (0.055871, -0.998438, 0.000763)), ((0.358778, -0.006894, 
            0.526667), (0.067641, -0.997709, 0.000722)), ((0.326572, -0.009126, 
            0.526667), (0.069854, -0.997557, 0.000701)), ((0.294229, -0.011387, 
            0.526667), (0.069674, -0.99757, 0.000704)), ((0.262144, -0.014437, 
            0.526667), (0.107257, -0.994231, -7e-05)), ((0.22997, -0.018241, 0.526667), 
            (0.122438, -0.992476, -0.000469)), ((0.199584, -0.019572, 0.526667), (0.0, 
            -0.999994, 0.003434)), ((0.120113, -0.019572, 0.526667), (0.0, -0.999994, 
            0.003434)), ((0.093739, -0.019541, 0.526667), (-0.001899, -0.999992, 
            0.003531)), ((0.07163, -0.018977, 0.526667), (-0.039628, -0.999199, 
            0.005608)), ((0.052595, -0.017797, 0.526667), (-0.074595, -0.997185, 
            0.007657)), ((0.036297, -0.01607, 0.526667), (-0.123513, -0.992286, 
            0.010673)), ((0.022907, -0.013794, 0.526667), (-0.194483, -0.980788, 
            0.015231)), ((0.012585, -0.010998, 0.526667), (-0.304348, -0.952295, 
            0.02251)), ((0.005611, -0.007712, 0.526667), (-0.509417, -0.85975, 
            0.036388)), ((0.002001, -0.002193, 0.526667), (-0.932433, -0.355412, 
            0.065202)), ((0.001368, 0.001398, 0.526667), (-0.987308, 0.143028, 
            0.069039)), ((0.002674, 0.006464, 0.526667), (-0.927883, 0.367154, 
            0.065055)), ((0.006139, 0.013386, 0.526667), (-0.817545, 0.572957, 
            0.057801)), ((0.011899, 0.020522, 0.526667), (-0.700158, 0.71222, 
            0.050206)), ((0.019945, 0.027638, 0.526667), (-0.591304, 0.805285, 
            0.043305)), ((0.030224, 0.034515, 0.526667), (-0.490483, 0.870662, 
            0.03707)), ((0.042702, 0.040948, 0.526667), (-0.398994, 0.91641, 
            0.031577)), ((0.057336, 0.046785, 0.526667), (-0.318094, 0.947678, 
            0.026888)), ((0.074047, 0.052213, 0.526667), (-0.292195, 0.95602, 
            0.025454)), ((0.094146, 0.058034, 0.526667), (-0.256786, 0.96618, 
            0.023606)), ((0.120113, 0.062311, 0.526667), (0.0, 0.99994, 0.010933)), ((
            0.199584, 0.062311, 0.526667), (0.0, 0.99994, 0.010933)), ((0.231923, 
            0.061597, 0.526667), (0.055922, 0.998393, 0.009145)), ((0.266864, 0.05918, 
            0.526667), (0.10351, 0.994597, 0.007916)), ((0.295581, 0.055952, 0.526667), 
            (0.127802, 0.991772, 0.007403)), ((0.324547, 0.052009, 0.526667), (
            0.148768, 0.988847, 0.007063)), ((0.353498, 0.047471, 0.526667), (0.167005, 
            0.985932, 0.006859)), ((0.393485, 0.040156, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.448357, 0.029429, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.471804, 0.024611, 0.526667), (0.205114, 0.978714, 
            0.006823)), ((0.495668, 0.019581, 0.526667), (0.206854, 0.978348, 
            0.006846)), ((0.517488, 0.01496, 0.526667), (0.207325, 0.978248, 
            0.006854)), ((0.536915, 0.010826, 0.526667), (0.208594, 0.977978, 
            0.006881)), ((0.553668, 0.007233, 0.526667), (0.210343, 0.977603, 
            0.006925)), ((0.567521, 0.004274, 0.526667), (0.208057, 0.978093, 
            0.006861)), ((0.583553, 0.001123, 0.526667), (0.187214, 0.982299, 
            0.006234)), ((0.000203, 0.001398, 0.476667), (-0.989669, 0.14337, 0.0)), ((
            0.00151, 0.006476, 0.476667), (-0.929852, 0.367933, 0.0)), ((0.004984, 
            0.013419, 0.476667), (-0.818914, 0.573916, 0.0)), ((0.010758, 0.020575, 
            0.476667), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.476667), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.476667), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.476667), (-0.399193, 0.916867, 
            0.0)), ((0.056319, 0.046917, 0.476667), (-0.318209, 0.94802, 0.0)), ((
            0.073077, 0.052361, 0.476667), (-0.29229, 0.95633, 0.0)), ((0.093229, 
            0.058198, 0.476667), (-0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 
            0.476667), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.476667), (0.0, 1.0, 
            0.0)), ((0.231395, 0.061779, 0.476667), (0.055924, 0.998435, 0.0)), ((
            0.266448, 0.059356, 0.476667), (0.103513, 0.994628, 0.0)), ((0.295248, 
            0.056119, 0.476667), (0.127805, 0.991799, 0.0)), ((0.324299, 0.052166, 
            0.476667), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 0.476667), (
            0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.476667), (0.191844, 
            0.981425, 0.0)), ((0.44176, 0.030833, 0.476667), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.476667), (0.205119, 0.978737, 0.0)), ((
            0.488195, 0.021277, 0.476667), (0.206858, 0.978371, 0.0)), ((0.5108, 
            0.016494, 0.476667), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.476667), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.476667), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.476667), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.476667), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.476667), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.476667), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.476667), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.476667), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.476667), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.476667), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.476667), (0.055871, -0.998438, 
            0.0)), ((0.416659, -0.003417, 0.476667), (0.055871, -0.998438, 0.0)), ((
            0.369303, -0.006193, 0.476667), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.476667), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.476667), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.476667), (
            0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.476667), (0.122438, 
            -0.992476, 0.0)), ((0.208773, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((
            0.129143, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 
            0.476667), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.476667), 
            (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.476667), (-0.074597, 
            -0.997214, 0.0)), ((0.040354, -0.016754, 0.476667), (-0.12352, -0.992342, 
            0.0)), ((0.025938, -0.014654, 0.476667), (-0.194505, -0.980901, 0.0)), ((
            0.014545, -0.012018, 0.476667), (-0.304425, -0.952536, 0.0)), ((0.006398, 
            -0.008883, 0.476667), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 
            0.476667), (-0.934421, -0.35617, 0.0)), ((0.583833, 0.000321, 0.42), (
            -0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 0.42), (-0.022329, 
            -0.999751, 0.0)), ((0.55165, 0.001302, 0.42), (-0.004709, -0.999989, 0.0)), 
            ((0.533034, 0.001276, 0.42), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.42), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.42), 
            (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.42), (0.055871, 
            -0.998438, 0.0)), ((0.416659, -0.003417, 0.42), (0.055871, -0.998438, 
            0.0)), ((0.369303, -0.006193, 0.42), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.42), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.42), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.42), (0.107257, 
            -0.994231, 0.0)), ((0.240235, -0.016966, 0.42), (0.122438, -0.992476, 
            0.0)), ((0.208773, -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.129143, 
            -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.42), (
            -0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.42), (-0.039629, 
            -0.999214, 0.0)), ((0.057635, -0.018302, 0.42), (-0.074597, -0.997214, 
            0.0)), ((0.040354, -0.016754, 0.42), (-0.12352, -0.992342, 0.0)), ((
            0.025938, -0.014654, 0.42), (-0.194505, -0.980901, 0.0)), ((0.014545, 
            -0.012018, 0.42), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 
            0.42), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.42), (
            -0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 0.42), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.42), (-0.929852, 0.367933, 0.0)), ((
            0.004984, 0.013419, 0.42), (-0.818914, 0.573916, 0.0)), ((0.010758, 
            0.020575, 0.42), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.42), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.42), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.42), (-0.399193, 0.916867, 0.0)), 
            ((0.056319, 0.046917, 0.42), (-0.318209, 0.94802, 0.0)), ((0.073077, 
            0.052361, 0.42), (-0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.42), (
            -0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 0.42), (0.0, 1.0, 0.0)), 
            ((0.198897, 0.062493, 0.42), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.42), 
            (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.42), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.42), (0.127805, 0.991799, 0.0)), (
            (0.324299, 0.052166, 0.42), (0.148772, 0.988872, 0.0)), ((0.353335, 
            0.047615, 0.42), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.42), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.42), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.42), (0.205119, 0.978737, 0.0)), ((0.488195, 
            0.021277, 0.42), (0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.42), (
            0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 0.42), (0.208599, 0.978001, 
            0.0)), ((0.548842, 0.00839, 0.42), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.42), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.42), (
            0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.333333), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.333333), (-0.929852, 0.367933, 
            0.0)), ((0.004984, 0.013419, 0.333333), (-0.818914, 0.573916, 0.0)), ((
            0.010758, 0.020575, 0.333333), (-0.701042, 0.71312, 0.0)), ((0.018825, 
            0.027712, 0.333333), (-0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 
            0.333333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 0.333333), (
            -0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.333333), (-0.318209, 
            0.94802, 0.0)), ((0.073077, 0.052361, 0.333333), (-0.29229, 0.95633, 0.0)), 
            ((0.093229, 0.058198, 0.333333), (-0.256858, 0.966449, 0.0)), ((0.119267, 
            0.062493, 0.333333), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.333333), (
            0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.333333), (0.055924, 0.998435, 
            0.0)), ((0.266448, 0.059356, 0.333333), (0.103513, 0.994628, 0.0)), ((
            0.295248, 0.056119, 0.333333), (0.127805, 0.991799, 0.0)), ((0.324299, 
            0.052166, 0.333333), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 
            0.333333), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.333333), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.333333), (0.191844, 
            0.981425, 0.0)), ((0.463659, 0.026434, 0.333333), (0.205119, 0.978737, 
            0.0)), ((0.488195, 0.021277, 0.333333), (0.206858, 0.978371, 0.0)), ((
            0.5108, 0.016494, 0.333333), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.333333), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.333333), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.333333), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.333333), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.333333), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.333333), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.333333), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.333333), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.333333), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.333333), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.333333), (0.055871, -0.998438, 
            0.0)), ((0.416659, -0.003417, 0.333333), (0.055871, -0.998438, 0.0)), ((
            0.369303, -0.006193, 0.333333), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.333333), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.333333), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.333333), (
            0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.333333), (0.122438, 
            -0.992476, 0.0)), ((0.208773, -0.019629, 0.333333), (0.0, -1.0, 0.0)), ((
            0.129143, -0.019629, 0.333333), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 
            0.333333), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.333333), 
            (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.333333), (-0.074597, 
            -0.997214, 0.0)), ((0.040354, -0.016754, 0.333333), (-0.12352, -0.992342, 
            0.0)), ((0.025938, -0.014654, 0.333333), (-0.194505, -0.980901, 0.0)), ((
            0.014545, -0.012018, 0.333333), (-0.304425, -0.952536, 0.0)), ((0.006398, 
            -0.008883, 0.333333), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 
            0.333333), (-0.934421, -0.35617, 0.0)), ((0.583833, 0.000321, 0.233333), (
            -0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 0.233333), (-0.022329, 
            -0.999751, 0.0)), ((0.55165, 0.001302, 0.233333), (-0.004709, -0.999989, 
            0.0)), ((0.533034, 0.001276, 0.233333), (0.011647, -0.999932, 0.0)), ((
            0.511306, 0.000915, 0.233333), (0.025441, -0.999676, 0.0)), ((0.486823, 
            0.000187, 0.233333), (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 
            0.233333), (0.055871, -0.998438, 0.0)), ((0.416659, -0.003417, 0.233333), (
            0.055871, -0.998438, 0.0)), ((0.369303, -0.006193, 0.233333), (0.067641, 
            -0.99771, 0.0)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((
            0.33714, -0.008398, 0.233333), (0.069854, -0.997557, 0.0)), ((0.304702, 
            -0.010667, 0.233333), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 
            0.233333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.233333), (
            0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.233333), (0.0, -1.0, 
            0.0)), ((0.129143, -0.019629, 0.233333), (0.0, -1.0, 0.0)), ((0.101105, 
            -0.019614, 0.233333), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 
            0.233333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.233333), 
            (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.233333), (-0.12352, 
            -0.992342, 0.0)), ((0.025938, -0.014654, 0.233333), (-0.194505, -0.980901, 
            0.0)), ((0.014545, -0.012018, 0.233333), (-0.304425, -0.952536, 0.0)), ((
            0.006398, -0.008883, 0.233333), (-0.509755, -0.86032, 0.0)), ((0.001672, 
            -0.004386, 0.233333), (-0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 
            0.233333), (-0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.233333), (
            -0.929852, 0.367933, 0.0)), ((0.004984, 0.013419, 0.233333), (-0.818914, 
            0.573916, 0.0)), ((0.010758, 0.020575, 0.233333), (-0.701042, 0.71312, 
            0.0)), ((0.018825, 0.027712, 0.233333), (-0.591859, 0.806041, 0.0)), ((
            0.029132, 0.034609, 0.233333), (-0.49082, 0.871261, 0.0)), ((0.041644, 
            0.041062, 0.233333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 
            0.233333), (-0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.233333), (
            -0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.233333), (-0.256858, 
            0.966449, 0.0)), ((0.119267, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((
            0.198897, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 
            0.233333), (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.233333), (
            0.103513, 0.994628, 0.0)), ((0.295248, 0.056119, 0.233333), (0.127805, 
            0.991799, 0.0)), ((0.324299, 0.052166, 0.233333), (0.148772, 0.988872, 
            0.0)), ((0.353335, 0.047615, 0.233333), (0.167009, 0.985955, 0.0)), ((
            0.393406, 0.040285, 0.233333), (0.191844, 0.981425, 0.0)), ((0.44176, 
            0.030833, 0.233333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 
            0.233333), (0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.233333), (
            0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.233333), (0.20733, 
            0.978271, 0.0)), ((0.53112, 0.01218, 0.233333), (0.208599, 0.978001, 0.0)), 
            ((0.548842, 0.00839, 0.233333), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.233333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 
            0.233333), (0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.133333), (
            -0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.133333), (-0.929852, 
            0.367933, 0.0)), ((0.004984, 0.013419, 0.133333), (-0.818914, 0.573916, 
            0.0)), ((0.010758, 0.020575, 0.133333), (-0.701042, 0.71312, 0.0)), ((
            0.018825, 0.027712, 0.133333), (-0.591859, 0.806041, 0.0)), ((0.029132, 
            0.034609, 0.133333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 
            0.133333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.133333), (
            -0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.133333), (-0.29229, 
            0.95633, 0.0)), ((0.093229, 0.058198, 0.133333), (-0.256858, 0.966449, 
            0.0)), ((0.119267, 0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.198897, 
            0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.133333), (
            0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.133333), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.133333), (0.127805, 0.991799, 
            0.0)), ((0.324299, 0.052166, 0.133333), (0.148772, 0.988872, 0.0)), ((
            0.353335, 0.047615, 0.133333), (0.167009, 0.985955, 0.0)), ((0.393406, 
            0.040285, 0.133333), (0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 
            0.133333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 0.133333), (
            0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.133333), (0.206858, 
            0.978371, 0.0)), ((0.5108, 0.016494, 0.133333), (0.20733, 0.978271, 0.0)), 
            ((0.53112, 0.01218, 0.133333), (0.208599, 0.978001, 0.0)), ((0.548842, 
            0.00839, 0.133333), (0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 
            0.133333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.133333), (
            0.187218, 0.982318, 0.0)), ((0.583833, 0.000321, 0.133333), (-0.051907, 
            -0.998652, 0.0)), ((0.566831, 0.001066, 0.133333), (-0.022329, -0.999751, 
            0.0)), ((0.55165, 0.001302, 0.133333), (-0.004709, -0.999989, 0.0)), ((
            0.533034, 0.001276, 0.133333), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.133333), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 
            0.133333), (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.133333), (
            0.055871, -0.998438, 0.0)), ((0.416659, -0.003417, 0.133333), (0.055871, 
            -0.998438, 0.0)), ((0.369303, -0.006193, 0.133333), (0.067641, -0.99771, 
            0.0)), ((0.33714, -0.008398, 0.133333), (0.069854, -0.997557, 0.0)), ((
            0.304702, -0.010667, 0.133333), (0.069674, -0.99757, 0.0)), ((0.272393, 
            -0.01333, 0.133333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 
            0.133333), (0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.133333), (
            0.0, -1.0, 0.0)), ((0.129143, -0.019629, 0.133333), (0.0, -1.0, 0.0)), ((
            0.101105, -0.019614, 0.133333), (-0.001899, -0.999998, 0.0)), ((0.077593, 
            -0.019307, 0.133333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 
            0.133333), (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.133333), 
            (-0.12352, -0.992342, 0.0)), ((0.025938, -0.014654, 0.133333), (-0.194505, 
            -0.980901, 0.0)), ((0.014545, -0.012018, 0.133333), (-0.304425, -0.952536, 
            0.0)), ((0.006398, -0.008883, 0.133333), (-0.509755, -0.86032, 0.0)), ((
            0.001672, -0.004386, 0.133333), (-0.934421, -0.35617, 0.0)), ((0.583833, 
            0.000321, 0.033333), (-0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 
            0.033333), (-0.022329, -0.999751, 0.0)), ((0.55165, 0.001302, 0.033333), (
            -0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 0.033333), (0.011647, 
            -0.999932, 0.0)), ((0.511306, 0.000915, 0.033333), (0.025441, -0.999676, 
            0.0)), ((0.486823, 0.000187, 0.033333), (0.03742, -0.9993, 0.0)), ((
            0.457924, -0.001108, 0.033333), (0.055871, -0.998438, 0.0)), ((0.416659, 
            -0.003417, 0.033333), (0.055871, -0.998438, 0.0)), ((0.369303, -0.006193, 
            0.033333), (0.067641, -0.99771, 0.0)), ((0.33714, -0.008398, 0.033333), (
            0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 0.033333), (0.069674, 
            -0.99757, 0.0)), ((0.272393, -0.01333, 0.033333), (0.107257, -0.994231, 
            0.0)), ((0.240235, -0.016966, 0.033333), (0.122438, -0.992476, 0.0)), ((
            0.208773, -0.019629, 0.033333), (0.0, -1.0, 0.0)), ((0.129143, -0.019629, 
            0.033333), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.033333), (-0.001899, 
            -0.999998, 0.0)), ((0.077593, -0.019307, 0.033333), (-0.039629, -0.999214, 
            0.0)), ((0.057635, -0.018302, 0.033333), (-0.074597, -0.997214, 0.0)), ((
            0.040354, -0.016754, 0.033333), (-0.12352, -0.992342, 0.0)), ((0.025938, 
            -0.014654, 0.033333), (-0.194505, -0.980901, 0.0)), ((0.014545, -0.012018, 
            0.033333), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 0.033333), 
            (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.033333), (-0.934421, 
            -0.35617, 0.0)), ((0.000203, 0.001398, 0.033333), (-0.989669, 0.14337, 
            0.0)), ((0.00151, 0.006476, 0.033333), (-0.929852, 0.367933, 0.0)), ((
            0.004984, 0.013419, 0.033333), (-0.818914, 0.573916, 0.0)), ((0.010758, 
            0.020575, 0.033333), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 
            0.033333), (-0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.033333), (
            -0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 0.033333), (-0.399193, 
            0.916867, 0.0)), ((0.056319, 0.046917, 0.033333), (-0.318209, 0.94802, 
            0.0)), ((0.073077, 0.052361, 0.033333), (-0.29229, 0.95633, 0.0)), ((
            0.093229, 0.058198, 0.033333), (-0.256858, 0.966449, 0.0)), ((0.119267, 
            0.062493, 0.033333), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.033333), (
            0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.033333), (0.055924, 0.998435, 
            0.0)), ((0.266448, 0.059356, 0.033333), (0.103513, 0.994628, 0.0)), ((
            0.295248, 0.056119, 0.033333), (0.127805, 0.991799, 0.0)), ((0.324299, 
            0.052166, 0.033333), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 
            0.033333), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.033333), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.033333), (0.191844, 
            0.981425, 0.0)), ((0.463659, 0.026434, 0.033333), (0.205119, 0.978737, 
            0.0)), ((0.488195, 0.021277, 0.033333), (0.206858, 0.978371, 0.0)), ((
            0.5108, 0.016494, 0.033333), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.033333), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.033333), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.033333), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.033333), (0.187218, 0.982318, 
            0.0)), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='skin-tip', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.222296, -0.054011, 
            2.34181), (-0.921246, 0.06129, 0.384121)), ((0.222992, -0.051389, 2.34232), 
            (-0.87222, 0.272656, 0.406069)), ((0.224821, -0.047822, 2.343013), (
            -0.778269, 0.476044, 0.409487)), ((0.227854, -0.044147, 2.343728), (
            -0.674525, 0.621555, 0.398353)), ((0.232086, -0.040482, 2.34444), (
            -0.575154, 0.723888, 0.38103)), ((0.23749, -0.036942, 2.345128), (
            -0.480724, 0.799158, 0.360904)), ((0.244048, -0.033631, 2.345772), (
            -0.393325, 0.854192, 0.340076)), ((0.251736, -0.030628, 2.346356), (
            -0.314905, 0.893516, 0.3201)), ((0.260514, -0.027832, 2.346899), (
            -0.289604, 0.904327, 0.313563)), ((0.271089, -0.02483, 2.347483), (
            -0.254874, 0.917694, 0.304757)), ((0.280898, -0.022657, 2.347905), (0.0, 
            0.972107, 0.23454)), ((0.334478, -0.022657, 2.347905), (0.0, 0.972107, 
            0.23454)), ((0.343446, -0.023035, 2.347832), (0.055897, 0.973641, 
            0.22113)), ((0.361734, -0.024284, 2.347589), (0.10349, 0.972127, 
            0.210377)), ((0.376803, -0.02595, 2.347265), (0.127788, 0.970349, 
            0.205166)), ((0.392002, -0.027983, 2.34687), (0.14876, 0.968233, 
            0.200986)), ((0.407192, -0.030323, 2.346415), (0.167001, 0.965944, 
            0.197645)), ((0.424035, -0.033301, 2.345836), (0.191841, 0.962164, 
            0.193486)), ((0.44294, -0.036928, 2.345131), (0.191841, 0.962164, 
            0.193486)), ((0.46489, -0.041203, 2.3443), (0.205116, 0.95972, 0.192005)), 
            ((0.477717, -0.043849, 2.343786), (0.206856, 0.959379, 0.19184)), ((
            0.489531, -0.046303, 2.343309), (0.207327, 0.959285, 0.191803)), ((
            0.500148, -0.048516, 2.342878), (0.208596, 0.959025, 0.191724)), ((
            0.509405, -0.05046, 2.342501), (0.210345, 0.958661, 0.191637)), ((0.519289, 
            -0.052536, 2.342097), (0.208059, 0.959142, 0.191723)), ((0.527678, 
            -0.054156, 2.341782), (0.187215, 0.963321, 0.192257)), ((0.52753, 
            -0.054581, 2.3417), (-0.051907, -0.980043, -0.191889)), ((0.518659, 
            -0.054201, 2.341773), (-0.022329, -0.981195, -0.191722)), ((0.507824, 
            -0.054068, 2.341799), (-0.004709, -0.981439, -0.191717)), ((0.497539, 
            -0.054135, 2.341786), (0.011647, -0.981353, -0.19186)), ((0.485685, 
            -0.054381, 2.341739), (0.025441, -0.981037, -0.192143)), ((0.472452, 
            -0.054817, 2.341654), (0.03742, -0.980572, -0.192558)), ((0.445475, 
            -0.056214, 2.341382), (0.055871, -0.979511, -0.193484)), ((0.426566, 
            -0.057252, 2.34118), (0.055871, -0.979511, -0.193484)), ((0.409977, 
            -0.058291, 2.340978), (0.067641, -0.978511, -0.194782)), ((0.393081, 
            -0.05944, 2.340755), (0.069854, -0.978298, -0.195074)), ((0.37611, 
            -0.060605, 2.340529), (0.069673, -0.978316, -0.195046)), ((0.359274, 
            -0.06217, 2.340224), (0.107249, -0.973584, -0.201572)), ((0.342394, 
            -0.064127, 2.339844), (0.122423, -0.971186, -0.204476)), ((0.334478, 
            -0.064831, 2.339707), (0.0, -0.984214, -0.176982)), ((0.280898, -0.064831, 
            2.339707), (0.0, -0.984214, -0.176982)), ((0.270879, -0.064815, 2.33971), (
            -0.001899, -0.984319, -0.176387)), ((0.259259, -0.064528, 2.339766), (
            -0.039614, -0.985704, -0.163761)), ((0.249259, -0.063923, 2.339884), (
            -0.074539, -0.985673, -0.151302)), ((0.240694, -0.063036, 2.340056), (
            -0.123318, -0.983434, -0.132851)), ((0.233654, -0.061866, 2.340283), (
            -0.193825, -0.975453, -0.104511)), ((0.228223, -0.060429, 2.340563), (
            -0.302046, -0.951537, -0.05784)), ((0.224547, -0.058739, 2.340891), (
            -0.499339, -0.865632, 0.036618)), ((0.222638, -0.055911, 2.341441), (
            -0.876188, -0.394143, 0.277389)), ((0.529442, -0.037388, 2.253371), (
            -0.051907, -0.980043, -0.191889)), ((0.519487, -0.036963, 2.253453), (
            -0.022329, -0.981195, -0.191722)), ((0.507386, -0.036814, 2.253482), (
            -0.004709, -0.981439, -0.191717)), ((0.49584, -0.036888, 2.253468), (
            0.011647, -0.981353, -0.19186)), ((0.482531, -0.037164, 2.253414), (
            0.025441, -0.981037, -0.192143)), ((0.467672, -0.037653, 2.253319), (
            0.03742, -0.980572, -0.192558)), ((0.443431, -0.038888, 2.253079), (
            0.055871, -0.979511, -0.193484)), ((0.419603, -0.040196, 2.252825), (
            0.055871, -0.979511, -0.193484)), ((0.397508, -0.041553, 2.252561), (
            0.067641, -0.978511, -0.194782)), ((0.37853, -0.042843, 2.25231), (
            0.069854, -0.978298, -0.195074)), ((0.359467, -0.044151, 2.252056), (
            0.069673, -0.978316, -0.195046)), ((0.340554, -0.045907, 2.251715), (
            0.107249, -0.973584, -0.201572)), ((0.321594, -0.048104, 2.251288), (
            0.122423, -0.971186, -0.204476)), ((0.310697, -0.048904, 2.251132), (0.0, 
            -0.984214, -0.176982)), ((0.253593, -0.048904, 2.251132), (0.0, -0.984214, 
            -0.176982)), ((0.241246, -0.048886, 2.251136), (-0.001899, -0.984319, 
            -0.176387)), ((0.228186, -0.048565, 2.251198), (-0.039614, -0.985704, 
            -0.163761)), ((0.216947, -0.047887, 2.25133), (-0.074539, -0.985673, 
            -0.151302)), ((0.20732, -0.046891, 2.251523), (-0.123318, -0.983434, 
            -0.132851)), ((0.199406, -0.045579, 2.251779), (-0.193825, -0.975453, 
            -0.104511)), ((0.193299, -0.043965, 2.252092), (-0.302046, -0.951537, 
            -0.05784)), ((0.189162, -0.042067, 2.252461), (-0.499339, -0.865632, 
            0.036618)), ((0.18701, -0.038898, 2.253077), (-0.876188, -0.394143, 
            0.277389)), ((0.186622, -0.03674, 2.253497), (-0.921246, 0.06129, 
            0.384121)), ((0.187408, -0.03379, 2.25407), (-0.87222, 0.272656, 
            0.406069)), ((0.189468, -0.029783, 2.254849), (-0.778269, 0.476044, 
            0.409487)), ((0.192879, -0.025654, 2.255652), (-0.674525, 0.621555, 
            0.398353)), ((0.197639, -0.021538, 2.256452), (-0.575154, 0.723888, 
            0.38103)), ((0.203714, -0.017563, 2.257224), (-0.480724, 0.799158, 
            0.360904)), ((0.211084, -0.013845, 2.257947), (-0.393325, 0.854192, 
            0.340076)), ((0.219725, -0.010473, 2.258603), (-0.314905, 0.893516, 
            0.3201)), ((0.22959, -0.007331, 2.259213), (-0.289604, 0.904327, 
            0.313563)), ((0.241482, -0.003957, 2.259869), (-0.254874, 0.917694, 
            0.304757)), ((0.253593, -0.001531, 2.260341), (0.0, 0.972107, 0.23454)), ((
            0.310697, -0.001531, 2.260341), (0.0, 0.972107, 0.23454)), ((0.32279, 
            -0.001959, 2.260257), (0.055897, 0.973641, 0.22113)), ((0.343311, 
            -0.003364, 2.259984), (0.10349, 0.972127, 0.210377)), ((0.360238, 
            -0.005237, 2.25962), (0.127788, 0.970349, 0.205166)), ((0.377311, 
            -0.007523, 2.259176), (0.14876, 0.968233, 0.200986)), ((0.394373, 
            -0.010152, 2.258665), (0.167001, 0.965944, 0.197645)), ((0.416776, 
            -0.014167, 2.257884), (0.191841, 0.962164, 0.193486)), ((0.440553, 
            -0.018729, 2.256998), (0.191841, 0.962164, 0.193486)), ((0.459175, 
            -0.022373, 2.256289), (0.205116, 0.95972, 0.192005)), ((0.473579, 
            -0.025345, 2.255712), (0.206856, 0.959379, 0.19184)), ((0.486846, 
            -0.028101, 2.255176), (0.207327, 0.959285, 0.191803)), ((0.498766, 
            -0.030586, 2.254693), (0.208596, 0.959025, 0.191724)), ((0.509158, 
            -0.032768, 2.254269), (0.210345, 0.958661, 0.191637)), ((0.520198, 
            -0.035087, 2.253818), (0.208059, 0.959142, 0.191723)), ((0.52961, 
            -0.036906, 2.253464), (0.187215, 0.963321, 0.192257)), ((0.146234, 
            -0.017202, 2.153516), (-0.921246, 0.06129, 0.384121)), ((0.147115, 
            -0.013888, 2.15416), (-0.87222, 0.272656, 0.406069)), ((0.149429, 
            -0.009383, 2.155036), (-0.778269, 0.476044, 0.409487)), ((0.153262, 
            -0.004741, 2.155938), (-0.674525, 0.621555, 0.398353)), ((0.15861, 
            -0.000114, 2.156838), (-0.575154, 0.723888, 0.38103)), ((0.165438, 
            0.004356, 2.157707), (-0.480724, 0.799158, 0.360904)), ((0.173722, 
            0.008536, 2.158519), (-0.393325, 0.854192, 0.340076)), ((0.183435, 
            0.012328, 2.159256), (-0.314905, 0.893516, 0.3201)), ((0.194523, 0.015859, 
            2.159943), (-0.289604, 0.904327, 0.313563)), ((0.207886, 0.019652, 
            2.16068), (-0.254874, 0.917694, 0.304757)), ((0.222507, 0.022386, 
            2.161211), (0.0, 0.972107, 0.23454)), ((0.283599, 0.022386, 2.161211), (
            0.0, 0.972107, 0.23454)), ((0.299281, 0.021906, 2.161118), (0.055897, 
            0.973641, 0.22113)), ((0.322362, 0.020328, 2.160811), (0.10349, 0.972127, 
            0.210377)), ((0.341392, 0.018223, 2.160402), (0.127788, 0.970349, 
            0.205166)), ((0.360586, 0.015654, 2.159903), (0.14876, 0.968233, 
            0.200986)), ((0.379768, 0.012699, 2.159328), (0.167001, 0.965944, 
            0.197645)), ((0.407978, 0.007606, 2.158338), (0.191841, 0.962164, 
            0.193486)), ((0.438166, 0.001813, 2.157212), (0.191841, 0.962164, 
            0.193486)), ((0.452625, -0.001041, 2.156657), (0.205116, 0.95972, 
            0.192005)), ((0.468821, -0.004382, 2.156008), (0.206856, 0.959379, 
            0.19184)), ((0.483737, -0.007481, 2.155406), (0.207327, 0.959285, 
            0.191803)), ((0.497142, -0.010275, 2.154862), (0.208596, 0.959025, 
            0.191724)), ((0.508827, -0.012728, 2.154386), (0.210345, 0.958661, 
            0.191637)), ((0.521267, -0.015342, 2.153877), (0.208059, 0.959142, 
            0.191723)), ((0.531854, -0.017387, 2.15348), (0.187215, 0.963321, 
            0.192257)), ((0.531666, -0.017927, 2.153375), (-0.051907, -0.980043, 
            -0.191889)), ((0.52047, -0.017448, 2.153468), (-0.022329, -0.981195, 
            -0.191722)), ((0.506833, -0.017281, 2.153501), (-0.004709, -0.981439, 
            -0.191717)), ((0.493849, -0.017365, 2.153484), (0.011647, -0.981353, 
            -0.19186)), ((0.478884, -0.017674, 2.153424), (0.025441, -0.981037, 
            -0.192143)), ((0.459047, -0.01834, 2.153295), (0.03742, -0.980572, 
            -0.192558)), ((0.441387, -0.019258, 2.153116), (0.055871, -0.979511, 
            -0.193484)), ((0.411163, -0.020919, 2.152794), (0.055871, -0.979511, 
            -0.193484)), ((0.38329, -0.022611, 2.152465), (0.067641, -0.978511, 
            -0.194782)), ((0.361953, -0.024061, 2.152183), (0.069854, -0.978298, 
            -0.195074)), ((0.340522, -0.025532, 2.151897), (0.069673, -0.978316, 
            -0.195046)), ((0.31926, -0.027507, 2.151513), (0.107249, -0.973584, 
            -0.201572)), ((0.297944, -0.029977, 2.151033), (0.122423, -0.971186, 
            -0.204476)), ((0.283599, -0.030872, 2.150859), (0.0, -0.984214, 
            -0.176982)), ((0.222507, -0.030872, 2.150859), (0.0, -0.984214, 
            -0.176982)), ((0.207621, -0.030852, 2.150863), (-0.001899, -0.984319, 
            -0.176387)), ((0.192941, -0.030491, 2.150933), (-0.039614, -0.985704, 
            -0.163761)), ((0.180309, -0.029728, 2.151081), (-0.074539, -0.985673, 
            -0.151302)), ((0.169489, -0.028608, 2.151299), (-0.123318, -0.983434, 
            -0.132851)), ((0.160594, -0.027132, 2.151586), (-0.193825, -0.975453, 
            -0.104511)), ((0.153732, -0.025317, 2.151939), (-0.302046, -0.951537, 
            -0.05784)), ((0.149085, -0.023183, 2.152353), (-0.499339, -0.865632, 
            0.036618)), ((0.146668, -0.019617, 2.153047), (-0.876188, -0.394143, 
            0.277389)), ((0.533357, -0.003653, 2.08), (-0.051251, -0.980093, 
            -0.191811)), ((0.521236, -0.003128, 2.080042), (-0.022048, -0.981247, 
            -0.19149)), ((0.506371, -0.002944, 2.080059), (-0.00465, -0.981495, 
            -0.19143)), ((0.492322, -0.003037, 2.080049), (0.0115, -0.981406, 
            -0.191598)), ((0.476133, -0.003376, 2.080023), (0.025119, -0.981078, 
            -0.191977)), ((0.458062, -0.003979, 2.07998), (0.036944, -0.980588, 
            -0.192566)), ((0.44016, -0.004856, 2.079912), (0.055216, -0.979493, 
            -0.193767)), ((0.404411, -0.006839, 2.079763), (0.055109, -0.979381, 
            -0.194358)), ((0.372759, -0.008782, 2.079668), (0.066745, -0.978281, 
            -0.196243)), ((0.349691, -0.01037, 2.079569), (0.068922, -0.977971, 
            -0.197034)), ((0.326523, -0.01198, 2.07947), (0.068739, -0.977886, 
            -0.197517)), ((0.303539, -0.014143, 2.07932), (0.105703, -0.973066, 
            -0.204866)), ((0.280494, -0.016848, 2.079146), (0.120592, -0.970497, 
            -0.208791)), ((0.26355, -0.017825, 2.079145), (0.0, -0.983527, -0.180763)), 
            ((0.199533, -0.017825, 2.079145), (0.0, -0.983527, -0.180763)), ((0.182881, 
            -0.017803, 2.079147), (-0.001868, -0.983637, -0.18015)), ((0.167028, 
            -0.017406, 2.079183), (-0.039247, -0.985181, -0.166968)), ((0.153382, 
            -0.016569, 2.079242), (-0.073977, -0.985317, -0.153878)), ((0.141696, 
            -0.015343, 2.079326), (-0.122696, -0.983282, -0.134545)), ((0.132092, 
            -0.013724, 2.079433), (-0.193559, -0.975457, -0.104969)), ((0.124686, 
            -0.011734, 2.079563), (-0.303316, -0.95122, -0.05638)), ((0.119679, 
            -0.009393, 2.079713), (-0.505936, -0.861534, 0.042276)), ((0.117104, 
            -0.005474, 2.079991), (-0.887331, -0.358133, 0.290489)), ((0.116775, 
            -0.001931, 2.080174), (-0.914895, 0.109414, 0.388581)), ((0.117527, 
            0.000766, 2.080198), (-0.858275, 0.316949, 0.403619)), ((0.120023, 
            0.005703, 2.080496), (-0.760791, 0.511954, 0.398873)), ((0.124164, 
            0.010788, 2.080805), (-0.657481, 0.649404, 0.382091)), ((0.129942, 
            0.015855, 2.081116), (-0.560265, 0.745443, 0.361133)), ((0.137319, 
            0.020748, 2.08142), (-0.468525, 0.815872, 0.338875)), ((0.14627, 0.025323, 
            2.081707), (-0.383785, 0.867286, 0.317056)), ((0.156764, 0.029473, 
            2.081971), (-0.307708, 0.903978, 0.296882)), ((0.168743, 0.033339, 
            2.082206), (-0.283224, 0.91445, 0.289078)), ((0.183166, 0.037491, 
            2.082447), (-0.249498, 0.927212, 0.279334)), ((0.199533, 0.040483, 
            2.082723), (0.0, 0.97607, 0.217455)), ((0.26355, 0.040483, 2.082723), (0.0, 
            0.97607, 0.217455)), ((0.281918, 0.03996, 2.082722), (0.055107, 0.977006, 
            0.205969)), ((0.306906, 0.038233, 2.082628), (0.102133, 0.975058, 
            0.197057)), ((0.327477, 0.035929, 2.082498), (0.12617, 0.973036, 
            0.193084)), ((0.348226, 0.033117, 2.082335), (0.146928, 0.970706, 
            0.190108)), ((0.368963, 0.029882, 2.082144), (0.164987, 0.968222, 
            0.187948)), ((0.400935, 0.023997, 2.081929), (0.189443, 0.964183, 
            0.18564)), ((0.436733, 0.017072, 2.081346), (0.190042, 0.963722, 
            0.187414)), ((0.447739, 0.014841, 2.081227), (0.202669, 0.961281, 
            0.186717)), ((0.465254, 0.011182, 2.080997), (0.204365, 0.960759, 
            0.187555)), ((0.481388, 0.007789, 2.080782), (0.204806, 0.960488, 
            0.188459)), ((0.495888, 0.00473, 2.080587), (0.206038, 0.96007, 0.189246)), 
            ((0.508532, 0.002043, 2.080414), (0.207745, 0.959566, 0.189934)), ((
            0.519148, -0.000218, 2.080265), (0.205467, 0.95991, 0.190675)), ((0.533558, 
            -0.003064, 2.080001), (0.184843, 0.963823, 0.192039)), ((0.106671, 
            0.001029, 2.033333), (-0.987308, 0.143028, 0.069039)), ((0.107632, 
            0.004756, 2.033333), (-0.927883, 0.367154, 0.065055)), ((0.110181, 
            0.009848, 2.033333), (-0.817545, 0.572957, 0.057801)), ((0.114418, 
            0.015097, 2.033333), (-0.700158, 0.71222, 0.050206)), ((0.120337, 0.020332, 
            2.033333), (-0.591304, 0.805285, 0.043305)), ((0.127899, 0.02539, 
            2.033333), (-0.490483, 0.870662, 0.03707)), ((0.137078, 0.030123, 
            2.033333), (-0.398994, 0.91641, 0.031577)), ((0.147843, 0.034417, 
            2.033333), (-0.318094, 0.947678, 0.026888)), ((0.160137, 0.03841, 
            2.033333), (-0.292195, 0.95602, 0.025454)), ((0.174922, 0.042692, 
            2.033333), (-0.256786, 0.96618, 0.023606)), ((0.191823, 0.045837, 
            2.033333), (0.0, 0.99994, 0.010933)), ((0.256893, 0.045837, 2.033333), (
            0.0, 0.99994, 0.010933)), ((0.276275, 0.045312, 2.033333), (0.055922, 
            0.998393, 0.009145)), ((0.301978, 0.043534, 2.033333), (0.10351, 0.994597, 
            0.007916)), ((0.323103, 0.041159, 2.033333), (0.127802, 0.991772, 
            0.007403)), ((0.344411, 0.038259, 2.033333), (0.148768, 0.988847, 
            0.007063)), ((0.365708, 0.034921, 2.033333), (0.167005, 0.985932, 
            0.006859)), ((0.398236, 0.028931, 2.033333), (0.19184, 0.981403, 
            0.006707)), ((0.438597, 0.021041, 2.033333), (0.19184, 0.981403, 
            0.006707)), ((0.452733, 0.018105, 2.033333), (0.205114, 0.978714, 
            0.006823)), ((0.470289, 0.014404, 2.033333), (0.206854, 0.978348, 
            0.006846)), ((0.48634, 0.011006, 2.033333), (0.207325, 0.978248, 
            0.006854)), ((0.500631, 0.007964, 2.033333), (0.208594, 0.977978, 
            0.006881)), ((0.512955, 0.005321, 2.033333), (0.210343, 0.977603, 
            0.006925)), ((0.523146, 0.003144, 2.033333), (0.208057, 0.978093, 
            0.006861)), ((0.534939, 0.000826, 2.033333), (0.187214, 0.982299, 
            0.006234)), ((0.534736, 0.000236, 2.033333), (-0.051907, -0.99865, 
            -0.001728)), ((0.522269, 0.000782, 2.033333), (-0.022329, -0.99975, 
            -0.000839)), ((0.511132, 0.000955, 2.033333), (-0.004709, -0.999989, 
            -0.000353)), ((0.497475, 0.000936, 2.033333), (0.011647, -0.999932, 
            4.8e-05)), ((0.481536, 0.00067, 2.033333), (0.025441, -0.999676, 
            0.000337)), ((0.463576, 0.000136, 2.033333), (0.03742, -0.999299, 
            0.000538)), ((0.445484, -0.00064, 2.033333), (0.055871, -0.998438, 
            0.000763)), ((0.401834, -0.003082, 2.033333), (0.055871, -0.998438, 
            0.000763)), ((0.369593, -0.005072, 2.033333), (0.067641, -0.997709, 
            0.000722)), ((0.345901, -0.006713, 2.033333), (0.069854, -0.997557, 
            0.000701)), ((0.322108, -0.008376, 2.033333), (0.069674, -0.99757, 
            0.000704)), ((0.298506, -0.01062, 2.033333), (0.107257, -0.994231, 
            -7e-05)), ((0.274838, -0.013418, 2.033333), (0.122438, -0.992476, 
            -0.000469)), ((0.256893, -0.014398, 2.033333), (0.0, -0.999994, 0.003434)), 
            ((0.191823, -0.014398, 2.033333), (0.0, -0.999994, 0.003434)), ((0.174623, 
            -0.014375, 2.033333), (-0.001899, -0.999992, 0.003531)), ((0.158359, 
            -0.01396, 2.033333), (-0.039628, -0.999199, 0.005608)), ((0.144356, 
            -0.013092, 2.033333), (-0.074595, -0.997185, 0.007657)), ((0.132366, 
            -0.011821, 2.033333), (-0.123513, -0.992286, 0.010673)), ((0.122516, 
            -0.010147, 2.033333), (-0.194483, -0.980788, 0.015231)), ((0.114923, 
            -0.00809, 2.033333), (-0.304348, -0.952295, 0.02251)), ((0.109792, 
            -0.005673, 2.033333), (-0.509417, -0.85975, 0.036388)), ((0.107137, 
            -0.001614, 2.033333), (-0.932433, -0.355412, 0.065202)), ((0.537292, 
            0.000241, 1.953333), (-0.051907, -0.99865, -0.001728)), ((0.524593, 
            0.000797, 1.953333), (-0.022329, -0.99975, -0.000839)), ((0.513239, 
            0.000973, 1.953333), (-0.004709, -0.999989, -0.000353)), ((0.499317, 
            0.000954, 1.953333), (0.011647, -0.999932, 4.8e-05)), ((0.48307, 0.000682, 
            1.953333), (0.025441, -0.999676, 0.000337)), ((0.464765, 0.000138, 
            1.953333), (0.03742, -0.999299, 0.000538)), ((0.446053, -0.000669, 
            1.953333), (0.055871, -0.998438, 0.000763)), ((0.401673, -0.003152, 
            1.953333), (0.055871, -0.998438, 0.000763)), ((0.369055, -0.005166, 
            1.953333), (0.067641, -0.997709, 0.000722)), )+\
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.344912, -0.006839, 
            1.953333), (0.069854, -0.997557, 0.000701)), ((0.320666, -0.008534, 
            1.953333), (0.069674, -0.99757, 0.000704)), ((0.296612, -0.010819, 
            1.953333), (0.107257, -0.994231, -7e-05)), ((0.272493, -0.013669, 
            1.953333), (0.122438, -0.992476, -0.000469)), ((0.253913, -0.014672, 
            1.953333), (0.0, -0.999994, 0.003434)), ((0.188079, -0.014672, 1.953333), (
            0.0, -0.999994, 0.003434)), ((0.170357, -0.014649, 1.953333), (-0.001899, 
            -0.999992, 0.003531)), ((0.153778, -0.014227, 1.953333), (-0.039628, 
            -0.999199, 0.005608)), ((0.139505, -0.013343, 1.953333), (-0.074595, 
            -0.997185, 0.007657)), ((0.127283, -0.012049, 1.953333), (-0.123513, 
            -0.992286, 0.010673)), ((0.117241, -0.010344, 1.953333), (-0.194483, 
            -0.980788, 0.015231)), ((0.1095, -0.008248, 1.953333), (-0.304348, 
            -0.952295, 0.02251)), ((0.104267, -0.005785, 1.953333), (-0.509417, 
            -0.85975, 0.036388)), ((0.101557, -0.001652, 1.953333), (-0.932433, 
            -0.355412, 0.065202)), ((0.10108, 0.001053, 1.953333), (-0.987308, 
            0.143028, 0.069039)), ((0.102062, 0.004854, 1.953333), (-0.927883, 
            0.367154, 0.065055)), ((0.104662, 0.010044, 1.953333), (-0.817545, 
            0.572957, 0.057801)), ((0.108983, 0.015393, 1.953333), (-0.700158, 0.71222, 
            0.050206)), ((0.115017, 0.020728, 1.953333), (-0.591304, 0.805285, 
            0.043305)), ((0.122726, 0.025883, 1.953333), (-0.490483, 0.870662, 
            0.03707)), ((0.132083, 0.030705, 1.953333), (-0.398994, 0.91641, 
            0.031577)), ((0.143056, 0.03508, 1.953333), (-0.318094, 0.947678, 
            0.026888)), ((0.155586, 0.039149, 1.953333), (-0.292195, 0.95602, 
            0.025454)), ((0.170661, 0.043514, 1.953333), (-0.256786, 0.96618, 
            0.023606)), ((0.188079, 0.046712, 1.953333), (0.0, 0.99994, 0.010933)), ((
            0.253913, 0.046712, 1.953333), (0.0, 0.99994, 0.010933)), ((0.273964, 
            0.046174, 1.953333), (0.055922, 0.998393, 0.009145)), ((0.300147, 0.044361, 
            1.953333), (0.10351, 0.994597, 0.007916)), ((0.321675, 0.04194, 1.953333), 
            (0.127802, 0.991772, 0.007403)), ((0.34339, 0.038984, 1.953333), (0.148768, 
            0.988847, 0.007063)), ((0.365093, 0.035582, 1.953333), (0.167005, 0.985932, 
            0.006859)), ((0.398015, 0.029521, 1.953333), (0.19184, 0.981403, 
            0.006707)), ((0.439051, 0.021499, 1.953333), (0.19184, 0.981403, 
            0.006707)), ((0.453717, 0.018456, 1.953333), (0.205114, 0.978714, 
            0.006823)), ((0.471609, 0.014685, 1.953333), (0.206854, 0.978348, 
            0.006846)), ((0.487969, 0.011221, 1.953333), (0.207325, 0.978248, 
            0.006854)), ((0.502536, 0.008121, 1.953333), (0.208594, 0.977978, 
            0.006881)), ((0.515099, 0.005427, 1.953333), (0.210343, 0.977603, 
            0.006925)), ((0.525487, 0.003207, 1.953333), (0.208057, 0.978093, 
            0.006861)), ((0.5375, 0.000846, 1.953333), (0.187214, 0.982299, 0.006234)), 
            ((0.094091, 0.001078, 1.853333), (-0.987308, 0.143028, 0.069039)), ((
            0.095096, 0.004967, 1.853333), (-0.927883, 0.367154, 0.065055)), ((
            0.097757, 0.010279, 1.853333), (-0.817545, 0.572957, 0.057801)), ((
            0.102179, 0.015753, 1.853333), (-0.700158, 0.71222, 0.050206)), ((0.108354, 
            0.021213, 1.853333), (-0.591304, 0.805285, 0.043305)), ((0.116243, 
            0.026488, 1.853333), (-0.490483, 0.870662, 0.03707)), ((0.125819, 0.031423, 
            1.853333), (-0.398994, 0.91641, 0.031577)), ((0.137049, 0.035901, 
            1.853333), (-0.318094, 0.947678, 0.026888)), ((0.149872, 0.040065, 
            1.853333), (-0.292195, 0.95602, 0.025454)), ((0.1653, 0.044532, 1.853333), 
            (-0.256786, 0.96618, 0.023606)), ((0.183318, 0.047805, 1.853333), (0.0, 
            0.99994, 0.010933)), ((0.250109, 0.047805, 1.853333), (0.0, 0.99994, 
            0.010933)), ((0.27102, 0.047255, 1.853333), (0.055922, 0.998393, 
            0.009145)), ((0.297816, 0.0454, 1.853333), (0.10351, 0.994597, 0.007916)), 
            ((0.319848, 0.042922, 1.853333), (0.127802, 0.991772, 0.007403)), ((
            0.342071, 0.039897, 1.853333), (0.148768, 0.988847, 0.007063)), ((0.364283, 
            0.036415, 1.853333), (0.167005, 0.985932, 0.006859)), ((0.397699, 0.030266, 
            1.853333), (0.19184, 0.981403, 0.006707)), ((0.439699, 0.022056, 1.853333), 
            (0.19184, 0.981403, 0.006707)), ((0.454983, 0.018888, 1.853333), (0.205114, 
            0.978714, 0.006823)), ((0.473294, 0.015028, 1.853333), (0.206854, 0.978348, 
            0.006846)), ((0.490037, 0.011483, 1.853333), (0.207325, 0.978248, 
            0.006854)), ((0.504945, 0.008311, 1.853333), (0.208594, 0.977978, 
            0.006881)), ((0.517801, 0.005554, 1.853333), (0.210343, 0.977603, 
            0.006925)), ((0.528433, 0.003282, 1.853333), (0.208057, 0.978093, 
            0.006861)), ((0.540727, 0.000866, 1.853333), (0.187214, 0.982299, 
            0.006234)), ((0.540514, 0.000247, 1.853333), (-0.051907, -0.99865, 
            -0.001728)), ((0.527517, 0.000816, 1.853333), (-0.022329, -0.99975, 
            -0.000839)), ((0.515898, 0.000996, 1.853333), (-0.004709, -0.999989, 
            -0.000353)), ((0.50165, 0.000976, 1.853333), (0.011647, -0.999932, 
            4.8e-05)), ((0.485024, 0.000698, 1.853333), (0.025441, -0.999676, 
            0.000337)), ((0.46629, 0.000141, 1.853333), (0.03742, -0.999299, 
            0.000538)), ((0.446866, -0.0007, 1.853333), (0.055871, -0.998438, 
            0.000763)), ((0.401443, -0.003242, 1.853333), (0.055871, -0.998438, 
            0.000763)), ((0.368337, -0.005287, 1.853333), (0.067641, -0.997709, 
            0.000722)), ((0.343629, -0.006999, 1.853333), (0.069854, -0.997557, 
            0.000701)), ((0.318815, -0.008733, 1.853333), (0.069674, -0.99757, 
            0.000704)), ((0.294198, -0.011072, 1.853333), (0.107257, -0.994231, 
            -7e-05)), ((0.269515, -0.01399, 1.853333), (0.122438, -0.992476, 
            -0.000469)), ((0.250109, -0.015016, 1.853333), (0.0, -0.999994, 0.003434)), 
            ((0.183318, -0.015016, 1.853333), (0.0, -0.999994, 0.003434)), ((0.164988, 
            -0.014992, 1.853333), (-0.001899, -0.999992, 0.003531)), ((0.148021, 
            -0.01456, 1.853333), (-0.039628, -0.999199, 0.005608)), ((0.133414, 
            -0.013655, 1.853333), (-0.074595, -0.997185, 0.007657)), ((0.120907, 
            -0.012331, 1.853333), (-0.123513, -0.992286, 0.010673)), ((0.11063, 
            -0.010586, 1.853333), (-0.194483, -0.980788, 0.015231)), ((0.102707, 
            -0.008441, 1.853333), (-0.304348, -0.952295, 0.02251)), ((0.097352, 
            -0.005921, 1.853333), (-0.509417, -0.85975, 0.036388)), ((0.094579, 
            -0.00169, 1.853333), (-0.932433, -0.355412, 0.065202)), ((0.543736, 
            0.000253, 1.753333), (-0.051907, -0.99865, -0.001728)), ((0.530442, 
            0.000835, 1.753333), (-0.022329, -0.99975, -0.000839)), ((0.518557, 
            0.001019, 1.753333), (-0.004709, -0.999989, -0.000353)), ((0.503983, 
            0.000998, 1.753333), (0.011647, -0.999932, 4.8e-05)), ((0.486977, 0.000714, 
            1.753333), (0.025441, -0.999676, 0.000337)), ((0.467815, 0.000144, 
            1.753333), (0.03742, -0.999299, 0.000538)), ((0.447679, -0.000731, 
            1.753333), (0.055871, -0.998438, 0.000763)), ((0.401213, -0.003331, 
            1.753333), (0.055871, -0.998438, 0.000763)), ((0.367619, -0.005408, 
            1.753333), (0.067641, -0.997709, 0.000722)), ((0.342346, -0.007159, 
            1.753333), (0.069854, -0.997557, 0.000701)), ((0.316964, -0.008933, 
            1.753333), (0.069674, -0.99757, 0.000704)), ((0.291784, -0.011325, 
            1.753333), (0.107257, -0.994231, -7e-05)), ((0.266536, -0.01431, 1.753333), 
            (0.122438, -0.992476, -0.000469)), ((0.246304, -0.015359, 1.753333), (0.0, 
            -0.999994, 0.003434)), ((0.178558, -0.015359, 1.753333), (0.0, -0.999994, 
            0.003434)), ((0.15962, -0.015335, 1.753333), (-0.001899, -0.999992, 
            0.003531)), ((0.142265, -0.014893, 1.753333), (-0.039628, -0.999199, 
            0.005608)), ((0.127323, -0.013968, 1.753333), (-0.074595, -0.997185, 
            0.007657)), ((0.11453, -0.012613, 1.753333), (-0.123513, -0.992286, 
            0.010673)), ((0.104019, -0.010828, 1.753333), (-0.194483, -0.980788, 
            0.015231)), ((0.095915, -0.008634, 1.753333), (-0.304348, -0.952295, 
            0.02251)), ((0.090438, -0.006056, 1.753333), (-0.509417, -0.85975, 
            0.036388)), ((0.087601, -0.001729, 1.753333), (-0.932433, -0.355412, 
            0.065202)), ((0.087102, 0.001102, 1.753333), (-0.987308, 0.143028, 
            0.069039)), ((0.088129, 0.005081, 1.753333), (-0.927883, 0.367154, 
            0.065055)), ((0.090852, 0.010513, 1.753333), (-0.817545, 0.572957, 
            0.057801)), ((0.095374, 0.016113, 1.753333), (-0.700158, 0.71222, 
            0.050206)), ((0.101691, 0.021698, 1.753333), (-0.591304, 0.805285, 
            0.043305)), ((0.10976, 0.027094, 1.753333), (-0.490483, 0.870662, 
            0.03707)), ((0.119555, 0.032142, 1.753333), (-0.398994, 0.91641, 
            0.031577)), ((0.131042, 0.036721, 1.753333), (-0.318094, 0.947678, 
            0.026888)), ((0.144158, 0.040981, 1.753333), (-0.292195, 0.95602, 
            0.025454)), ((0.159938, 0.04555, 1.753333), (-0.256786, 0.96618, 
            0.023606)), ((0.178558, 0.048899, 1.753333), (0.0, 0.99994, 0.010933)), ((
            0.246304, 0.048899, 1.753333), (0.0, 0.99994, 0.010933)), ((0.268076, 
            0.048336, 1.753333), (0.055922, 0.998393, 0.009145)), ((0.295485, 0.046438, 
            1.753333), (0.10351, 0.994597, 0.007916)), ((0.318021, 0.043904, 1.753333), 
            (0.127802, 0.991772, 0.007403)), ((0.340753, 0.04081, 1.753333), (0.148768, 
            0.988847, 0.007063)), ((0.363472, 0.037248, 1.753333), (0.167005, 0.985932, 
            0.006859)), ((0.397384, 0.031011, 1.753333), (0.19184, 0.981403, 
            0.006707)), ((0.440348, 0.022612, 1.753333), (0.19184, 0.981403, 
            0.006707)), ((0.456249, 0.01932, 1.753333), (0.205114, 0.978714, 
            0.006823)), ((0.474979, 0.015372, 1.753333), (0.206854, 0.978348, 
            0.006846)), ((0.492104, 0.011746, 1.753333), (0.207325, 0.978248, 
            0.006854)), ((0.507353, 0.008501, 1.753333), (0.208594, 0.977978, 
            0.006881)), ((0.520503, 0.00568, 1.753333), (0.210343, 0.977603, 
            0.006925)), ((0.531378, 0.003357, 1.753333), (0.208057, 0.978093, 
            0.006861)), ((0.543953, 0.000885, 1.753333), (0.187214, 0.982299, 
            0.006234)), ((0.080113, 0.001127, 1.653333), (-0.987308, 0.143028, 
            0.069039)), ((0.081163, 0.005194, 1.653333), (-0.927883, 0.367154, 
            0.065055)), ((0.083946, 0.010748, 1.653333), (-0.817545, 0.572957, 
            0.057801)), ((0.08857, 0.016473, 1.653333), (-0.700158, 0.71222, 
            0.050206)), ((0.095028, 0.022182, 1.653333), (-0.591304, 0.805285, 
            0.043305)), ((0.103277, 0.027699, 1.653333), (-0.490483, 0.870662, 
            0.03707)), ((0.113291, 0.03286, 1.653333), (-0.398994, 0.91641, 0.031577)), 
            ((0.125034, 0.037542, 1.653333), (-0.318094, 0.947678, 0.026888)), ((
            0.138444, 0.041897, 1.653333), (-0.292195, 0.95602, 0.025454)), ((0.154577, 
            0.046568, 1.653333), (-0.256786, 0.96618, 0.023606)), ((0.173798, 0.049992, 
            1.653333), (0.0, 0.99994, 0.010933)), ((0.2425, 0.049992, 1.653333), (0.0, 
            0.99994, 0.010933)), ((0.265132, 0.049417, 1.653333), (0.055922, 0.998393, 
            0.009145)), ((0.293154, 0.047477, 1.653333), (0.10351, 0.994597, 
            0.007916)), ((0.316194, 0.044886, 1.653333), (0.127802, 0.991772, 
            0.007403)), ((0.339434, 0.041722, 1.653333), (0.148768, 0.988847, 
            0.007063)), ((0.362661, 0.038081, 1.653333), (0.167005, 0.985932, 
            0.006859)), ((0.397068, 0.031756, 1.653333), (0.19184, 0.981403, 
            0.006707)), ((0.440996, 0.023169, 1.653333), (0.19184, 0.981403, 
            0.006707)), ((0.457515, 0.019752, 1.653333), (0.205114, 0.978714, 
            0.006823)), ((0.476664, 0.015715, 1.653333), (0.206854, 0.978348, 
            0.006846)), ((0.494172, 0.012008, 1.653333), (0.207325, 0.978248, 
            0.006854)), ((0.509762, 0.008691, 1.653333), (0.208594, 0.977978, 
            0.006881)), ((0.523206, 0.005807, 1.653333), (0.210343, 0.977603, 
            0.006925)), ((0.534323, 0.003432, 1.653333), (0.208057, 0.978093, 
            0.006861)), ((0.54718, 0.000905, 1.653333), (0.187214, 0.982299, 
            0.006234)), ((0.546958, 0.000258, 1.653333), (-0.051907, -0.99865, 
            -0.001728)), ((0.533366, 0.000853, 1.653333), (-0.022329, -0.99975, 
            -0.000839)), ((0.521215, 0.001042, 1.653333), (-0.004709, -0.999989, 
            -0.000353)), ((0.506317, 0.001021, 1.653333), (0.011647, -0.999932, 
            4.8e-05)), ((0.48893, 0.00073, 1.653333), (0.025441, -0.999676, 0.000337)), 
            ((0.469339, 0.000148, 1.653333), (0.03742, -0.999299, 0.000538)), ((
            0.448492, -0.000762, 1.653333), (0.055871, -0.998438, 0.000763)), ((
            0.400984, -0.00342, 1.653333), (0.055871, -0.998438, 0.000763)), ((
            0.366901, -0.005529, 1.653333), (0.067641, -0.997709, 0.000722)), ((
            0.341062, -0.007319, 1.653333), (0.069854, -0.997557, 0.000701)), ((
            0.315113, -0.009133, 1.653333), (0.069674, -0.99757, 0.000704)), ((
            0.289371, -0.011579, 1.653333), (0.107257, -0.994231, -7e-05)), ((0.263558, 
            -0.01463, 1.653333), (0.122438, -0.992476, -0.000469)), ((0.2425, 
            -0.015703, 1.653333), (0.0, -0.999994, 0.003434)), ((0.173798, -0.015703, 
            1.653333), (0.0, -0.999994, 0.003434)), ((0.154251, -0.015678, 1.653333), (
            -0.001899, -0.999992, 0.003531)), ((0.136508, -0.015226, 1.653333), (
            -0.039628, -0.999199, 0.005608)), ((0.121233, -0.01428, 1.653333), (
            -0.074595, -0.997185, 0.007657)), ((0.108154, -0.012895, 1.653333), (
            -0.123513, -0.992286, 0.010673)), ((0.097407, -0.01107, 1.653333), (
            -0.194483, -0.980788, 0.015231)), ((0.089123, -0.008827, 1.653333), (
            -0.304348, -0.952295, 0.02251)), ((0.083523, -0.006191, 1.653333), (
            -0.509417, -0.85975, 0.036388)), ((0.080623, -0.001767, 1.653333), (
            -0.932433, -0.355412, 0.065202)), ((0.55018, 0.000264, 1.553333), (
            -0.051907, -0.99865, -0.001728)), ((0.536291, 0.000872, 1.553333), (
            -0.022329, -0.99975, -0.000839)), ((0.523874, 0.001064, 1.553333), (
            -0.004709, -0.999989, -0.000353)), ((0.50865, 0.001043, 1.553333), (
            0.011647, -0.999932, 4.8e-05)), ((0.490883, 0.000746, 1.553333), (0.025441, 
            -0.999676, 0.000337)), ((0.470864, 0.000151, 1.553333), (0.03742, 
            -0.999299, 0.000538)), ((0.449305, -0.000793, 1.553333), (0.055871, 
            -0.998438, 0.000763)), ((0.400754, -0.00351, 1.553333), (0.055871, 
            -0.998438, 0.000763)), ((0.366183, -0.00565, 1.553333), (0.067641, 
            -0.997709, 0.000722)), ((0.339779, -0.007479, 1.553333), (0.069854, 
            -0.997557, 0.000701)), ((0.313262, -0.009333, 1.553333), (0.069674, 
            -0.99757, 0.000704)), ((0.286957, -0.011832, 1.553333), (0.107257, 
            -0.994231, -7e-05)), ((0.26058, -0.01495, 1.553333), (0.122438, -0.992476, 
            -0.000469)), ((0.238695, -0.016046, 1.553333), (0.0, -0.999994, 0.003434)), 
            ((0.169038, -0.016046, 1.553333), (0.0, -0.999994, 0.003434)), ((0.148882, 
            -0.016021, 1.553333), (-0.001899, -0.999992, 0.003531)), ((0.130751, 
            -0.015559, 1.553333), (-0.039628, -0.999199, 0.005608)), ((0.115142, 
            -0.014592, 1.553333), (-0.074595, -0.997185, 0.007657)), ((0.101777, 
            -0.013177, 1.553333), (-0.123513, -0.992286, 0.010673)), ((0.090796, 
            -0.011312, 1.553333), (-0.194483, -0.980788, 0.015231)), ((0.08233, 
            -0.00902, 1.553333), (-0.304348, -0.952295, 0.02251)), ((0.076608, 
            -0.006326, 1.553333), (-0.509417, -0.85975, 0.036388)), ((0.073645, 
            -0.001805, 1.553333), (-0.932433, -0.355412, 0.065202)), ((0.073124, 
            0.001151, 1.553333), (-0.987308, 0.143028, 0.069039)), ((0.074197, 
            0.005307, 1.553333), (-0.927883, 0.367154, 0.065055)), ((0.077041, 
            0.010983, 1.553333), (-0.817545, 0.572957, 0.057801)), ((0.081765, 
            0.016833, 1.553333), (-0.700158, 0.71222, 0.050206)), ((0.088364, 0.022667, 
            1.553333), (-0.591304, 0.805285, 0.043305)), ((0.096794, 0.028305, 
            1.553333), (-0.490483, 0.870662, 0.03707)), ((0.107027, 0.033578, 
            1.553333), (-0.398994, 0.91641, 0.031577)), ((0.119027, 0.038363, 
            1.553333), (-0.318094, 0.947678, 0.026888)), ((0.13273, 0.042813, 
            1.553333), (-0.292195, 0.95602, 0.025454)), ((0.149215, 0.047587, 
            1.553333), (-0.256786, 0.96618, 0.023606)), ((0.169038, 0.051085, 
            1.553333), (0.0, 0.99994, 0.010933)), ((0.238695, 0.051085, 1.553333), (
            0.0, 0.99994, 0.010933)), ((0.262187, 0.050498, 1.553333), (0.055922, 
            0.998393, 0.009145)), ((0.290823, 0.048515, 1.553333), (0.10351, 0.994597, 
            0.007916)), ((0.314367, 0.045868, 1.553333), (0.127802, 0.991772, 
            0.007403)), ((0.338115, 0.042635, 1.553333), (0.148768, 0.988847, 
            0.007063)), ((0.36185, 0.038914, 1.553333), (0.167005, 0.985932, 
            0.006859)), ((0.396752, 0.032501, 1.553333), (0.19184, 0.981403, 
            0.006707)), ((0.441645, 0.023725, 1.553333), (0.19184, 0.981403, 
            0.006707)), ((0.458781, 0.020183, 1.553333), (0.205114, 0.978714, 
            0.006823)), ((0.478348, 0.016059, 1.553333), (0.206854, 0.978348, 
            0.006846)), ((0.49624, 0.01227, 1.553333), (0.207325, 0.978248, 0.006854)), 
            ((0.51217, 0.00888, 1.553333), (0.208594, 0.977978, 0.006881)), ((0.525908, 
            0.005934, 1.553333), (0.210343, 0.977603, 0.006925)), ((0.537269, 0.003507, 
            1.553333), (0.208057, 0.978093, 0.006861)), ((0.550407, 0.000925, 
            1.553333), (0.187214, 0.982299, 0.006234)), ((0.066135, 0.001176, 
            1.453333), (-0.987308, 0.143028, 0.069039)), ((0.067231, 0.005421, 
            1.453333), (-0.927883, 0.367154, 0.065055)), ((0.070135, 0.011218, 
            1.453333), (-0.817545, 0.572957, 0.057801)), ((0.074961, 0.017193, 
            1.453333), (-0.700158, 0.71222, 0.050206)), ((0.081701, 0.023152, 
            1.453333), (-0.591304, 0.805285, 0.043305)), ((0.090311, 0.02891, 
            1.453333), (-0.490483, 0.870662, 0.03707)), ((0.100762, 0.034297, 
            1.453333), (-0.398994, 0.91641, 0.031577)), ((0.11302, 0.039184, 1.453333), 
            (-0.318094, 0.947678, 0.026888)), ((0.127016, 0.043729, 1.453333), (
            -0.292195, 0.95602, 0.025454)), ((0.143854, 0.048605, 1.453333), (
            -0.256786, 0.96618, 0.023606)), ((0.164278, 0.052179, 1.453333), (0.0, 
            0.99994, 0.010933)), ((0.234891, 0.052179, 1.453333), (0.0, 0.99994, 
            0.010933)), ((0.259243, 0.051578, 1.453333), (0.055922, 0.998393, 
            0.009145)), ((0.288492, 0.049554, 1.453333), (0.10351, 0.994597, 
            0.007916)), ((0.312539, 0.046849, 1.453333), (0.127802, 0.991772, 
            0.007403)), ((0.336796, 0.043548, 1.453333), (0.148768, 0.988847, 
            0.007063)), ((0.36104, 0.039747, 1.453333), (0.167005, 0.985932, 
            0.006859)), ((0.396437, 0.033246, 1.453333), (0.19184, 0.981403, 
            0.006707)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.442293, 
            0.024282, 1.453333), (0.19184, 0.981403, 0.006707)), ((0.460047, 0.020615, 
            1.453333), (0.205114, 0.978714, 0.006823)), ((0.480033, 0.016402, 
            1.453333), (0.206854, 0.978348, 0.006846)), ((0.498307, 0.012533, 
            1.453333), (0.207325, 0.978248, 0.006854)), ((0.514578, 0.00907, 1.453333), 
            (0.208594, 0.977978, 0.006881)), ((0.528611, 0.006061, 1.453333), (
            0.210343, 0.977603, 0.006925)), ((0.540214, 0.003582, 1.453333), (0.208057, 
            0.978093, 0.006861)), ((0.553634, 0.000944, 1.453333), (0.187214, 0.982299, 
            0.006234)), ((0.553402, 0.00027, 1.453333), (-0.051907, -0.99865, 
            -0.001728)), ((0.539215, 0.000891, 1.453333), (-0.022329, -0.99975, 
            -0.000839)), ((0.526533, 0.001087, 1.453333), (-0.004709, -0.999989, 
            -0.000353)), ((0.510983, 0.001065, 1.453333), (0.011647, -0.999932, 
            4.8e-05)), ((0.492836, 0.000762, 1.453333), (0.025441, -0.999676, 
            0.000337)), ((0.472389, 0.000154, 1.453333), (0.03742, -0.999299, 
            0.000538)), ((0.450118, -0.000824, 1.453333), (0.055871, -0.998438, 
            0.000763)), ((0.400524, -0.003599, 1.453333), (0.055871, -0.998438, 
            0.000763)), ((0.365465, -0.005771, 1.453333), (0.067641, -0.997709, 
            0.000722)), ((0.338496, -0.007639, 1.453333), (0.069854, -0.997557, 
            0.000701)), ((0.311411, -0.009533, 1.453333), (0.069674, -0.99757, 
            0.000704)), ((0.284543, -0.012086, 1.453333), (0.107257, -0.994231, 
            -7e-05)), ((0.257601, -0.01527, 1.453333), (0.122438, -0.992476, 
            -0.000469)), ((0.234891, -0.01639, 1.453333), (0.0, -0.999994, 0.003434)), 
            ((0.164278, -0.01639, 1.453333), (0.0, -0.999994, 0.003434)), ((0.143513, 
            -0.016363, 1.453333), (-0.001899, -0.999992, 0.003531)), ((0.124995, 
            -0.015892, 1.453333), (-0.039628, -0.999199, 0.005608)), ((0.109052, 
            -0.014905, 1.453333), (-0.074595, -0.997185, 0.007657)), ((0.095401, 
            -0.013459, 1.453333), (-0.123513, -0.992286, 0.010673)), ((0.084184, 
            -0.011554, 1.453333), (-0.194483, -0.980788, 0.015231)), ((0.075538, 
            -0.009213, 1.453333), (-0.304348, -0.952295, 0.02251)), ((0.069693, 
            -0.006462, 1.453333), (-0.509417, -0.85975, 0.036388)), ((0.066667, 
            -0.001844, 1.453333), (-0.932433, -0.355412, 0.065202)), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='Skin web-after boom', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.307997, 
            -0.03714, 2.34509), ), ((0.282886, -0.017976, 2.257144), ), ((0.254459, 
            0.008379, 2.158489), ), ((0.233911, 0.020538, 2.082645), ), ((0.226827, 
            0.025695, 2.033333), ), ((0.223529, 0.026091, 1.953333), ), ((0.219406, 
            0.026705, 1.853333), ), ((0.215283, 0.02732, 1.753333), ), ((0.21116, 
            0.027934, 1.653333), ), ((0.207037, 0.028549, 1.553333), ), ((0.202914, 
            0.029163, 1.453333), ), ((0.198791, 0.029777, 1.353333), ), ((0.194668, 
            0.030392, 1.253333), ), ((0.190545, 0.031006, 1.153333), ), ((0.186422, 
            0.031621, 1.053333), ), ((0.182299, 0.032235, 0.953333), ), ((0.178176, 
            0.03285, 0.853333), ), ((0.174053, 0.033464, 0.753333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar cap bottom mid', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.207124, 
            -0.016733, 1.353333), ), ((0.182124, -0.016733, 1.353333), ), ((0.203001, 
            -0.017076, 1.253333), ), ((0.178001, -0.017076, 1.253333), ), ((0.198878, 
            -0.01742, 1.153333), ), ((0.173878, -0.01742, 1.153333), ), ((0.194755, 
            -0.017763, 1.053333), ), ((0.169755, -0.017763, 1.053333), ), ((0.190632, 
            -0.018107, 0.953333), ), ((0.165632, -0.018107, 0.953333), ), ((0.186509, 
            -0.01845, 0.853333), ), ((0.161509, -0.01845, 0.853333), ), ((0.182386, 
            -0.018794, 0.753333), ), ((0.157386, -0.018794, 0.753333), ), ((0.154637, 
            -0.019023, 0.686667), ), ((0.178263, -0.019137, 0.653333), ), ((0.187971, 
            -0.019023, 0.686667), ), ((0.161596, -0.019137, 0.653333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar cap bottom root', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.175239, 
            -0.019389, 0.58), ), ((0.150239, -0.019389, 0.58), ), ((0.17304, -0.019572, 
            0.526667), ), ((0.14804, -0.019572, 0.526667), ), ((0.180687, -0.019629, 
            0.476667), ), ((0.155687, -0.019629, 0.476667), ), ((0.180687, -0.019629, 
            0.42), ), ((0.155687, -0.019629, 0.42), ), ((0.180687, -0.019629, 
            0.333333), ), ((0.155687, -0.019629, 0.333333), ), ((0.180687, -0.019629, 
            0.233333), ), ((0.155687, -0.019629, 0.233333), ), ((0.180687, -0.019629, 
            0.133333), ), ((0.155687, -0.019629, 0.133333), ), ((0.180687, -0.019629, 
            0.033333), ), ((0.155687, -0.019629, 0.033333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar cap bottom tip', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.31633, 
            -0.064831, 2.339707), ), ((0.29133, -0.064831, 2.339707), ), ((0.291219, 
            -0.048904, 2.251132), ), ((0.266219, -0.048904, 2.251132), ), ((0.262792, 
            -0.030872, 2.150859), ), ((0.237792, -0.030872, 2.150859), ), ((0.241945, 
            -0.017825, 2.079145), ), ((0.216945, -0.017825, 2.079145), ), ((0.23516, 
            -0.014398, 2.033333), ), ((0.21016, -0.014398, 2.033333), ), ((0.231862, 
            -0.014672, 1.953333), ), ((0.206862, -0.014672, 1.953333), ), ((0.227739, 
            -0.015016, 1.853333), ), ((0.202739, -0.015016, 1.853333), ), ((0.223616, 
            -0.015359, 1.753333), ), ((0.198616, -0.015359, 1.753333), ), ((0.219493, 
            -0.015703, 1.653333), ), ((0.194493, -0.015703, 1.653333), ), ((0.21537, 
            -0.016046, 1.553333), ), ((0.19037, -0.016046, 1.553333), ), ((0.211247, 
            -0.01639, 1.453333), ), ((0.186247, -0.01639, 1.453333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar cap top mid', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.182124, 
            0.053272, 1.353333), ), ((0.207124, 0.053272, 1.353333), ), ((0.178001, 
            0.054365, 1.253333), ), ((0.203001, 0.054365, 1.253333), ), ((0.173878, 
            0.055459, 1.153333), ), ((0.198878, 0.055459, 1.153333), ), ((0.169755, 
            0.056552, 1.053333), ), ((0.194755, 0.056552, 1.053333), ), ((0.165632, 
            0.057646, 0.953333), ), ((0.190632, 0.057646, 0.953333), ), ((0.161509, 
            0.058739, 0.853333), ), ((0.186509, 0.058739, 0.853333), ), ((0.157386, 
            0.059832, 0.753333), ), ((0.182386, 0.059832, 0.753333), ), ((0.161596, 
            0.060926, 0.653333), ), ((0.187971, 0.060561, 0.686667), ), ((0.154637, 
            0.060561, 0.686667), ), ((0.178263, 0.060926, 0.653333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar cap top root', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.150239, 
            0.061727, 0.58), ), ((0.175239, 0.061727, 0.58), ), ((0.14804, 0.062311, 
            0.526667), ), ((0.17304, 0.062311, 0.526667), ), ((0.147353, 0.062493, 
            0.476667), ), ((0.172353, 0.062493, 0.476667), ), ((0.147353, 0.062493, 
            0.42), ), ((0.172353, 0.062493, 0.42), ), ((0.147353, 0.062493, 0.333333), 
            ), ((0.172353, 0.062493, 0.333333), ), ((0.147353, 0.062493, 0.233333), ), 
            ((0.172353, 0.062493, 0.233333), ), ((0.147353, 0.062493, 0.133333), ), ((
            0.172353, 0.062493, 0.133333), ), ((0.147353, 0.062493, 0.033333), ), ((
            0.172353, 0.062493, 0.033333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar cap top tip', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.29133, 
            -0.022657, 2.347905), ), ((0.31633, -0.022657, 2.347905), ), ((0.266219, 
            -0.001531, 2.260341), ), ((0.291219, -0.001531, 2.260341), ), ((0.237792, 
            0.022386, 2.161211), ), ((0.262792, 0.022386, 2.161211), ), ((0.216945, 
            0.040483, 2.082723), ), ((0.241945, 0.040483, 2.082723), ), ((0.21016, 
            0.045837, 2.033333), ), ((0.23516, 0.045837, 2.033333), ), ((0.206862, 
            0.046712, 1.953333), ), ((0.231862, 0.046712, 1.953333), ), ((0.202739, 
            0.047805, 1.853333), ), ((0.227739, 0.047805, 1.853333), ), ((0.198616, 
            0.048899, 1.753333), ), ((0.223616, 0.048899, 1.753333), ), ((0.194493, 
            0.049992, 1.653333), ), ((0.219493, 0.049992, 1.653333), ), ((0.19037, 
            0.051085, 1.553333), ), ((0.21537, 0.051085, 1.553333), ), ((0.186247, 
            0.052179, 1.453333), ), ((0.211247, 0.052179, 1.453333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='spar web before boom', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.186596, 
            0.034078, 0.653333), ), ((0.153263, 0.034078, 0.653333), ), ((0.191906, 
            0.034593, 0.58), ), ((0.189707, 0.034937, 0.526667), ), ((0.18902, 
            0.007745, 0.476667), ), ((0.18902, 0.007745, 0.42), ), ((0.18902, 0.007745, 
            0.333333), ), ((0.18902, 0.007745, 0.233333), ), ((0.18902, 0.007745, 
            0.133333), ), ((0.18902, 0.007745, 0.033333), ), ((0.141906, 0.034593, 
            0.58), ), ((0.139707, 0.034937, 0.526667), ), ((0.13902, 0.007745, 
            0.476667), ), ((0.13902, 0.007745, 0.42), ), ((0.13902, 0.007745, 
            0.333333), ), ((0.13902, 0.007745, 0.233333), ), ((0.13902, 0.007745, 
            0.133333), ), ((0.13902, 0.007745, 0.033333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='top flange', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.29133, -0.022657, 
            2.347905), ), ((0.31633, -0.022657, 2.347905), ), ((0.266219, -0.001531, 
            2.260341), ), ((0.291219, -0.001531, 2.260341), ), ((0.237792, 0.022386, 
            2.161211), ), ((0.262792, 0.022386, 2.161211), ), ((0.216945, 0.040483, 
            2.082723), ), ((0.241945, 0.040483, 2.082723), ), ((0.21016, 0.045837, 
            2.033333), ), ((0.23516, 0.045837, 2.033333), ), ((0.206862, 0.046712, 
            1.953333), ), ((0.231862, 0.046712, 1.953333), ), ((0.202739, 0.047805, 
            1.853333), ), ((0.227739, 0.047805, 1.853333), ), ((0.198616, 0.048899, 
            1.753333), ), ((0.223616, 0.048899, 1.753333), ), ((0.194493, 0.049992, 
            1.653333), ), ((0.219493, 0.049992, 1.653333), ), ((0.19037, 0.051085, 
            1.553333), ), ((0.21537, 0.051085, 1.553333), ), ((0.186247, 0.052179, 
            1.453333), ), ((0.211247, 0.052179, 1.453333), ), ((0.182124, 0.053272, 
            1.353333), ), ((0.207124, 0.053272, 1.353333), ), ((0.178001, 0.054365, 
            1.253333), ), ((0.203001, 0.054365, 1.253333), ), ((0.173878, 0.055459, 
            1.153333), ), ((0.198878, 0.055459, 1.153333), ), ((0.169755, 0.056552, 
            1.053333), ), ((0.194755, 0.056552, 1.053333), ), ((0.165632, 0.057646, 
            0.953333), ), ((0.190632, 0.057646, 0.953333), ), ((0.161509, 0.058739, 
            0.853333), ), ((0.186509, 0.058739, 0.853333), ), ((0.157386, 0.059832, 
            0.753333), ), ((0.182386, 0.059832, 0.753333), ), ((0.161596, 0.060926, 
            0.653333), ), ((0.187971, 0.060561, 0.686667), ), ((0.154637, 0.060561, 
            0.686667), ), ((0.178263, 0.060926, 0.653333), ), ((0.150239, 0.061727, 
            0.58), ), ((0.175239, 0.061727, 0.58), ), ((0.14804, 0.062311, 0.526667), 
            ), ((0.17304, 0.062311, 0.526667), ), ((0.147353, 0.062493, 0.476667), ), (
            (0.172353, 0.062493, 0.476667), ), ((0.147353, 0.062493, 0.42), ), ((
            0.172353, 0.062493, 0.42), ), ((0.147353, 0.062493, 0.333333), ), ((
            0.172353, 0.062493, 0.333333), ), ((0.147353, 0.062493, 0.233333), ), ((
            0.172353, 0.062493, 0.233333), ), ((0.147353, 0.062493, 0.133333), ), ((
            0.172353, 0.062493, 0.133333), ), ((0.147353, 0.062493, 0.033333), ), ((
            0.172353, 0.062493, 0.033333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='bottom flange', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.31633, -0.064831, 
            2.339707), ), ((0.29133, -0.064831, 2.339707), ), ((0.291219, -0.048904, 
            2.251132), ), ((0.266219, -0.048904, 2.251132), ), ((0.262792, -0.030872, 
            2.150859), ), ((0.237792, -0.030872, 2.150859), ), ((0.241945, -0.017825, 
            2.079145), ), ((0.216945, -0.017825, 2.079145), ), ((0.23516, -0.014398, 
            2.033333), ), ((0.21016, -0.014398, 2.033333), ), ((0.231862, -0.014672, 
            1.953333), ), ((0.206862, -0.014672, 1.953333), ), ((0.227739, -0.015016, 
            1.853333), ), ((0.202739, -0.015016, 1.853333), ), ((0.223616, -0.015359, 
            1.753333), ), ((0.198616, -0.015359, 1.753333), ), ((0.219493, -0.015703, 
            1.653333), ), ((0.194493, -0.015703, 1.653333), ), ((0.21537, -0.016046, 
            1.553333), ), ((0.19037, -0.016046, 1.553333), ), ((0.211247, -0.01639, 
            1.453333), ), ((0.186247, -0.01639, 1.453333), ), ((0.207124, -0.016733, 
            1.353333), ), ((0.182124, -0.016733, 1.353333), ), ((0.203001, -0.017076, 
            1.253333), ), ((0.178001, -0.017076, 1.253333), ), ((0.198878, -0.01742, 
            1.153333), ), ((0.173878, -0.01742, 1.153333), ), ((0.194755, -0.017763, 
            1.053333), ), ((0.169755, -0.017763, 1.053333), ), ((0.190632, -0.018107, 
            0.953333), ), ((0.165632, -0.018107, 0.953333), ), ((0.186509, -0.01845, 
            0.853333), ), ((0.161509, -0.01845, 0.853333), ), ((0.182386, -0.018794, 
            0.753333), ), ((0.157386, -0.018794, 0.753333), ), ((0.154637, -0.019023, 
            0.686667), ), ((0.178263, -0.019137, 0.653333), ), ((0.187971, -0.019023, 
            0.686667), ), ((0.161596, -0.019137, 0.653333), ), ((0.175239, -0.019389, 
            0.58), ), ((0.150239, -0.019389, 0.58), ), ((0.17304, -0.019572, 0.526667), 
            ), ((0.14804, -0.019572, 0.526667), ), ((0.180687, -0.019629, 0.476667), ), 
            ((0.155687, -0.019629, 0.476667), ), ((0.180687, -0.019629, 0.42), ), ((
            0.155687, -0.019629, 0.42), ), ((0.180687, -0.019629, 0.333333), ), ((
            0.155687, -0.019629, 0.333333), ), ((0.180687, -0.019629, 0.233333), ), ((
            0.155687, -0.019629, 0.233333), ), ((0.180687, -0.019629, 0.133333), ), ((
            0.155687, -0.019629, 0.133333), ), ((0.180687, -0.019629, 0.033333), ), ((
            0.155687, -0.019629, 0.033333), ), ))        

        mdb.models['Model-1'].parts['Wing'].Surface(name='base rib', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.208773, 0.035119, 
            0.0), ), ((0.10131, 0.032971, 0.0), ), ((0.229442, 0.009076, 0.0), ), ((
            0.119267, 0.007745, 0.0), ), ((0.147353, 0.007745, 0.0), ), ((0.172353, 
            0.007745, 0.0), ), ((0.44176, 0.019318, 0.0), ), ))

        #singularity sets        
        mdb.models['Model-1'].parts['Wing'].DatumPlaneByPrincipalPlane(offset=0.025, 
            principalPlane=XYPLANE)
        mdb.models['Model-1'].parts['Wing'].PartitionFaceByDatumPlane(datumPlane=
            mdb.models['Model-1'].parts['Wing'].datums[136], faces=
            mdb.models['Model-1'].parts['Wing'].faces)
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.006651, 0.015798, 
            0.016667), ), ((0.013196, 0.022973, 0.016667), ), ((0.022015, 0.030055, 
            0.016667), ), ((0.03306, 0.036822, 0.016667), ), ((0.046301, 0.04309, 
            0.016667), ), ((0.061681, 0.048716, 0.016667), ), ((0.079113, 0.054205, 
            0.016667), ), ((0.10131, 0.060345, 0.016667), ), ((0.129143, 0.062493, 
            0.016667), ), ((0.208773, 0.062493, 0.016667), ), ((0.244141, 0.061065, 
            0.016667), ), ((0.27601, 0.058361, 0.016667), ), ((0.304924, 0.054872, 
            0.016667), ), ((0.333999, 0.050707, 0.016667), ), ((0.362972, 0.045983, 
            0.016667), ), ((0.414203, 0.03622, 0.016667), ), ((0.472039, 0.024678, 
            0.016667), ), ((0.495972, 0.019633, 0.016667), ), ((0.517853, 0.015, 
            0.016667), ), ((0.537335, 0.010854, 0.016667), ), ((0.554134, 0.007251, 
            0.016667), ), ((0.568024, 0.004283, 0.016667), ), ((0.584108, 0.001123, 
            0.016667), ), ((0.577665, 0.000641, 0.016667), ), ((0.562164, 0.00117, 
            0.016667), ), ((0.545803, 0.001329, 0.016667), ), ((0.526111, 0.001196, 
            0.016667), ), ((0.503423, 0.000714, 0.016667), ), ((0.478105, -0.00014, 
            0.016667), ), ((0.446462, -0.001749, 0.016667), ), ((0.398318, -0.004443, 
            0.016667), ), ((0.358628, -0.006917, 0.016667), ), ((0.326327, -0.009155, 
            0.016667), ), ((0.293889, -0.011422, 0.016667), ), ((0.26171, -0.014483, 
            0.016667), ), ((0.229442, -0.018298, 0.016667), ), ((0.198897, -0.019629, 
            0.016667), ), ((0.119267, -0.019629, 0.016667), ), ((0.09282, -0.019598, 
            0.016667), ), ((0.070651, -0.019031, 0.016667), ), ((0.051562, -0.017847, 
            0.016667), ), ((0.035219, -0.016115, 0.016667), ), ((0.021793, -0.013832, 
            0.016667), ), ((0.011444, -0.011027, 0.016667), ), ((0.004453, -0.007731, 
            0.016667), ), ((0.000836, -0.002193, 0.016667), ), ((0.000405, 0.002797, 
            0.016667), ), ((0.002413, 0.008758, 0.016667), ), ((0.44852, 0.029512, 
            0.016667), ), ), name='Singularity-skin')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.00151, 0.006476, 
            0.05), (-0.929852, 0.367933, 0.0)), ((0.000203, 0.001398, 0.05), (
            -0.989669, 0.14337, 0.0)), ((0.001672, -0.004386, 0.05), (-0.934421, 
            -0.35617, 0.0)), ((0.006398, -0.008883, 0.05), (-0.509755, -0.86032, 0.0)), 
            ((0.014545, -0.012018, 0.05), (-0.304425, -0.952536, 0.0)), ((0.025938, 
            -0.014654, 0.05), (-0.194505, -0.980901, 0.0)), ((0.040354, -0.016754, 
            0.05), (-0.12352, -0.992342, 0.0)), ((0.057635, -0.018302, 0.05), (
            -0.074597, -0.997214, 0.0)), ((0.077593, -0.019307, 0.05), (-0.039629, 
            -0.999214, 0.0)), ((0.101105, -0.019614, 0.05), (-0.001899, -0.999998, 
            0.0)), ((0.129143, -0.019629, 0.05), (0.0, -1.0, 0.0)), ((0.208773, 
            -0.019629, 0.05), (0.0, -1.0, 0.0)), ((0.240235, -0.016966, 0.05), (
            0.122438, -0.992476, 0.0)), ((0.272393, -0.01333, 0.05), (0.107257, 
            -0.994231, 0.0)), ((0.304702, -0.010667, 0.05), (0.069674, -0.99757, 0.0)), 
            ((0.33714, -0.008398, 0.05), (0.069854, -0.997557, 0.0)), ((0.369303, 
            -0.006193, 0.05), (0.067641, -0.99771, 0.0)), ((0.44176, 0.030833, 0.05), (
            0.191844, 0.981425, 0.0)), ((0.416659, -0.003417, 0.05), (0.055871, 
            -0.998438, 0.0)), ((0.457924, -0.001108, 0.05), (0.055871, -0.998438, 
            0.0)), ((0.486823, 0.000187, 0.05), (0.03742, -0.9993, 0.0)), ((0.511306, 
            0.000915, 0.05), (0.025441, -0.999676, 0.0)), ((0.533034, 0.001276, 0.05), 
            (0.011647, -0.999932, 0.0)), ((0.55165, 0.001302, 0.05), (-0.004709, 
            -0.999989, 0.0)), ((0.566831, 0.001066, 0.05), (-0.022329, -0.999751, 
            0.0)), ((0.005096, 0.001386, 0.58), (-0.987308, 0.143028, 0.069039)), ((
            0.00639, 0.006405, 0.58), (-0.927883, 0.367154, 0.065055)), ((0.009823, 
            0.013263, 0.58), (-0.817545, 0.572957, 0.057801)), ((0.015529, 0.020332, 
            0.58), (-0.700158, 0.71222, 0.050206)), ((0.0235, 0.027381, 0.58), (
            -0.591304, 0.805285, 0.043305)), ((0.033684, 0.034193, 0.58), (-0.490483, 
            0.870662, 0.03707)), ((0.046045, 0.040566, 0.58), (-0.398994, 0.91641, 
            0.031577)), ((0.060543, 0.046349, 0.58), (-0.318094, 0.947678, 0.026888)), 
            ((0.077098, 0.051726, 0.58), (-0.292195, 0.95602, 0.025454)), ((0.097011, 
            0.057492, 0.58), (-0.256786, 0.96618, 0.023606)), ((0.122663, 0.061727, 
            0.58), (0.0, 0.99994, 0.010933)), ((0.201623, 0.061727, 0.58), (0.0, 
            0.99994, 0.010933)), ((0.233501, 0.06102, 0.58), (0.055922, 0.998393, 
            0.009145)), ((0.268113, 0.058626, 0.58), (0.10351, 0.994597, 0.007916)), ((
            0.296561, 0.055427, 0.58), (0.127802, 0.991772, 0.007403)), ((0.325256, 
            0.051522, 0.58), (0.148768, 0.988847, 0.007063)), ((0.353937, 0.047026, 
            0.58), (0.167005, 0.985932, 0.006859)), ((0.393658, 0.039757, 0.58), (
            0.19184, 0.981403, 0.006707)), ((0.448001, 0.029135, 0.58), (0.19184, 
            0.981403, 0.006707)), ((0.471124, 0.024382, 0.58), (0.205114, 0.978714, 
            0.006823)), ((0.494765, 0.019399, 0.58), (0.206854, 0.978348, 0.006846)), (
            (0.516381, 0.014821, 0.58), (0.207325, 0.978248, 0.006854)), ((0.535627, 
            0.010726, 0.58), (0.208594, 0.977978, 0.006881)), ((0.552224, 0.007166, 
            0.58), (0.210343, 0.977603, 0.006925)), ((0.565948, 0.004234, 0.58), (
            0.208057, 0.978093, 0.006861)), ((0.581828, 0.001113, 0.58), (0.187214, 
            0.982299, 0.006234)), ((0.581556, 0.000318, 0.58), (-0.051907, -0.99865, 
            -0.001728)), ((0.564767, 0.001053, 0.58), (-0.022329, -0.99975, 
            -0.000839)), ((0.549768, 0.001286, 0.58), (-0.004709, -0.999989, 
            -0.000353)), ((0.531376, 0.00126, 0.58), (0.011647, -0.999932, 4.8e-05)), (
            (0.509911, 0.000903, 0.58), (0.025441, -0.999676, 0.000337)), ((0.485725, 
            0.000183, 0.58), (0.03742, -0.999299, 0.000538)), ((0.457274, -0.00109, 
            0.58), (0.055871, -0.998438, 0.000763)), ((0.398502, -0.004379, 0.58), (
            0.055871, -0.998438, 0.000763)), ((0.359168, -0.00683, 0.58), (0.067641, 
            -0.997709, 0.000722)), ((0.327263, -0.00904, 0.58), (0.069854, -0.997557, 
            0.000701)), ((0.295222, -0.01128, 0.58), (0.069674, -0.99757, 0.000704)), (
            (0.263437, -0.014301, 0.58), (0.107257, -0.994231, -7e-05)), ((0.231565, 
            -0.018069, 0.58), (0.122438, -0.992476, -0.000469)), ((0.201623, -0.019389, 
            0.58), (0.0, -0.999994, 0.003434)), ((0.122663, -0.019389, 0.58), (0.0, 
            -0.999994, 0.003434)), ((0.096607, -0.019358, 0.58), (-0.001899, -0.999992, 
            0.003531)), ((0.074704, -0.018799, 0.58), (-0.039628, -0.999199, 
            0.005608)), ((0.055847, -0.017631, 0.58), (-0.074595, -0.997185, 
            0.007657)), ((0.039701, -0.01592, 0.58), (-0.123513, -0.992286, 0.010673)), 
            ((0.026435, -0.013665, 0.58), (-0.194483, -0.980788, 0.015231)), ((0.01621, 
            -0.010896, 0.58), (-0.304348, -0.952295, 0.02251)), ((0.0093, -0.00764, 
            0.58), (-0.509417, -0.85975, 0.036388)), ((0.005723, -0.002174, 0.58), (
            -0.932433, -0.355412, 0.065202)), ((0.583278, 0.000321, 0.526667), (
            -0.051907, -0.99865, -0.001728)), ((0.56633, 0.001063, 0.526667), (
            -0.022329, -0.99975, -0.000839)), ((0.55119, 0.001298, 0.526667), (
            -0.004709, -0.999989, -0.000353)), ((0.532624, 0.001272, 0.526667), (
            0.011647, -0.999932, 4.8e-05)), ((0.510957, 0.000911, 0.526667), (0.025441, 
            -0.999676, 0.000337)), ((0.486543, 0.000185, 0.526667), (0.03742, 
            -0.999299, 0.000538)), ((0.457721, -0.001106, 0.526667), (0.055871, 
            -0.998438, 0.000763)), ((0.398376, -0.004427, 0.526667), (0.055871, 
            -0.998438, 0.000763)), ((0.358778, -0.006894, 0.526667), (0.067641, 
            -0.997709, 0.000722)), ((0.326572, -0.009126, 0.526667), (0.069854, 
            -0.997557, 0.000701)), ((0.294229, -0.011387, 0.526667), (0.069674, 
            -0.99757, 0.000704)), ((0.262144, -0.014437, 0.526667), (0.107257, 
            -0.994231, -7e-05)), ((0.22997, -0.018241, 0.526667), (0.122438, -0.992476, 
            -0.000469)), ((0.199584, -0.019572, 0.526667), (0.0, -0.999994, 0.003434)), 
            ((0.120113, -0.019572, 0.526667), (0.0, -0.999994, 0.003434)), ((0.093739, 
            -0.019541, 0.526667), (-0.001899, -0.999992, 0.003531)), ((0.07163, 
            -0.018977, 0.526667), (-0.039628, -0.999199, 0.005608)), ((0.052595, 
            -0.017797, 0.526667), (-0.074595, -0.997185, 0.007657)), ((0.036297, 
            -0.01607, 0.526667), (-0.123513, -0.992286, 0.010673)), ((0.022907, 
            -0.013794, 0.526667), (-0.194483, -0.980788, 0.015231)), ((0.012585, 
            -0.010998, 0.526667), (-0.304348, -0.952295, 0.02251)), ((0.005611, 
            -0.007712, 0.526667), (-0.509417, -0.85975, 0.036388)), ((0.002001, 
            -0.002193, 0.526667), (-0.932433, -0.355412, 0.065202)), ((0.001368, 
            0.001398, 0.526667), (-0.987308, 0.143028, 0.069039)), ((0.002674, 
            0.006464, 0.526667), (-0.927883, 0.367154, 0.065055)), ((0.006139, 
            0.013386, 0.526667), (-0.817545, 0.572957, 0.057801)), ((0.011899, 
            0.020522, 0.526667), (-0.700158, 0.71222, 0.050206)), ((0.019945, 0.027638, 
            0.526667), (-0.591304, 0.805285, 0.043305)), ((0.030224, 0.034515, 
            0.526667), (-0.490483, 0.870662, 0.03707)), ((0.042702, 0.040948, 
            0.526667), (-0.398994, 0.91641, 0.031577)), ((0.057336, 0.046785, 
            0.526667), (-0.318094, 0.947678, 0.026888)), ((0.074047, 0.052213, 
            0.526667), (-0.292195, 0.95602, 0.025454)), ((0.094146, 0.058034, 
            0.526667), (-0.256786, 0.96618, 0.023606)), ((0.120113, 0.062311, 
            0.526667), (0.0, 0.99994, 0.010933)), ((0.199584, 0.062311, 0.526667), (
            0.0, 0.99994, 0.010933)), ((0.231923, 0.061597, 0.526667), (0.055922, 
            0.998393, 0.009145)), ((0.266864, 0.05918, 0.526667), (0.10351, 0.994597, 
            0.007916)), ((0.295581, 0.055952, 0.526667), (0.127802, 0.991772, 
            0.007403)), ((0.324547, 0.052009, 0.526667), (0.148768, 0.988847, 
            0.007063)), ((0.353498, 0.047471, 0.526667), (0.167005, 0.985932, 
            0.006859)), ((0.393485, 0.040156, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.448357, 0.029429, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.471804, 0.024611, 0.526667), (0.205114, 0.978714, 
            0.006823)), ((0.495668, 0.019581, 0.526667), (0.206854, 0.978348, 
            0.006846)), ((0.517488, 0.01496, 0.526667), (0.207325, 0.978248, 
            0.006854)), ((0.536915, 0.010826, 0.526667), (0.208594, 0.977978, 
            0.006881)), ((0.553668, 0.007233, 0.526667), (0.210343, 0.977603, 
            0.006925)), ((0.567521, 0.004274, 0.526667), (0.208057, 0.978093, 
            0.006861)), ((0.583553, 0.001123, 0.526667), (0.187214, 0.982299, 
            0.006234)), ((0.000203, 0.001398, 0.476667), (-0.989669, 0.14337, 0.0)), ((
            0.00151, 0.006476, 0.476667), (-0.929852, 0.367933, 0.0)), ((0.004984, 
            0.013419, 0.476667), (-0.818914, 0.573916, 0.0)), ((0.010758, 0.020575, 
            0.476667), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.476667), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.476667), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.476667), (-0.399193, 0.916867, 
            0.0)), ((0.056319, 0.046917, 0.476667), (-0.318209, 0.94802, 0.0)), ((
            0.073077, 0.052361, 0.476667), (-0.29229, 0.95633, 0.0)), ((0.093229, 
            0.058198, 0.476667), (-0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 
            0.476667), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.476667), (0.0, 1.0, 
            0.0)), ((0.231395, 0.061779, 0.476667), (0.055924, 0.998435, 0.0)), ((
            0.266448, 0.059356, 0.476667), (0.103513, 0.994628, 0.0)), ((0.295248, 
            0.056119, 0.476667), (0.127805, 0.991799, 0.0)), ((0.324299, 0.052166, 
            0.476667), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 0.476667), (
            0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.476667), (0.191844, 
            0.981425, 0.0)), ((0.44176, 0.030833, 0.476667), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.476667), (0.205119, 0.978737, 0.0)), ((
            0.488195, 0.021277, 0.476667), (0.206858, 0.978371, 0.0)), ((0.5108, 
            0.016494, 0.476667), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.476667), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.476667), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.476667), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.476667), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.476667), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.476667), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.476667), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.476667), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.476667), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.476667), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.476667), (0.055871, -0.998438, 
            0.0)), ((0.416659, -0.003417, 0.476667), (0.055871, -0.998438, 0.0)), ((
            0.369303, -0.006193, 0.476667), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.476667), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.476667), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.476667), (
            0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.476667), (0.122438, 
            -0.992476, 0.0)), ((0.208773, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((
            0.129143, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 
            0.476667), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.476667), 
            (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.476667), (-0.074597, 
            -0.997214, 0.0)), ((0.040354, -0.016754, 0.476667), (-0.12352, -0.992342, 
            0.0)), ((0.025938, -0.014654, 0.476667), (-0.194505, -0.980901, 0.0)), ((
            0.014545, -0.012018, 0.476667), (-0.304425, -0.952536, 0.0)), ((0.006398, 
            -0.008883, 0.476667), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 
            0.476667), (-0.934421, -0.35617, 0.0)), ((0.583833, 0.000321, 0.42), (
            -0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 0.42), (-0.022329, 
            -0.999751, 0.0)), ((0.55165, 0.001302, 0.42), (-0.004709, -0.999989, 0.0)), 
            ((0.533034, 0.001276, 0.42), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.42), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.42), 
            (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.42), (0.055871, 
            -0.998438, 0.0)), ((0.416659, -0.003417, 0.42), (0.055871, -0.998438, 
            0.0)), ((0.369303, -0.006193, 0.42), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.42), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.42), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.42), (0.107257, 
            -0.994231, 0.0)), ((0.240235, -0.016966, 0.42), (0.122438, -0.992476, 
            0.0)), ((0.208773, -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.129143, 
            -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.42), (
            -0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.42), (-0.039629, 
            -0.999214, 0.0)), ((0.057635, -0.018302, 0.42), (-0.074597, -0.997214, 
            0.0)), ((0.040354, -0.016754, 0.42), (-0.12352, -0.992342, 0.0)), ((
            0.025938, -0.014654, 0.42), (-0.194505, -0.980901, 0.0)), ((0.014545, 
            -0.012018, 0.42), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 
            0.42), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.42), (
            -0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 0.42), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.42), (-0.929852, 0.367933, 0.0)), ((
            0.004984, 0.013419, 0.42), (-0.818914, 0.573916, 0.0)), ((0.010758, 
            0.020575, 0.42), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.42), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.42), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.42), (-0.399193, 0.916867, 0.0)), 
            ((0.056319, 0.046917, 0.42), (-0.318209, 0.94802, 0.0)), ((0.073077, 
            0.052361, 0.42), (-0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.42), (
            -0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 0.42), (0.0, 1.0, 0.0)), 
            ((0.198897, 0.062493, 0.42), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.42), 
            (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.42), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.42), (0.127805, 0.991799, 0.0)), (
            (0.324299, 0.052166, 0.42), (0.148772, 0.988872, 0.0)), ((0.353335, 
            0.047615, 0.42), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.42), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.42), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.42), (0.205119, 0.978737, 0.0)), ((0.488195, 
            0.021277, 0.42), (0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.42), (
            0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 0.42), (0.208599, 0.978001, 
            0.0)), ((0.548842, 0.00839, 0.42), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.42), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.42), (
            0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.333333), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.333333), (-0.929852, 0.367933, 
            0.0)), ((0.004984, 0.013419, 0.333333), (-0.818914, 0.573916, 0.0)), ((
            0.010758, 0.020575, 0.333333), (-0.701042, 0.71312, 0.0)), ((0.018825, 
            0.027712, 0.333333), (-0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 
            0.333333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 0.333333), (
            -0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.333333), (-0.318209, 
            0.94802, 0.0)), ((0.073077, 0.052361, 0.333333), (-0.29229, 0.95633, 0.0)), 
            ((0.093229, 0.058198, 0.333333), (-0.256858, 0.966449, 0.0)), ((0.119267, 
            0.062493, 0.333333), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.333333), (
            0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.333333), (0.055924, 0.998435, 
            0.0)), ((0.266448, 0.059356, 0.333333), (0.103513, 0.994628, 0.0)), ((
            0.295248, 0.056119, 0.333333), (0.127805, 0.991799, 0.0)), ((0.324299, 
            0.052166, 0.333333), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 
            0.333333), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.333333), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.333333), (0.191844, 
            0.981425, 0.0)), ((0.463659, 0.026434, 0.333333), (0.205119, 0.978737, 
            0.0)), ((0.488195, 0.021277, 0.333333), (0.206858, 0.978371, 0.0)), ((
            0.5108, 0.016494, 0.333333), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.333333), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.333333), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.333333), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.333333), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.333333), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.333333), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.333333), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.333333), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.333333), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.333333), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.333333), (0.055871, -0.998438, 
            0.0)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.416659, 
            -0.003417, 0.333333), (0.055871, -0.998438, 0.0)), ((0.369303, -0.006193, 
            0.333333), (0.067641, -0.99771, 0.0)), ((0.33714, -0.008398, 0.333333), (
            0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 0.333333), (0.069674, 
            -0.99757, 0.0)), ((0.272393, -0.01333, 0.333333), (0.107257, -0.994231, 
            0.0)), ((0.240235, -0.016966, 0.333333), (0.122438, -0.992476, 0.0)), ((
            0.208773, -0.019629, 0.333333), (0.0, -1.0, 0.0)), ((0.129143, -0.019629, 
            0.333333), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.333333), (-0.001899, 
            -0.999998, 0.0)), ((0.077593, -0.019307, 0.333333), (-0.039629, -0.999214, 
            0.0)), ((0.057635, -0.018302, 0.333333), (-0.074597, -0.997214, 0.0)), ((
            0.040354, -0.016754, 0.333333), (-0.12352, -0.992342, 0.0)), ((0.025938, 
            -0.014654, 0.333333), (-0.194505, -0.980901, 0.0)), ((0.014545, -0.012018, 
            0.333333), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 0.333333), 
            (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.333333), (-0.934421, 
            -0.35617, 0.0)), ((0.583833, 0.000321, 0.233333), (-0.051907, -0.998652, 
            0.0)), ((0.566831, 0.001066, 0.233333), (-0.022329, -0.999751, 0.0)), ((
            0.55165, 0.001302, 0.233333), (-0.004709, -0.999989, 0.0)), ((0.533034, 
            0.001276, 0.233333), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 
            0.233333), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.233333), (
            0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.233333), (0.055871, 
            -0.998438, 0.0)), ((0.416659, -0.003417, 0.233333), (0.055871, -0.998438, 
            0.0)), ((0.369303, -0.006193, 0.233333), (0.067641, -0.99771, 0.0)), ((
            0.33714, -0.008398, 0.233333), (0.069854, -0.997557, 0.0)), ((0.304702, 
            -0.010667, 0.233333), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 
            0.233333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.233333), (
            0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.233333), (0.0, -1.0, 
            0.0)), ((0.129143, -0.019629, 0.233333), (0.0, -1.0, 0.0)), ((0.101105, 
            -0.019614, 0.233333), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 
            0.233333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.233333), 
            (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.233333), (-0.12352, 
            -0.992342, 0.0)), ((0.025938, -0.014654, 0.233333), (-0.194505, -0.980901, 
            0.0)), ((0.014545, -0.012018, 0.233333), (-0.304425, -0.952536, 0.0)), ((
            0.006398, -0.008883, 0.233333), (-0.509755, -0.86032, 0.0)), ((0.001672, 
            -0.004386, 0.233333), (-0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 
            0.233333), (-0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.233333), (
            -0.929852, 0.367933, 0.0)), ((0.004984, 0.013419, 0.233333), (-0.818914, 
            0.573916, 0.0)), ((0.010758, 0.020575, 0.233333), (-0.701042, 0.71312, 
            0.0)), ((0.018825, 0.027712, 0.233333), (-0.591859, 0.806041, 0.0)), ((
            0.029132, 0.034609, 0.233333), (-0.49082, 0.871261, 0.0)), ((0.041644, 
            0.041062, 0.233333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 
            0.233333), (-0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.233333), (
            -0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.233333), (-0.256858, 
            0.966449, 0.0)), ((0.119267, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((
            0.198897, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 
            0.233333), (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.233333), (
            0.103513, 0.994628, 0.0)), ((0.295248, 0.056119, 0.233333), (0.127805, 
            0.991799, 0.0)), ((0.324299, 0.052166, 0.233333), (0.148772, 0.988872, 
            0.0)), ((0.353335, 0.047615, 0.233333), (0.167009, 0.985955, 0.0)), ((
            0.393406, 0.040285, 0.233333), (0.191844, 0.981425, 0.0)), ((0.44176, 
            0.030833, 0.233333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 
            0.233333), (0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.233333), (
            0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.233333), (0.20733, 
            0.978271, 0.0)), ((0.53112, 0.01218, 0.233333), (0.208599, 0.978001, 0.0)), 
            ((0.548842, 0.00839, 0.233333), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.233333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 
            0.233333), (0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.133333), (
            -0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.133333), (-0.929852, 
            0.367933, 0.0)), ((0.004984, 0.013419, 0.133333), (-0.818914, 0.573916, 
            0.0)), ((0.010758, 0.020575, 0.133333), (-0.701042, 0.71312, 0.0)), ((
            0.018825, 0.027712, 0.133333), (-0.591859, 0.806041, 0.0)), ((0.029132, 
            0.034609, 0.133333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 
            0.133333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.133333), (
            -0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.133333), (-0.29229, 
            0.95633, 0.0)), ((0.093229, 0.058198, 0.133333), (-0.256858, 0.966449, 
            0.0)), ((0.119267, 0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.198897, 
            0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.133333), (
            0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.133333), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.133333), (0.127805, 0.991799, 
            0.0)), ((0.324299, 0.052166, 0.133333), (0.148772, 0.988872, 0.0)), ((
            0.353335, 0.047615, 0.133333), (0.167009, 0.985955, 0.0)), ((0.393406, 
            0.040285, 0.133333), (0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 
            0.133333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 0.133333), (
            0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.133333), (0.206858, 
            0.978371, 0.0)), ((0.5108, 0.016494, 0.133333), (0.20733, 0.978271, 0.0)), 
            ((0.53112, 0.01218, 0.133333), (0.208599, 0.978001, 0.0)), ((0.548842, 
            0.00839, 0.133333), (0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 
            0.133333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.133333), (
            0.187218, 0.982318, 0.0)), ((0.583833, 0.000321, 0.133333), (-0.051907, 
            -0.998652, 0.0)), ((0.566831, 0.001066, 0.133333), (-0.022329, -0.999751, 
            0.0)), ((0.55165, 0.001302, 0.133333), (-0.004709, -0.999989, 0.0)), ((
            0.533034, 0.001276, 0.133333), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.133333), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 
            0.133333), (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.133333), (
            0.055871, -0.998438, 0.0)), ((0.416659, -0.003417, 0.133333), (0.055871, 
            -0.998438, 0.0)), ((0.369303, -0.006193, 0.133333), (0.067641, -0.99771, 
            0.0)), ((0.33714, -0.008398, 0.133333), (0.069854, -0.997557, 0.0)), ((
            0.304702, -0.010667, 0.133333), (0.069674, -0.99757, 0.0)), ((0.272393, 
            -0.01333, 0.133333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 
            0.133333), (0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.133333), (
            0.0, -1.0, 0.0)), ((0.129143, -0.019629, 0.133333), (0.0, -1.0, 0.0)), ((
            0.101105, -0.019614, 0.133333), (-0.001899, -0.999998, 0.0)), ((0.077593, 
            -0.019307, 0.133333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 
            0.133333), (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.133333), 
            (-0.12352, -0.992342, 0.0)), ((0.025938, -0.014654, 0.133333), (-0.194505, 
            -0.980901, 0.0)), ((0.014545, -0.012018, 0.133333), (-0.304425, -0.952536, 
            0.0)), ((0.006398, -0.008883, 0.133333), (-0.509755, -0.86032, 0.0)), ((
            0.001672, -0.004386, 0.133333), (-0.934421, -0.35617, 0.0)), ((0.583833, 
            0.000321, 0.05), (-0.051907, -0.998652, 0.0)), ((0.004984, 0.013419, 0.05), 
            (-0.818914, 0.573916, 0.0)), ((0.010758, 0.020575, 0.05), (-0.701042, 
            0.71312, 0.0)), ((0.018825, 0.027712, 0.05), (-0.591859, 0.806041, 0.0)), (
            (0.029132, 0.034609, 0.05), (-0.49082, 0.871261, 0.0)), ((0.041644, 
            0.041062, 0.05), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.05), 
            (-0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.05), (-0.29229, 
            0.95633, 0.0)), ((0.093229, 0.058198, 0.05), (-0.256858, 0.966449, 0.0)), (
            (0.119267, 0.062493, 0.05), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.05), 
            (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.05), (0.055924, 0.998435, 0.0)), 
            ((0.266448, 0.059356, 0.05), (0.103513, 0.994628, 0.0)), ((0.295248, 
            0.056119, 0.05), (0.127805, 0.991799, 0.0)), ((0.324299, 0.052166, 0.05), (
            0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 0.05), (0.167009, 
            0.985955, 0.0)), ((0.393406, 0.040285, 0.05), (0.191844, 0.981425, 0.0)), (
            (0.463659, 0.026434, 0.05), (0.205119, 0.978737, 0.0)), ((0.488195, 
            0.021277, 0.05), (0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.05), (
            0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 0.05), (0.208599, 0.978001, 
            0.0)), ((0.548842, 0.00839, 0.05), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.05), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.05), (
            0.187218, 0.982318, 0.0)), ), name='skin-root-NS')    
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.172353, -0.019629, 
            0.016667), ), ((0.147353, -0.019629, 0.016667), ), ), name=
            'SCB-Singularity')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.155687, -0.019629, 
            0.05), ), ((0.180687, -0.019629, 0.05), ), ((0.175239, -0.019389, 0.58), ), 
            ((0.150239, -0.019389, 0.58), ), ((0.17304, -0.019572, 0.526667), ), ((
            0.14804, -0.019572, 0.526667), ), ((0.180687, -0.019629, 0.476667), ), ((
            0.155687, -0.019629, 0.476667), ), ((0.180687, -0.019629, 0.42), ), ((
            0.155687, -0.019629, 0.42), ), ((0.180687, -0.019629, 0.333333), ), ((
            0.155687, -0.019629, 0.333333), ), ((0.180687, -0.019629, 0.233333), ), ((
            0.155687, -0.019629, 0.233333), ), ((0.180687, -0.019629, 0.133333), ), ((
            0.155687, -0.019629, 0.133333), ), ), name='SCB-NS')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.155687, 0.062493, 
            0.016667), ), ((0.180687, 0.062493, 0.016667), ), ), name=
            'SCT-Singularity')
        mdb.models['Model-1'].parts['Wing'].Set(faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.150239, 0.061727, 
            0.58), ), ((0.175239, 0.061727, 0.58), ), ((0.14804, 0.062311, 0.526667), 
            ), ((0.17304, 0.062311, 0.526667), ), ((0.147353, 0.062493, 0.476667), ), (
            (0.172353, 0.062493, 0.476667), ), ((0.147353, 0.062493, 0.42), ), ((
            0.172353, 0.062493, 0.42), ), ((0.147353, 0.062493, 0.333333), ), ((
            0.172353, 0.062493, 0.333333), ), ((0.147353, 0.062493, 0.233333), ), ((
            0.172353, 0.062493, 0.233333), ), ((0.147353, 0.062493, 0.133333), ), ((
            0.172353, 0.062493, 0.133333), ), ((0.147353, 0.062493, 0.05), ), ((
            0.172353, 0.062493, 0.05), ), ), name='SCT-NS')    

        mdb.models['Model-1'].parts['Wing'].Surface(name='SCB-NS', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.155687, -0.019629, 
            0.05), ), ((0.180687, -0.019629, 0.05), ), ((0.175239, -0.019389, 0.58), ), 
            ((0.150239, -0.019389, 0.58), ), ((0.17304, -0.019572, 0.526667), ), ((
            0.14804, -0.019572, 0.526667), ), ((0.180687, -0.019629, 0.476667), ), ((
            0.155687, -0.019629, 0.476667), ), ((0.180687, -0.019629, 0.42), ), ((
            0.155687, -0.019629, 0.42), ), ((0.180687, -0.019629, 0.333333), ), ((
            0.155687, -0.019629, 0.333333), ), ((0.180687, -0.019629, 0.233333), ), ((
            0.155687, -0.019629, 0.233333), ), ((0.180687, -0.019629, 0.133333), ), ((
            0.155687, -0.019629, 0.133333), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='SCB-Singularity', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.172353, -0.019629, 
            0.016667), ), ((0.147353, -0.019629, 0.016667), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='SCT-NS', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.150239, 0.061727, 
            0.58), ), ((0.175239, 0.061727, 0.58), ), ((0.14804, 0.062311, 0.526667), 
            ), ((0.17304, 0.062311, 0.526667), ), ((0.147353, 0.062493, 0.476667), ), (
            (0.172353, 0.062493, 0.476667), ), ((0.147353, 0.062493, 0.42), ), ((
            0.172353, 0.062493, 0.42), ), ((0.147353, 0.062493, 0.333333), ), ((
            0.172353, 0.062493, 0.333333), ), ((0.147353, 0.062493, 0.233333), ), ((
            0.172353, 0.062493, 0.233333), ), ((0.147353, 0.062493, 0.133333), ), ((
            0.172353, 0.062493, 0.133333), ), ((0.147353, 0.062493, 0.05), ), ((
            0.172353, 0.062493, 0.05), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='SCT-Singularity', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.155687, 0.062493, 
            0.016667), ), ((0.180687, 0.062493, 0.016667), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='Singularity-skin', 
            side1Faces=mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.006651, 
            0.015798, 0.016667), ), ((0.013196, 0.022973, 0.016667), ), ((0.022015, 
            0.030055, 0.016667), ), ((0.03306, 0.036822, 0.016667), ), ((0.046301, 
            0.04309, 0.016667), ), ((0.061681, 0.048716, 0.016667), ), ((0.079113, 
            0.054205, 0.016667), ), ((0.10131, 0.060345, 0.016667), ), ((0.129143, 
            0.062493, 0.016667), ), ((0.208773, 0.062493, 0.016667), ), ((0.244141, 
            0.061065, 0.016667), ), ((0.27601, 0.058361, 0.016667), ), ((0.304924, 
            0.054872, 0.016667), ), ((0.333999, 0.050707, 0.016667), ), ((0.362972, 
            0.045983, 0.016667), ), ((0.414203, 0.03622, 0.016667), ), ((0.472039, 
            0.024678, 0.016667), ), ((0.495972, 0.019633, 0.016667), ), ((0.517853, 
            0.015, 0.016667), ), ((0.537335, 0.010854, 0.016667), ), ((0.554134, 
            0.007251, 0.016667), ), ((0.568024, 0.004283, 0.016667), ), ((0.584108, 
            0.001123, 0.016667), ), ((0.577665, 0.000641, 0.016667), ), ((0.562164, 
            0.00117, 0.016667), ), ((0.545803, 0.001329, 0.016667), ), ((0.526111, 
            0.001196, 0.016667), ), ((0.503423, 0.000714, 0.016667), ), ((0.478105, 
            -0.00014, 0.016667), ), ((0.446462, -0.001749, 0.016667), ), ((0.398318, 
            -0.004443, 0.016667), ), ((0.358628, -0.006917, 0.016667), ), ((0.326327, 
            -0.009155, 0.016667), ), ((0.293889, -0.011422, 0.016667), ), ((0.26171, 
            -0.014483, 0.016667), ), ((0.229442, -0.018298, 0.016667), ), ((0.198897, 
            -0.019629, 0.016667), ), ((0.119267, -0.019629, 0.016667), ), ((0.09282, 
            -0.019598, 0.016667), ), ((0.070651, -0.019031, 0.016667), ), ((0.051562, 
            -0.017847, 0.016667), ), ((0.035219, -0.016115, 0.016667), ), ((0.021793, 
            -0.013832, 0.016667), ), ((0.011444, -0.011027, 0.016667), ), ((0.004453, 
            -0.007731, 0.016667), ), ((0.000836, -0.002193, 0.016667), ), ((0.000405, 
            0.002797, 0.016667), ), ((0.002413, 0.008758, 0.016667), ), ((0.44852, 
            0.029512, 0.016667), ), ))
        mdb.models['Model-1'].parts['Wing'].Surface(name='skin-root-NS', side1Faces=
            mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.00151, 0.006476, 
            0.05), (-0.929852, 0.367933, 0.0)), ((0.000203, 0.001398, 0.05), (
            -0.989669, 0.14337, 0.0)), ((0.001672, -0.004386, 0.05), (-0.934421, 
            -0.35617, 0.0)), ((0.006398, -0.008883, 0.05), (-0.509755, -0.86032, 0.0)), 
            ((0.014545, -0.012018, 0.05), (-0.304425, -0.952536, 0.0)), ((0.025938, 
            -0.014654, 0.05), (-0.194505, -0.980901, 0.0)), ((0.040354, -0.016754, 
            0.05), (-0.12352, -0.992342, 0.0)), ((0.057635, -0.018302, 0.05), (
            -0.074597, -0.997214, 0.0)), ((0.077593, -0.019307, 0.05), (-0.039629, 
            -0.999214, 0.0)), ((0.101105, -0.019614, 0.05), (-0.001899, -0.999998, 
            0.0)), ((0.129143, -0.019629, 0.05), (0.0, -1.0, 0.0)), ((0.208773, 
            -0.019629, 0.05), (0.0, -1.0, 0.0)), ((0.240235, -0.016966, 0.05), (
            0.122438, -0.992476, 0.0)), ((0.272393, -0.01333, 0.05), (0.107257, 
            -0.994231, 0.0)), ((0.304702, -0.010667, 0.05), (0.069674, -0.99757, 0.0)), 
            ((0.33714, -0.008398, 0.05), (0.069854, -0.997557, 0.0)), ((0.369303, 
            -0.006193, 0.05), (0.067641, -0.99771, 0.0)), ((0.44176, 0.030833, 0.05), (
            0.191844, 0.981425, 0.0)), ((0.416659, -0.003417, 0.05), (0.055871, 
            -0.998438, 0.0)), ((0.457924, -0.001108, 0.05), (0.055871, -0.998438, 
            0.0)), ((0.486823, 0.000187, 0.05), (0.03742, -0.9993, 0.0)), ((0.511306, 
            0.000915, 0.05), (0.025441, -0.999676, 0.0)), ((0.533034, 0.001276, 0.05), 
            (0.011647, -0.999932, 0.0)), ((0.55165, 0.001302, 0.05), (-0.004709, 
            -0.999989, 0.0)), ((0.566831, 0.001066, 0.05), (-0.022329, -0.999751, 
            0.0)), ((0.005096, 0.001386, 0.58), (-0.987308, 0.143028, 0.069039)), ((
            0.00639, 0.006405, 0.58), (-0.927883, 0.367154, 0.065055)), ((0.009823, 
            0.013263, 0.58), (-0.817545, 0.572957, 0.057801)), ((0.015529, 0.020332, 
            0.58), (-0.700158, 0.71222, 0.050206)), ((0.0235, 0.027381, 0.58), (
            -0.591304, 0.805285, 0.043305)), ((0.033684, 0.034193, 0.58), (-0.490483, 
            0.870662, 0.03707)), ((0.046045, 0.040566, 0.58), (-0.398994, 0.91641, 
            0.031577)), ((0.060543, 0.046349, 0.58), (-0.318094, 0.947678, 0.026888)), 
            ((0.077098, 0.051726, 0.58), (-0.292195, 0.95602, 0.025454)), ((0.097011, 
            0.057492, 0.58), (-0.256786, 0.96618, 0.023606)), ((0.122663, 0.061727, 
            0.58), (0.0, 0.99994, 0.010933)), ((0.201623, 0.061727, 0.58), (0.0, 
            0.99994, 0.010933)), ((0.233501, 0.06102, 0.58), (0.055922, 0.998393, 
            0.009145)), ((0.268113, 0.058626, 0.58), (0.10351, 0.994597, 0.007916)), ((
            0.296561, 0.055427, 0.58), (0.127802, 0.991772, 0.007403)), ((0.325256, 
            0.051522, 0.58), (0.148768, 0.988847, 0.007063)), ((0.353937, 0.047026, 
            0.58), (0.167005, 0.985932, 0.006859)), ((0.393658, 0.039757, 0.58), (
            0.19184, 0.981403, 0.006707)), ((0.448001, 0.029135, 0.58), (0.19184, 
            0.981403, 0.006707)), ((0.471124, 0.024382, 0.58), (0.205114, 0.978714, 
            0.006823)), ((0.494765, 0.019399, 0.58), (0.206854, 0.978348, 0.006846)), (
            (0.516381, 0.014821, 0.58), (0.207325, 0.978248, 0.006854)), ((0.535627, 
            0.010726, 0.58), (0.208594, 0.977978, 0.006881)), ((0.552224, 0.007166, 
            0.58), (0.210343, 0.977603, 0.006925)), ((0.565948, 0.004234, 0.58), (
            0.208057, 0.978093, 0.006861)), ((0.581828, 0.001113, 0.58), (0.187214, 
            0.982299, 0.006234)), ((0.581556, 0.000318, 0.58), (-0.051907, -0.99865, 
            -0.001728)), ((0.564767, 0.001053, 0.58), (-0.022329, -0.99975, 
            -0.000839)), ((0.549768, 0.001286, 0.58), (-0.004709, -0.999989, 
            -0.000353)), ((0.531376, 0.00126, 0.58), (0.011647, -0.999932, 4.8e-05)), (
            (0.509911, 0.000903, 0.58), (0.025441, -0.999676, 0.000337)), ((0.485725, 
            0.000183, 0.58), (0.03742, -0.999299, 0.000538)), ((0.457274, -0.00109, 
            0.58), (0.055871, -0.998438, 0.000763)), ((0.398502, -0.004379, 0.58), (
            0.055871, -0.998438, 0.000763)), ((0.359168, -0.00683, 0.58), (0.067641, 
            -0.997709, 0.000722)), ((0.327263, -0.00904, 0.58), (0.069854, -0.997557, 
            0.000701)), ((0.295222, -0.01128, 0.58), (0.069674, -0.99757, 0.000704)), (
            (0.263437, -0.014301, 0.58), (0.107257, -0.994231, -7e-05)), ((0.231565, 
            -0.018069, 0.58), (0.122438, -0.992476, -0.000469)), ((0.201623, -0.019389, 
            0.58), (0.0, -0.999994, 0.003434)), ((0.122663, -0.019389, 0.58), (0.0, 
            -0.999994, 0.003434)), ((0.096607, -0.019358, 0.58), (-0.001899, -0.999992, 
            0.003531)), ((0.074704, -0.018799, 0.58), (-0.039628, -0.999199, 
            0.005608)), ((0.055847, -0.017631, 0.58), (-0.074595, -0.997185, 
            0.007657)), ((0.039701, -0.01592, 0.58), (-0.123513, -0.992286, 0.010673)), 
            ((0.026435, -0.013665, 0.58), (-0.194483, -0.980788, 0.015231)), ((0.01621, 
            -0.010896, 0.58), (-0.304348, -0.952295, 0.02251)), ((0.0093, -0.00764, 
            0.58), (-0.509417, -0.85975, 0.036388)), ((0.005723, -0.002174, 0.58), (
            -0.932433, -0.355412, 0.065202)), ((0.583278, 0.000321, 0.526667), (
            -0.051907, -0.99865, -0.001728)), ((0.56633, 0.001063, 0.526667), (
            -0.022329, -0.99975, -0.000839)), ((0.55119, 0.001298, 0.526667), (
            -0.004709, -0.999989, -0.000353)), ((0.532624, 0.001272, 0.526667), (
            0.011647, -0.999932, 4.8e-05)), ((0.510957, 0.000911, 0.526667), (0.025441, 
            -0.999676, 0.000337)), ((0.486543, 0.000185, 0.526667), (0.03742, 
            -0.999299, 0.000538)), ((0.457721, -0.001106, 0.526667), (0.055871, 
            -0.998438, 0.000763)), ((0.398376, -0.004427, 0.526667), (0.055871, 
            -0.998438, 0.000763)), ((0.358778, -0.006894, 0.526667), (0.067641, 
            -0.997709, 0.000722)), ((0.326572, -0.009126, 0.526667), (0.069854, 
            -0.997557, 0.000701)), ((0.294229, -0.011387, 0.526667), (0.069674, 
            -0.99757, 0.000704)), ((0.262144, -0.014437, 0.526667), (0.107257, 
            -0.994231, -7e-05)), ((0.22997, -0.018241, 0.526667), (0.122438, -0.992476, 
            -0.000469)), ((0.199584, -0.019572, 0.526667), (0.0, -0.999994, 0.003434)), 
            ((0.120113, -0.019572, 0.526667), (0.0, -0.999994, 0.003434)), ((0.093739, 
            -0.019541, 0.526667), (-0.001899, -0.999992, 0.003531)), ((0.07163, 
            -0.018977, 0.526667), (-0.039628, -0.999199, 0.005608)), ((0.052595, 
            -0.017797, 0.526667), (-0.074595, -0.997185, 0.007657)), ((0.036297, 
            -0.01607, 0.526667), (-0.123513, -0.992286, 0.010673)), ((0.022907, 
            -0.013794, 0.526667), (-0.194483, -0.980788, 0.015231)), ((0.012585, 
            -0.010998, 0.526667), (-0.304348, -0.952295, 0.02251)), ((0.005611, 
            -0.007712, 0.526667), (-0.509417, -0.85975, 0.036388)), ((0.002001, 
            -0.002193, 0.526667), (-0.932433, -0.355412, 0.065202)), ((0.001368, 
            0.001398, 0.526667), (-0.987308, 0.143028, 0.069039)), ((0.002674, 
            0.006464, 0.526667), (-0.927883, 0.367154, 0.065055)), ((0.006139, 
            0.013386, 0.526667), (-0.817545, 0.572957, 0.057801)), ((0.011899, 
            0.020522, 0.526667), (-0.700158, 0.71222, 0.050206)), ((0.019945, 0.027638, 
            0.526667), (-0.591304, 0.805285, 0.043305)), ((0.030224, 0.034515, 
            0.526667), (-0.490483, 0.870662, 0.03707)), ((0.042702, 0.040948, 
            0.526667), (-0.398994, 0.91641, 0.031577)), ((0.057336, 0.046785, 
            0.526667), (-0.318094, 0.947678, 0.026888)), ((0.074047, 0.052213, 
            0.526667), (-0.292195, 0.95602, 0.025454)), ((0.094146, 0.058034, 
            0.526667), (-0.256786, 0.96618, 0.023606)), ((0.120113, 0.062311, 
            0.526667), (0.0, 0.99994, 0.010933)), ((0.199584, 0.062311, 0.526667), (
            0.0, 0.99994, 0.010933)), ((0.231923, 0.061597, 0.526667), (0.055922, 
            0.998393, 0.009145)), ((0.266864, 0.05918, 0.526667), (0.10351, 0.994597, 
            0.007916)), ((0.295581, 0.055952, 0.526667), (0.127802, 0.991772, 
            0.007403)), ((0.324547, 0.052009, 0.526667), (0.148768, 0.988847, 
            0.007063)), ((0.353498, 0.047471, 0.526667), (0.167005, 0.985932, 
            0.006859)), ((0.393485, 0.040156, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.448357, 0.029429, 0.526667), (0.19184, 0.981403, 
            0.006707)), ((0.471804, 0.024611, 0.526667), (0.205114, 0.978714, 
            0.006823)), ((0.495668, 0.019581, 0.526667), (0.206854, 0.978348, 
            0.006846)), ((0.517488, 0.01496, 0.526667), (0.207325, 0.978248, 
            0.006854)), ((0.536915, 0.010826, 0.526667), (0.208594, 0.977978, 
            0.006881)), ((0.553668, 0.007233, 0.526667), (0.210343, 0.977603, 
            0.006925)), ((0.567521, 0.004274, 0.526667), (0.208057, 0.978093, 
            0.006861)), ((0.583553, 0.001123, 0.526667), (0.187214, 0.982299, 
            0.006234)), ((0.000203, 0.001398, 0.476667), (-0.989669, 0.14337, 0.0)), ((
            0.00151, 0.006476, 0.476667), (-0.929852, 0.367933, 0.0)), ((0.004984, 
            0.013419, 0.476667), (-0.818914, 0.573916, 0.0)), ((0.010758, 0.020575, 
            0.476667), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.476667), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.476667), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.476667), (-0.399193, 0.916867, 
            0.0)), ((0.056319, 0.046917, 0.476667), (-0.318209, 0.94802, 0.0)), ((
            0.073077, 0.052361, 0.476667), (-0.29229, 0.95633, 0.0)), ((0.093229, 
            0.058198, 0.476667), (-0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 
            0.476667), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.476667), (0.0, 1.0, 
            0.0)), ((0.231395, 0.061779, 0.476667), (0.055924, 0.998435, 0.0)), ((
            0.266448, 0.059356, 0.476667), (0.103513, 0.994628, 0.0)), ((0.295248, 
            0.056119, 0.476667), (0.127805, 0.991799, 0.0)), ((0.324299, 0.052166, 
            0.476667), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 0.476667), (
            0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.476667), (0.191844, 
            0.981425, 0.0)), ((0.44176, 0.030833, 0.476667), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.476667), (0.205119, 0.978737, 0.0)), ((
            0.488195, 0.021277, 0.476667), (0.206858, 0.978371, 0.0)), ((0.5108, 
            0.016494, 0.476667), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.476667), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.476667), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.476667), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.476667), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.476667), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.476667), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.476667), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.476667), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.476667), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.476667), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.476667), (0.055871, -0.998438, 
            0.0)), ((0.416659, -0.003417, 0.476667), (0.055871, -0.998438, 0.0)), ((
            0.369303, -0.006193, 0.476667), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.476667), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.476667), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.476667), (
            0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.476667), (0.122438, 
            -0.992476, 0.0)), ((0.208773, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((
            0.129143, -0.019629, 0.476667), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 
            0.476667), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.476667), 
            (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.476667), (-0.074597, 
            -0.997214, 0.0)), ((0.040354, -0.016754, 0.476667), (-0.12352, -0.992342, 
            0.0)), ((0.025938, -0.014654, 0.476667), (-0.194505, -0.980901, 0.0)), ((
            0.014545, -0.012018, 0.476667), (-0.304425, -0.952536, 0.0)), ((0.006398, 
            -0.008883, 0.476667), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 
            0.476667), (-0.934421, -0.35617, 0.0)), ((0.583833, 0.000321, 0.42), (
            -0.051907, -0.998652, 0.0)), ((0.566831, 0.001066, 0.42), (-0.022329, 
            -0.999751, 0.0)), ((0.55165, 0.001302, 0.42), (-0.004709, -0.999989, 0.0)), 
            ((0.533034, 0.001276, 0.42), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.42), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.42), 
            (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.42), (0.055871, 
            -0.998438, 0.0)), ((0.416659, -0.003417, 0.42), (0.055871, -0.998438, 
            0.0)), ((0.369303, -0.006193, 0.42), (0.067641, -0.99771, 0.0)), ((0.33714, 
            -0.008398, 0.42), (0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 
            0.42), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 0.42), (0.107257, 
            -0.994231, 0.0)), ((0.240235, -0.016966, 0.42), (0.122438, -0.992476, 
            0.0)), ((0.208773, -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.129143, 
            -0.019629, 0.42), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.42), (
            -0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 0.42), (-0.039629, 
            -0.999214, 0.0)), ((0.057635, -0.018302, 0.42), (-0.074597, -0.997214, 
            0.0)), ((0.040354, -0.016754, 0.42), (-0.12352, -0.992342, 0.0)), ((
            0.025938, -0.014654, 0.42), (-0.194505, -0.980901, 0.0)), ((0.014545, 
            -0.012018, 0.42), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 
            0.42), (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.42), (
            -0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 0.42), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.42), (-0.929852, 0.367933, 0.0)), ((
            0.004984, 0.013419, 0.42), (-0.818914, 0.573916, 0.0)), ((0.010758, 
            0.020575, 0.42), (-0.701042, 0.71312, 0.0)), ((0.018825, 0.027712, 0.42), (
            -0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 0.42), (-0.49082, 
            0.871261, 0.0)), ((0.041644, 0.041062, 0.42), (-0.399193, 0.916867, 0.0)), 
            ((0.056319, 0.046917, 0.42), (-0.318209, 0.94802, 0.0)), ((0.073077, 
            0.052361, 0.42), (-0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.42), (
            -0.256858, 0.966449, 0.0)), ((0.119267, 0.062493, 0.42), (0.0, 1.0, 0.0)), 
            ((0.198897, 0.062493, 0.42), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.42), 
            (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.42), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.42), (0.127805, 0.991799, 0.0)), (
            (0.324299, 0.052166, 0.42), (0.148772, 0.988872, 0.0)), ((0.353335, 
            0.047615, 0.42), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.42), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.42), (0.191844, 0.981425, 
            0.0)), ((0.463659, 0.026434, 0.42), (0.205119, 0.978737, 0.0)), ((0.488195, 
            0.021277, 0.42), (0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.42), (
            0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 0.42), (0.208599, 0.978001, 
            0.0)), ((0.548842, 0.00839, 0.42), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.42), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.42), (
            0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.333333), (-0.989669, 
            0.14337, 0.0)), ((0.00151, 0.006476, 0.333333), (-0.929852, 0.367933, 
            0.0)), ((0.004984, 0.013419, 0.333333), (-0.818914, 0.573916, 0.0)), ((
            0.010758, 0.020575, 0.333333), (-0.701042, 0.71312, 0.0)), ((0.018825, 
            0.027712, 0.333333), (-0.591859, 0.806041, 0.0)), ((0.029132, 0.034609, 
            0.333333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 0.333333), (
            -0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.333333), (-0.318209, 
            0.94802, 0.0)), ((0.073077, 0.052361, 0.333333), (-0.29229, 0.95633, 0.0)), 
            ((0.093229, 0.058198, 0.333333), (-0.256858, 0.966449, 0.0)), ((0.119267, 
            0.062493, 0.333333), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.333333), (
            0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.333333), (0.055924, 0.998435, 
            0.0)), ((0.266448, 0.059356, 0.333333), (0.103513, 0.994628, 0.0)), ((
            0.295248, 0.056119, 0.333333), (0.127805, 0.991799, 0.0)), ((0.324299, 
            0.052166, 0.333333), (0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 
            0.333333), (0.167009, 0.985955, 0.0)), ((0.393406, 0.040285, 0.333333), (
            0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 0.333333), (0.191844, 
            0.981425, 0.0)), ((0.463659, 0.026434, 0.333333), (0.205119, 0.978737, 
            0.0)), ((0.488195, 0.021277, 0.333333), (0.206858, 0.978371, 0.0)), ((
            0.5108, 0.016494, 0.333333), (0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 
            0.333333), (0.208599, 0.978001, 0.0)), ((0.548842, 0.00839, 0.333333), (
            0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 0.333333), (0.208062, 
            0.978116, 0.0)), ((0.578216, 0.002246, 0.333333), (0.187218, 0.982318, 
            0.0)), ((0.583833, 0.000321, 0.333333), (-0.051907, -0.998652, 0.0)), ((
            0.566831, 0.001066, 0.333333), (-0.022329, -0.999751, 0.0)), ((0.55165, 
            0.001302, 0.333333), (-0.004709, -0.999989, 0.0)), ((0.533034, 0.001276, 
            0.333333), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 0.333333), (
            0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.333333), (0.03742, 
            -0.9993, 0.0)), ((0.457924, -0.001108, 0.333333), (0.055871, -0.998438, 
            0.0)), )+mdb.models['Model-1'].parts['Wing'].faces.findAt(((0.416659, 
            -0.003417, 0.333333), (0.055871, -0.998438, 0.0)), ((0.369303, -0.006193, 
            0.333333), (0.067641, -0.99771, 0.0)), ((0.33714, -0.008398, 0.333333), (
            0.069854, -0.997557, 0.0)), ((0.304702, -0.010667, 0.333333), (0.069674, 
            -0.99757, 0.0)), ((0.272393, -0.01333, 0.333333), (0.107257, -0.994231, 
            0.0)), ((0.240235, -0.016966, 0.333333), (0.122438, -0.992476, 0.0)), ((
            0.208773, -0.019629, 0.333333), (0.0, -1.0, 0.0)), ((0.129143, -0.019629, 
            0.333333), (0.0, -1.0, 0.0)), ((0.101105, -0.019614, 0.333333), (-0.001899, 
            -0.999998, 0.0)), ((0.077593, -0.019307, 0.333333), (-0.039629, -0.999214, 
            0.0)), ((0.057635, -0.018302, 0.333333), (-0.074597, -0.997214, 0.0)), ((
            0.040354, -0.016754, 0.333333), (-0.12352, -0.992342, 0.0)), ((0.025938, 
            -0.014654, 0.333333), (-0.194505, -0.980901, 0.0)), ((0.014545, -0.012018, 
            0.333333), (-0.304425, -0.952536, 0.0)), ((0.006398, -0.008883, 0.333333), 
            (-0.509755, -0.86032, 0.0)), ((0.001672, -0.004386, 0.333333), (-0.934421, 
            -0.35617, 0.0)), ((0.583833, 0.000321, 0.233333), (-0.051907, -0.998652, 
            0.0)), ((0.566831, 0.001066, 0.233333), (-0.022329, -0.999751, 0.0)), ((
            0.55165, 0.001302, 0.233333), (-0.004709, -0.999989, 0.0)), ((0.533034, 
            0.001276, 0.233333), (0.011647, -0.999932, 0.0)), ((0.511306, 0.000915, 
            0.233333), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 0.233333), (
            0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.233333), (0.055871, 
            -0.998438, 0.0)), ((0.416659, -0.003417, 0.233333), (0.055871, -0.998438, 
            0.0)), ((0.369303, -0.006193, 0.233333), (0.067641, -0.99771, 0.0)), ((
            0.33714, -0.008398, 0.233333), (0.069854, -0.997557, 0.0)), ((0.304702, 
            -0.010667, 0.233333), (0.069674, -0.99757, 0.0)), ((0.272393, -0.01333, 
            0.233333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 0.233333), (
            0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.233333), (0.0, -1.0, 
            0.0)), ((0.129143, -0.019629, 0.233333), (0.0, -1.0, 0.0)), ((0.101105, 
            -0.019614, 0.233333), (-0.001899, -0.999998, 0.0)), ((0.077593, -0.019307, 
            0.233333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 0.233333), 
            (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.233333), (-0.12352, 
            -0.992342, 0.0)), ((0.025938, -0.014654, 0.233333), (-0.194505, -0.980901, 
            0.0)), ((0.014545, -0.012018, 0.233333), (-0.304425, -0.952536, 0.0)), ((
            0.006398, -0.008883, 0.233333), (-0.509755, -0.86032, 0.0)), ((0.001672, 
            -0.004386, 0.233333), (-0.934421, -0.35617, 0.0)), ((0.000203, 0.001398, 
            0.233333), (-0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.233333), (
            -0.929852, 0.367933, 0.0)), ((0.004984, 0.013419, 0.233333), (-0.818914, 
            0.573916, 0.0)), ((0.010758, 0.020575, 0.233333), (-0.701042, 0.71312, 
            0.0)), ((0.018825, 0.027712, 0.233333), (-0.591859, 0.806041, 0.0)), ((
            0.029132, 0.034609, 0.233333), (-0.49082, 0.871261, 0.0)), ((0.041644, 
            0.041062, 0.233333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 
            0.233333), (-0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.233333), (
            -0.29229, 0.95633, 0.0)), ((0.093229, 0.058198, 0.233333), (-0.256858, 
            0.966449, 0.0)), ((0.119267, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((
            0.198897, 0.062493, 0.233333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 
            0.233333), (0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.233333), (
            0.103513, 0.994628, 0.0)), ((0.295248, 0.056119, 0.233333), (0.127805, 
            0.991799, 0.0)), ((0.324299, 0.052166, 0.233333), (0.148772, 0.988872, 
            0.0)), ((0.353335, 0.047615, 0.233333), (0.167009, 0.985955, 0.0)), ((
            0.393406, 0.040285, 0.233333), (0.191844, 0.981425, 0.0)), ((0.44176, 
            0.030833, 0.233333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 
            0.233333), (0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.233333), (
            0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.233333), (0.20733, 
            0.978271, 0.0)), ((0.53112, 0.01218, 0.233333), (0.208599, 0.978001, 0.0)), 
            ((0.548842, 0.00839, 0.233333), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.233333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 
            0.233333), (0.187218, 0.982318, 0.0)), ((0.000203, 0.001398, 0.133333), (
            -0.989669, 0.14337, 0.0)), ((0.00151, 0.006476, 0.133333), (-0.929852, 
            0.367933, 0.0)), ((0.004984, 0.013419, 0.133333), (-0.818914, 0.573916, 
            0.0)), ((0.010758, 0.020575, 0.133333), (-0.701042, 0.71312, 0.0)), ((
            0.018825, 0.027712, 0.133333), (-0.591859, 0.806041, 0.0)), ((0.029132, 
            0.034609, 0.133333), (-0.49082, 0.871261, 0.0)), ((0.041644, 0.041062, 
            0.133333), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.133333), (
            -0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.133333), (-0.29229, 
            0.95633, 0.0)), ((0.093229, 0.058198, 0.133333), (-0.256858, 0.966449, 
            0.0)), ((0.119267, 0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.198897, 
            0.062493, 0.133333), (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.133333), (
            0.055924, 0.998435, 0.0)), ((0.266448, 0.059356, 0.133333), (0.103513, 
            0.994628, 0.0)), ((0.295248, 0.056119, 0.133333), (0.127805, 0.991799, 
            0.0)), ((0.324299, 0.052166, 0.133333), (0.148772, 0.988872, 0.0)), ((
            0.353335, 0.047615, 0.133333), (0.167009, 0.985955, 0.0)), ((0.393406, 
            0.040285, 0.133333), (0.191844, 0.981425, 0.0)), ((0.44176, 0.030833, 
            0.133333), (0.191844, 0.981425, 0.0)), ((0.463659, 0.026434, 0.133333), (
            0.205119, 0.978737, 0.0)), ((0.488195, 0.021277, 0.133333), (0.206858, 
            0.978371, 0.0)), ((0.5108, 0.016494, 0.133333), (0.20733, 0.978271, 0.0)), 
            ((0.53112, 0.01218, 0.133333), (0.208599, 0.978001, 0.0)), ((0.548842, 
            0.00839, 0.133333), (0.210348, 0.977627, 0.0)), ((0.563725, 0.005198, 
            0.133333), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.133333), (
            0.187218, 0.982318, 0.0)), ((0.583833, 0.000321, 0.133333), (-0.051907, 
            -0.998652, 0.0)), ((0.566831, 0.001066, 0.133333), (-0.022329, -0.999751, 
            0.0)), ((0.55165, 0.001302, 0.133333), (-0.004709, -0.999989, 0.0)), ((
            0.533034, 0.001276, 0.133333), (0.011647, -0.999932, 0.0)), ((0.511306, 
            0.000915, 0.133333), (0.025441, -0.999676, 0.0)), ((0.486823, 0.000187, 
            0.133333), (0.03742, -0.9993, 0.0)), ((0.457924, -0.001108, 0.133333), (
            0.055871, -0.998438, 0.0)), ((0.416659, -0.003417, 0.133333), (0.055871, 
            -0.998438, 0.0)), ((0.369303, -0.006193, 0.133333), (0.067641, -0.99771, 
            0.0)), ((0.33714, -0.008398, 0.133333), (0.069854, -0.997557, 0.0)), ((
            0.304702, -0.010667, 0.133333), (0.069674, -0.99757, 0.0)), ((0.272393, 
            -0.01333, 0.133333), (0.107257, -0.994231, 0.0)), ((0.240235, -0.016966, 
            0.133333), (0.122438, -0.992476, 0.0)), ((0.208773, -0.019629, 0.133333), (
            0.0, -1.0, 0.0)), ((0.129143, -0.019629, 0.133333), (0.0, -1.0, 0.0)), ((
            0.101105, -0.019614, 0.133333), (-0.001899, -0.999998, 0.0)), ((0.077593, 
            -0.019307, 0.133333), (-0.039629, -0.999214, 0.0)), ((0.057635, -0.018302, 
            0.133333), (-0.074597, -0.997214, 0.0)), ((0.040354, -0.016754, 0.133333), 
            (-0.12352, -0.992342, 0.0)), ((0.025938, -0.014654, 0.133333), (-0.194505, 
            -0.980901, 0.0)), ((0.014545, -0.012018, 0.133333), (-0.304425, -0.952536, 
            0.0)), ((0.006398, -0.008883, 0.133333), (-0.509755, -0.86032, 0.0)), ((
            0.001672, -0.004386, 0.133333), (-0.934421, -0.35617, 0.0)), ((0.583833, 
            0.000321, 0.05), (-0.051907, -0.998652, 0.0)), ((0.004984, 0.013419, 0.05), 
            (-0.818914, 0.573916, 0.0)), ((0.010758, 0.020575, 0.05), (-0.701042, 
            0.71312, 0.0)), ((0.018825, 0.027712, 0.05), (-0.591859, 0.806041, 0.0)), (
            (0.029132, 0.034609, 0.05), (-0.49082, 0.871261, 0.0)), ((0.041644, 
            0.041062, 0.05), (-0.399193, 0.916867, 0.0)), ((0.056319, 0.046917, 0.05), 
            (-0.318209, 0.94802, 0.0)), ((0.073077, 0.052361, 0.05), (-0.29229, 
            0.95633, 0.0)), ((0.093229, 0.058198, 0.05), (-0.256858, 0.966449, 0.0)), (
            (0.119267, 0.062493, 0.05), (0.0, 1.0, 0.0)), ((0.198897, 0.062493, 0.05), 
            (0.0, 1.0, 0.0)), ((0.231395, 0.061779, 0.05), (0.055924, 0.998435, 0.0)), 
            ((0.266448, 0.059356, 0.05), (0.103513, 0.994628, 0.0)), ((0.295248, 
            0.056119, 0.05), (0.127805, 0.991799, 0.0)), ((0.324299, 0.052166, 0.05), (
            0.148772, 0.988872, 0.0)), ((0.353335, 0.047615, 0.05), (0.167009, 
            0.985955, 0.0)), ((0.393406, 0.040285, 0.05), (0.191844, 0.981425, 0.0)), (
            (0.463659, 0.026434, 0.05), (0.205119, 0.978737, 0.0)), ((0.488195, 
            0.021277, 0.05), (0.206858, 0.978371, 0.0)), ((0.5108, 0.016494, 0.05), (
            0.20733, 0.978271, 0.0)), ((0.53112, 0.01218, 0.05), (0.208599, 0.978001, 
            0.0)), ((0.548842, 0.00839, 0.05), (0.210348, 0.977627, 0.0)), ((0.563725, 
            0.005198, 0.05), (0.208062, 0.978116, 0.0)), ((0.578216, 0.002246, 0.05), (
            0.187218, 0.982318, 0.0)), ))

    createsurfaces()


    #step
    mdb.models['Model-1'].StaticStep(initialInc=1, maxInc=1, maxNumInc=1, 
        minInc=0.00001, name='Step-1', previous='Initial')

    #boundary condition
    mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
        name='BC-1', region=
        mdb.models['Model-1'].rootAssembly.instances['Wing-1'].sets['Base Rib'])
    #load
    load_list_test=[('load-1',((3898.30508474576)/4.9)*Load_Factor),
        ('load-2',(6322.03389830509/4.9)*Load_Factor),
        ('load-3',(4403.71176537977/4.9)*Load_Factor),
        ('load-4',(3983.51315310947/4.9)*Load_Factor),
        ('load-5',(3950.31535590166/4.9)*Load_Factor),
        ('load-6',(3827.68555357388/4.9)*Load_Factor),
        ('load-7',(3733.29416948157/4.9)*Load_Factor),
        ('load-8',(3593.11282283643/4.9)*Load_Factor)]

    for i in load_list_test:
        # mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
            # distributionType=UNIFORM, field='', magnitude=i[1], name=i[0], region=
            # mdb.models['Model-1'].rootAssembly.instances['Wing-1'].surfaces[i[0]])
        mdb.models['Model-1'].SurfaceTraction(createStepName='Step-1', directionVector=
            ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), distributionType=UNIFORM, field='', 
            follower=OFF, localCsys=None, magnitude=i[1], name=i[0], region=
            mdb.models['Model-1'].rootAssembly.instances['Wing-1'].surfaces[i[0]], 
            resultant=ON, traction=GENERAL)

    #mesh
    mdb.models['Model-1'].parts['Wing'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=global_mesh_size)
    mdb.models['Model-1'].parts['Wing'].setMeshControls(algorithm=MEDIAL_AXIS, regions=mdb.models['Model-1'].parts['Wing'].faces)
    mdb.models['Model-1'].parts['Wing'].generateMesh()

    #table of material properties

    mdb.models['Model-1'].Material(name='Airex C70 Foam')
    mdb.models['Model-1'].materials['Airex C70 Foam'].Elastic(table=((105000000.0, 0.29), ))
    mdb.models['Model-1'].materials['Airex C70 Foam'].Density(table=((130.0, ), ))

    Composite_properties=[('Glass Fibre-86GSM',21390000000.0,21390000000.0,0.133,1570000000.0,1570000000.0,1570000000.0,0.000083,1714.0,260000000.0,130000000.0,260000000.0,130000000.0,59000000.0,59000000.0,),
        ('Glass Fibre-106GSM',21140000000.0,21140000000.0,0.088,1280000000.0,1280000000.0,1280000000.0,0.00009,1794.0,347000000.0,173500000.0,347000000.0,173500000.0,64000000.0,64000000.0,),
        ('Carbon Fibre-200GSM',48900000000.0,48900000000.0,0.16,6500000000.0,6500000000.0,6500000000.0,0.00032,1364.0,464200000.0,232100000.0,464200000.0,232100000.0,125000000.0,125000000.0,),
        ('Carbon UD Tape',82100000000.0,32200000000.0,0.28,1310000000.0,1310000000.0,1310000000.0,0.0003,1553.0,1073400000.0,590370000.0,561000000.0,561000000.0,57000000.0,57000000.0,)]

 
    #create material properties function
    def Createcomposites(mat_name, E_1, E_2, v, G_12, G_13, G_23, thickness, density, long_ten_strength, long_comp_strength, tran_ten_strength, tran_comp_strength, long_shear_strength, tran_shear_strength):
        mdb.models['Model-1'].Material(name=mat_name)
        mdb.models['Model-1'].materials[mat_name].Elastic(table=((
            E_1, E_2, v, G_12, G_13, G_23), ), type=LAMINA)
        mdb.models['Model-1'].materials[mat_name].Density(table=((density, ), 
            ))
        mdb.models['Model-1'].materials[mat_name].HashinDamageInitiation(
            table=((long_ten_strength, long_comp_strength, tran_ten_strength, tran_comp_strength, long_shear_strength, tran_shear_strength), ))
        mdb.models['Model-1'].materials[mat_name].elastic.FailStress(table=(
            (long_ten_strength, long_comp_strength, tran_ten_strength, tran_comp_strength, long_shear_strength, -0.5, 0.0), 
            ))
    #
    for i in Composite_properties:
            Createcomposites(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], i[12], i[13], i[14])
    #######
    #######
        
    def assigncompositelayup1():
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='ribs', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['ribs'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['ribs'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName='Rib-Ply-1', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Ribs'], suppressed=False, 
            thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['ribs'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName='Rib-Ply-2', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Ribs'], suppressed=False, 
            thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['ribs'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName='Rib-Ply-3', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Ribs'], suppressed=False, 
            thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['ribs'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName='Rib-Ply-4', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Ribs'], suppressed=False, 
            thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['ribs'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=mdb.models['Model-1'].parts['Wing'].surfaces['Ribs'], 
            orientationType=DISCRETE, primaryAxisDefinition=VECTOR, 
            primaryAxisDirection=AXIS_2, primaryAxisVector=(1.0, 0.0, 0.0), 
            stackDirection=STACK_3)

        #base rib layup
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='base rib', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'base rib Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'base rib Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'base rib Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Airex C70 Foam', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='base rib Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.0015, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'base rib Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'base rib Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'base rib Ply-7', region=
            mdb.models['Model-1'].parts['Wing'].sets['Base Rib'], suppressed=False, 
            thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['base rib'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=mdb.models['Model-1'].parts['Wing'].surfaces['base rib'], 
            orientationType=DISCRETE, primaryAxisDefinition=VECTOR, 
            primaryAxisDirection=AXIS_1, primaryAxisVector=(1.0, 0.0, 0.0), 
            stackDirection=STACK_3)
            #rear spar AB
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.435, -0.034317, 
            2.340017), ), ((0.435, -0.01599, 2.248882), ), ((0.435, 0.004094, 
            2.149007), ), ((0.435, 0.01845, 2.07601), ), ((0.435, 0.021767, 2.03), ), (
            (0.435, 0.022348, 1.945), ), ((0.435, 0.023031, 1.845), ), ((0.435, 
            0.023715, 1.745), ), ((0.435, 0.024398, 1.645), ), ((0.435, 0.025081, 
            1.545), ), ((0.435, 0.025765, 1.445), ), ((0.435, 0.026448, 1.345), ), ((
            0.435, 0.027131, 1.245), ), ((0.435, 0.027815, 1.145), ), ((0.435, 
            0.028498, 1.045), ), ((0.435, 0.029182, 0.945), ), ((0.435, 0.029865, 
            0.845), ), ((0.435, 0.030548, 0.745), ), ((0.435, 0.031232, 0.645), ), ((
            0.435, 0.03171, 0.575), ), ((0.435, 0.032069, 0.5225), ), ((0.435, 
            0.032154, 0.4725), ), ((0.435, 0.032154, 0.415), ), ((0.435, 0.032154, 
            0.325), ), ((0.435, 0.032154, 0.225), ), ((0.435, 0.032154, 0.125), ), ((
            0.435, 0.032154, 0.025), ), ), name='Set-18')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='rear spar after boom', offsetType=MIDDLE_SURFACE, symmetric=
            False, thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar after boom'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'rear spar-AB Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-after boom'], 
            suppressed=False, thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'rear spar-AB Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-after boom'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'rear spar-AB Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-after boom'], 
            suppressed=False, thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar after boom'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['Skin web-after boom'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-18'], stackDirection=STACK_3)
            #rear spar-BB
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.435, -0.034317, 
            2.340017), ), ((0.435, -0.01599, 2.248882), ), ((0.435, 0.004094, 
            2.149007), ), ((0.435, 0.01845, 2.07601), ), ((0.435, 0.021767, 2.03), ), (
            (0.435, 0.022348, 1.945), ), ((0.435, 0.023031, 1.845), ), ((0.435, 
            0.023715, 1.745), ), ((0.435, 0.024398, 1.645), ), ((0.435, 0.025081, 
            1.545), ), ((0.435, 0.025765, 1.445), ), ((0.435, 0.026448, 1.345), ), ((
            0.435, 0.027131, 1.245), ), ((0.435, 0.027815, 1.145), ), ((0.435, 
            0.028498, 1.045), ), ((0.435, 0.029182, 0.945), ), ((0.435, 0.029865, 
            0.845), ), ((0.435, 0.030548, 0.745), ), ((0.435, 0.031232, 0.645), ), ((
            0.435, 0.03171, 0.575), ), ((0.435, 0.032069, 0.5225), ), ((0.435, 
            0.032154, 0.4725), ), ((0.435, 0.032154, 0.415), ), ((0.435, 0.032154, 
            0.325), ), ((0.435, 0.032154, 0.225), ), ((0.435, 0.032154, 0.125), ), ((
            0.435, 0.032154, 0.025), ), ), name='Set-19')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='rear spar before boom', offsetType=MIDDLE_SURFACE, symmetric=
            False, thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar before boom'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'rear spar-BB-Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-before boom'], 
            suppressed=False, thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'rear spar-BB-Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-before boom'], 
            suppressed=False, thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'rear spar-BB-Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-before boom'], 
            suppressed=False, thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-106GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'rear spar-BB-Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Rear Spar-before boom'], 
            suppressed=False, thickness=9e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['rear spar before boom'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['Rear Spar-before boom'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-19'], stackDirection=STACK_3)


        
          #spar web before boom  
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.331457, -0.063854, 
            2.334276), ), ((0.306457, -0.063854, 2.334276), ), ((0.281457, -0.063854, 
            2.334276), ), ((0.305517, -0.047401, 2.242776), ), ((0.280517, -0.047401, 
            2.242776), ), ((0.255517, -0.047401, 2.242776), ), ((0.27709, -0.02937, 
            2.142503), ), ((0.25209, -0.02937, 2.142503), ), ((0.22709, -0.02937, 
            2.142503), ), ((0.257191, -0.016945, 2.074359), ), ((0.232191, -0.016945, 
            2.074359), ), ((0.207191, -0.016945, 2.074359), ), ((0.25169, -0.014409, 
            2.03), ), ((0.22669, -0.014409, 2.03), ), ((0.20169, -0.014409, 2.03), ), (
            (0.248185, -0.014701, 1.945), ), ((0.223185, -0.014701, 1.945), ), ((
            0.198185, -0.014701, 1.945), ), ((0.244062, -0.015044, 1.845), ), ((
            0.219062, -0.015044, 1.845), ), ((0.194062, -0.015044, 1.845), ), ((
            0.239939, -0.015388, 1.745), ), ((0.214939, -0.015388, 1.745), ), ((
            0.189939, -0.015388, 1.745), ), ((0.235816, -0.015731, 1.645), ), ((
            0.210816, -0.015731, 1.645), ), ((0.185816, -0.015731, 1.645), ), ((
            0.231693, -0.016075, 1.545), ), ((0.206693, -0.016075, 1.545), ), ((
            0.181693, -0.016075, 1.545), ), ((0.22757, -0.016418, 1.445), ), ((0.20257, 
            -0.016418, 1.445), ), ((0.17757, -0.016418, 1.445), ), ((0.223447, 
            -0.016762, 1.345), ), ((0.198447, -0.016762, 1.345), ), ((0.173447, 
            -0.016762, 1.345), ), ((0.219324, -0.017105, 1.245), ), ((0.194324, 
            -0.017105, 1.245), ), ((0.169324, -0.017105, 1.245), ), ((0.215201, 
            -0.017449, 1.145), ), ((0.190201, -0.017449, 1.145), ), ((0.165201, 
            -0.017449, 1.145), ), ((0.211078, -0.017792, 1.045), ), ((0.186078, 
            -0.017792, 1.045), ), ((0.161078, -0.017792, 1.045), ), ((0.206955, 
            -0.018135, 0.945), ), ((0.181955, -0.018135, 0.945), ), ((0.156955, 
            -0.018135, 0.945), ), ((0.202832, -0.018479, 0.845), ), ((0.177832, 
            -0.018479, 0.845), ), ((0.152832, -0.018479, 0.845), ), ((0.198709, 
            -0.018822, 0.745), ), ((0.173709, -0.018822, 0.745), ), ((0.148709, 
            -0.018822, 0.745), ), ((0.150836, -0.019166, 0.645), ), ((0.144586, 
            -0.019166, 0.645), ), ((0.188336, -0.019166, 0.645), ), ((0.169586, 
            -0.019166, 0.645), ), ((0.194586, -0.019166, 0.645), ), ((0.1917, 
            -0.019406, 0.575), ), ((0.1667, -0.019406, 0.575), ), ((0.1417, -0.019406, 
            0.575), ), ((0.189535, -0.019586, 0.5225), ), ((0.164535, -0.019586, 
            0.5225), ), ((0.139535, -0.019586, 0.5225), ), ((0.18902, -0.019629, 
            0.4725), ), ((0.16402, -0.019629, 0.4725), ), ((0.13902, -0.019629, 
            0.4725), ), ((0.18902, -0.019629, 0.415), ), ((0.16402, -0.019629, 0.415), 
            ), ((0.13902, -0.019629, 0.415), ), ((0.18902, -0.019629, 0.325), ), ((
            0.16402, -0.019629, 0.325), ), ((0.13902, -0.019629, 0.325), ), ((0.18902, 
            -0.019629, 0.225), ), ((0.16402, -0.019629, 0.225), ), ((0.13902, 
            -0.019629, 0.225), ), ((0.18902, -0.019629, 0.125), ), ((0.16402, 
            -0.019629, 0.125), ), ((0.13902, -0.019629, 0.125), ), ((0.18902, 
            -0.019629, 0.025), ), ((0.16402, -0.019629, 0.025), ), ((0.13902, 
            -0.019629, 0.025), ), ), name='Set-23')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='spar web before boom', offsetType=MIDDLE_SURFACE, symmetric=
            False, thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-BB Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar web-before boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-BB Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar web-before boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Airex C70 Foam', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='spar web-BB Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar web-before boom'], 
            suppressed=False, thickness=0.0015, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-BB Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar web-before boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-BB Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar web-before boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web before boom'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['spar web before boom'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-23'], stackDirection=STACK_3)
            #spar web after boom
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.281457, -0.021362, 
            2.342535), ), ((0.306457, -0.021362, 2.342535), ), ((0.331457, -0.021362, 
            2.342535), ), ((0.255517, 0.000463, 2.25208), ), ((0.280517, 0.000463, 
            2.25208), ), ((0.305517, 0.000463, 2.25208), ), ((0.22709, 0.024379, 
            2.15295), ), ((0.25209, 0.024379, 2.15295), ), ((0.27709, 0.024379, 
            2.15295), ), ((0.207191, 0.041749, 2.077042), ), ((0.232191, 0.041749, 
            2.077042), ), ((0.257191, 0.041749, 2.077042), ), ((0.20169, 0.045874, 
            2.03), ), ((0.22669, 0.045874, 2.03), ), ((0.25169, 0.045874, 2.03), ), ((
            0.198185, 0.046803, 1.945), ), ((0.223185, 0.046803, 1.945), ), ((0.248185, 
            0.046803, 1.945), ), ((0.194062, 0.047896, 1.845), ), ((0.219062, 0.047896, 
            1.845), ), ((0.244062, 0.047896, 1.845), ), ((0.189939, 0.04899, 1.745), ), 
            ((0.214939, 0.04899, 1.745), ), ((0.239939, 0.04899, 1.745), ), ((0.185816, 
            0.050083, 1.645), ), ((0.210816, 0.050083, 1.645), ), ((0.235816, 0.050083, 
            1.645), ), ((0.181693, 0.051176, 1.545), ), ((0.206693, 0.051176, 1.545), 
            ), ((0.231693, 0.051176, 1.545), ), ((0.17757, 0.05227, 1.445), ), ((
            0.20257, 0.05227, 1.445), ), ((0.22757, 0.05227, 1.445), ), ((0.173447, 
            0.053363, 1.345), ), ((0.198447, 0.053363, 1.345), ), ((0.223447, 0.053363, 
            1.345), ), ((0.169324, 0.054457, 1.245), ), ((0.194324, 0.054457, 1.245), 
            ), ((0.219324, 0.054457, 1.245), ), ((0.165201, 0.05555, 1.145), ), ((
            0.190201, 0.05555, 1.145), ), ((0.215201, 0.05555, 1.145), ), ((0.161078, 
            0.056643, 1.045), ), ((0.186078, 0.056643, 1.045), ), ((0.211078, 0.056643, 
            1.045), ), ((0.156955, 0.057737, 0.945), ), ((0.181955, 0.057737, 0.945), 
            ), ((0.206955, 0.057737, 0.945), ), ((0.152832, 0.05883, 0.845), ), ((
            0.177832, 0.05883, 0.845), ), ((0.202832, 0.05883, 0.845), ), ((0.148709, 
            0.059923, 0.745), ), ((0.173709, 0.059923, 0.745), ), ((0.198709, 0.059923, 
            0.745), ), ((0.150836, 0.061017, 0.645), ), ((0.169586, 0.061017, 0.645), 
            ), ((0.188336, 0.061017, 0.645), ), ((0.194586, 0.061017, 0.645), ), ((
            0.144586, 0.061017, 0.645), ), ((0.1417, 0.061782, 0.575), ), ((0.1667, 
            0.061782, 0.575), ), ((0.1917, 0.061782, 0.575), ), ((0.139535, 0.062356, 
            0.5225), ), ((0.164535, 0.062356, 0.5225), ), ((0.189535, 0.062356, 
            0.5225), ), ((0.13902, 0.062493, 0.4725), ), ((0.16402, 0.062493, 0.4725), 
            ), ((0.18902, 0.062493, 0.4725), ), ((0.13902, 0.062493, 0.415), ), ((
            0.16402, 0.062493, 0.415), ), ((0.18902, 0.062493, 0.415), ), ((0.13902, 
            0.062493, 0.325), ), ((0.16402, 0.062493, 0.325), ), ((0.18902, 0.062493, 
            0.325), ), ((0.13902, 0.062493, 0.225), ), ((0.16402, 0.062493, 0.225), ), 
            ((0.18902, 0.062493, 0.225), ), ((0.13902, 0.062493, 0.125), ), ((0.16402, 
            0.062493, 0.125), ), ((0.18902, 0.062493, 0.125), ), ((0.13902, 0.062493, 
            0.025), ), ((0.16402, 0.062493, 0.025), ), ((0.18902, 0.062493, 0.025), ), 
            ), name='Set-24')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='spar web after boom', offsetType=MIDDLE_SURFACE, symmetric=
            False, thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web after boom'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-AB Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar Web-after boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-AB Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar Web-after boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-AB Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar Web-after boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web after boom'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar web-AB Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar Web-after boom'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar web after boom'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['Skin web-after boom'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-24'], stackDirection=STACK_3)


        #skin mid
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.562195, 0.0, 1.345), 
            ), ((0.546396, 0.000821, 1.345), ), ((0.534441, 0.001088, 1.345), ), ((
            0.519462, 0.001159, 1.345), ), ((0.501728, 0.000952, 1.345), ), ((0.481536, 
            0.000438, 1.345), ), ((0.459202, -0.000398, 1.345), ), ((0.435, -0.001752, 
            1.345), ), ((0.382855, -0.00467, 1.345), ), ((0.355509, -0.006524, 1.345), 
            ), ((0.327809, -0.008464, 1.345), ), ((0.30011, -0.010399, 1.345), ), ((
            0.272743, -0.013351, 1.345), ), ((0.245096, -0.016762, 1.345), ), ((
            0.223447, -0.016762, 1.345), ), ((0.173447, -0.016762, 1.345), ), ((
            0.151798, -0.016762, 1.345), ), ((0.130574, -0.016721, 1.345), ), ((
            0.11279, -0.016016, 1.345), ), ((0.097232, -0.014852, 1.345), ), ((
            0.084078, -0.013215, 1.345), ), ((0.073458, -0.011109, 1.345), ), ((
            0.065513, -0.00857, 1.345), ), ((0.06053, -0.005617, 1.345), ), ((0.058389, 
            0.0, 1.345), ), ((0.058908, 0.003582, 1.345), ), ((0.06122, 0.009426, 
            1.345), ), ((0.065493, 0.015522, 1.345), ), ((0.07174, 0.021664, 1.345), ), 
            ((0.079911, 0.027664, 1.345), ), ((0.089973, 0.033332, 1.345), ), ((
            0.101903, 0.038526, 1.345), ), ((0.115636, 0.043136, 1.345), ), ((0.131098, 
            0.047862, 1.345), ), ((0.151798, 0.053363, 1.345), ), ((0.173447, 0.053363, 
            1.345), ), ((0.223447, 0.053363, 1.345), ), ((0.245096, 0.053363, 1.345), 
            ), ((0.277746, 0.051534, 1.345), ), ((0.302241, 0.048985, 1.345), ), ((
            0.327029, 0.045791, 1.345), ), ((0.351876, 0.042053, 1.345), ), ((0.376563, 
            0.037871, 1.345), ), ((0.435, 0.026448, 1.345), ), ((0.447156, 0.024072, 
            1.345), ), ((0.468623, 0.019573, 1.345), ), ((0.488544, 0.015361, 1.345), 
            ), ((0.50661, 0.011532, 1.345), ), ((0.522531, 0.008136, 1.345), ), ((
            0.536088, 0.005219, 1.345), ), ((0.547101, 0.002877, 1.345), ), ((0.051396, 
            0.0, 1.245), ), ((0.051926, 0.003655, 1.245), ), ((0.054286, 0.009619, 
            1.245), ), ((0.058645, 0.01584, 1.245), ), ((0.065021, 0.022108, 1.245), ), 
            ((0.07336, 0.028231, 1.245), ), ((0.083627, 0.034015, 1.245), ), ((
            0.095802, 0.039315, 1.245), ), ((0.109817, 0.04402, 1.245), ), ((0.125595, 
            0.048842, 1.245), ), ((0.14672, 0.054457, 1.245), ), ((0.169324, 0.054457, 
            1.245), ), ((0.219324, 0.054457, 1.245), ), ((0.241929, 0.054457, 1.245), 
            ), ((0.275248, 0.05259, 1.245), ), ((0.300245, 0.049989, 1.245), ), ((
            0.32554, 0.046729, 1.245), ), ((0.350897, 0.042914, 1.245), ), ((0.376089, 
            0.038647, 1.245), ), ((0.435, 0.027131, 1.245), ), ((0.448129, 0.024565, 
            1.245), ), ((0.470036, 0.019974, 1.245), ), ((0.490365, 0.015676, 1.245), 
            ), ((0.508801, 0.011768, 1.245), ), ((0.525048, 0.008303, 1.245), ), ((
            0.538883, 0.005326, 1.245), ), ((0.550122, 0.002936, 1.245), ), ((0.565525, 
            0.0, 1.245), ), ((0.549402, 0.000838, 1.245), ), ((0.537202, 0.001111, 
            1.245), ), ((0.521917, 0.001182, 1.245), ), ((0.503819, 0.000972, 1.245), 
            ), ((0.483213, 0.000447, 1.245), ), ((0.460422, -0.000406, 1.245), ), ((
            0.435, -0.001829, 1.245), ), ((0.382511, -0.004766, 1.245), ), ((0.354604, 
            -0.006658, 1.245), ), ((0.326337, -0.008637, 1.245), ), ((0.29807, 
            -0.010612, 1.245), ), ((0.270143, -0.013624, 1.245), ), ((0.241929, 
            -0.017105, 1.245), ), ((0.219324, -0.017105, 1.245), ), ((0.169324, 
            -0.017105, 1.245), ), ((0.14672, -0.017105, 1.245), ), ((0.125061, 
            -0.017064, 1.245), ), ((0.106912, -0.016344, 1.245), ), ((0.091036, 
            -0.015157, 1.245), ), ((0.077612, -0.013486, 1.245), ), ((0.066774, 
            -0.011337, 1.245), ), ((0.058666, -0.008745, 1.245), ), ((0.053581, 
            -0.005733, 1.245), ), ((0.568855, 0.0, 1.145), ), ((0.552408, 0.000855, 
            1.145), ), ((0.539963, 0.001133, 1.145), ), ((0.524371, 0.001206, 1.145), 
            ), ((0.50591, 0.000991, 1.145), ), ((0.48489, 0.000456, 1.145), ), ((
            0.461641, -0.000414, 1.145), ), ((0.435, -0.001905, 1.145), ), ((0.382166, 
            -0.004862, 1.145), ), ((0.353699, -0.006792, 1.145), ), ((0.324865, 
            -0.008811, 1.145), ), ((0.29603, -0.010825, 1.145), ), ((0.267542, 
            -0.013898, 1.145), ), ((0.238761, -0.017449, 1.145), ), ((0.215201, 
            -0.017449, 1.145), ), ((0.165201, -0.017449, 1.145), ), ((0.141641, 
            -0.017449, 1.145), ), ((0.119547, -0.017407, 1.145), ), ((0.101034, 
            -0.016672, 1.145), ), ((0.084839, -0.015461, 1.145), ), ((0.071145, 
            -0.013756, 1.145), ), ((0.06009, -0.011564, 1.145), ), ((0.051819, 
            -0.008921, 1.145), ), ((0.046632, -0.005848, 1.145), ), ((0.044404, 0.0, 
            1.145), ), ((0.044944, 0.003729, 1.145), ), ((0.047351, 0.009812, 1.145), 
            ), ((0.051798, 0.016158, 1.145), ), ((0.058301, 0.022551, 1.145), ), ((
            0.066808, 0.028798, 1.145), ), ((0.077281, 0.034698, 1.145), ), ((0.0897, 
            0.040105, 1.145), ), ((0.103997, 0.044904, 1.145), ), ((0.120092, 0.049823, 
            1.145), ), ((0.141641, 0.05555, 1.145), ), ((0.165201, 0.05555, 1.145), ), 
            ((0.215201, 0.05555, 1.145), ), ((0.238761, 0.05555, 1.145), ), ((0.27275, 
            0.053646, 1.145), ), ((0.298249, 0.050992, 1.145), ), ((0.324052, 0.047667, 
            1.145), ), ((0.349918, 0.043776, 1.145), ), ((0.375616, 0.039423, 1.145), 
            ), ((0.435, 0.027815, 1.145), ), ((0.449102, 0.025058, 1.145), ), ((
            0.471449, 0.020375, 1.145), ), ((0.492186, 0.015991, 1.145), ), ((0.510992, 
            0.012005, 1.145), ), ((0.527565, 0.00847, 1.145), ), ((0.541678, 0.005433, 
            1.145), ), ((0.553143, 0.002995, 1.145), ), ((0.037411, 0.0, 1.045), ), ((
            0.037962, 0.003802, 1.045), ), ((0.040416, 0.010006, 1.045), ), ((0.044951, 
            0.016476, 1.045), ), ((0.051582, 0.022995, 1.045), ), ((0.060256, 0.029364, 
            1.045), ), ((0.070936, 0.035381, 1.045), ), ((0.083599, 0.040894, 1.045), 
            ), ((0.098177, 0.045787, 1.045), ), ((0.114589, 0.050804, 1.045), ), ((
            0.136562, 0.056643, 1.045), ), ((0.161078, 0.056643, 1.045), ), ((0.211078, 
            0.056643, 1.045), ), ((0.235594, 0.056643, 1.045), ), ((0.270252, 0.054702, 
            1.045), ), ((0.296252, 0.051996, 1.045), ), ((0.322563, 0.048606, 1.045), 
            ), ((0.348938, 0.044638, 1.045), ), ((0.375142, 0.040199, 1.045), ), ((
            0.435, 0.028498, 1.045), ), ((0.450075, 0.025552, 1.045), ), ((0.472861, 
            0.020776, 1.045), ), ((0.494006, 0.016305, 1.045), ), ((0.513183, 0.012241, 
            1.045), ), ((0.530082, 0.008637, 1.045), ), ((0.544473, 0.00554, 1.045), ), 
            ((0.556163, 0.003054, 1.045), ), ((0.572185, 0.0, 1.045), ), ((0.555415, 
            0.000872, 1.045), ), ((0.542724, 0.001155, 1.045), ), ((0.526825, 0.00123, 
            1.045), ), ((0.508001, 0.001011, 1.045), ), ((0.486568, 0.000465, 1.045), 
            ), ((0.462861, -0.000422, 1.045), ), ((0.435, -0.001982, 1.045), ), ((
            0.381821, -0.004957, 1.045), ), ((0.352794, -0.006925, 1.045), ), ((
            0.323392, -0.008984, 1.045), ), ((0.29399, -0.011038, 1.045), ), ((
            0.264941, -0.014172, 1.045), ), ((0.235594, -0.017792, 1.045), ), ((
            0.211078, -0.017792, 1.045), ), ((0.161078, -0.017792, 1.045), ), ((
            0.136562, -0.017792, 1.045), ), ((0.114033, -0.017749, 1.045), ), ((
            0.095156, -0.017, 1.045), ), ((0.078642, -0.015765, 1.045), ), ((0.064679, 
            -0.014027, 1.045), ), ((0.053406, -0.011792, 1.045), ), ((0.044973, 
            -0.009097, 1.045), ), ((0.039684, -0.005963, 1.045), ), ((0.575515, 0.0, 
            0.945), ), ((0.558421, 0.000889, 0.945), ), ((0.545486, 0.001177, 0.945), 
            ), ((0.52928, 0.001254, 0.945), ), ((0.510092, 0.00103, 0.945), ), ((
            0.488245, 0.000474, 0.945), ), ((0.464081, -0.000431, 0.945), ), ((0.435, 
            -0.002058, 0.945), ), ((0.381477, -0.005053, 0.945), ), ((0.351889, 
            -0.007059, 0.945), ), ((0.32192, -0.009158, 0.945), ), ((0.29195, 
            -0.011251, 0.945), ), ((0.26234, -0.014445, 0.945), ), ((0.232427, 
            -0.018135, 0.945), ), ((0.206955, -0.018135, 0.945), ), ((0.156955, 
            -0.018135, 0.945), ), ((0.131483, -0.018135, 0.945), ), ((0.10852, 
            -0.018092, 0.945), ), ((0.089278, -0.017329, 0.945), ), ((0.072445, 
            -0.016069, 0.945), ), ((0.058213, -0.014298, 0.945), ), ((0.046722, 
            -0.012019, 0.945), ), ((0.038126, -0.009272, 0.945), ), ((0.032735, 
            -0.006078, 0.945), ), ((0.030418, 0.0, 0.945), ), ((0.03098, 0.003876, 
            0.945), ), ((0.033482, 0.010199, 0.945), ), ((0.038104, 0.016794, 0.945), 
            ), ((0.044863, 0.023439, 0.945), ), ((0.053705, 0.029931, 0.945), ), ((
            0.06459, 0.036064, 0.945), ), ((0.077498, 0.041684, 0.945), ), ((0.092358, 
            0.046671, 0.945), ), ((0.109087, 0.051784, 0.945), ), ((0.131483, 0.057737, 
            0.945), ), ((0.156955, 0.057737, 0.945), ), ((0.206955, 0.057737, 0.945), 
            ), ((0.232427, 0.057737, 0.945), ), ((0.267753, 0.055758, 0.945), ), ((
            0.294256, 0.053, 0.945), ), ((0.321075, 0.049544, 0.945), ), ((0.347959, 
            0.045499, 0.945), ), ((0.374669, 0.040975, 0.945), ), ((0.435, 0.029182, 
            0.945), ), ((0.451048, 0.026045, 0.945), ), ((0.474274, 0.021177, 0.945), 
            ), ((0.495827, 0.01662, 0.945), ), ((0.515374, 0.012477, 0.945), ), ((
            0.532599, 0.008803, 0.945), ), ((0.547268, 0.005647, 0.945), ), )+\
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.559184, 0.003113, 
            0.945), ), ((0.023425, 0.0, 0.845), ), ((0.023998, 0.003949, 0.845), ), ((
            0.026547, 0.010392, 0.845), ), ((0.031257, 0.017112, 0.845), ), ((0.038144, 
            0.023883, 0.845), ), ((0.047153, 0.030498, 0.845), ), ((0.058245, 0.036747, 
            0.845), ), ((0.071397, 0.042473, 0.845), ), ((0.086538, 0.047555, 0.845), 
            ), ((0.103584, 0.052765, 0.845), ), ((0.126404, 0.05883, 0.845), ), ((
            0.152832, 0.05883, 0.845), ), ((0.202832, 0.05883, 0.845), ), ((0.22926, 
            0.05883, 0.845), ), ((0.265255, 0.056814, 0.845), ), ((0.29226, 0.054003, 
            0.845), ), ((0.319586, 0.050482, 0.845), ), ((0.346979, 0.046361, 0.845), 
            ), ((0.374195, 0.041751, 0.845), ), ((0.435, 0.029865, 0.845), ), ((
            0.45202, 0.026538, 0.845), ), ((0.475687, 0.021578, 0.845), ), ((0.497648, 
            0.016935, 0.845), ), ((0.517565, 0.012714, 0.845), ), ((0.535117, 0.00897, 
            0.845), ), ((0.550063, 0.005754, 0.845), ), ((0.562204, 0.003171, 0.845), 
            ), ((0.578845, 0.0, 0.845), ), ((0.561427, 0.000905, 0.845), ), ((0.548247, 
            0.0012, 0.845), ), ((0.531734, 0.001277, 0.845), ), ((0.512183, 0.00105, 
            0.845), ), ((0.489922, 0.000483, 0.845), ), ((0.4653, -0.000439, 0.845), ), 
            ((0.435, -0.002134, 0.845), ), ((0.381132, -0.005149, 0.845), ), ((
            0.350984, -0.007193, 0.845), ), ((0.320447, -0.009331, 0.845), ), ((
            0.28991, -0.011464, 0.845), ), ((0.25974, -0.014719, 0.845), ), ((0.22926, 
            -0.018479, 0.845), ), ((0.202832, -0.018479, 0.845), ), ((0.152832, 
            -0.018479, 0.845), ), ((0.126404, -0.018479, 0.845), ), ((0.103006, 
            -0.018434, 0.845), ), ((0.0834, -0.017657, 0.845), ), ((0.066248, 
            -0.016374, 0.845), ), ((0.051746, -0.014569, 0.845), ), ((0.040038, 
            -0.012247, 0.845), ), ((0.031279, -0.009448, 0.845), ), ((0.025786, 
            -0.006193, 0.845), ), ((0.582175, 0.0, 0.745), ), ((0.564433, 0.000922, 
            0.745), ), ((0.551008, 0.001222, 0.745), ), ((0.534189, 0.001301, 0.745), 
            ), ((0.514274, 0.001069, 0.745), ), ((0.491599, 0.000492, 0.745), ), ((
            0.46652, -0.000447, 0.745), ), ((0.435, -0.002211, 0.745), ), ((0.380788, 
            -0.005244, 0.745), ), ((0.350079, -0.007326, 0.745), ), ((0.318975, 
            -0.009504, 0.745), ), ((0.28787, -0.011677, 0.745), ), ((0.257139, 
            -0.014992, 0.745), ), ((0.226093, -0.018822, 0.745), ), ((0.198709, 
            -0.018822, 0.745), ), ((0.148709, -0.018822, 0.745), ), ((0.121326, 
            -0.018822, 0.745), ), ((0.097492, -0.018777, 0.745), ), ((0.077522, 
            -0.017985, 0.745), ), ((0.060052, -0.016678, 0.745), ), ((0.04528, 
            -0.014839, 0.745), ), ((0.033354, -0.012475, 0.745), ), ((0.024432, 
            -0.009623, 0.745), ), ((0.018837, -0.006308, 0.745), ), ((0.016433, 0.0, 
            0.745), ), ((0.017016, 0.004022, 0.745), ), ((0.019612, 0.010585, 0.745), 
            ), ((0.02441, 0.017431, 0.745), ), ((0.031425, 0.024327, 0.745), ), ((
            0.040601, 0.031065, 0.745), ), ((0.051899, 0.037429, 0.745), ), ((0.065296, 
            0.043262, 0.745), ), ((0.080718, 0.048439, 0.745), ), ((0.098081, 0.053745, 
            0.745), ), ((0.121326, 0.059923, 0.745), ), ((0.148709, 0.059923, 0.745), 
            ), ((0.198709, 0.059923, 0.745), ), ((0.226093, 0.059923, 0.745), ), ((
            0.262757, 0.05787, 0.745), ), ((0.290263, 0.055007, 0.745), ), ((0.318098, 
            0.05142, 0.745), ), ((0.346, 0.047222, 0.745), ), ((0.373721, 0.042527, 
            0.745), ), ((0.435, 0.030548, 0.745), ), ((0.452993, 0.027031, 0.745), ), (
            (0.477099, 0.021979, 0.745), ), ((0.499469, 0.017249, 0.745), ), ((
            0.519756, 0.01295, 0.745), ), ((0.537634, 0.009137, 0.745), ), ((0.552858, 
            0.005861, 0.745), ), ((0.565225, 0.00323, 0.745), ), ((0.194586, 0.061017, 
            0.645), ), ((0.144586, -0.019166, 0.645), ), ((0.585505, 0.0, 0.645), ), ((
            0.567439, 0.000939, 0.645), ), ((0.553769, 0.001244, 0.645), ), ((0.536643, 
            0.001325, 0.645), ), ((0.516365, 0.001089, 0.645), ), ((0.493277, 0.000501, 
            0.645), ), ((0.46774, -0.000455, 0.645), ), ((0.435, -0.002287, 0.645), ), 
            ((0.380443, -0.00534, 0.645), ), ((0.349174, -0.00746, 0.645), ), ((
            0.317502, -0.009678, 0.645), ), ((0.28583, -0.01189, 0.645), ), ((0.254538, 
            -0.015266, 0.645), ), ((0.222925, -0.019166, 0.645), ), ((0.194586, 
            -0.019166, 0.645), ), ((0.116247, -0.019166, 0.645), ), ((0.091979, 
            -0.01912, 0.645), ), ((0.071644, -0.018313, 0.645), ), ((0.053855, 
            -0.016982, 0.645), ), ((0.038814, -0.01511, 0.645), ), ((0.02667, 
            -0.012702, 0.645), ), ((0.017586, -0.009799, 0.645), ), ((0.011888, 
            -0.006423, 0.645), ), ((0.00944, 0.0, 0.645), ), ((0.010033, 0.004096, 
            0.645), ), ((0.012678, 0.010778, 0.645), ), ((0.017563, 0.017749, 0.645), 
            ), ((0.024706, 0.024771, 0.645), ), ((0.03405, 0.031632, 0.645), ), ((
            0.045554, 0.038112, 0.645), ), ((0.059195, 0.044052, 0.645), ), ((0.074898, 
            0.049323, 0.645), ), ((0.092578, 0.054726, 0.645), ), ((0.116247, 0.061017, 
            0.645), ), ((0.144586, 0.061017, 0.645), ), ((0.222925, 0.061017, 0.645), 
            ), ((0.260259, 0.058926, 0.645), ), ((0.288267, 0.056011, 0.645), ), ((
            0.316609, 0.052359, 0.645), ), ((0.345021, 0.048084, 0.645), ), ((0.373248, 
            0.043303, 0.645), ), ((0.435, 0.031232, 0.645), ), ((0.453966, 0.027524, 
            0.645), ), ((0.478512, 0.02238, 0.645), ), ((0.50129, 0.017564, 0.645), ), 
            ((0.521947, 0.013186, 0.645), ), ((0.540151, 0.009303, 0.645), ), ((
            0.555653, 0.005968, 0.645), ), ((0.568246, 0.003289, 0.645), ), ), name=
            'Set-20')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='skin mid', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin mid Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-mid'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'skin mid Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-mid'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin mid Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-mid'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'skin mid Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-mid'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin mid Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-mid'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin mid'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=mdb.models['Model-1'].parts['Wing'].surfaces['skin-mid'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-20'], stackDirection=STACK_3)
            
            #skin tip
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.219999, -0.053698, 
            2.33625), ), ((0.220319, -0.051527, 2.336672), ), ((0.221747, -0.047986, 
            2.33736), ), ((0.224384, -0.044292, 2.338078), ), ((0.22824, -0.04057, 
            2.338802), ), ((0.233285, -0.036934, 2.339508), ), ((0.239495, -0.0335, 
            2.340176), ), ((0.24686, -0.030352, 2.340788), ), ((0.255338, -0.027559, 
            2.341331), ), ((0.264882, -0.024695, 2.341887), ), ((0.277661, -0.021362, 
            2.342535), ), ((0.281457, -0.021362, 2.342535), ), ((0.331457, -0.021362, 
            2.342535), ), ((0.335253, -0.021362, 2.342535), ), ((0.355408, -0.02247, 
            2.34232), ), ((0.370529, -0.024015, 2.34202), ), ((0.38583, -0.02595, 
            2.341644), ), ((0.401169, -0.028215, 2.341203), ), ((0.416408, -0.030749, 
            2.340711), ), ((0.435, -0.034317, 2.340017), ), ((0.459985, -0.039111, 
            2.339085), ), ((0.473237, -0.041837, 2.338555), ), ((0.485534, -0.044389, 
            2.338059), ), ((0.496686, -0.04671, 2.337608), ), ((0.506514, -0.048767, 
            2.337208), ), ((0.514883, -0.050535, 2.336865), ), ((0.521681, -0.051954, 
            2.336589), ), ((0.530999, -0.053698, 2.33625), ), ((0.521246, -0.0532, 
            2.336347), ), ((0.513866, -0.053038, 2.336378), ), ((0.50462, -0.052995, 
            2.336386), ), ((0.493673, -0.053121, 2.336362), ), ((0.481208, -0.053432, 
            2.336302), ), ((0.467421, -0.053939, 2.336203), ), ((0.435, -0.05572, 
            2.335857), ), ((0.420292, -0.056528, 2.3357), ), ((0.403411, -0.057651, 
            2.335482), ), ((0.386312, -0.058826, 2.335253), ), ((0.369213, -0.059999, 
            2.335025), ), ((0.35232, -0.061788, 2.334677), ), ((0.335253, -0.063854, 
            2.334276), ), ((0.331457, -0.063854, 2.334276), ), ((0.281457, -0.063854, 
            2.334276), ), ((0.277661, -0.063854, 2.334276), ), ((0.264559, -0.06383, 
            2.33428), ), ((0.253581, -0.063403, 2.334364), ), ((0.243977, -0.062697, 
            2.334501), ), ((0.235857, -0.061705, 2.334693), ), ((0.229301, -0.060429, 
            2.334942), ), ((0.224396, -0.05889, 2.335241), ), ((0.221321, -0.057101, 
            2.335588), ), ((0.533439, -0.03596, 2.245), ), ((0.522453, -0.0354, 
            2.245109), ), ((0.514141, -0.035218, 2.245144), ), ((0.503726, -0.035169, 
            2.245154), ), ((0.491395, -0.03531, 2.245126), ), ((0.477355, -0.035661, 
            2.245058), ), ((0.461826, -0.036232, 2.244947), ), ((0.435, -0.037706, 
            2.244661), ), ((0.40874, -0.039148, 2.24438), ), ((0.389725, -0.040413, 
            2.244134), ), ((0.370465, -0.041737, 2.243877), ), ((0.351206, -0.043058, 
            2.24362), ), ((0.332177, -0.045073, 2.243229), ), ((0.312953, -0.047401, 
            2.242776), ), ((0.305517, -0.047401, 2.242776), ), ((0.255517, -0.047401, 
            2.242776), ), ((0.248081, -0.047401, 2.242776), ), ((0.233323, -0.047373, 
            2.242782), ), ((0.220958, -0.046892, 2.242875), ), ((0.21014, -0.046098, 
            2.24303), ), ((0.200994, -0.04498, 2.243247), ), ((0.193609, -0.043543, 
            2.243526), ), ((0.188085, -0.04181, 2.243863), ), ((0.18462, -0.039795, 
            2.244255), ), ((0.183131, -0.03596, 2.245), ), ((0.183492, -0.033515, 
            2.245475), ), ((0.1851, -0.029527, 2.246251), ), ((0.188071, -0.025366, 
            2.247059), ), ((0.192415, -0.021174, 2.247874), ), ((0.198097, -0.017078, 
            2.24867), ), ((0.205092, -0.01321, 2.249422), ), ((0.213387, -0.009664, 
            2.250111), ), ((0.222937, -0.006518, 2.250723), ), ((0.233688, -0.003293, 
            2.25135), ), ((0.248081, 0.000463, 2.25208), ), ((0.255517, 0.000463, 
            2.25208), ), ((0.305517, 0.000463, 2.25208), ), ((0.312953, 0.000463, 
            2.25208), ), ((0.335655, -0.000786, 2.251837), ), ((0.352687, -0.002526, 
            2.251499), ), ((0.369922, -0.004706, 2.251075), ), ((0.3872, -0.007257, 
            2.250579), ), ((0.404365, -0.010112, 2.250025), ), ((0.435, -0.01599, 
            2.248882), ), ((0.45345, -0.01953, 2.248194), ), ((0.468376, -0.022601, 
            2.247597), ), ((0.482228, -0.025476, 2.247038), ), ((0.49479, -0.028089, 
            2.24653), ), ((0.505859, -0.030407, 2.246079), ), ((0.515286, -0.032398, 
            2.245692), ), ((0.522944, -0.033997, 2.245382), ), ((0.142729, -0.016522, 
            2.145), ), ((0.143134, -0.013777, 2.145534), ), ((0.14494, -0.009297, 
            2.146404), ), ((0.148276, -0.004625, 2.147313), ), ((0.153153, 8.2e-05, 
            2.148228), ), ((0.159534, 0.004682, 2.149122), ), ((0.16739, 0.009026, 
            2.149966), ), ((0.176705, 0.013007, 2.15074), ), ((0.187429, 0.01654, 
            2.151427), ), ((0.199502, 0.020163, 2.152131), ), ((0.215665, 0.024379, 
            2.15295), ), ((0.22709, 0.024379, 2.15295), ), ((0.27709, 0.024379, 
            2.15295), ), ((0.288514, 0.024379, 2.15295), ), ((0.314008, 0.022978, 
            2.152678), ), ((0.333135, 0.021024, 2.152298), ), ((0.352489, 0.018575, 
            2.151822), ), ((0.371891, 0.01571, 2.151265), ), ((0.391167, 0.012505, 
            2.150642), ), ((0.435, 0.004094, 2.149007), ), ((0.446288, 0.001928, 
            2.148586), ), ((0.46305, -0.00152, 2.147916), ), ((0.478605, -0.004748, 
            2.147289), ), ((0.492711, -0.007683, 2.146718), ), ((0.505142, -0.010286, 
            2.146212), ), ((0.515728, -0.012522, 2.145778), ), ((0.524328, -0.014317, 
            2.145429), ), ((0.536113, -0.016522, 2.145), ), ((0.523777, -0.015893, 
            2.145122), ), ((0.514442, -0.015688, 2.145162), ), ((0.502747, -0.015634, 
            2.145173), ), ((0.488899, -0.015792, 2.145142), ), ((0.473133, -0.016186, 
            2.145065), ), ((0.455694, -0.016827, 2.144941), ), ((0.435, -0.017964, 
            2.14472), ), ((0.39608, -0.020102, 2.144304), ), ((0.374727, -0.021523, 
            2.144028), ), ((0.353099, -0.02301, 2.143739), ), ((0.331471, -0.024493, 
            2.143451), ), ((0.310102, -0.026755, 2.143011), ), ((0.288514, -0.02937, 
            2.142503), ), ((0.27709, -0.02937, 2.142503), ), ((0.22709, -0.02937, 
            2.142503), ), ((0.215665, -0.02937, 2.142503), ), ((0.199093, -0.029339, 
            2.142509), ), ((0.185206, -0.028798, 2.142614), ), ((0.173059, -0.027906, 
            2.142787), ), ((0.162787, -0.026651, 2.143031), ), ((0.154495, -0.025037, 
            2.143345), ), ((0.148291, -0.023091, 2.143723), ), ((0.144401, -0.020828, 
            2.144163), ), ((0.537985, -0.002916, 2.075), ), ((0.524703, -0.002228, 
            2.075031), ), ((0.514653, -0.002005, 2.075042), ), ((0.502061, -0.001946, 
            2.075044), ), ((0.487152, -0.002119, 2.075036), ), ((0.470177, -0.002549, 
            2.075017), ), ((0.451401, -0.003249, 2.074985), ), ((0.435, -0.004162, 
            2.074934), ), ((0.387218, -0.006825, 2.074821), ), ((0.364229, -0.008376, 
            2.07475), ), ((0.340943, -0.01, 2.074676), ), ((0.317656, -0.011619, 
            2.074602), ), ((0.29465, -0.01409, 2.074489), ), ((0.271407, -0.016945, 
            2.074359), ), ((0.257191, -0.016945, 2.074359), ), ((0.207191, -0.016945, 
            2.074359), ), ((0.192974, -0.016945, 2.074359), ), ((0.175132, -0.016911, 
            2.07436), ), ((0.160181, -0.016321, 2.074387), ), ((0.147102, -0.015347, 
            2.074432), ), ((0.136043, -0.013976, 2.074494), ), ((0.127115, -0.012214, 
            2.074575), ), ((0.120436, -0.010089, 2.074672), ), ((0.116247, -0.007617, 
            2.074785), ), ((0.114447, -0.002916, 2.075), ), ((0.114883, 8.2e-05, 
            2.075137), ), ((0.116827, 0.004974, 2.075361), ), ((0.120419, 0.010076, 
            2.075594), ), ((0.125671, 0.015217, 2.075829), ), ((0.132541, 0.020239, 
            2.076059), ), ((0.140999, 0.024983, 2.076276), ), ((0.151028, 0.02933, 
            2.076474), ), ((0.162574, 0.033189, 2.076651), ), ((0.175572, 0.037144, 
            2.076832), ), ((0.192974, 0.041749, 2.077042), ), ((0.207191, 0.041749, 
            2.077042), ), ((0.257191, 0.041749, 2.077042), ), ((0.271407, 0.041749, 
            2.077042), ), ((0.298856, 0.040218, 2.076972), ), ((0.319448, 0.038084, 
            2.076875), ), ((0.340286, 0.035411, 2.076752), ), ((0.361175, 0.032282, 
            2.076609), ), ((0.381928, 0.028782, 2.076449), ), ((0.435, 0.01845, 
            2.07601), ), ((0.441275, 0.017232, 2.075921), ), ((0.459322, 0.013467, 
            2.075749), ), ((0.476068, 0.009941, 2.075588), ), ((0.491256, 0.006737, 
            2.075441), ), ((0.50464, 0.003894, 2.075311), ), ((0.516038, 0.001453, 
            2.0752), ), ((0.525296, -0.000508, 2.07511), ), ((0.106289, 0.0, 2.03), ), 
            ((0.106735, 0.003079, 2.03), ), ((0.108723, 0.008103, 2.03), ), ((0.112395, 
            0.013344, 2.03), ), ((0.117766, 0.018623, 2.03), ), ((0.124791, 0.023781, 
            2.03), ), ((0.13344, 0.028654, 2.03), ), ((0.143695, 0.033119, 2.03), ), ((
            0.155502, 0.037082, 2.03), ), ((0.168793, 0.041144, 2.03), ), ((0.186588, 
            0.045874, 2.03), ), ((0.20169, 0.045874, 2.03), ), ((0.25169, 0.045874, 
            2.03), ), ((0.266791, 0.045874, 2.03), ), ((0.294859, 0.044301, 2.03), ), (
            (0.315916, 0.04211, 2.03), ), ((0.337225, 0.039364, 2.03), ), ((0.358585, 
            0.036151, 2.03), ), ((0.379807, 0.032556, 2.03), ), ((0.435, 0.021767, 
            2.03), ), ((0.440492, 0.020693, 2.03), ), ((0.458946, 0.016826, 2.03), ), (
            (0.476071, 0.013205, 2.03), ), ((0.491602, 0.009914, 2.03), ), ((0.505288, 
            0.006995, 2.03), ), ((0.516942, 0.004487, 2.03), ), ((0.52641, 0.002473, 
            2.03), ), ((0.539386, 0.0, 2.03), ), ((0.525804, 0.000706, 2.03), ), ((
            0.515526, 0.000935, 2.03), ), ((0.50265, 0.000996, 2.03), ), ((0.487405, 
            0.000819, 2.03), ), ((0.470047, 0.000377, 2.03), ), ((0.450848, -0.000342, 
            2.03), ), ((0.435, -0.001229, 2.03), ), ((0.385216, -0.004015, 2.03), ), ((
            0.361708, -0.005609, 2.03), ), ((0.337896, -0.007276, 2.03), ), ((0.314084, 
            -0.008939, 2.03), ), ((0.290558, -0.011477, 2.03), ), ((0.266791, 
            -0.014409, 2.03), ), ((0.25169, -0.014409, 2.03), ), ((0.20169, -0.014409, 
            2.03), ), ((0.186588, -0.014409, 2.03), ), ((0.168343, -0.014374, 2.03), ), 
            ((0.153055, -0.013768, 2.03), ), ((0.139681, -0.012768, 2.03), ), ((
            0.128372, -0.01136, 2.03), ), ((0.119243, -0.00955, 2.03), ), ((0.112413, 
            -0.007367, 2.03), ), )+mdb.models['Model-1'].parts['Wing'].edges.findAt(((
            0.108129, -0.004829, 2.03), ), ((0.542216, 0.0, 1.945), ), ((0.528359, 
            0.00072, 1.945), ), ((0.517873, 0.000954, 1.945), ), ((0.504736, 0.001016, 
            1.945), ), ((0.489183, 0.000835, 1.945), ), ((0.471472, 0.000384, 1.945), 
            ), ((0.451884, -0.000349, 1.945), ), ((0.435, -0.001294, 1.945), ), ((
            0.384923, -0.004096, 1.945), ), ((0.360938, -0.005722, 1.945), ), ((
            0.336644, -0.007423, 1.945), ), ((0.31235, -0.00912, 1.945), ), ((0.288348, 
            -0.01171, 1.945), ), ((0.264099, -0.014701, 1.945), ), ((0.248185, 
            -0.014701, 1.945), ), ((0.198185, -0.014701, 1.945), ), ((0.182271, 
            -0.014701, 1.945), ), ((0.163656, -0.014666, 1.945), ), ((0.148058, 
            -0.014047, 1.945), ), ((0.134413, -0.013026, 1.945), ), ((0.122876, 
            -0.01159, 1.945), ), ((0.113561, -0.009743, 1.945), ), ((0.106593, 
            -0.007516, 1.945), ), ((0.102223, -0.004927, 1.945), ), ((0.100345, 0.0, 
            1.945), ), ((0.1008, 0.003142, 1.945), ), ((0.102828, 0.008267, 1.945), ), 
            ((0.106575, 0.013614, 1.945), ), ((0.112055, 0.019, 1.945), ), ((0.119222, 
            0.024263, 1.945), ), ((0.128046, 0.029234, 1.945), ), ((0.138509, 0.03379, 
            1.945), ), ((0.150555, 0.037833, 1.945), ), ((0.164116, 0.041978, 1.945), 
            ), ((0.182271, 0.046803, 1.945), ), ((0.198185, 0.046803, 1.945), ), ((
            0.248185, 0.046803, 1.945), ), ((0.264099, 0.046803, 1.945), ), ((0.292736, 
            0.045199, 1.945), ), ((0.314219, 0.042963, 1.945), ), ((0.335959, 0.040162, 
            1.945), ), ((0.357752, 0.036883, 1.945), ), ((0.379404, 0.033215, 1.945), 
            ), ((0.435, 0.022348, 1.945), ), ((0.441319, 0.021113, 1.945), ), ((
            0.460147, 0.017167, 1.945), ), ((0.477619, 0.013473, 1.945), ), ((0.493464, 
            0.010114, 1.945), ), ((0.507427, 0.007136, 1.945), ), ((0.519318, 0.004578, 
            1.945), ), ((0.528977, 0.002523, 1.945), ), ((0.093352, 0.0, 1.845), ), ((
            0.093818, 0.003215, 1.845), ), ((0.095894, 0.008461, 1.845), ), ((0.099728, 
            0.013932, 1.845), ), ((0.105335, 0.019444, 1.845), ), ((0.11267, 0.02483, 
            1.845), ), ((0.1217, 0.029917, 1.845), ), ((0.132408, 0.034579, 1.845), ), 
            ((0.144735, 0.038717, 1.845), ), ((0.158613, 0.042958, 1.845), ), ((
            0.177192, 0.047896, 1.845), ), ((0.194062, 0.047896, 1.845), ), ((0.244062, 
            0.047896, 1.845), ), ((0.260932, 0.047896, 1.845), ), ((0.290237, 0.046255, 
            1.845), ), ((0.312223, 0.043967, 1.845), ), ((0.334471, 0.0411, 1.845), ), 
            ((0.356773, 0.037745, 1.845), ), ((0.378931, 0.033991, 1.845), ), ((0.435, 
            0.023031, 1.845), ), ((0.442292, 0.021606, 1.845), ), ((0.46156, 0.017568, 
            1.845), ), ((0.47944, 0.013787, 1.845), ), ((0.495655, 0.010351, 1.845), ), 
            ((0.509945, 0.007303, 1.845), ), ((0.522113, 0.004685, 1.845), ), ((
            0.531998, 0.002582, 1.845), ), ((0.545546, 0.0, 1.845), ), ((0.531365, 
            0.000737, 1.845), ), ((0.520634, 0.000977, 1.845), ), ((0.507191, 0.00104, 
            1.845), ), ((0.491274, 0.000855, 1.845), ), ((0.47315, 0.000393, 1.845), ), 
            ((0.453104, -0.000357, 1.845), ), ((0.435, -0.00137, 1.845), ), ((0.384579, 
            -0.004192, 1.845), ), ((0.360033, -0.005856, 1.845), ), ((0.335172, 
            -0.007597, 1.845), ), ((0.31031, -0.009333, 1.845), ), ((0.285747, 
            -0.011983, 1.845), ), ((0.260932, -0.015044, 1.845), ), ((0.244062, 
            -0.015044, 1.845), ), ((0.194062, -0.015044, 1.845), ), ((0.177192, 
            -0.015044, 1.845), ), ((0.158143, -0.015008, 1.845), ), ((0.14218, 
            -0.014375, 1.845), ), ((0.128216, -0.013331, 1.845), ), ((0.11641, 
            -0.011861, 1.845), ), ((0.106877, -0.009971, 1.845), ), ((0.099746, 
            -0.007692, 1.845), ), ((0.095274, -0.005042, 1.845), ), ((0.548876, 0.0, 
            1.745), ), ((0.534371, 0.000754, 1.745), ), ((0.523396, 0.000999, 1.745), 
            ), ((0.509645, 0.001064, 1.745), ), ((0.493365, 0.000874, 1.745), ), ((
            0.474827, 0.000402, 1.745), ), ((0.454324, -0.000365, 1.745), ), ((0.435, 
            -0.001447, 1.745), ), ((0.384234, -0.004288, 1.745), ), ((0.359128, 
            -0.00599, 1.745), ), ((0.333699, -0.00777, 1.745), ), ((0.30827, -0.009546, 
            1.745), ), ((0.283146, -0.012257, 1.745), ), ((0.257765, -0.015388, 1.745), 
            ), ((0.239939, -0.015388, 1.745), ), ((0.189939, -0.015388, 1.745), ), ((
            0.172114, -0.015388, 1.745), ), ((0.152629, -0.015351, 1.745), ), ((
            0.136302, -0.014703, 1.745), ), ((0.12202, -0.013635, 1.745), ), ((
            0.109943, -0.012132, 1.745), ), ((0.100193, -0.010198, 1.745), ), ((0.0929, 
            -0.007867, 1.745), ), ((0.088325, -0.005157, 1.745), ), ((0.08636, 0.0, 
            1.745), ), ((0.086836, 0.003288, 1.745), ), ((0.088959, 0.008654, 1.745), 
            ), ((0.092881, 0.01425, 1.745), ), ((0.098616, 0.019888, 1.745), ), ((
            0.106118, 0.025397, 1.745), ), ((0.115355, 0.0306, 1.745), ), ((0.126307, 
            0.035369, 1.745), ), ((0.138915, 0.039601, 1.745), ), ((0.15311, 0.043939, 
            1.745), ), ((0.172114, 0.04899, 1.745), ), ((0.189939, 0.04899, 1.745), ), 
            ((0.239939, 0.04899, 1.745), ), ((0.257765, 0.04899, 1.745), ), ((0.287739, 
            0.047311, 1.745), ), ((0.310227, 0.04497, 1.745), ), ((0.332982, 0.042038, 
            1.745), ), ((0.355794, 0.038606, 1.745), ), ((0.378457, 0.034767, 1.745), 
            ), ((0.435, 0.023715, 1.745), ), ((0.443265, 0.022099, 1.745), ), ((
            0.462973, 0.017969, 1.745), ), ((0.481261, 0.014102, 1.745), ), ((0.497846, 
            0.010587, 1.745), ), ((0.512462, 0.00747, 1.745), ), ((0.524908, 0.004792, 
            1.745), ), ((0.535019, 0.002641, 1.745), ), ((0.079367, 0.0, 1.645), ), ((
            0.079854, 0.003362, 1.645), ), ((0.082024, 0.008847, 1.645), ), ((0.086034, 
            0.014568, 1.645), ), ((0.091897, 0.020332, 1.645), ), ((0.099567, 0.025964, 
            1.645), ), ((0.109009, 0.031283, 1.645), ), ((0.120206, 0.036158, 1.645), 
            ), ((0.133096, 0.040484, 1.645), ), ((0.147607, 0.04492, 1.645), ), ((
            0.167035, 0.050083, 1.645), ), ((0.185816, 0.050083, 1.645), ), ((0.235816, 
            0.050083, 1.645), ), ((0.254597, 0.050083, 1.645), ), ((0.285241, 0.048367, 
            1.645), ), ((0.30823, 0.045974, 1.645), ), ((0.331494, 0.042976, 1.645), ), 
            ((0.354814, 0.039468, 1.645), ), ((0.377983, 0.035543, 1.645), ), ((0.435, 
            0.024398, 1.645), ), ((0.444238, 0.022592, 1.645), ), ((0.464385, 0.01837, 
            1.645), ), ((0.483081, 0.014417, 1.645), ), ((0.500037, 0.010823, 1.645), 
            ), ((0.514979, 0.007636, 1.645), ), ((0.527703, 0.004899, 1.645), ), ((
            0.538039, 0.0027, 1.645), ), ((0.552206, 0.0, 1.645), ), ((0.537377, 
            0.000771, 1.645), ), ((0.526157, 0.001021, 1.645), ), ((0.512099, 0.001088, 
            1.645), ), ((0.495456, 0.000894, 1.645), ), ((0.476504, 0.000411, 1.645), 
            ), ((0.455543, -0.000374, 1.645), ), ((0.435, -0.001523, 1.645), ), ((
            0.383889, -0.004383, 1.645), ), ((0.358224, -0.006123, 1.645), ), ((
            0.332227, -0.007944, 1.645), ), ((0.30623, -0.009759, 1.645), ), ((
            0.280546, -0.01253, 1.645), ), ((0.254597, -0.015731, 1.645), ), ((
            0.235816, -0.015731, 1.645), ), ((0.185816, -0.015731, 1.645), ), ((
            0.167035, -0.015731, 1.645), ), ((0.147115, -0.015694, 1.645), ), ((
            0.130424, -0.015032, 1.645), ), ((0.115823, -0.013939, 1.645), ), ((
            0.103477, -0.012403, 1.645), ), ((0.09351, -0.010426, 1.645), ), ((
            0.086053, -0.008043, 1.645), ), ((0.081376, -0.005272, 1.645), ), ((
            0.555536, 0.0, 1.545), ), ((0.540384, 0.000788, 1.545), ), ((0.528918, 
            0.001044, 1.545), ), ((0.514554, 0.001111, 1.545), ), ((0.497547, 0.000913, 
            1.545), ), ((0.478181, 0.00042, 1.545), ), ((0.456763, -0.000382, 1.545), 
            ), ((0.435, -0.0016, 1.545), ), ((0.383545, -0.004479, 1.545), ), ((
            0.357319, -0.006257, 1.545), ), ((0.330754, -0.008117, 1.545), ), ((
            0.30419, -0.009972, 1.545), ), ((0.277945, -0.012804, 1.545), ), ((0.25143, 
            -0.016075, 1.545), ), ((0.231693, -0.016075, 1.545), ), ((0.181693, 
            -0.016075, 1.545), ), ((0.161956, -0.016075, 1.545), ), ((0.141602, 
            -0.016036, 1.545), ), ((0.124546, -0.01536, 1.545), ), ((0.109626, 
            -0.014244, 1.545), ), ((0.097011, -0.012673, 1.545), ), ((0.086826, 
            -0.010654, 1.545), ), ((0.079206, -0.008219, 1.545), ), ((0.074428, 
            -0.005387, 1.545), ), ((0.072374, 0.0, 1.545), ), ((0.072872, 0.003435, 
            1.545), ), ((0.07509, 0.00904, 1.545), ), ((0.079187, 0.014886, 1.545), ), 
            ((0.085178, 0.020776, 1.545), ), ((0.093015, 0.02653, 1.545), ), ((
            0.102664, 0.031966, 1.545), ), ((0.114105, 0.036947, 1.545), ), ((0.127276, 
            0.041368, 1.545), ), ((0.142104, 0.0459, 1.545), ), ((0.161956, 0.051176, 
            1.545), ), ((0.181693, 0.051176, 1.545), ), ((0.231693, 0.051176, 1.545), 
            ), ((0.25143, 0.051176, 1.545), ), ((0.282743, 0.049423, 1.545), ), ((
            0.306234, 0.046978, 1.545), ), ((0.330006, 0.043915, 1.545), ), ((0.353835, 
            0.040329, 1.545), ), ((0.37751, 0.036319, 1.545), ), ((0.435, 0.025081, 
            1.545), ), ((0.44521, 0.023085, 1.545), ), ((0.465798, 0.018771, 1.545), ), 
            ((0.484902, 0.014732, 1.545), ), ((0.502228, 0.01106, 1.545), ), ((
            0.517496, 0.007803, 1.545), ), )+\
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.530498, 0.005006, 
            1.545), ), ((0.54106, 0.002759, 1.545), ), ((0.065382, 0.0, 1.445), ), ((
            0.06589, 0.003509, 1.445), ), ((0.068155, 0.009233, 1.445), ), ((0.07234, 
            0.015204, 1.445), ), ((0.078459, 0.02122, 1.445), ), ((0.086463, 0.027097, 
            1.445), ), ((0.096318, 0.032649, 1.445), ), ((0.108004, 0.037737, 1.445), 
            ), ((0.121456, 0.042252, 1.445), ), ((0.136601, 0.046881, 1.445), ), ((
            0.156877, 0.05227, 1.445), ), ((0.17757, 0.05227, 1.445), ), ((0.22757, 
            0.05227, 1.445), ), ((0.248263, 0.05227, 1.445), ), ((0.280244, 0.050478, 
            1.445), ), ((0.304238, 0.047981, 1.445), ), ((0.328517, 0.044853, 1.445), 
            ), ((0.352856, 0.041191, 1.445), ), ((0.377036, 0.037095, 1.445), ), ((
            0.435, 0.025765, 1.445), ), ((0.446183, 0.023579, 1.445), ), ((0.467211, 
            0.019172, 1.445), ), ((0.486723, 0.015046, 1.445), ), ((0.504419, 0.011296, 
            1.445), ), ((0.520013, 0.00797, 1.445), ), ((0.533293, 0.005112, 1.445), ), 
            ((0.544081, 0.002818, 1.445), ), ((0.558865, 0.0, 1.445), ), ((0.54339, 
            0.000804, 1.445), ), ((0.531679, 0.001066, 1.445), ), ((0.517008, 0.001135, 
            1.445), ), ((0.499638, 0.000933, 1.445), ), ((0.479859, 0.000429, 1.445), 
            ), ((0.457983, -0.00039, 1.445), ), ((0.435, -0.001676, 1.445), ), ((
            0.3832, -0.004575, 1.445), ), ((0.356414, -0.006391, 1.445), ), ((0.329282, 
            -0.008291, 1.445), ), ((0.30215, -0.010186, 1.445), ), ((0.275344, 
            -0.013077, 1.445), ), ((0.248263, -0.016418, 1.445), ), ((0.22757, 
            -0.016418, 1.445), ), ((0.17757, -0.016418, 1.445), ), ((0.156877, 
            -0.016418, 1.445), ), ((0.136088, -0.016379, 1.445), ), ((0.118668, 
            -0.015688, 1.445), ), ((0.103429, -0.014548, 1.445), ), ((0.090544, 
            -0.012944, 1.445), ), ((0.080142, -0.010881, 1.445), ), ((0.072359, 
            -0.008394, 1.445), ), ((0.067479, -0.005502, 1.445), ), ), name='Set-22')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='skin tip', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin tip'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin tip Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-tip'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'skin tip Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-tip'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'skin tip Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-tip'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin tip Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Skin-tip'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin tip'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=mdb.models['Model-1'].parts['Wing'].surfaces['skin-tip'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-22'], stackDirection=STACK_3)




        #spar cap bottom mid

        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.331457, -0.063854, 
            2.334276), ), ((0.306457, -0.063854, 2.334276), ), ((0.281457, -0.063854, 
            2.334276), ), ((0.305517, -0.047401, 2.242776), ), ((0.280517, -0.047401, 
            2.242776), ), ((0.255517, -0.047401, 2.242776), ), ((0.27709, -0.02937, 
            2.142503), ), ((0.25209, -0.02937, 2.142503), ), ((0.22709, -0.02937, 
            2.142503), ), ((0.257191, -0.016945, 2.074359), ), ((0.232191, -0.016945, 
            2.074359), ), ((0.207191, -0.016945, 2.074359), ), ((0.25169, -0.014409, 
            2.03), ), ((0.22669, -0.014409, 2.03), ), ((0.20169, -0.014409, 2.03), ), (
            (0.248185, -0.014701, 1.945), ), ((0.223185, -0.014701, 1.945), ), ((
            0.198185, -0.014701, 1.945), ), ((0.244062, -0.015044, 1.845), ), ((
            0.219062, -0.015044, 1.845), ), ((0.194062, -0.015044, 1.845), ), ((
            0.239939, -0.015388, 1.745), ), ((0.214939, -0.015388, 1.745), ), ((
            0.189939, -0.015388, 1.745), ), ((0.235816, -0.015731, 1.645), ), ((
            0.210816, -0.015731, 1.645), ), ((0.185816, -0.015731, 1.645), ), ((
            0.231693, -0.016075, 1.545), ), ((0.206693, -0.016075, 1.545), ), ((
            0.181693, -0.016075, 1.545), ), ((0.22757, -0.016418, 1.445), ), ((0.20257, 
            -0.016418, 1.445), ), ((0.17757, -0.016418, 1.445), ), ((0.223447, 
            -0.016762, 1.345), ), ((0.198447, -0.016762, 1.345), ), ((0.173447, 
            -0.016762, 1.345), ), ((0.219324, -0.017105, 1.245), ), ((0.194324, 
            -0.017105, 1.245), ), ((0.169324, -0.017105, 1.245), ), ((0.215201, 
            -0.017449, 1.145), ), ((0.190201, -0.017449, 1.145), ), ((0.165201, 
            -0.017449, 1.145), ), ((0.211078, -0.017792, 1.045), ), ((0.186078, 
            -0.017792, 1.045), ), ((0.161078, -0.017792, 1.045), ), ((0.206955, 
            -0.018135, 0.945), ), ((0.181955, -0.018135, 0.945), ), ((0.156955, 
            -0.018135, 0.945), ), ((0.202832, -0.018479, 0.845), ), ((0.177832, 
            -0.018479, 0.845), ), ((0.152832, -0.018479, 0.845), ), ((0.198709, 
            -0.018822, 0.745), ), ((0.173709, -0.018822, 0.745), ), ((0.148709, 
            -0.018822, 0.745), ), ((0.144586, -0.019166, 0.645), ), ((0.169586, 
            -0.019166, 0.645), ), ((0.194586, -0.019166, 0.645), ), ((0.1917, 
            -0.019406, 0.575), ), ((0.1667, -0.019406, 0.575), ), ((0.1417, -0.019406, 
            0.575), ), ((0.189535, -0.019586, 0.5225), ), ((0.164535, -0.019586, 
            0.5225), ), ((0.139535, -0.019586, 0.5225), ), ((0.18902, -0.019629, 
            0.4725), ), ((0.16402, -0.019629, 0.4725), ), ((0.13902, -0.019629, 
            0.4725), ), ((0.18902, -0.019629, 0.415), ), ((0.16402, -0.019629, 0.415), 
            ), ((0.13902, -0.019629, 0.415), ), ((0.18902, -0.019629, 0.325), ), ((
            0.16402, -0.019629, 0.325), ), ((0.13902, -0.019629, 0.325), ), ((0.18902, 
            -0.019629, 0.225), ), ((0.16402, -0.019629, 0.225), ), ((0.13902, 
            -0.019629, 0.225), ), ((0.18902, -0.019629, 0.125), ), ((0.16402, 
            -0.019629, 0.125), ), ((0.13902, -0.019629, 0.125), ), ((0.18902, 
            -0.019629, 0.025), ), ((0.16402, -0.019629, 0.025), ), ((0.13902, 
            -0.019629, 0.025), ), ), name='Set-27')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='spar cap bottom mid', offsetType=MIDDLE_SURFACE, symmetric=
            False, thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BM Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BM Ply-1a', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='spar cap BM Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BM Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap BM Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BM Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap BM Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BM Ply-7', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap bottom mid'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom mid'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['spar cap bottom mid'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-27'], stackDirection=STACK_3)
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.281457, -0.021362, 
            2.342535), ), ((0.306457, -0.021362, 2.342535), ), ((0.331457, -0.021362, 
            2.342535), ), ((0.255517, 0.000463, 2.25208), ), ((0.280517, 0.000463, 
            2.25208), ), ((0.305517, 0.000463, 2.25208), ), ((0.22709, 0.024379, 
            2.15295), ), ((0.25209, 0.024379, 2.15295), ), ((0.27709, 0.024379, 
            2.15295), ), ((0.207191, 0.041749, 2.077042), ), ((0.232191, 0.041749, 
            2.077042), ), ((0.257191, 0.041749, 2.077042), ), ((0.20169, 0.045874, 
            2.03), ), ((0.22669, 0.045874, 2.03), ), ((0.25169, 0.045874, 2.03), ), ((
            0.198185, 0.046803, 1.945), ), ((0.223185, 0.046803, 1.945), ), ((0.248185, 
            0.046803, 1.945), ), ((0.194062, 0.047896, 1.845), ), ((0.219062, 0.047896, 
            1.845), ), ((0.244062, 0.047896, 1.845), ), ((0.189939, 0.04899, 1.745), ), 
            ((0.214939, 0.04899, 1.745), ), ((0.239939, 0.04899, 1.745), ), ((0.185816, 
            0.050083, 1.645), ), ((0.210816, 0.050083, 1.645), ), ((0.235816, 0.050083, 
            1.645), ), ((0.181693, 0.051176, 1.545), ), ((0.206693, 0.051176, 1.545), 
            ), ((0.231693, 0.051176, 1.545), ), ((0.17757, 0.05227, 1.445), ), ((
            0.20257, 0.05227, 1.445), ), ((0.22757, 0.05227, 1.445), ), ((0.173447, 
            0.053363, 1.345), ), ((0.198447, 0.053363, 1.345), ), ((0.223447, 0.053363, 
            1.345), ), ((0.169324, 0.054457, 1.245), ), ((0.194324, 0.054457, 1.245), 
            ), ((0.219324, 0.054457, 1.245), ), ((0.165201, 0.05555, 1.145), ), ((
            0.190201, 0.05555, 1.145), ), ((0.215201, 0.05555, 1.145), ), ((0.161078, 
            0.056643, 1.045), ), ((0.186078, 0.056643, 1.045), ), ((0.211078, 0.056643, 
            1.045), ), ((0.156955, 0.057737, 0.945), ), ((0.181955, 0.057737, 0.945), 
            ), ((0.206955, 0.057737, 0.945), ), ((0.152832, 0.05883, 0.845), ), ((
            0.177832, 0.05883, 0.845), ), ((0.202832, 0.05883, 0.845), ), ((0.148709, 
            0.059923, 0.745), ), ((0.173709, 0.059923, 0.745), ), ((0.198709, 0.059923, 
            0.745), ), ((0.169586, 0.061017, 0.645), ), ((0.194586, 0.061017, 0.645), 
            ), ((0.144586, 0.061017, 0.645), ), ((0.1417, 0.061782, 0.575), ), ((
            0.1667, 0.061782, 0.575), ), ((0.1917, 0.061782, 0.575), ), ((0.139535, 
            0.062356, 0.5225), ), ((0.164535, 0.062356, 0.5225), ), ((0.189535, 
            0.062356, 0.5225), ), ((0.13902, 0.062493, 0.4725), ), ((0.16402, 0.062493, 
            0.4725), ), ((0.18902, 0.062493, 0.4725), ), ((0.13902, 0.062493, 0.415), 
            ), ((0.16402, 0.062493, 0.415), ), ((0.18902, 0.062493, 0.415), ), ((
            0.13902, 0.062493, 0.325), ), ((0.16402, 0.062493, 0.325), ), ((0.18902, 
            0.062493, 0.325), ), ((0.13902, 0.062493, 0.225), ), ((0.16402, 0.062493, 
            0.225), ), ((0.18902, 0.062493, 0.225), ), ((0.13902, 0.062493, 0.125), ), 
            ((0.16402, 0.062493, 0.125), ), ((0.18902, 0.062493, 0.125), ), ((0.13902, 
            0.062493, 0.025), ), ((0.16402, 0.062493, 0.025), ), ((0.18902, 0.062493, 
            0.025), ), ), name='Set-28')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='spar cap top mid', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap top mid Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap top mid Ply-1a', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='spar cap top mid Ply-2', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], 
            suppressed=False, thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap top mid Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap top mid Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap top mid Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap top mid Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap top mid Ply-7', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top mid'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top mid'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['spar cap top mid'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-28'], stackDirection=STACK_3)

        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.331457, -0.063854, 
            2.334276), ), ((0.306457, -0.063854, 2.334276), ), ((0.281457, -0.063854, 
            2.334276), ), ((0.305517, -0.047401, 2.242776), ), ((0.280517, -0.047401, 
            2.242776), ), ((0.255517, -0.047401, 2.242776), ), ((0.27709, -0.02937, 
            2.142503), ), ((0.25209, -0.02937, 2.142503), ), ((0.22709, -0.02937, 
            2.142503), ), ((0.257191, -0.016945, 2.074359), ), ((0.232191, -0.016945, 
            2.074359), ), ((0.207191, -0.016945, 2.074359), ), ((0.25169, -0.014409, 
            2.03), ), ((0.22669, -0.014409, 2.03), ), ((0.20169, -0.014409, 2.03), ), (
            (0.248185, -0.014701, 1.945), ), ((0.223185, -0.014701, 1.945), ), ((
            0.198185, -0.014701, 1.945), ), ((0.244062, -0.015044, 1.845), ), ((
            0.219062, -0.015044, 1.845), ), ((0.194062, -0.015044, 1.845), ), ((
            0.239939, -0.015388, 1.745), ), ((0.214939, -0.015388, 1.745), ), ((
            0.189939, -0.015388, 1.745), ), ((0.235816, -0.015731, 1.645), ), ((
            0.210816, -0.015731, 1.645), ), ((0.185816, -0.015731, 1.645), ), ((
            0.231693, -0.016075, 1.545), ), ((0.206693, -0.016075, 1.545), ), ((
            0.181693, -0.016075, 1.545), ), ((0.22757, -0.016418, 1.445), ), ((0.20257, 
            -0.016418, 1.445), ), ((0.17757, -0.016418, 1.445), ), ((0.223447, 
            -0.016762, 1.345), ), ((0.198447, -0.016762, 1.345), ), ((0.173447, 
            -0.016762, 1.345), ), ((0.219324, -0.017105, 1.245), ), ((0.194324, 
            -0.017105, 1.245), ), ((0.169324, -0.017105, 1.245), ), ((0.215201, 
            -0.017449, 1.145), ), ((0.190201, -0.017449, 1.145), ), ((0.165201, 
            -0.017449, 1.145), ), ((0.211078, -0.017792, 1.045), ), ((0.186078, 
            -0.017792, 1.045), ), ((0.161078, -0.017792, 1.045), ), ((0.206955, 
            -0.018135, 0.945), ), ((0.181955, -0.018135, 0.945), ), ((0.156955, 
            -0.018135, 0.945), ), ((0.202832, -0.018479, 0.845), ), ((0.177832, 
            -0.018479, 0.845), ), ((0.152832, -0.018479, 0.845), ), ((0.198709, 
            -0.018822, 0.745), ), ((0.173709, -0.018822, 0.745), ), ((0.148709, 
            -0.018822, 0.745), ), ((0.150836, -0.019166, 0.645), ), ((0.144586, 
            -0.019166, 0.645), ), ((0.188336, -0.019166, 0.645), ), ((0.169586, 
            -0.019166, 0.645), ), ((0.194586, -0.019166, 0.645), ), ((0.1917, 
            -0.019406, 0.575), ), ((0.1667, -0.019406, 0.575), ), ((0.1417, -0.019406, 
            0.575), ), ((0.189535, -0.019586, 0.5225), ), ((0.164535, -0.019586, 
            0.5225), ), ((0.139535, -0.019586, 0.5225), ), ((0.18902, -0.019629, 
            0.4725), ), ((0.16402, -0.019629, 0.4725), ), ((0.13902, -0.019629, 
            0.4725), ), ((0.18902, -0.019629, 0.415), ), ((0.16402, -0.019629, 0.415), 
            ), ((0.13902, -0.019629, 0.415), ), ((0.18902, -0.019629, 0.325), ), ((
            0.16402, -0.019629, 0.325), ), ((0.13902, -0.019629, 0.325), ), ((0.18902, 
            -0.019629, 0.225), ), ((0.16402, -0.019629, 0.225), ), ((0.13902, 
            -0.019629, 0.225), ), ((0.18902, -0.019629, 0.125), ), ((0.16402, 
            -0.019629, 0.125), ), ((0.13902, -0.019629, 0.125), ), ((0.18902, 
            -0.019629, 0.025), ), ((0.16402, -0.019629, 0.025), ), ((0.13902, 
            -0.019629, 0.025), ), ), name='Set-29')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='spar cap bottom tip', offsetType=MIDDLE_SURFACE, symmetric=
            False, thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BT Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BT Ply-1a', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='spar cap BT Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BT Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap BT Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap BT Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap BT Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-bottom tip'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap bottom tip'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['spar cap bottom tip'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-29'], stackDirection=STACK_3)
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.281457, -0.021362, 
            2.342535), ), ((0.306457, -0.021362, 2.342535), ), ((0.331457, -0.021362, 
            2.342535), ), ((0.255517, 0.000463, 2.25208), ), ((0.280517, 0.000463, 
            2.25208), ), ((0.305517, 0.000463, 2.25208), ), ((0.22709, 0.024379, 
            2.15295), ), ((0.25209, 0.024379, 2.15295), ), ((0.27709, 0.024379, 
            2.15295), ), ((0.207191, 0.041749, 2.077042), ), ((0.232191, 0.041749, 
            2.077042), ), ((0.257191, 0.041749, 2.077042), ), ((0.20169, 0.045874, 
            2.03), ), ((0.22669, 0.045874, 2.03), ), ((0.25169, 0.045874, 2.03), ), ((
            0.198185, 0.046803, 1.945), ), ((0.223185, 0.046803, 1.945), ), ((0.248185, 
            0.046803, 1.945), ), ((0.194062, 0.047896, 1.845), ), ((0.219062, 0.047896, 
            1.845), ), ((0.244062, 0.047896, 1.845), ), ((0.189939, 0.04899, 1.745), ), 
            ((0.214939, 0.04899, 1.745), ), ((0.239939, 0.04899, 1.745), ), ((0.185816, 
            0.050083, 1.645), ), ((0.210816, 0.050083, 1.645), ), ((0.235816, 0.050083, 
            1.645), ), ((0.181693, 0.051176, 1.545), ), ((0.206693, 0.051176, 1.545), 
            ), ((0.231693, 0.051176, 1.545), ), ((0.17757, 0.05227, 1.445), ), ((
            0.20257, 0.05227, 1.445), ), ((0.22757, 0.05227, 1.445), ), ((0.173447, 
            0.053363, 1.345), ), ((0.198447, 0.053363, 1.345), ), ((0.223447, 0.053363, 
            1.345), ), ((0.169324, 0.054457, 1.245), ), ((0.194324, 0.054457, 1.245), 
            ), ((0.219324, 0.054457, 1.245), ), ((0.165201, 0.05555, 1.145), ), ((
            0.190201, 0.05555, 1.145), ), ((0.215201, 0.05555, 1.145), ), ((0.161078, 
            0.056643, 1.045), ), ((0.186078, 0.056643, 1.045), ), ((0.211078, 0.056643, 
            1.045), ), ((0.156955, 0.057737, 0.945), ), ((0.181955, 0.057737, 0.945), 
            ), ((0.206955, 0.057737, 0.945), ), ((0.152832, 0.05883, 0.845), ), ((
            0.177832, 0.05883, 0.845), ), ((0.202832, 0.05883, 0.845), ), ((0.148709, 
            0.059923, 0.745), ), ((0.173709, 0.059923, 0.745), ), ((0.198709, 0.059923, 
            0.745), ), ((0.150836, 0.061017, 0.645), ), ((0.169586, 0.061017, 0.645), 
            ), ((0.188336, 0.061017, 0.645), ), ((0.194586, 0.061017, 0.645), ), ((
            0.144586, 0.061017, 0.645), ), ((0.1417, 0.061782, 0.575), ), ((0.1667, 
            0.061782, 0.575), ), ((0.1917, 0.061782, 0.575), ), ((0.139535, 0.062356, 
            0.5225), ), ((0.164535, 0.062356, 0.5225), ), ((0.189535, 0.062356, 
            0.5225), ), ((0.13902, 0.062493, 0.4725), ), ((0.16402, 0.062493, 0.4725), 
            ), ((0.18902, 0.062493, 0.4725), ), ((0.13902, 0.062493, 0.415), ), ((
            0.16402, 0.062493, 0.415), ), ((0.18902, 0.062493, 0.415), ), ((0.13902, 
            0.062493, 0.325), ), ((0.16402, 0.062493, 0.325), ), ((0.18902, 0.062493, 
            0.325), ), ((0.13902, 0.062493, 0.225), ), ((0.16402, 0.062493, 0.225), ), 
            ((0.18902, 0.062493, 0.225), ), ((0.13902, 0.062493, 0.125), ), ((0.16402, 
            0.062493, 0.125), ), ((0.18902, 0.062493, 0.125), ), ((0.13902, 0.062493, 
            0.025), ), ((0.16402, 0.062493, 0.025), ), ((0.18902, 0.062493, 0.025), ), 
            ), name='Set-30')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='spar cap top tip', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap TT Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap TT Ply-1a', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='spar cap TT Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap TT Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap TT Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'spar cap TT Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'spar cap TT Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['Spar cap-top tip'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['spar cap top tip'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['spar cap top tip'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-30'], stackDirection=STACK_3)


        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.000608, 0.004195, 
            0.04375), ), ((0.003316, 0.011039, 0.04375), ), ((0.0, 0.0, 0.04375), ), ((
            0.002508, -0.006579, 0.04375), ), ((0.008343, -0.010036, 0.04375), ), ((
            0.017647, -0.01301, 0.04375), ), ((0.030084, -0.015476, 0.04375), ), ((
            0.045489, -0.017393, 0.04375), ), ((0.063708, -0.018756, 0.04375), ), ((
            0.084535, -0.019582, 0.04375), ), ((0.10939, -0.019629, 0.04375), ), ((
            0.21865, -0.019629, 0.04375), ), ((0.251027, -0.015635, 0.04375), ), ((
            0.283076, -0.012178, 0.04375), ), ((0.315514, -0.009912, 0.04375), ), ((
            0.347952, -0.007641, 0.04375), ), ((0.379978, -0.005469, 0.04375), ), ((
            0.455279, 0.02819, 0.04375), ), ((0.469386, -0.000466, 0.04375), ), ((
            0.495541, 0.000513, 0.04375), ), ((0.519188, 0.001115, 0.04375), ), ((
            0.539956, 0.001357, 0.04375), ), ((0.557497, 0.001274, 0.04375), ), ((
            0.571498, 0.000962, 0.04375), ), ((0.004545, 0.0, 0.575), ), ((0.005146, 
            0.004147, 0.575), ), ((0.007823, 0.010913, 0.575), ), ((0.01277, 0.017971, 
            0.575), ), ((0.020002, 0.025081, 0.575), ), ((0.029463, 0.032028, 0.575), 
            ), ((0.041112, 0.03859, 0.575), ), ((0.054924, 0.044604, 0.575), ), ((
            0.070825, 0.049941, 0.575), ), ((0.088726, 0.055413, 0.575), ), ((0.112692, 
            0.061782, 0.575), ), ((0.220708, 0.061782, 0.575), ), ((0.25851, 0.059665, 
            0.575), ), ((0.286869, 0.056713, 0.575), ), ((0.315567, 0.053015, 0.575), 
            ), ((0.344335, 0.048687, 0.575), ), ((0.372916, 0.043846, 0.575), ), ((
            0.454647, 0.02787, 0.575), ), ((0.479501, 0.022661, 0.575), ), ((0.502564, 
            0.017785, 0.575), ), ((0.523481, 0.013352, 0.575), ), ((0.541913, 0.00942, 
            0.575), ), ((0.557609, 0.006043, 0.575), ), ((0.57036, 0.003331, 0.575), ), 
            ((0.587836, 0.0, 0.575), ), ((0.569544, 0.000951, 0.575), ), ((0.555702, 
            0.00126, 0.575), ), ((0.538361, 0.001342, 0.575), ), ((0.517829, 0.001102, 
            0.575), ), ((0.494451, 0.000507, 0.575), ), ((0.468594, -0.000461, 0.575), 
            ), ((0.380202, -0.005407, 0.575), ), ((0.348541, -0.007554, 0.575), ), ((
            0.316471, -0.009799, 0.575), ), ((0.284402, -0.012039, 0.575), ), ((
            0.252718, -0.015457, 0.575), ), ((0.220708, -0.019406, 0.575), ), ((
            0.112692, -0.019406, 0.575), ), ((0.088119, -0.019359, 0.575), ), ((
            0.067529, -0.018543, 0.575), ), ((0.049517, -0.017195, 0.575), ), ((
            0.034287, -0.0153, 0.575), ), ((0.021991, -0.012862, 0.575), ), ((0.012793, 
            -0.009922, 0.575), ), ((0.007024, -0.006504, 0.575), ), ((0.589584, 0.0, 
            0.5225), ), ((0.571122, 0.00096, 0.5225), ), ((0.557152, 0.001272, 0.5225), 
            ), ((0.539649, 0.001354, 0.5225), ), ((0.518927, 0.001113, 0.5225), ), ((
            0.495331, 0.000512, 0.5225), ), ((0.469234, -0.000465, 0.5225), ), ((
            0.380021, -0.005457, 0.5225), ), ((0.348066, -0.007624, 0.5225), ), ((
            0.315698, -0.00989, 0.5225), ), ((0.283331, -0.012151, 0.5225), ), ((
            0.251352, -0.015601, 0.5225), ), ((0.219046, -0.019586, 0.5225), ), ((
            0.110025, -0.019586, 0.5225), ), ((0.085224, -0.019539, 0.5225), ), ((
            0.064443, -0.018715, 0.5225), ), ((0.046264, -0.017355, 0.5225), ), ((
            0.030892, -0.015442, 0.5225), ), ((0.018482, -0.012981, 0.5225), ), ((
            0.009198, -0.010014, 0.5225), ), ((0.003376, -0.006564, 0.5225), ), ((
            0.000874, 0.0, 0.5225), ), ((0.00148, 0.004186, 0.5225), ), ((0.004183, 
            0.011015, 0.5225), ), ((0.009175, 0.018138, 0.5225), ), ((0.016475, 
            0.025315, 0.5225), ), ((0.026024, 0.032326, 0.5225), ), ((0.03778, 
            0.038949, 0.5225), ), ((0.051721, 0.045019, 0.5225), ), ((0.067769, 
            0.050405, 0.5225), ), ((0.085837, 0.055927, 0.5225), ), ((0.110025, 
            0.062356, 0.5225), ), ((0.219046, 0.062356, 0.5225), ), ((0.257198, 
            0.060219, 0.5225), ), ((0.285821, 0.05724, 0.5225), ), ((0.314786, 
            0.053508, 0.5225), ), ((0.343821, 0.04914, 0.5225), ), ((0.372668, 
            0.044253, 0.5225), ), ((0.455158, 0.028129, 0.5225), ), ((0.480243, 
            0.022871, 0.5225), ), ((0.50352, 0.01795, 0.5225), ), ((0.524631, 0.013476, 
            0.5225), ), ((0.543235, 0.009508, 0.5225), ), ((0.559077, 0.006099, 
            0.5225), ), ((0.571946, 0.003362, 0.5225), ), ((0.0, 0.0, 0.4725), ), ((
            0.000608, 0.004195, 0.4725), ), ((0.003316, 0.011039, 0.4725), ), ((
            0.008319, 0.018178, 0.4725), ), ((0.015635, 0.02537, 0.4725), ), ((
            0.025205, 0.032397, 0.4725), ), ((0.036987, 0.039034, 0.4725), ), ((
            0.050958, 0.045117, 0.4725), ), ((0.067042, 0.050516, 0.4725), ), ((
            0.085149, 0.05605, 0.4725), ), ((0.10939, 0.062493, 0.4725), ), ((0.21865, 
            0.062493, 0.4725), ), ((0.256886, 0.060351, 0.4725), ), ((0.285572, 
            0.057366, 0.4725), ), ((0.3146, 0.053625, 0.4725), ), ((0.343699, 0.049247, 
            0.4725), ), ((0.372609, 0.04435, 0.4725), ), ((0.455279, 0.02819, 0.4725), 
            ), ((0.480419, 0.022922, 0.4725), ), ((0.503748, 0.017989, 0.4725), ), ((
            0.524905, 0.013505, 0.4725), ), ((0.543549, 0.009528, 0.4725), ), ((
            0.559426, 0.006112, 0.4725), ), ((0.572324, 0.003369, 0.4725), ), ((0.59, 
            0.0, 0.4725), ), ((0.571498, 0.000962, 0.4725), ), ((0.557497, 0.001274, 
            0.4725), ), ((0.539956, 0.001357, 0.4725), ), ((0.519188, 0.001115, 
            0.4725), ), ((0.495541, 0.000513, 0.4725), ), ((0.469386, -0.000466, 
            0.4725), ), ((0.379978, -0.005469, 0.4725), ), ((0.347952, -0.007641, 
            0.4725), ), ((0.315514, -0.009912, 0.4725), ), ((0.283076, -0.012178, 
            0.4725), ), ((0.251027, -0.015635, 0.4725), ), ((0.21865, -0.019629, 
            0.4725), ), ((0.10939, -0.019629, 0.4725), ), ((0.084535, -0.019582, 
            0.4725), ), ((0.063708, -0.018756, 0.4725), ), ((0.045489, -0.017393, 
            0.4725), ), ((0.030084, -0.015476, 0.4725), ), ((0.017647, -0.01301, 
            0.4725), ), ((0.008343, -0.010036, 0.4725), ), ((0.002508, -0.006579, 
            0.4725), ), ((0.59, 0.0, 0.415), ), ((0.571498, 0.000962, 0.415), ), ((
            0.557497, 0.001274, 0.415), ), ((0.539956, 0.001357, 0.415), ), ((0.519188, 
            0.001115, 0.415), ), ((0.495541, 0.000513, 0.415), ), ((0.469386, 
            -0.000466, 0.415), ), ((0.379978, -0.005469, 0.415), ), ((0.347952, 
            -0.007641, 0.415), ), ((0.315514, -0.009912, 0.415), ), ((0.283076, 
            -0.012178, 0.415), ), ((0.251027, -0.015635, 0.415), ), ((0.21865, 
            -0.019629, 0.415), ), ((0.10939, -0.019629, 0.415), ), ((0.084535, 
            -0.019582, 0.415), ), ((0.063708, -0.018756, 0.415), ), ((0.045489, 
            -0.017393, 0.415), ), ((0.030084, -0.015476, 0.415), ), ((0.017647, 
            -0.01301, 0.415), ), ((0.008343, -0.010036, 0.415), ), ((0.002508, 
            -0.006579, 0.415), ), ((0.0, 0.0, 0.415), ), ((0.000608, 0.004195, 0.415), 
            ), ((0.003316, 0.011039, 0.415), ), ((0.008319, 0.018178, 0.415), ), ((
            0.015635, 0.02537, 0.415), ), ((0.025205, 0.032397, 0.415), ), ((0.036987, 
            0.039034, 0.415), ), ((0.050958, 0.045117, 0.415), ), ((0.067042, 0.050516, 
            0.415), ), ((0.085149, 0.05605, 0.415), ), ((0.10939, 0.062493, 0.415), ), 
            ((0.21865, 0.062493, 0.415), ), ((0.256886, 0.060351, 0.415), ), ((
            0.285572, 0.057366, 0.415), ), ((0.3146, 0.053625, 0.415), ), ((0.343699, 
            0.049247, 0.415), ), ((0.372609, 0.04435, 0.415), ), ((0.455279, 0.02819, 
            0.415), ), ((0.480419, 0.022922, 0.415), ), ((0.503748, 0.017989, 0.415), 
            ), ((0.524905, 0.013505, 0.415), ), ((0.543549, 0.009528, 0.415), ), ((
            0.559426, 0.006112, 0.415), ), ((0.572324, 0.003369, 0.415), ), ((0.0, 0.0, 
            0.325), ), ((0.000608, 0.004195, 0.325), ), ((0.003316, 0.011039, 0.325), 
            ), ((0.008319, 0.018178, 0.325), ), ((0.015635, 0.02537, 0.325), ), ((
            0.025205, 0.032397, 0.325), ), ((0.036987, 0.039034, 0.325), ), ((0.050958, 
            0.045117, 0.325), ), ((0.067042, 0.050516, 0.325), ), ((0.085149, 0.05605, 
            0.325), ), ((0.10939, 0.062493, 0.325), ), ((0.21865, 0.062493, 0.325), ), 
            ((0.256886, 0.060351, 0.325), ), ((0.285572, 0.057366, 0.325), ), ((0.3146, 
            0.053625, 0.325), ), ((0.343699, 0.049247, 0.325), ), ((0.372609, 0.04435, 
            0.325), ), ((0.455279, 0.02819, 0.325), ), ((0.480419, 0.022922, 0.325), ), 
            ((0.503748, 0.017989, 0.325), ), ((0.524905, 0.013505, 0.325), ), ((
            0.543549, 0.009528, 0.325), ), ((0.559426, 0.006112, 0.325), ), ((0.572324, 
            0.003369, 0.325), ), ((0.59, 0.0, 0.325), ), ((0.571498, 0.000962, 0.325), 
            ), ((0.557497, 0.001274, 0.325), ), ((0.539956, 0.001357, 0.325), ), ((
            0.519188, 0.001115, 0.325), ), ((0.495541, 0.000513, 0.325), ), ((0.469386, 
            -0.000466, 0.325), ), ((0.379978, -0.005469, 0.325), ), ((0.347952, 
            -0.007641, 0.325), ), ((0.315514, -0.009912, 0.325), ), ((0.283076, 
            -0.012178, 0.325), ), ((0.251027, -0.015635, 0.325), ), ((0.21865, 
            -0.019629, 0.325), ), ((0.10939, -0.019629, 0.325), ), ((0.084535, 
            -0.019582, 0.325), ), ((0.063708, -0.018756, 0.325), ), ((0.045489, 
            -0.017393, 0.325), ), ((0.030084, -0.015476, 0.325), ), ((0.017647, 
            -0.01301, 0.325), ), ((0.008343, -0.010036, 0.325), ), ((0.002508, 
            -0.006579, 0.325), ), ((0.59, 0.0, 0.225), ), ((0.571498, 0.000962, 0.225), 
            ), ((0.557497, 0.001274, 0.225), ), ((0.539956, 0.001357, 0.225), ), ((
            0.519188, 0.001115, 0.225), ), )+\
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.495541, 0.000513, 
            0.225), ), ((0.469386, -0.000466, 0.225), ), ((0.379978, -0.005469, 0.225), 
            ), ((0.347952, -0.007641, 0.225), ), ((0.315514, -0.009912, 0.225), ), ((
            0.283076, -0.012178, 0.225), ), ((0.251027, -0.015635, 0.225), ), ((
            0.21865, -0.019629, 0.225), ), ((0.10939, -0.019629, 0.225), ), ((0.084535, 
            -0.019582, 0.225), ), ((0.063708, -0.018756, 0.225), ), ((0.045489, 
            -0.017393, 0.225), ), ((0.030084, -0.015476, 0.225), ), ((0.017647, 
            -0.01301, 0.225), ), ((0.008343, -0.010036, 0.225), ), ((0.002508, 
            -0.006579, 0.225), ), ((0.0, 0.0, 0.225), ), ((0.000608, 0.004195, 0.225), 
            ), ((0.003316, 0.011039, 0.225), ), ((0.008319, 0.018178, 0.225), ), ((
            0.015635, 0.02537, 0.225), ), ((0.025205, 0.032397, 0.225), ), ((0.036987, 
            0.039034, 0.225), ), ((0.050958, 0.045117, 0.225), ), ((0.067042, 0.050516, 
            0.225), ), ((0.085149, 0.05605, 0.225), ), ((0.10939, 0.062493, 0.225), ), 
            ((0.21865, 0.062493, 0.225), ), ((0.256886, 0.060351, 0.225), ), ((
            0.285572, 0.057366, 0.225), ), ((0.3146, 0.053625, 0.225), ), ((0.343699, 
            0.049247, 0.225), ), ((0.372609, 0.04435, 0.225), ), ((0.455279, 0.02819, 
            0.225), ), ((0.480419, 0.022922, 0.225), ), ((0.503748, 0.017989, 0.225), 
            ), ((0.524905, 0.013505, 0.225), ), ((0.543549, 0.009528, 0.225), ), ((
            0.559426, 0.006112, 0.225), ), ((0.572324, 0.003369, 0.225), ), ((0.0, 0.0, 
            0.125), ), ((0.000608, 0.004195, 0.125), ), ((0.003316, 0.011039, 0.125), 
            ), ((0.008319, 0.018178, 0.125), ), ((0.015635, 0.02537, 0.125), ), ((
            0.025205, 0.032397, 0.125), ), ((0.036987, 0.039034, 0.125), ), ((0.050958, 
            0.045117, 0.125), ), ((0.067042, 0.050516, 0.125), ), ((0.085149, 0.05605, 
            0.125), ), ((0.10939, 0.062493, 0.125), ), ((0.21865, 0.062493, 0.125), ), 
            ((0.256886, 0.060351, 0.125), ), ((0.285572, 0.057366, 0.125), ), ((0.3146, 
            0.053625, 0.125), ), ((0.343699, 0.049247, 0.125), ), ((0.372609, 0.04435, 
            0.125), ), ((0.455279, 0.02819, 0.125), ), ((0.480419, 0.022922, 0.125), ), 
            ((0.503748, 0.017989, 0.125), ), ((0.524905, 0.013505, 0.125), ), ((
            0.543549, 0.009528, 0.125), ), ((0.559426, 0.006112, 0.125), ), ((0.572324, 
            0.003369, 0.125), ), ((0.59, 0.0, 0.125), ), ((0.571498, 0.000962, 0.125), 
            ), ((0.557497, 0.001274, 0.125), ), ((0.539956, 0.001357, 0.125), ), ((
            0.519188, 0.001115, 0.125), ), ((0.495541, 0.000513, 0.125), ), ((0.469386, 
            -0.000466, 0.125), ), ((0.379978, -0.005469, 0.125), ), ((0.347952, 
            -0.007641, 0.125), ), ((0.315514, -0.009912, 0.125), ), ((0.283076, 
            -0.012178, 0.125), ), ((0.251027, -0.015635, 0.125), ), ((0.21865, 
            -0.019629, 0.125), ), ((0.10939, -0.019629, 0.125), ), ((0.084535, 
            -0.019582, 0.125), ), ((0.063708, -0.018756, 0.125), ), ((0.045489, 
            -0.017393, 0.125), ), ((0.030084, -0.015476, 0.125), ), ((0.017647, 
            -0.01301, 0.125), ), ((0.008343, -0.010036, 0.125), ), ((0.002508, 
            -0.006579, 0.125), ), ((0.59, 0.0, 0.04375), ), ((0.008319, 0.018178, 
            0.04375), ), ((0.015635, 0.02537, 0.04375), ), ((0.025205, 0.032397, 
            0.04375), ), ((0.036987, 0.039034, 0.04375), ), ((0.050958, 0.045117, 
            0.04375), ), ((0.067042, 0.050516, 0.04375), ), ((0.085149, 0.05605, 
            0.04375), ), ((0.10939, 0.062493, 0.04375), ), ((0.21865, 0.062493, 
            0.04375), ), ((0.256886, 0.060351, 0.04375), ), ((0.285572, 0.057366, 
            0.04375), ), ((0.3146, 0.053625, 0.04375), ), ((0.343699, 0.049247, 
            0.04375), ), ((0.372609, 0.04435, 0.04375), ), ((0.480419, 0.022922, 
            0.04375), ), ((0.503748, 0.017989, 0.04375), ), ((0.524905, 0.013505, 
            0.04375), ), ((0.543549, 0.009528, 0.04375), ), ((0.559426, 0.006112, 
            0.04375), ), ((0.572324, 0.003369, 0.04375), ), ), name='Set-34')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='skin-root-NS', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin-root-NS-Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['skin-root-NS'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'skin-root-NS-Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['skin-root-NS'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin-root-NS-Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['skin-root-NS'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin-root-NS-Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['skin-root-NS'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'skin-root-NS-Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['skin-root-NS'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'skin-root-NS-Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['skin-root-NS'], suppressed=False, 
            thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['skin-root-NS'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['skin-root-NS'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-34'], stackDirection=STACK_3)

        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.16402, -0.019629, 
            0.04375), ), ((0.1667, -0.019406, 0.575), ), ((0.164535, -0.019586, 
            0.5225), ), ((0.16402, -0.019629, 0.4725), ), ((0.16402, -0.019629, 0.415), 
            ), ((0.16402, -0.019629, 0.325), ), ((0.16402, -0.019629, 0.225), ), ((
            0.16402, -0.019629, 0.125), ), ), name='Set-35')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='SCB-NS', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-NS Ply-1', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-NS Ply-2', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='SCB-NS Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], suppressed=False, 
            thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-NS Ply-4', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCB-NS Ply-5', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-NS Ply-6', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-NS Ply-7', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCB-NS Ply-8', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-NS Ply-9', region=mdb.models['Model-1'].parts['Wing'].sets['SCB-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-NS'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=mdb.models['Model-1'].parts['Wing'].surfaces['SCB-NS'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-35'], stackDirection=STACK_3)

        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.1667, 0.061782, 
            0.575), ), ((0.164535, 0.062356, 0.5225), ), ((0.16402, 0.062493, 0.4725), 
            ), ((0.16402, 0.062493, 0.415), ), ((0.16402, 0.062493, 0.325), ), ((
            0.16402, 0.062493, 0.225), ), ((0.16402, 0.062493, 0.125), ), ((0.16402, 
            0.062493, 0.04375), ), ), name='Set-36')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='SCT-NS', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-NS Ply-1', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-NS Ply-2', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='SCT-NS Ply-3', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], suppressed=False, 
            thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-NS Ply-4', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCT-NS Ply-5', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-NS Ply-6', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-NS Ply-7', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCT-NS Ply-8', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-NS Ply-9', region=mdb.models['Model-1'].parts['Wing'].sets['SCT-NS'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-NS'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=mdb.models['Model-1'].parts['Wing'].surfaces['SCT-NS'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-36'], stackDirection=STACK_3)

        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.16402, 0.062493, 
            0.00625), )), name='Set-37')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='SCT-Singularity', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-Singularity Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-Singularity Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='SCT-Singularity Ply-3', 
            region=mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], 
            suppressed=False, thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-Singularity Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCT-Singularity Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-Singularity Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-Singularity Ply-7', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCT-Singularity Ply-8', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCT-Singularity Ply-9', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCT-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCT-Singularity'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['SCT-Singularity'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-37'], stackDirection=STACK_3)
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.18902, -0.019629, 
            0.00625), ), ((0.16402, -0.019629, 0.00625), ), ((0.13902, -0.019629, 
            0.00625), ), ), name='Set-38')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='SCB-Singularity', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-Singularity Ply-1', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon Fibre-200GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-Singularity Ply-2', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=0.00032, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Carbon UD Tape', numIntPoints=3, orientationType=
            SPECIFY_ORIENT, orientationValue=0.0, plyName='SCB-Singularity Ply-3', 
            region=mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], 
            suppressed=False, thickness=0.0003, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-Singularity Ply-4', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCB-Singularity Ply-5', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-Singularity Ply-6', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-Singularity Ply-7', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName=
            'SCB-Singularity Ply-8', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName=
            'SCB-Singularity Ply-9', region=
            mdb.models['Model-1'].parts['Wing'].sets['SCB-Singularity'], suppressed=
            False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['SCB-Singularity'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['SCB-Singularity'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-38'], stackDirection=STACK_3)
        mdb.models['Model-1'].parts['Wing'].Set(edges=
            mdb.models['Model-1'].parts['Wing'].edges.findAt(((0.008319, 0.018178, 
            0.00625), ), ((0.003316, 0.011039, 0.00625), ), ((0.015635, 0.02537, 
            0.00625), ), ((0.025205, 0.032397, 0.00625), ), ((0.036987, 0.039034, 
            0.00625), ), ((0.050958, 0.045117, 0.00625), ), ((0.067042, 0.050516, 
            0.00625), ), ((0.085149, 0.05605, 0.00625), ), ((0.10939, 0.062493, 
            0.00625), ), ((0.13902, 0.062493, 0.00625), ), ((0.18902, 0.062493, 
            0.00625), ), ((0.21865, 0.062493, 0.00625), ), ((0.256886, 0.060351, 
            0.00625), ), ((0.285572, 0.057366, 0.00625), ), ((0.3146, 0.053625, 
            0.00625), ), ((0.343699, 0.049247, 0.00625), ), ((0.372609, 0.04435, 
            0.00625), ), ((0.435, 0.032154, 0.00625), ), ((0.480419, 0.022922, 
            0.00625), ), ((0.455279, 0.02819, 0.00625), ), ((0.503748, 0.017989, 
            0.00625), ), ((0.524905, 0.013505, 0.00625), ), ((0.543549, 0.009528, 
            0.00625), ), ((0.559426, 0.006112, 0.00625), ), ((0.572324, 0.003369, 
            0.00625), ), ((0.59, 0.0, 0.00625), ), ((0.571498, 0.000962, 0.00625), ), (
            (0.557497, 0.001274, 0.00625), ), ((0.539956, 0.001357, 0.00625), ), ((
            0.519188, 0.001115, 0.00625), ), ((0.495541, 0.000513, 0.00625), ), ((
            0.469386, -0.000466, 0.00625), ), ((0.435, -0.00239, 0.00625), ), ((
            0.379978, -0.005469, 0.00625), ), ((0.347952, -0.007641, 0.00625), ), ((
            0.315514, -0.009912, 0.00625), ), ((0.283076, -0.012178, 0.00625), ), ((
            0.251027, -0.015635, 0.00625), ), ((0.21865, -0.019629, 0.00625), ), ((
            0.18902, -0.019629, 0.00625), ), ((0.13902, -0.019629, 0.00625), ), ((
            0.10939, -0.019629, 0.00625), ), ((0.084535, -0.019582, 0.00625), ), ((
            0.063708, -0.018756, 0.00625), ), ((0.045489, -0.017393, 0.00625), ), ((
            0.030084, -0.015476, 0.00625), ), ((0.017647, -0.01301, 0.00625), ), ((
            0.008343, -0.010036, 0.00625), ), ((0.002508, -0.006579, 0.00625), ), ((
            0.0, 0.0, 0.00625), ), ((0.000608, 0.004195, 0.00625), ), ), name='Set-39')
        mdb.models['Model-1'].parts['Wing'].CompositeLayup(description='', elementType=
            SHELL, name='Singularity-skin', offsetType=MIDDLE_SURFACE, symmetric=False, 
            thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].Section(
            integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
            temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName='Ply-1', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Singularity-skin'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName='Ply-2', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Singularity-skin'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName='Ply-3', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Singularity-skin'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName='Ply-4', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Singularity-skin'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=0.0, plyName='Ply-5', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Singularity-skin'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].CompositePly(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, material='Glass Fibre-86GSM', numIntPoints=3, 
            orientationType=SPECIFY_ORIENT, orientationValue=45.0, plyName='Ply-6', 
            region=mdb.models['Model-1'].parts['Wing'].sets['Singularity-skin'], 
            suppressed=False, thickness=8.3e-05, thicknessType=SPECIFY_THICKNESS)
        mdb.models['Model-1'].parts['Wing'].compositeLayups['Singularity-skin'].ReferenceOrientation(
            additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
            , axis=AXIS_3, flipNormalDirection=False, flipPrimaryDirection=False, 
            localCsys=None, normalAxisDefinition=SURFACE, normalAxisDirection=AXIS_3, 
            normalAxisRegion=
            mdb.models['Model-1'].parts['Wing'].surfaces['Singularity-skin'], 
            orientationType=DISCRETE, primaryAxisDefinition=EDGE, primaryAxisDirection=
            AXIS_1, primaryAxisRegion=
            mdb.models['Model-1'].parts['Wing'].sets['Set-39'], stackDirection=STACK_3)

                   
    assigncompositelayup1()


    #code to delete plies





    field_output_list=[('base rib','Wing-1.base rib'),
        ('rear spar after boom','Wing-1.rear spar after boom'),
        ('rear spar before boom','Wing-1.rear spar before boom'),
        ('ribs','Wing-1.ribs'),
        ('skin mid','Wing-1.skin mid'),
        ('skin tip','Wing-1.skin tip'),
        ('spar cap bottom mid','Wing-1.spar cap bottom mid'),
        ('spar cap bottom tip','Wing-1.spar cap bottom tip'),
        ('spar cap top mid','Wing-1.spar cap top mid'),
        ('spar cap top tip','Wing-1.spar cap top tip'),
        ('spar web after boom','Wing-1.spar web after boom'),
        ('spar web before boom','Wing-1.spar web before boom'),
        ('SCB-NS','Wing-1.SCB-NS'),
        ('SCT-NS','Wing-1.SCT-NS'),
        ('skin-root-NS','Wing-1.skin-root-NS')]


    def fieldoutputs(field_output_name,lay_up_name):
        mdb.models['Model-1'].FieldOutputRequest(createStepName='Step-1', 
            layupLocationMethod=SPECIFIED, layupNames=(lay_up_name, ), name=
            field_output_name, outputAtPlyBottom=False, outputAtPlyMid=True, outputAtPlyTop=
            False, rebar=EXCLUDE, variables=('S', 'E', 'HSNFTCRT', 'HSNFCCRT', 
            'HSNMTCRT', 'HSNMCCRT','CFAILURE'))
    for i in field_output_list:
        fieldoutputs(i[0], i[1])
CreateStandardModel()
