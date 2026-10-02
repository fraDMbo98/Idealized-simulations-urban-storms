import math
import numpy as np
import globalVariables as GV
global modelParams, TrackerParams, modelData

def boxBounds(startValue, minimum, maximum, boxSize):
    if (startValue < minimum + boxSize):
        minValue = minimum
    else:
        minValue = startValue - boxSize
    if (startValue > maximum - boxSize):
        maxValue = maximum
    else:
        maxValue = startValue + boxSize
    return [minValue, maxValue]

#======================================================================
# found2Dupdraft
#
# Function that takes the point of a found 2D updraft and returns a 
# dictionary of calculations pertaining to that point.
#
# This function should be updated to reflect the user's interests

def found2Dupdraft(x, y, z):
    data = {}
    #print(GV.__dict__)
    #print(GV.modelData['U'])
    #print(dir(GV['modelData']))
    #print(GV['modelData']['U'])
    #print(GV.U)
    nearUMean = (GV.modelData['U'][x-1, y, z] + GV.modelData['U'][x+1, y, z] + GV.modelData['U'][x, y-1, z] + GV.modelData['U'][x, y+1, z] + \
                 GV.modelData['U'][x-1, y-1, z] + GV.modelData['U'][x-1, y+1, z] + GV.modelData['U'][x+1, y-1, z] + GV.modelData['U'][x+1, y+1, z]) / 8.0
    #data['nearUMean'] = nearUMean

    nearVMean = (GV.modelData['V'][x-1, y, z] + GV.modelData['V'][x+1, y, z] + GV.modelData['V'][x, y-1, z] + GV.modelData['V'][x, y+1, z] + \
                 GV.modelData['V'][x-1, y-1, z] + GV.modelData['V'][x-1, y+1, z] + GV.modelData['V'][x+1, y-1, z] + GV.modelData['V'][x+1, y+1, z]) / 8.0
    angle=np.pi/3
    #print(nearUMean,nearVMean)
    data['nearVMean'] = nearVMean#*np.cos(angle)+nearUMean*np.sin(angle)
    data['nearUMean'] = nearUMean#*np.cos(angle)-nearVMean*np.sin(angle)
    data['WInside'] = GV.modelData['W'][x, y, z]
    #print('rotation',data['nearVMean'],data['nearUMean'])
    #------------------------------------------
    # 3-box calculation
    #------------------------------------------
    xRange = boxBounds(x, 0, GV.modelParams['NX']-1, 3)
    minX = xRange[0]
    maxX = xRange[1]

    yRange = boxBounds(y, 0, GV.modelParams['NY']-1, 3)
    minY = yRange[0]
    maxY = yRange[1]

    totalPoints = 2 * (maxX - minX) + 2 * (maxY - minY)

    data['ThetaBox3'] =  (sum(GV.modelData['THETA'][minX:maxX+1, minY, z]) + \
                          sum(GV.modelData['THETA'][minX:maxX+1, maxY, z]) + \
                          sum(GV.modelData['THETA'][minX, minY:maxY+1, z]) + \
                          sum(GV.modelData['THETA'][maxX, minY:maxY+1, z])) / totalPoints

    #---------------------------------------
    # Other data calculations
    #---------------------------------------

    data['NearCoreIntQV'] = np.mean(GV.modelData['QV'][minX:maxX+1, minY:maxY+1, z])


    return data
