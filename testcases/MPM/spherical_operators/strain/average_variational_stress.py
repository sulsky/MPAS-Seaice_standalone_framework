from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def average_variational_stress():

    gridSizes = [2562, 10242, 40962, 163842]

    operatorMethods = ["mpmvar"]

    for operatorMethod in operatorMethods:

        print("  Operator Method: ", operatorMethod)

        for gridSize in gridSizes:

            print("    Gridsize: ", gridSize)

            filenameModify = "./output_%s_%i/output.2000.nc" %(operatorMethod, gridSize)
            fileModify = Dataset(filenameModify,"a")

            nVertices = len(fileModify.dimensions["nVertices"])
            nCells = len(fileModify.dimensions["nCells"])
            vertexDegree = len(fileModify.dimensions["vertexDegree"])
            maxEdges = len(fileModify.dimensions["maxEdges"])
            nTimes = len(fileModify.dimensions["Time"])

            areaCell = fileModify.variables["areaCell"][:]

            nEdgesOnCell = fileModify.variables["nEdgesOnCell"][:]
            cellsOnVertex = fileModify.variables["cellsOnVertex"][:]
            cellVerticesAtVertex = fileModify.variables["cellVerticesAtVertex"][:]
            verticesOnCell = fileModify.variables["verticesOnCell"][:]

            cellsOnVertex[:] -= 1
            cellVerticesAtVertex[:] -= 1
            verticesOnCell[:] -= 1

            stress11var = fileModify.variables["stress11var"][:]
            stress22var = fileModify.variables["stress22var"][:]
            stress12var = fileModify.variables["stress12var"][:]

            stress11varAvgVertex = np.zeros((nTimes,nVertices))
            stress22varAvgVertex = np.zeros((nTimes,nVertices))
            stress12varAvgVertex = np.zeros((nTimes,nVertices))

            stress11varAvg = np.zeros((nTimes,nCells,maxEdges))
            stress22varAvg = np.zeros((nTimes,nCells,maxEdges))
            stress12varAvg = np.zeros((nTimes,nCells,maxEdges))

            for iTime in range(0, nTimes):

                for iVertex in range(0, nVertices):

                    stress11avg = 0.0
                    stress22avg = 0.0
                    stress12avg = 0.0
                    denominator = 0.0

                    for iVertexDegree in range(0,vertexDegree):

                        iCell = cellsOnVertex[iVertex,iVertexDegree]

                        if (iCell <= nCells-1):

                            iVertexOnCell = cellVerticesAtVertex[iVertex,iVertexDegree]

                            stress11avg = stress11avg + areaCell[iCell] * stress11var[iTime,iCell,iVertexOnCell]
                            stress22avg = stress22avg + areaCell[iCell] * stress22var[iTime,iCell,iVertexOnCell]
                            stress12avg = stress12avg + areaCell[iCell] * stress12var[iTime,iCell,iVertexOnCell]
                            denominator = denominator + areaCell[iCell]

                    stress11varAvgVertex[iTime,iVertex] = stress11avg / denominator
                    stress22varAvgVertex[iTime,iVertex] = stress22avg / denominator
                    stress12varAvgVertex[iTime,iVertex] = stress12avg / denominator

            for iCell in range(0,nCells):
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    iVertex = verticesOnCell[iCell,iVertexOnCell]
                    stress11varAvg[iTime,iCell,iVertexOnCell] = stress11varAvgVertex[iTime,iVertex]
                    stress22varAvg[iTime,iCell,iVertexOnCell] = stress22varAvgVertex[iTime,iVertex]
                    stress12varAvg[iTime,iCell,iVertexOnCell] = stress12varAvgVertex[iTime,iVertex]


            try:
                var = fileModify.createVariable("stress11varAvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["stress11varAvgVertex"]
            var[:] = stress11varAvgVertex[:]

            try:
                var = fileModify.createVariable("stress22varAvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["stress22varAvgVertex"]
            var[:] = stress22varAvgVertex[:]

            try:
                var = fileModify.createVariable("stress12varAvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["stress12varAvgVertex"]
            var[:] = stress12varAvgVertex[:]

            try:
                var = fileModify.createVariable("stress11varAvg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["stress11varAvg"]
            var[:] = stress11varAvg[:]

            try:
                var = fileModify.createVariable("stress22varAvg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["stress22varAvg"]
            var[:] = stress22varAvg[:]

            try:
                var = fileModify.createVariable("stress12varAvg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["stress12varAvg"]
            var[:] = stress12varAvg[:]

            fileModify.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    average_variational_stress()
