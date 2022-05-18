from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def average_variational_strains():

    gridSizes = [2562, 10242, 40962, 163842]

    operatorMethods = ["wachspress","pwl","weakwachs"]

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

            strain11var = fileModify.variables["strain11var"][:]
            strain22var = fileModify.variables["strain22var"][:]
            strain12var = fileModify.variables["strain12var"][:]

            strain11varAvgVertex = np.zeros((nTimes,nVertices))
            strain22varAvgVertex = np.zeros((nTimes,nVertices))
            strain12varAvgVertex = np.zeros((nTimes,nVertices))

            strain11varAvg = np.zeros((nTimes,nCells,maxEdges))
            strain22varAvg = np.zeros((nTimes,nCells,maxEdges))
            strain12varAvg = np.zeros((nTimes,nCells,maxEdges))

            for iTime in range(0, nTimes):

                for iVertex in range(0, nVertices):

                    strain11avg = 0.0
                    strain22avg = 0.0
                    strain12avg = 0.0
                    denominator = 0.0

                    for iVertexDegree in range(0,vertexDegree):

                        iCell = cellsOnVertex[iVertex,iVertexDegree]

                        if (iCell <= nCells-1):

                            iVertexOnCell = cellVerticesAtVertex[iVertex,iVertexDegree]

                            strain11avg = strain11avg + areaCell[iCell] * strain11var[iTime,iCell,iVertexOnCell]
                            strain22avg = strain22avg + areaCell[iCell] * strain22var[iTime,iCell,iVertexOnCell]
                            strain12avg = strain12avg + areaCell[iCell] * strain12var[iTime,iCell,iVertexOnCell]
                            denominator = denominator + areaCell[iCell]

                    strain11varAvgVertex[iTime,iVertex] = strain11avg / denominator
                    strain22varAvgVertex[iTime,iVertex] = strain22avg / denominator
                    strain12varAvgVertex[iTime,iVertex] = strain12avg / denominator

            for iCell in range(0,nCells):
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    iVertex = verticesOnCell[iCell,iVertexOnCell]
                    strain11varAvg[iTime,iCell,iVertexOnCell] = strain11varAvgVertex[iTime,iVertex]
                    strain22varAvg[iTime,iCell,iVertexOnCell] = strain22varAvgVertex[iTime,iVertex]
                    strain12varAvg[iTime,iCell,iVertexOnCell] = strain12varAvgVertex[iTime,iVertex]


            try:
                var = fileModify.createVariable("strain11varAvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain11varAvgVertex"]
            var[:] = strain11varAvgVertex[:]

            try:
                var = fileModify.createVariable("strain22varAvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain22varAvgVertex"]
            var[:] = strain22varAvgVertex[:]

            try:
                var = fileModify.createVariable("strain12varAvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain12varAvgVertex"]
            var[:] = strain12varAvgVertex[:]

            try:
                var = fileModify.createVariable("strain11varAvg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["strain11varAvg"]
            var[:] = strain11varAvg[:]

            try:
                var = fileModify.createVariable("strain22varAvg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["strain22varAvg"]
            var[:] = strain22varAvg[:]

            try:
                var = fileModify.createVariable("strain12varAvg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["strain12varAvg"]
            var[:] = strain12varAvg[:]

            fileModify.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    average_variational_strains()
