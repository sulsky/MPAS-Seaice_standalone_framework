from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def average_weak_strains():

    operatorMethod = "weak"

    gridSizes = [2562, 10242, 40962, 163842]

    for gridSize in gridSizes:

        print("    GridSize: ", gridSize)

        filenameModify = "./output_%s_%i/output.2000.nc" %(operatorMethod, gridSize)
        fileModify = Dataset(filenameModify,"a")

        nVertices = len(fileModify.dimensions["nVertices"])
        nCells = len(fileModify.dimensions["nCells"])
        maxEdges = len(fileModify.dimensions["maxEdges"])
        vertexDegree = len(fileModify.dimensions["vertexDegree"])
        nTimes = len(fileModify.dimensions["Time"])

        nEdgesOnCell = fileModify.variables["nEdgesOnCell"][:]

        areaCell = fileModify.variables["areaCell"][:]

        cellsOnVertex = fileModify.variables["cellsOnVertex"][:]
        verticesOnCell = fileModify.variables["verticesOnCell"][:]

        cellsOnVertex[:] -= 1
        verticesOnCell[:] -= 1

        strain11weak = fileModify.variables["strain11weak"][:]
        strain22weak = fileModify.variables["strain22weak"][:]
        strain12weak = fileModify.variables["strain12weak"][:]

        strain11weakAvgVertex = np.zeros((nTimes,nVertices))
        strain22weakAvgVertex = np.zeros((nTimes,nVertices))
        strain12weakAvgVertex = np.zeros((nTimes,nVertices))

        strain11weakAvg = np.zeros((nTimes,nCells,maxEdges))
        strain22weakAvg = np.zeros((nTimes,nCells,maxEdges))
        strain12weakAvg = np.zeros((nTimes,nCells,maxEdges))

        for iTime in range(0, nTimes):

            for iVertex in range(0, nVertices):

                strain11avg = 0.0
                strain22avg = 0.0
                strain12avg = 0.0
                denominator = 0.0

                for iCellOnVertex in range(0,vertexDegree):
                    iCell = cellsOnVertex[iVertex,iCellOnVertex]
                    if (iCell < nCells):
                        strain11avg += areaCell[iCell] * strain11weak[iTime,iCell]
                        strain22avg += areaCell[iCell] * strain22weak[iTime,iCell]
                        strain12avg += areaCell[iCell] * strain12weak[iTime,iCell]
                        denominator += areaCell[iCell]

                if (denominator > 0.0):
                    strain11weakAvgVertex[iTime,iVertex] = strain11avg / denominator
                    strain22weakAvgVertex[iTime,iVertex] = strain22avg / denominator
                    strain12weakAvgVertex[iTime,iVertex] = strain12avg / denominator

            for iCell in range(0,nCells):
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    iVertex = verticesOnCell[iCell,iVertexOnCell]
                    strain11weakAvg[iTime,iCell,iVertexOnCell] = strain11weakAvgVertex[iTime,iVertex]
                    strain22weakAvg[iTime,iCell,iVertexOnCell] = strain22weakAvgVertex[iTime,iVertex]
                    strain12weakAvg[iTime,iCell,iVertexOnCell] = strain12weakAvgVertex[iTime,iVertex]

        try:
            var = fileModify.createVariable("strain11weakAvgVertex","d",dimensions=["Time","nVertices"])
        except:
            var = fileModify.variables["strain11weakAvgVertex"]
        var[:] = strain11weakAvgVertex[:]

        try:
            var = fileModify.createVariable("strain22weakAvgVertex","d",dimensions=["Time","nVertices"])
        except:
            var = fileModify.variables["strain22weakAvgVertex"]
        var[:] = strain22weakAvgVertex[:]

        try:
            var = fileModify.createVariable("strain12weakAvgVertex","d",dimensions=["Time","nVertices"])
        except:
            var = fileModify.variables["strain12weakAvgVertex"]
        var[:] = strain12weakAvgVertex[:]

        try:
            var = fileModify.createVariable("strain11weakAvg","d",dimensions=["Time","nCells","maxEdges"])
        except:
            var = fileModify.variables["strain11weakAvg"]
        var[:] = strain11weakAvg[:]

        try:
            var = fileModify.createVariable("strain22weakAvg","d",dimensions=["Time","nCells","maxEdges"])
        except:
            var = fileModify.variables["strain22weakAvg"]
        var[:] = strain22weakAvg[:]

        try:
            var = fileModify.createVariable("strain12weakAvg","d",dimensions=["Time","nCells","maxEdges"])
        except:
            var = fileModify.variables["strain12weakAvg"]
        var[:] = strain12weakAvg[:]

        fileModify.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    average_weak_strains()
