from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def average_variational_strains():

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

            strain11 = fileModify.variables["strain11"][:]
            strain22 = fileModify.variables["strain22"][:]
            strain12 = fileModify.variables["strain12"][:]

            strain11AvgVertex = np.zeros((nTimes,nVertices))
            strain22AvgVertex = np.zeros((nTimes,nVertices))
            strain12AvgVertex = np.zeros((nTimes,nVertices))

            strain11Avg = np.zeros((nTimes,nCells,maxEdges))
            strain22Avg = np.zeros((nTimes,nCells,maxEdges))
            strain12Avg = np.zeros((nTimes,nCells,maxEdges))

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

                            strain11avg = strain11avg + areaCell[iCell] * strain11[iTime,iCell,iVertexOnCell]
                            strain22avg = strain22avg + areaCell[iCell] * strain22[iTime,iCell,iVertexOnCell]
                            strain12avg = strain12avg + areaCell[iCell] * strain12[iTime,iCell,iVertexOnCell]
                            denominator = denominator + areaCell[iCell]

                    strain11AvgVertex[iTime,iVertex] = strain11avg / denominator
                    strain22AvgVertex[iTime,iVertex] = strain22avg / denominator
                    strain12AvgVertex[iTime,iVertex] = strain12avg / denominator

            for iCell in range(0,nCells):
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    iVertex = verticesOnCell[iCell,iVertexOnCell]
                    strain11Avg[iTime,iCell,iVertexOnCell] = strain11AvgVertex[iTime,iVertex]
                    strain22Avg[iTime,iCell,iVertexOnCell] = strain22AvgVertex[iTime,iVertex]
                    strain12Avg[iTime,iCell,iVertexOnCell] = strain12AvgVertex[iTime,iVertex]


            try:
                var = fileModify.createVariable("strain11AvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain11AvgVertex"]
            var[:] = strain11AvgVertex[:]

            try:
                var = fileModify.createVariable("strain22AvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain22AvgVertex"]
            var[:] = strain22AvgVertex[:]

            try:
                var = fileModify.createVariable("strain12AvgVertex","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain12AvgVertex"]
            var[:] = strain12AvgVertex[:]

            try:
                var = fileModify.createVariable("strain11Avg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["strain11Avg"]
            var[:] = strain11Avg[:]

            try:
                var = fileModify.createVariable("strain22Avg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["strain22Avg"]
            var[:] = strain22Avg[:]

            try:
                var = fileModify.createVariable("strain12Avg","d",dimensions=["Time","nCells","maxEdges"])
            except:
                var = fileModify.variables["strain12Avg"]
            var[:] = strain12Avg[:]

            fileModify.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    average_strains()
