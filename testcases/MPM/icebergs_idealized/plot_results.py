import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def plot_figure(testType, uColor, vColor):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testType)))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    latIcebergs = []
    lonIcebergs = []

    uVelocityIcebergs = []
    vVelocityIcebergs = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        latIceberg = filein.variables["latIceberg"][0,0]
        lonIceberg = filein.variables["lonIceberg"][0,0]

        uVelocityIceberg = filein.variables["uVelocityIceberg"][0,0]
        vVelocityIceberg = filein.variables["vVelocityIceberg"][0,0]

        filein.close()

        latIcebergs.append(latIceberg)
        lonIcebergs.append(lonIceberg)

        uVelocityIcebergs.append(uVelocityIceberg)
        vVelocityIcebergs.append(vVelocityIceberg)

    fig, axis = plt.subplots()
    axis2 = axis.twinx()

    axis.plot(times, lonIcebergs, color=uColor, linestyle="dashed", label="lonIceberg")
    axis.plot(times, latIcebergs, color=vColor, linestyle="dashed", label="latIceberg")

    axis2.plot(times, uVelocityIcebergs, color=uColor, linestyle="solid", label="uVelocityIceberg")
    axis2.plot(times, vVelocityIcebergs, color=vColor, linestyle="solid", label="vVelocityIceberg")

    axis.legend(loc="upper left")
    axis2.legend(loc="lower right")

    axis.set_title("Position/velocity %s" %(testType))

    axis.set_xlabel("Time (days)")

    axis.set_ylabel("Position (rads)")
    axis2.set_ylabel("Velocity (m/s)")

    plt.tight_layout()

    plt.savefig("results_%s.png" %(testType),dpi=300)

#-------------------------------------------------------------------------------

def plot_coriolis():

    filenames = sorted(glob.glob("./output_coriolis/icebergs_output.*"))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    xs = []
    ys = []
    zs = []
    vs = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        posnIBGeo = filein.variables["posnIBGeo"][0,0,:]

        uVelocityIceberg = filein.variables["uVelocityIceberg"][0,0]
        vVelocityIceberg = filein.variables["vVelocityIceberg"][0,0]
        v = sqrt(pow(uVelocityIceberg,2)+
                 pow(vVelocityIceberg,2))
        vs.append(v)

        filein.close()

        xs.append(posnIBGeo[0])
        ys.append(posnIBGeo[1])
        zs.append(posnIBGeo[2])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    vs = np.array(vs)

    # Create sphere parameterization
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Sphere of radius R centered at (cx, cy, cz)
    R = sqrt(xs[0]*xs[0]+ys[0]*ys[0]+zs[0]*zs[0])
    cx, cy, cz = 0.0, 0.0, 0.0

    x = cx + R * np.outer(np.cos(u), np.sin(v))
    y = cy + R * np.outer(np.sin(u), np.sin(v))
    z = cz + R * np.outer(np.ones_like(u), np.cos(v))

    fig, axis = plt.subplots(subplot_kw={"projection": "3d"})

    # Plot transparent surface
    axis.plot_surface(
        x, y, z,
        alpha=0.3,          # transparency (0=fully transparent, 1=opaque)
        linewidth=0,
        antialiased=True
    )

    lats = np.deg2rad(np.arange(-60, 61, 30))  # degrees → radians
    phi = np.linspace(0, 2*np.pi, 200)

    for lat in lats:
        z = cz + R * np.sin(lat)
        r_xy = R * np.cos(lat)
        x = cx + r_xy * np.cos(phi)
        y = cy + r_xy * np.sin(phi)
        axis.plot(x, y, z, color="grey", linewidth=0.4)

    lons = np.deg2rad(np.arange(0, 360, 30))
    theta = np.linspace(0, np.pi, 200)

    for lon in lons:
        x = cx + R * np.sin(theta) * np.cos(lon)
        y = cy + R * np.sin(theta) * np.sin(lon)
        z = cz + R * np.cos(theta)
        axis.plot(x, y, z, color="grey", linewidth=0.4)

    axis.plot(xs, ys, zs, linewidth=0.5)
    axis.view_init(elev=20, azim=45)

    axis.set_title("Coriolis trajectory")

    plt.tight_layout()

    plt.savefig("results_coriolis.png",dpi=300)

    fig, axis = plt.subplots()

    axis.plot(times, vs)

    axis.set_xlabel("Time (days)")

    axis.set_title("Coriolis speed")

    plt.tight_layout()

    plt.savefig("results_coriolis_v.png",dpi=300)

#-------------------------------------------------------------------------------

def plot_seaice_forces(testcase):

    filenames = sorted(glob.glob("./output_%s/icebergs_output.*" %(testcase)))

    filein = Dataset(filenames[0],"r")
    time1 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    filein = Dataset(filenames[1],"r")
    time2 = filein.variables["daysSinceStartOfSim"][0]
    filein.close()

    dt = time2-time1
    times = np.arange(0,len(filenames)) * dt

    fius = []
    fivs = []
    fsus = []
    fsvs = []

    for filename in filenames:

        filein = Dataset(filename,"r")

        nVertices = len(filein.dimensions["nVertices"])

        uSeaiceForceIceberg = filein.variables["uSeaiceForceIceberg"][0,0]
        vSeaiceForceIceberg = filein.variables["vSeaiceForceIceberg"][0,0]

        icebergForceU = filein.variables["icebergForceU"][0,:]
        icebergForceV = filein.variables["icebergForceV"][0,:]

        areaTriangle = filein.variables["areaTriangle"][:]

        filein.close()

        icebergForceUTotal = 0.0
        icebergForceVTotal = 0.0
        for iVertex in range(0,nVertices):
            icebergForceUTotal += icebergForceU[iVertex] * areaTriangle[iVertex]
            icebergForceVTotal += icebergForceV[iVertex] * areaTriangle[iVertex]

        fius.append(uSeaiceForceIceberg)
        fivs.append(vSeaiceForceIceberg)
        fsus.append(icebergForceUTotal)
        fsvs.append(icebergForceVTotal)

    fig, axis = plt.subplots()

    axis.plot(times, fius, label="u iceberg", marker="x")
    axis.plot(times, fivs, label="v iceberg", marker="+")
    axis.plot(times, fsus, label="u seaice", marker="x")
    axis.plot(times, fsvs, label="v seaice", marker="+")

    axis.set_title("Sea ice force")

    axis.legend()
    axis.set_xlabel("Time (days)")
    axis.set_ylabel("Force (N)")

    plt.tight_layout()

    plt.savefig("results_seaice_force_%s.png" %(testcase),dpi=300)

#-------------------------------------------------------------------------------

def plot_testcase(testcase):

    if (testcase == "air_drag"):
        plot_figure("air_drag","green","red")
    elif (testcase == "ocean_drag"):
        plot_figure("ocean_drag","green","red")
    elif (testcase == "surface_tilt"):
        plot_figure("surface_tilt","red","green")
    elif (testcase == "wave_radiation"):
        plot_figure("wave_radiation","green","red")
    elif (testcase == "coriolis"):
        plot_coriolis()
    elif (testcase == "seaice_force_low_conc"):
        plot_figure("seaice_force_low_conc","green","red")
        plot_seaice_forces("seaice_force_low_conc")
    elif (testcase == "seaice_force_mid_conc"):
        plot_figure("seaice_force_mid_conc","green","red")
        plot_seaice_forces("seaice_force_mid_conc")
    elif (testcase == "seaice_force_high_conc"):
        plot_figure("seaice_force_high_conc","green","red")
        plot_seaice_forces("seaice_force_high_conc")

#-------------------------------------------------------------------------------

def plot_results(testcase):

    if (testcase is not None):

        plot_testcase(testcase)

    else:

        testcases = ["air_drag",
                     "ocean_drag",
                     "surface_tilt",
                     "wave_radiation",
                     "coriolis",
                     "seaice_force_low_conc",
                     "seaice_force_mid_conc",
                     "seaice_force_high_conc"]

        for testcase in testcases:
            plot_testcase(testcase)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='testcase', default=None)

    args = parser.parse_args()

    plot_results(args.testcase)
