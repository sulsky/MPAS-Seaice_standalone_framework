from stress_divergence_scaling import get_use_vertex, get_norm_integral_triangle, get_norm, get_resolution, scaling_lines
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def stress_divergence_scaling_per_random():

    gridType = "hex"
    operatorMethod = "pwl"


    grids = ["0082x0094",
             "0164x0188",
             "0328x0376",
             "0656x0752"]

    meshScales = [0.0,0.01,0.02,0.05,0.1]


    # scaling lines
    xMin = 2e-3
    xMax = 3.5e-3
    yMin = {"L2":1.5e-3,
            "Linf":5e-1}

    # plot
    cm = 1/2.54  # centimeters in inches
    plt.rc('font', family='Times New Roman', size=8)
    plt.rc('mathtext',fontset="stix")
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


    fig, axes = plt.subplots(2, 2, figsize=(15*cm,15*cm))

    scaling_lines(axes[0,0], xMin, xMax, yMin["L2"])
    scaling_lines(axes[0,1], xMin, xMax, yMin["L2"])
    scaling_lines(axes[1,0], xMin, xMax, yMin["Linf"])
    scaling_lines(axes[1,1], xMin, xMax, yMin["Linf"])

    for meshScale in meshScales:

        print("Mesh scale: ", meshScale)

        testName = "random_%f" %(meshScale)

        x = []
        y_L2U = []
        y_L2V = []
        y_LinfU = []
        y_LinfV = []

        for grid in grids:

            print("  Grid: ", grid)

            filenameIC = "ic_%s_%s_%s.nc" %(testName, gridType, grid)
            filename = "output_%s_%s_%s_%s/output.2000.nc" %(testName, gridType, operatorMethod, grid)
            print("    ", filename, filenameIC)

            useVertex = get_use_vertex(filename)

            if (gridType == "hex"):
                #L2_normU, L2_normV, Linf_normU, Linf_normV = get_norm_integral_triangle(filenameIC, filename, useVertex)
                L2_normU, L2_normV, Linf_normU, Linf_normV = get_norm(filenameIC, filename, useVertex)
            elif (gridType == "quad"):
                #L2_normU, L2_normV, Linf_normU, Linf_normV = get_norm_integral_square(filenameIC, filename, useVertex)
                L2_normU, L2_normV, Linf_normU, Linf_normV = get_norm(filenameIC, filename, useVertex)

            x.append(get_resolution(filename, useVertex))
            y_L2U.append(L2_normU)
            y_L2V.append(L2_normV)
            y_LinfU.append(Linf_normU)
            y_LinfV.append(Linf_normV)

        axes[0,0].loglog(x, y_L2U, markersize=5.0, label="%f" %(meshScale))
        axes[0,1].loglog(x, y_L2V, markersize=5.0, label="%f" %(meshScale))
        axes[1,0].loglog(x, y_LinfU, markersize=5.0, label="%f" %(meshScale))
        axes[1,1].loglog(x, y_LinfV, markersize=5.0, label="%f" %(meshScale))

    axes[0,0].legend(frameon=False, loc=2, fontsize=8, handlelength=2)
    axes[0,0].set_xlabel(None)
    axes[0,0].set_ylabel(r"$L_2$ error norm")
    axes[0,0].set_title(r"(a) $(\nabla \cdot \sigma)_u$",loc="left")
    axes[0,0].set_xticks(ticks=[3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3],minor=True)
    axes[0,0].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[0,0].set_xticks(ticks=[2e-3,1e-2],minor=False)
    axes[0,0].set_xticklabels(labels=[r'$2\times 10^{-3}$',r'$10^{-2}$'],minor=False)

    axes[0,1].legend(frameon=False, loc=2, fontsize=8, handlelength=2)
    axes[0,1].set_xlabel(None)
    axes[0,1].set_ylabel(None)
    axes[0,1].set_title(r"(b) $(\nabla \cdot \sigma)_v$",loc="left")
    axes[0,1].set_xticks(ticks=[3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3],minor=True)
    axes[0,1].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[0,1].set_xticks(ticks=[2e-3,1e-2],minor=False)
    axes[0,1].set_xticklabels(labels=[r'$2\times 10^{-3}$',r'$10^{-2}$'],minor=False)

    axes[1,0].legend(frameon=False, loc=2, fontsize=8, handlelength=2)
    axes[1,0].set_xlabel("Grid resolution")
    axes[1,0].set_ylabel(r"$L_\infty$ error norm")
    axes[1,0].set_title(r"(c) $(\nabla \cdot \sigma)_u$",loc="left")
    axes[1,0].set_xticks(ticks=[3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3],minor=True)
    axes[1,0].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[1,0].set_xticks(ticks=[2e-3,1e-2],minor=False)
    axes[1,0].set_xticklabels(labels=[r'$2\times 10^{-3}$',r'$10^{-2}$'],minor=False)

    axes[1,1].legend(frameon=False, loc=2, fontsize=8, handlelength=2)
    axes[1,1].set_xlabel("Grid resolution")
    axes[1,1].set_ylabel(None)
    axes[1,1].set_title(r"(d) $(\nabla \cdot \sigma)_v$",loc="left")
    axes[1,1].set_xticks(ticks=[3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3],minor=True)
    axes[1,1].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[1,1].set_xticks(ticks=[2e-3,1e-2],minor=False)
    axes[1,1].set_xticklabels(labels=[r'$2\times 10^{-3}$',r'$10^{-2}$'],minor=False)


    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    filenameOut = "stress_divergence_scaling_per_random.png"
    plt.savefig(filenameOut,dpi=300)
    filenameOut = "stress_divergence_scaling_per_random.eps"
    plt.savefig(filenameOut)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_scaling_per_random()
