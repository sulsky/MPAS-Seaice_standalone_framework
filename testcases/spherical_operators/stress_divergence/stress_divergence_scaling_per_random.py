from stress_divergence_scaling import scaling_lines, get_norm, get_resolution
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def stress_divergence_scaling_per_random():

    meshScales = [0.0,0.01,0.02,0.04]
    method = "pwl_alt"
    resolutions = [2562,10242,40962,163842]
    latitudeLimit = 20.0

    # plot options
    legendLabels = {"wachspress":"Wachs.",
                    "pwl":"PWL",
                    "weak":"Weak",
                    "wachspress_alt":"Wachs. alt",
                    "pwl_alt":"PWL alt."}

    lineColors = {"wachspress":"black",
                  "pwl":"grey",
                  "weak":"red",
                  "wachspress_alt":"black",
                  "pwl_alt":"grey"}

    lineStyles = {"wachspress":"solid",
                  "pwl":"solid",
                  "weak":"solid",
                  "wachspress_alt":"dashed",
                  "pwl_alt":"dashed"}

    markers = {"wachspress":"+",
               "pwl":"x",
               "weak":"^",
               "wachspress_alt":"+",
               "pwl_alt":"x"}

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


    # scaling lines
    xMin = 2.5e-2
    xMax = 5e-2
    yMin_L2 = 2.5e-3
    yMin_Linf = 2e-1

    fig, axes = plt.subplots(2,2,figsize=(15*cm,15*cm))

    scaling_lines(axes[0,0], xMin, xMax, yMin_L2)
    scaling_lines(axes[0,1], xMin, xMax, yMin_L2)
    scaling_lines(axes[1,0], xMin, xMax, yMin_Linf)
    scaling_lines(axes[1,1], xMin, xMax, yMin_Linf)


    for meshScale in meshScales:

        print("Mesh scale: ", meshScale)

        testName = "random_%f" %(meshScale)

        x = []
        y_L2U = []
        y_L2V = []
        y_LinfU = []
        y_LinfV = []

        for resolution in resolutions:

            filename = "./output_%s_%s_%i/output.2000.nc" %(testName,method,resolution)
            filenameIC = "./ic_%s.%i.nc" %(testName,resolution)

            print(filename, filenameIC)

            L2_normU, L2_normV, Linf_normU, Linf_normV = get_norm(filenameIC, filename, latitudeLimit)

            x.append(get_resolution(filename, latitudeLimit))
            y_L2U.append(L2_normU)
            y_L2V.append(L2_normV)
            y_LinfU.append(Linf_normU)
            y_LinfV.append(Linf_normV)

        axes[0,0].loglog(x,y_L2U, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])
        axes[0,1].loglog(x,y_L2V, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])
        axes[1,0].loglog(x,y_LinfU, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])
        axes[1,1].loglog(x,y_LinfV, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])

    axes[0,0].legend(frameon=False, loc=2, fontsize=8, handlelength=4)
    axes[0,0].set_xlabel("Grid resolution")
    axes[0,0].set_ylabel(r"$L_2$ error norm")
    axes[0,0].set_title(r'$(\nabla \cdot \sigma)_u$')
    axes[0,0].set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes[0,0].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[0,0].set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes[0,0].set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)

    axes[0,1].legend(frameon=False, loc=2, fontsize=8, handlelength=4)
    axes[0,1].set_xlabel("Grid resolution")
    axes[0,1].set_ylabel(r"$L_2$ error norm")
    axes[0,1].set_title(r'$(\nabla \cdot \sigma)_v$')
    axes[0,1].set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes[0,1].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[0,1].set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes[0,1].set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)

    axes[1,0].legend(frameon=False, loc=2, fontsize=8, handlelength=4)
    axes[1,0].set_xlabel("Grid resolution")
    axes[1,0].set_ylabel(r"$L_\infty$ error norm")
    axes[1,0].set_title(r'$(\nabla \cdot \sigma)_u$')
    axes[1,0].set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes[1,0].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[1,0].set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes[1,0].set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)

    axes[1,1].legend(frameon=False, loc=2, fontsize=8, handlelength=4)
    axes[1,1].set_xlabel("Grid resolution")
    axes[1,1].set_ylabel(r"$L_\infty$ error norm")
    axes[1,1].set_title(r'$(\nabla \cdot \sigma)_v$')
    axes[1,1].set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes[1,1].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes[1,1].set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes[1,1].set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_divergence_scaling_per_random_%s.png" %(method),dpi=400)



#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_scaling_per_random()
