import matplotlib.pyplot as plt
def plt_metric(metrics):
    
    diams = []
    for metric in metrics:
        for m in metric: diams.append(m)

    weights_water = np.ones_like(diams)/float(len(diams))
    water_hist, water_edges, water_p = plt.hist(diams, weights=weights_water, bins=40, alpha=0.5, color='midnightblue',label='N1888+')
    plt.legend()
    plt.xlabel("Diameter Length")
    plt.ylabel("Probability")

    plt.savefig("diam_h3o_water.png")
