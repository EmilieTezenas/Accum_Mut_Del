import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import math
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

#Parameters used to plot the heatmap
mu = 0.05
r = 0.05
s0 = 0.9
s = 0.1
h = 0.25
shom = s
shet = h*s
params1 = [mu, r, s0, shom, shet]


#Computing the matrix that contains either the value of the extinction probability or the extinction time, depending on the criticity of the branching process
def matrice_finale(n, params, R):

    matriceE = np.zeros((n+1,n+1))
    matriceP = np.zeros((n+1,n+1))
    nb_E = 0
    nb_P = 0

    mu = params[0]
    r = params[1]
    s0 = params[2]
    shom = params[3]
    shet = params[4]

    Echelle_axes = np.linspace(0,1,n+1)

    for k in range(n+1):
        omega = Echelle_axes[n-k]
        for j in range(k):
            xhet = Echelle_axes[j]

            #Computing the death rate to text whether the process is super or sub-critical
            d = s0 + shom*omega + shet*xhet
            taux = mu + R*0.75*r + d

            if taux <= 1:
                P = (mu + d)/(1-R*0.5*r)
                matriceP[k][j] = P
                nb_P = nb_P +1
            else:
                E = 1/(R*0.75*r-1)*np.log(1 - (1-R*0.75*r)/(mu+d))
                matriceE[k][j] = E
                nb_E = nb_E +1

    if nb_E == 0:
        minE = 0
    else:
        minE = np.min(matriceE[matriceE>0])

    if nb_P == 0:
        minP = 0
    else:
        minP = np.min(matriceP[matriceP>0])



    return(matriceE, matriceP, minE, minP)

#Function to rescale the matrix in order to plot colors only for the good range of values
def rescale_matrix(A, min, max):
    (n,m) = np.shape(A)

    A_rescale = np.zeros((n,m))

    for i in range(n):
        for j in range(m):
            if A[i][j] != 0:
                A_rescale[i][j] = (A[i][j] - 0.95*min)/(max-0.95*min)

    return A_rescale



#Plots the heatmap
def heatmap_E_P(n, params):

    #Computing the quantities for the recombining and non-recombining case
    matriceE0, matriceP0, minE0, minP0 = matrice_finale(n, params, 0)
    matriceE1, matriceP1, minE1, minP1 = matrice_finale(n, params, 1)

    #Computing the min and max, for the summary matrix and for the indications on the colorbar
    maxP = max(np.max(matriceP0), np.max(matriceP1))
    maxE = max(np.max(matriceE0), np.max(matriceE1))


    #Building a unique matrice for the two cases, with rescaling with the same value
    #Treating different cases to rescale values. MaxP and maxE can't be null at the same time
    if maxP == 0:
        minE = min(minE0, minE1)
        matrice_finaleR0 = rescale_matrix(matriceE0, minE, maxE) - matriceP0
        matrice_finaleR1 = rescale_matrix(matriceE1, minE, maxE) - matriceP1
        #Positionning indications on colorbar
        minP = 0
        pos_min_P = -0.05
        pos_min_E = 0.05
    elif maxE == 0:
        minP = min(minP0, minP1)
        matrice_finaleR0 = matriceE0 - rescale_matrix(matriceP0, minP, maxP)
        matrice_finaleR1 = matriceE1 - rescale_matrix(matriceP1, minP, maxP)
        #Positionning indications on colorbar
        pos_min_P = -0.05
        pos_min_E = 0
        minE = 0
    else:
        minP = min(minP0, minP1)
        minE = min(minE0, minE1)
        matrice_finaleR0 = rescale_matrix(matriceE0, minE, maxE) - rescale_matrix(matriceP0, minP, maxP)
        matrice_finaleR1 = rescale_matrix(matriceE1, minE, maxE) - rescale_matrix(matriceP1, minP, maxP)
        #Positionning indications on colorbar
        pos_min_P = -0.05
        pos_min_E = 0.05

    #Computing the min and max for the colorbar (min can be either -1 or 0, max can be either 1 or 0)
    mini = min(np.min(matrice_finaleR0), np.min(matrice_finaleR1))
    maxi = max(np.max(matrice_finaleR0), np.max(matrice_finaleR1))


    Echelle_axes = np.linspace(0,1,n+1)

    #Editing datas so that it can be understood by seaborn
    data_E_P0 = pd.DataFrame(matrice_finaleR0, columns = np.around(Echelle_axes, decimals=3), index = np.around(Echelle_axes[::-1], decimals=3))
    data_E_P1 = pd.DataFrame(matrice_finaleR1, columns = np.around(Echelle_axes, decimals=3), index = np.around(Echelle_axes[::-1], decimals=3))


    #Editing figure to have two panels and a single colorbar
    fig, axs = plt.subplots(1,2, figsize=(18, 10), sharex=False, sharey=True)
    fig.subplots_adjust(wspace=0.65)
    #For the add.axes : left, bottom, width, height
    cbar_ax = fig.add_axes([.47, .15, .03, .7])



    #Plotting the two heatmaps
    sns.heatmap(data_E_P0, vmin=-1, vmax=1, cmap="PuOr_r", center = 0, cbar=True, square=False, xticklabels=200, yticklabels=200, ax=axs[0], cbar_kws={"ticks":[-1,pos_min_P, pos_min_E,1]}, cbar_ax=cbar_ax)
    sns.heatmap(data_E_P1, vmin=-1, vmax=1, cmap="PuOr_r", center = 0, cbar=True, square=False, xticklabels=200, yticklabels=200, ax=axs[1], cbar_kws={"ticks":[-1,pos_min_P, pos_min_E,1]}, cbar_ax=cbar_ax)
    axs[1].tick_params( which='both', top=False, bottom=True, left=False, right=False)

    axs[0].set_xlabel(r'$x_{het}$', fontsize=20)
    axs[1].set_xlabel(r'$x_{het}$', fontsize=20)

    cbar = axs[0].collections[0].colorbar
    cbar.set_ticklabels([round(maxP, 1), round(minP,2), round(minE, 1), round(maxE, 1)])

    #Common xlabels and ylabels
    ax = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

    ax.xaxis.set_label_coords(.5, -.11)
    ax.yaxis.set_label_coords(-.04, .5)
    plt.xlabel(r'$x_{het}$', fontsize=20)
    plt.ylabel(r'$\omega$', fontsize=20)
    ax.set_title(label =r'$\mu$ = %.2f, $r$=%.2f, $s_0$= %.2f, $s_{hom}$=%.2f, $s_{het}$=%.3f, $h$=%.2f' %(params[0], params[1], params[2], params[3], params[4], params[4]/params[3]), y=1.12, fontsize=20)
    #For scientific writing, .1E
    plt.text(0.15, 1.05, r'$R=0$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'semibold', fontsize = 20)
    plt.text(0.8, 1.05, r'$R=1$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'semibold', fontsize = 20)

    #To indicate the meaning of the colors alongside the colorbar
    plt.arrow(0.55, 0.55, 0, 0.35, clip_on=False, color='dimgrey', head_width=0.01)
    plt.text(0.56, 0.75, r'$\mathbb{E}[T_{ext}]$', fontdict=None, clip_on=False, rotation = 0, color='dimgrey', fontweight = 'semibold', fontsize = 15)
    plt.arrow(0.55, 0.45, 0, -0.35, clip_on=False, color='dimgrey', head_width=0.01)
    plt.text(0.57, 0.25, r'$\mathbb{P}_{ext}$', fontdict=None, clip_on=False, rotation = 0, color='dimgrey', fontweight = 'semibold', fontsize = 15)

    #To indicate the super/sub-criticality
    plt.arrow(x=0.41, y=0.5, dx=0.08, dy=0, color='dimgrey', linestyle='dotted', clip_on=False, head_width=0)
    plt.arrow(x=0.41, y=0.505, dx=0.08, dy=0, color='dimgrey', linestyle='dotted', clip_on=False, head_width=0)
    plt.text(0.415, 0.75, 'Sub-criticality', verticalalignment = 'center', fontdict=None, clip_on=False, rotation = 90, color='dimgrey', fontweight = 'semibold', fontsize = 15)
    plt.text(0.415, 0.25, 'Super-criticality', verticalalignment = 'center',fontdict=None, clip_on=False, rotation = 90, color='dimgrey', fontweight = 'semibold', fontsize = 15)





    plt.show()
    plt.savefig(f"Criticality_UnitypeBP_rescale_mu{params[0]}_r{params[1]}_s0{params[2]}_s{params[3]}_h{params[4]/params[3]}.png", bbox_inches="tight")



#Command to launch:
#heatmap_E_P(1000, params1)



