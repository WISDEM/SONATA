#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:49:51 2018

@author: gu32kij
"""
# Third party modules
import matplotlib.pyplot as plt
import numpy as np

# First party modules

# from matplotlib2tikz import save as tikz_save

def plot_beam_properties(data, sigma=None, ref=None, x_offset=0, description=True):
    """
    generates a plot of the beamproperties using the 
    
    Parameters
    ----------
    Data : np.ndarray 
        Massterms(6), Stiffness(21), damping(1), and coordinate(1) stacked 
        horizontally for the use in DYMORE/PYMORE/MARC  
    sigma: np.ndarray
        standart deviation of Data
    
    ref: np.ndarray
        of reference to select = 'all', 'Massterms', 'Stiffness'
    
    x_offset : float 
        x offset if curvilinear coordinate doen't start at the rotational 
        center of the rotordefault UH-60A value = 0.8178698
            
    """

    # Inertia Properties:
    ylabel = [(r"m_{00}", "kg/m"), (r"X_{m2}", "m"), (r"X_{m3}", "m"), (r"m_{33}", "kg m"), (r"m_{23}", "kg m"), (r"m_{22}", "kg m")]
    # x = data[:,-1] + x_offset
    x = data[:, -1] + x_offset
    x = x / x[-1]
    data = np.copy(data)
    data[:, 1] = data[:, 1] / data[:, 0]
    data[:, 2] = data[:, 2] / data[:, 0]
    if isinstance(ref, (np.ndarray)):
        ref = np.copy(ref)
        ref_x = ref[:, -1] + x_offset
        ref_x = ref_x / ref_x[-1]
        ref[:, 1] = ref[:, 1] / ref[:, 0]
        ref[:, 2] = ref[:, 2] / ref[:, 0]
    fig1, ax1 = plt.subplots(3, 2, sharex=True)
    fig1.subplots_adjust(wspace=0.3, hspace=0.3)
    # fig1.suptitle('Inertia Properties', fontsize=14)

    c = 0
    for i in range(3):
        for j in range(2):
            if isinstance(ref, (np.ndarray)):
                (ref_line,) = ax1[i][j].plot(ref_x, ref[:, c], "-.b")
            if isinstance(sigma, (np.ndarray)):
                ax1[i][j].fill_between(
                    x, data[:, c] - sigma[:, c], data[:, c] + sigma[:, c], alpha=0.6, edgecolor="r", linestyle=":", facecolor="r", antialiased=True,
                )
                # ax1[i][j].errorbar(x, data[:,c], yerr=sigma[:,c])
            (line,) = ax1[i][j].plot(x, data[:, c], "--k.", markersize=2)
            ax1[i][j].ticklabel_format(axis="y", style="sci")
            tmp = r"$%s \quad [%s]$" % (ylabel[c][0], ylabel[c][1])
            ax1[i][j].set_ylabel(tmp)

            if i == 2:
                ax1[i][j].set_xlabel(r"radius $\quad [1/R]$")

            c += 1

    # ax1[2,1].legend([ref_line, line], ['UH-60A Reference', 'Inertial Properties'])
    ax1[2, 1].set_ylim(0, 0.15)

    # 6x6 Stiffness Properties:
    TSunits = np.array(
        [
            ["N", "N", "N", "Nm", "Nm", "Nm"],
            ["N", "N", "N", "Nm", "Nm", "Nm"],
            ["N", "N", "N", "Nm", "Nm", "Nm"],
            ["Nm", "Nm", "Nm", "Nm^2", "Nm^2", "Nm^2"],
            ["Nm", "Nm", "Nm", "Nm^2", "Nm^2", "Nm^2"],
            ["Nm", "Nm", "Nm", "Nm^2", "Nm^2", "Nm^2"],
        ]
    )

    ut = np.zeros((6, 6))
    ut[np.triu_indices(6)] = 1

    fig2, ax2 = plt.subplots(6, 6)
    fig2.subplots_adjust(wspace=0.5, hspace=0.5)
    fig2.suptitle("6x6 Stiffness Matrix", fontsize=14)
    c = 6
    for j in range(6):
        for i in range(6):
            if ut[i, j] == 1:
                if isinstance(ref, (np.ndarray)):
                    ax2[i][j].plot(ref_x, ref[:, c], "--b")
                if isinstance(sigma, (np.ndarray)):
                    ax2[i][j].fill_between(
                        x, data[:, c] - sigma[:, c], data[:, c] + sigma[:, c], alpha=0.6, edgecolor="r", linestyle=":", facecolor="r", antialiased=True,
                    )
                ax2[i][j].plot(x, data[:, c], "--k.")
                ylabel = r"$k_{%s%s} \quad [%s]$" % (i + 1, j + 1, TSunits[i, j])
                ax2[i][j].set_ylabel(ylabel)
                # ax2[i][j].set_ylim(0)
                # ax2[i][j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

                if i == j:
                    ax2[i][j].set_xlabel(r"radius $\quad [1/R]$")

                c += 1
            else:
                fig2.delaxes(ax2[i][j])
                # fig2.text(0,0,r'$A=B\cdotC$')
                # fig2.text(0,0,r'$\left(\begin{matrix}k11&k12&k13&k14&k15&k16\\k12&k22&k23&k24&k25&k26\\k13&k23&k33&k34&k35&k36\\k14&k24&k34&k44&k45&k46\\k15&k25&k35&k45&k55&k56\\k16&k26&k36&k46&k56&k66\end{matrix}\right)$')

    desc = (
        r"\begin{minipage}[b]{10cm} "
        r"\underline{\textbf{Description:}} \\"
        r"The 6x6 sectional stiffness matrix, TS (Timoshenko Stiffness Matrix) "
        r"(1-extension; 2,3-shear, 4-twist; 5,6-bending) relates the sectional axial "
        r"strain, $\epsilon_1$, transverse shearing strains, $\epsilon_2$ and $\epsilon_3$, "
        r"twisting curvatures, $\kappa_1$ and two bending curvatures, $\kappa_2$ and $\kappa_3$, "
        r"to the axial force, $F_1$, transverse shear forces, $F_2$ and $F_3$, twisting "
        r"moment, $M_1$, and two bending moments, $M_2$ and $M_3$. The relationship between "
        r"these sectional strains and sectional stress resultants takes the form of a"
        r"symmetric, 6x6 matrix: \\"
        r"$$ \left( \begin{matrix} F_{1} \\ F_{2} \\ F_{3} \\ M_{1} \\ M_{2} \\ M_{3} \end{matrix} \right) = "
        r"\left( \begin{matrix} k_{11} & k_{12} & k_{13} & k_{14} & k_{15} & k_{16} \\ "
        r"k_{12} & k_{22} & k_{23} & k_{24} & k_{25} & k_{26} \\ "
        r"k_{13} & k_{23} & k_{33} & k_{34} & k_{35} & k_{36} \\ "
        r"k_{14} & k_{24} & k_{34} & k_{44} & k_{45} & k_{46} \\ "
        r"k_{15} & k_{25} & k_{35} & k_{45} & k_{55} & k_{56} \\ "
        r"k_{16} & k_{26} & k_{36} & k_{46} & k_{56} & k_{66} \end{matrix} \right) "
        r"\cdot \left( \begin{matrix} \epsilon_{1} \\ \epsilon_{2} \\ \epsilon_{3} \\ \kappa_{1} \\ \kappa_{2} \\ \kappa_{3} \end{matrix} \right) $$ \end{minipage} "
    )

    if description == True:
        from matplotlib import rcParams

        plt.figtext(0.05, 0.05, desc, usetex=True, wrap=True, bbox=dict(ec=(1.0, 0.5, 0.5), fc=(1.0, 0.8, 0.8)))

    # ============ ONLY FOR PAPER!
    ut = np.zeros((6, 6))
    ut[np.triu_indices(6)] = 1

    fig3, ax3 = plt.subplots(3, 2)
    fig3.subplots_adjust(wspace=0.5, hspace=0.3)
    # fig3.suptitle('6x6 Stiffness Matrix', fontsize=14)
    c = 6
    p = 0
    q = 0
    for j in range(6):
        for i in range(6):
            if ut[i, j] == 1:
                if i == j:
                    if p == 3:
                        p = 0
                        q = 1

                    if isinstance(ref, (np.ndarray)):
                        ax3[p][q].plot(ref_x, ref[:, c], "-.b")
                    if isinstance(sigma, (np.ndarray)):
                        ax3[p][q].fill_between(
                            x, data[:, c] - sigma[:, c], data[:, c] + sigma[:, c], alpha=0.6, edgecolor="r", linestyle=":", facecolor="r", antialiased=True,
                        )
                    ax3[p][q].plot(x, data[:, c], "--k.", markersize=2)
                    ylabel = r"$k_{%s%s} \quad [%s]$" % (i + 1, j + 1, TSunits[i, j])
                    ax3[p][q].set_ylabel(ylabel)
                    # ax2[i][j].set_ylim(0)
                    # ax2[i][j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                    ax3[p][q].set_ylim(0)
                    if p == 2:
                        ax3[p][q].set_xlabel(r"radius $\quad [1/R]$")
                    p += 1
                #                else:
                #                    fig3.delaxes(ax3[i][j])

                c += 1

    # ax3[2,1].legend([ref_line, line], ['UH-60A Reference', 'Inertial Properties'], loc='upper center')
    ax3[2, 1].set_ylim(0, 1.4e7)

    plt.show()
