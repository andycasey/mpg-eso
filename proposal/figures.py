import matplotlib
from matplotlib.ticker import MaxNLocator
import numpy as np

from astropy import coordinates as coord
from astropy import units as u
from astropy.table import Table

from colormaps import (magma, inferno, plasma, viridis)

def filter_response():
    """
    Show the normalised response of the VPHAS g, WFI 871 and DDO51 filters.
    """

    fig, ax = plt.subplots()

    ddo51 = np.loadtxt("filter-responses/DDO51_filter.dat")
    ddo51[:, 0] /= 10. # Convert A to nm
    ddo51[:, 1] /= ddo51[:, 1].max()

    ax.plot(ddo51[:, 0], ddo51[:, 1], linestyle="--", c="#91AA9D", lw=2, label=r"KPNO DDO51")

    wfi_871 = np.loadtxt("filter-responses/wfi_filter_871.dat")
    wfi_871[:, 1] /= wfi_871[:, 1].max()

    ax.plot(wfi_871[:, 0], wfi_871[:, 1], c="#3E606F", lw=2, label=r"WFI 871",
        zorder=10)

    sdss_g = np.loadtxt("filter-responses/sdss_filter_g.dat")
    sdss_g[:, 1] /= sdss_g[:, 1].max()
    ax.plot(sdss_g[:, 0]/10., sdss_g[:, 1], c="#0874D4", lw=2, label=r"VPHAS g")

    ax.set_ylim(0, 1.05)
    ax.set_xlim(350, 550)

    ax.legend(loc="upper left", frameon=False)

    ax.set_xlabel(r"Wavelength [nm]")
    ax.set_ylabel(r"Normalised response")

    fig.tight_layout()

    return fig


def giant_selection(ax=None):

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    kepler = Table.read("data.fits.gz")

    # Scatter points of all stars coloured by logg
    x = kepler["kic_gmag_01"] - kepler["kic_kmag_01"]
    y = kepler["kic_gmag_01"] - kepler["kic_d51mag"]
    scat = ax.scatter(x, y, c=kepler["kic_logg_01"],
        vmin=0, vmax=5, cmap="terrain")


    xlims = [-1, 7]
    ylims = [-0.15, 0.65]
    xlims, ylims = map(np.array, (xlims, ylims))

    # Selection lines.
    x1 = 0.64/0.18
    x2 = 0.36/0.18
    x = np.array([x1, x2])
    ax.plot(x, -0.08 * x + 0.46, c="k", lw=2)

    x = np.array([x1, xlims[1]])
    y1 = +0.10 * x - 0.18
    ax.plot(x, y1, c="k", lw=2)

    x = np.array([x2, xlims[1]])
    y2 = +0.10 * x + 0.10
    ax.plot(x, y2, c="k", lw=2)

    # Fill between?
    x = np.array([x2, xlims[1]])
    ax.fill(
        [x2, x1, xlims[1], xlims[1], (ylims[1] - 0.1)/0.1, x2], 
        [
            +0.10 * x2 + 0.10,
            +0.10 * x1 - 0.18,
            +0.10 * xlims[1] - 0.18,
            ylims[1],
            ylims[1],
            +0.10 * x2 + 0.10],
        facecolor='#D2D9B8',
        zorder=-100)

    # Colorbar.
    cbar = plt.colorbar(scat)
    cbar.set_label(r"$\log{g}$")
    cbar.locator = MaxNLocator(6)
    cbar.update_ticks()

    # Ticks, limits and labels.
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    
    ax.set_xlabel(r"$g - K_{s}$")
    ax.set_ylabel(r"$g - {\rm DDO}_{51}$")

    fig.tight_layout()

    return fig





def vphas_dr2_g():

    fig, ax = plt.subplots(figsize=(12, 4))

    fov_edge = 1.
    tile_points = lambda ra, dec: np.array([
        [ra - fov_edge/2, dec - fov_edge/2],
        [ra + fov_edge/2, dec - fov_edge/2],
        [ra + fov_edge/2, dec + fov_edge/2],
        [ra - fov_edge/2, dec + fov_edge/2],
        [ra - fov_edge/2, dec - fov_edge/2]
    ])

    # Draw the planned regions of the VPHAS Survey.
    ax.fill(
        [40, 10, 10, -10, -10, -200, -200, -10, -10,  10, 10, 40],
        [ 5,  5, 10,  10,   5,    5,   -5,  -5, -10, -10, -5, -5],
        facecolor="#cccccc", edgecolor="#cccccc", zorder=-1,
        label="Planned survey area")

    # Load the dr2 exposures.
    vphas = Table.read("vphas-dr2-ugr-exposures.csv")
    vphas_g = vphas[vphas["filter"] == "g_SDSS"]

    # Draw each tile
    for i, field in enumerate(vphas_g):
        print("doing {0}/{1}".format(i, len(vphas_g)))
        field_center = [field["RA_1"], field["Dec_1"]]

        # Get tile points.
        points = tile_points(*field_center)

        # Convert to galactic.
        c = coord.ICRS(ra=points[:, 0], dec=points[:, 1], unit=(u.deg, u.deg))

        l, b = c.galactic.l.value, c.galactic.b.value
        l[l > 180] -= 360

        kwd = {} if i > 0 else {"label": r"DR2 ($g$-band)"}
        ax.fill(l, b, facecolor="#0874D4", lw=0.5, **kwd)
        
    ax.set_xlim(40, -40)
    ax.set_ylim(-10, 10)

    ax.set_xticks(np.arange(-40, 50, 10)[::-1]) # [40 to -40]
    ax.set_yticks([-10, -5, 0, 5, 10])

    ax.set_xlabel(r"Galactic longitude, $l$ $[^\circ]$")
    ax.set_ylabel(r"Galactic latitude, $b$ $[^\circ]$")
    ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()

    return fig


if __name__ == "__main__":

    figures = {
        "filter-response.pdf": filter_response,
        "giant-selection.pdf": giant_selection,
        "vphas-g.pdf": vphas_dr2_g
    }

    for filename, function in figures.items():
        fig = function()
        fig.savefig(filename)

