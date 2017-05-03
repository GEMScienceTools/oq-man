import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from openquake.hazardlib.source import SimpleFaultSource
from openquake.hazardlib.source import CharacteristicFaultSource
from openquake.hazardlib.geo.surface import SimpleFaultSurface


def plt_mmax_area(model):
    """
    :parameter model:
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    mmaxs = get_mmax(model)
    areas = get_areas(model)
    plt.plot(mmaxs, areas, 'o')
    plt.xlabel(r'Maximum magnitude')
    plt.ylabel(r'Fault surface area [km2]')
    plt.yscale('log')
    plt.grid(which='both', linestyle='--')
    text = plt.text(0.025, 0.95,
                    'maximum magnitude {0:.2f}'.format(max(mmaxs)),
                    transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=2,
                                                  foreground="white")])
    return fig


def get_areas(model):
    """
    Get maximum magnitudes from MFDs

    :parameter model:
        A list of hazardlib source instances
    :returns:
        A list with the areas of each fault
    """
    msp = 2
    areas = []
    for src in model:
        if isinstance(src, SimpleFaultSource):
            fault_surface = SimpleFaultSurface.from_fault_data(
                                                src.fault_trace,
                                                src.upper_seismogenic_depth,
                                                src.lower_seismogenic_depth,
                                                src.dip,
                                                msp)
            area = fault_surface.get_area()
            areas.append(area)
    return areas


def plt_mmax_length(model):
    """
    :parameter model:
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    mmaxs = get_mmax(model)
    lngths = get_lengths(model)
    plt.plot(mmaxs, lngths, 'o')
    plt.xlabel(r'Maximum magnitude')
    plt.ylabel(r'Fault length [km]')
    plt.yscale('log')
    plt.grid(which='both', linestyle='--')
    text = plt.text(0.025, 0.95,
                    'maximum magnitude {0:.2f}'.format(max(mmaxs)),
                    transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=2,
                                                  foreground="yellow")])
    return fig


def get_mmax(model):
    """
    Get maximum magnitudes from MFDs

    :parameter model:
        A list of hazardlib source instances
    :returns:
        A list with the maximum magnitude values
    """
    mmaxs = []
    for src in model:
        if isinstance(src, SimpleFaultSource):
            min_mag, max_mag = src.mfd.get_min_max_mag()
            mmaxs.append(max_mag)
    return mmaxs


def _add_mmax_histogram(mmaxs):
    binw = 0.5
    mmin = np.floor(min(mmaxs)/binw)*binw
    mmax = np.ceil(max(mmaxs)/binw)*binw
    edges = np.arange(mmin, mmax+0.01*binw, step=binw)
    hist, bin_edges = np.histogram(mmaxs, bins=edges)
    plt.hist(mmaxs, bins=bin_edges, log=True, edgecolor='white')


def plt_mmax_histogram(mmaxs):
    """
    :parameter mmaxs:
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    _add_mmax_histogram(mmaxs)
    plt.xlabel(r'Maximum magnitude')
    plt.grid(which='both', linestyle='--')
    text = plt.text(0.7, 0.95, '{0:d} sources'.format(len(mmaxs)),
                    transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3,
                                                  foreground="y")])
    return fig


def get_lengths(model, trt='.*'):
    """
    Compute the length of all the simple fault sources

    :parameter model:
        A list of hazardlib source instances
    :returns:
        A list with the length [km] for all the simple fault sources included
        in the model.
    """
    lngths = {}
    for src in model:
        if (isinstance(src, SimpleFaultSource) and
                re.search(trt, src.tectonic_region_type)):
            lngths[src.source_id] = src.fault_trace.get_length()
        else:
            pass
            # print ('unsupported fault type', type(src).__name__)
    return lngths


def _add_length_histogram(flens):
    hist, bin_edges = np.histogram(flens)
    plt.hist(flens, bins=bin_edges, log=True, edgecolor='white')


def plt_length_histogram(flens):
    """
    :parameter model:
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    _add_length_histogram(flens)
    plt.xlabel(r'Fault length [km]')
    plt.grid(which='both', linestyle='--')
    plt.text(0.7, 0.95, '# of sources: {0:d}'.format(len(flens)),
             transform=ax.transAxes)
    return fig
