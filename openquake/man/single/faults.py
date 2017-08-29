import re
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from openquake.hazardlib.source import SimpleFaultSource
from openquake.hazardlib.source import CharacteristicFaultSource
from openquake.hazardlib.geo.surface import SimpleFaultSurface
from openquake.hazardlib.geo.geodetic import distance, min_geodetic_distance


def _get_bb_mesh_from_char_fs(fault, mesh_spacing):
    """
    :parameter fault:
        An instance of the
        :class:`openquake.hazardlib.source.CharacteristicFaultSource`
    :parameter float mesh_spacing:
        The spacing [km] of the grid used to represent the fault surface
    """
    srfA = fault.surface
    bbox = srfA.get_bounding_box()
    mesh = srfA.get_mesh()
    return mesh, bbox

def _get_bb_mesh_from_simple_fs(fault, mesh_spacing):
    """
    :parameter fault:
        An instance of the
        :class:`openquake.hazardlib.source.CharacteristicFaultSource`
    :parameter float mesh_spacing:
        The spacing [km] of the grid used to represent the fault surface
    """
    srfA = SimpleFaultSurface.from_fault_data(fault.fault_trace,
                                    fault.upper_seismogenic_depth,
                                    fault.lower_seismogenic_depth,
                                    fault.dip,
                                    mesh_spacing)
    bbox = srfA.get_bounding_box()
    mesh = srfA.get_mesh()
    return mesh, bbox

def fault_surface_distance(srcs, mesh_spacing, types=None):
    """
    :parameter srcs:
    :parameter float mesh_spacing:
        The spacing [km] of the grid used to represent the fault surface
    :parameter types:
    """
    #
    # Fault distance array
    dmax = 40
    lnkg = np.ones((len(srcs), len(srcs)))*dmax
    #
    # Check input information
    if types is None:
        types = (SimpleFaultSource, CharacteristicFaultSource)
    #
    # Process srcs
    for idxA in range(0, len(srcs)):
        srcA = srcs[idxA]
        if isinstance(srcA, types):
            if isinstance(srcA, SimpleFaultSource):
                meshA, bboxA = _get_bb_mesh_from_simple_fs(srcA, mesh_spacing)
            elif isinstance(srcA, CharacteristicFaultSource):
                meshA, bboxA = _get_bb_mesh_from_char_fs(srcA, mesh_spacing)
            else:
                raise ValueError('Unsupported fault type')
            #
            # Second loop
            for idxB in range(idxA, len(srcs)):
                srcB = srcs[idxB]
                if isinstance(srcB, types):
                    if isinstance(srcB, SimpleFaultSource):
                        meshB, bboxB = _get_bb_mesh_from_simple_fs(srcB,
                                                                   mesh_spacing)
                    elif isinstance(srcB, CharacteristicFaultSource):
                        meshB, bboxB = _get_bb_mesh_from_char_fs(srcB,
                                                                 mesh_spacing)
                    else:
                        raise ValueError('Unsupported fault type')
                    #
                    # Calculate the distance between the two bounding boxes
                    tmpd = min_geodetic_distance(
                        np.array([bboxA[0], bboxA[2], bboxA[0], bboxA[2]]),
                        np.array([bboxA[1], bboxA[3], bboxA[3], bboxA[1]]),
                        np.array([bboxB[0], bboxB[2], bboxB[0], bboxB[2]]),
                        np.array([bboxB[1], bboxB[3], bboxB[3], bboxB[1]]))
                    #
                    # Calculate the distance between the two fault surfaces
                    # (if needed)
                    if (np.amin(tmpd) > dmax):
                        mindst = np.amin(tmpd)
                    else:
                        mindst = 1e10
                        for lon, lat, dep in zip(meshB.lons.flatten(),
                                                 meshB.lats.flatten(),
                                                 meshB.depths.flatten()):
                            tmpd = distance(meshA.lons, meshA.lats,
                                            meshA.depths,
                                            lon, lat, dep)
                            mindst = min(mindst, np.amin(tmpd))
                    #
                    # Update the array
                    lnkg[idxA, idxB] = np.amin(tmpd)
    assert len(lnkg) == len(srcs)
    return lnkg


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
        A list of :class:`openquake.hazardlib.source` instances
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
        A dictionary with the length [km] for all the simple fault sources
        included in the model.
    """
    lngths = {}
    for src in model:
        if (isinstance(src, SimpleFaultSource) and
                re.search(trt, src.tectonic_region_type)):
            lngths[src.source_id] = src.fault_trace.get_length()
        else:
            ogging.warning('unsupported fault type', type(src).__name__)
    return lngths


def _add_length_histogram(flens):
    """
    This adds an histogram to

    :parameter flens:
        A list of values representing the lengths of fault traces
    """
    hist, bin_edges = np.histogram(flens)
    plt.hist(flens, bins=bin_edges, log=True, edgecolor='white')


def plt_length_histogram(flens):
    """
    :parameter flens:
        A list of values representing the lengths of fault traces
    :returns:
        A matplotlib figure
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    _add_length_histogram(flens)
    plt.xlabel(r'Fault length [km]')
    plt.grid(which='both', linestyle='--')
    plt.text(0.7, 0.95, '# of sources: {0:d}'.format(len(flens)),
             transform=ax.transAxes)
    return fig
