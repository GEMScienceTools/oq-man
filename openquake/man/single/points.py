""" :mod:`openquake.man.single.point` module. This module contains functions
for computing general characteristics of gridded seismicity sources.
Assumptions: (1) the nodes of the grid represent centroids of cells (2) the
grid notes have a constant spacing (either in distance or long/lat). The
spacing can be different along longitude and along latitude """

import re
import logging
import numpy as np

from rtree import index
from pyproj import Proj, transform

from openquake.man.mfd import get_rates_within_m_range
from openquake.hazardlib.geo.geodetic import azimuth
from openquake.hazardlib.geo.geodetic import geodetic_distance


def get_rates_density(model, sidx, mmint=-11.0, mmaxt=11.0, trt='*', area=None):
    """
    :parameter model:
        A list of openquake source point instances
    :parameter sidx:
        A rtree spatial index
    :parameter mmint:
        Minimum magnitude
    :parameter mmaxt:
        Minimum magnitude
    :parameter trt:
        Tectonic region type keyword
    :returns:
        A (key, value) dictionary, where key is the source ID and value
        corresponds to density of the rate of occurrence [eqks/(yr*km2)]
    """
    dens = {}
    #
    # compute the area of each cell of the grid
    if area is None:
        area = get_cell_areas(model, sidx)
    #
    #
    for cnt, src in enumerate(model):
        if trt == src.tectonic_region_type:
            trates = get_rates_within_m_range(src.mfd, mmint, mmaxt)
            mmin, mmax = src.mfd.get_min_max_mag()
            dens[src.source_id] = trates / area[src.source_id]
            if abs(np.ceil(cnt/10000)*10000)-cnt < 1:
                print ('WARNING:', cnt, len(model))
    #
    #
    return dens, area


def get_cell_areas(points, sidx):
    """
    :parameter points:
    :parameter sidx:

    TODO: we need to add a test for the international date line since this is
        currently not supported. This might be simply solved by replacing
        the geographic coordinates in the spatial index with projected
        coordinates
    """
    #
    # preparing the spatial index
    sidx = index.Index()
    #
    # creating the list of coordinates and spatial index
    coo = []
    for idx, p in enumerate(points):
        lo = p.location.longitude
        la = p.location.latitude
        sidx.insert(idx, (lo, la, lo, la))
        coo.append((lo, la))
    #
    # compute the area of each grid cell
    areas = []
    coloc = []
    for lop, lap, in coo:
        #
        # get the 5 nearest neighbors
        nnidx = list(sidx.nearest((lop, lap, lop, lap), 20))
        #
        # get the coordinates of vertexes of the polygon
        area, clc = _get_cell_area(lop, lap, coo, nnidx)
        areas.append(areas)
        coloc.append(clc)
    #
    #
    return areas, coloc


def _get_cell_area(rlo, rla, coo, nnidx):
    """
    :parameter rlo:
    :parameter rla:
    :parameter coo:
    :parameter nnidx:
    """
    alo = [coo[idx][0] for idx in nnidx]
    ala = [coo[idx][1] for idx in nnidx]
    #
    # Computing azimuths and distances
    azis = azimuth(rlo, rla, alo, ala)
    dsts = geodetic_distance(rlo, rla, alo, ala)
    #
    # Processing the selected nodes
    delta = 5.0
    colocated = 0
    nearest_nodes = {}
    for azi, dst, idx in zip(azis, dsts, nnidx):
        if dst < 0.5:
            if (abs(rlo - coo[idx][0]) < 0.005 and
                abs(rla - coo[idx][1]) < 0.005):
                colocated += 1
        # East
        if abs(azi-90) < delta:
            if 90 in nearest_nodes:
                if dst < nearest_nodes[90][0]:
                    nearest_nodes[90] = (dst, idx)
            else:
                nearest_nodes[90] = (dst, idx)
        # South
        elif abs(azi-180) < delta:
            if 180 in nearest_nodes:
                if dst < nearest_nodes[180][0]:
                    nearest_nodes[180] = (dst, idx)
            else:
                nearest_nodes[180] = (dst, idx)
        # West
        elif abs(azi-270) < delta:
            if 270 in nearest_nodes:
                if dst < nearest_nodes[270][0]:
                    nearest_nodes[270] = (dst, idx)
            else:
                nearest_nodes[270] = (dst, idx)
        # North
        elif abs(azi-360) < delta or azi < delta:
            if 0 in nearest_nodes:
                if dst < nearest_nodes[0][0]:
                    nearest_nodes[0] = (dst, idx)
            else:
                nearest_nodes[0] = (dst, idx)
        else:
            pass
    #
    # fix missing information
    out = np.nan
    try:
        fdsts = _get_final_dsts(nearest_nodes)
        out = (fdsts[0]+fdsts[2])/2*(fdsts[1]+fdsts[3])/2
    except:
        pass
        print('Node:', rlo, rla)
        print (nearest_nodes)
        for idx in nnidx:
            print ('  ', coo[idx][0], coo[idx][1])
    return out, colocated


def _get_final_dsts(nno):
    """
    :parameter nno:
        A dictionary where keys are angles (0, 90, 180 and 270) and values are
        tuples containing a distance and one index
    """
    #
    # array containing the final values of distance along the 4 main directions
    fd = []
    #
    # north
    if 0 in nno:
        fd.append(nno[0][0])
    elif 180 in nno:
        fd.append(nno[180][0])
    else:
        raise ValueError('Cannot define distance toward north')
    #
    # east
    if 90 in nno:
        fd.append(nno[90][0])
    elif 270 in nno:
        fd.append(nno[270][0])
    else:
        raise ValueError('Cannot define distance toward east')
    #
    # south
    if 180 in nno:
        fd.append(nno[180][0])
    elif 0 in nno:
        fd.append(nno[0][0])
    else:
        raise ValueError('Cannot define distance toward south')
    #
    # west
    if 270 in nno:
        fd.append(nno[270][0])
    elif 90 in nno:
        fd.append(nno[90][0])
    else:
        raise ValueError('Cannot define distance toward west')
    #
    # final result
    return fd
