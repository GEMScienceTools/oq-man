import re
import logging
import numpy as np
from pyproj import Proj, transform

from openquake.man.mfd import get_rates_within_m_range
from openquake.hazardlib.geo.geodetic import azimuth
from openquake.hazardlib.geo.geodetic import geodetic_distance

"""
def get_coo(lon, lat, lon0, lat0):
    dst = geodetic_distance(lon0, lat0, lon, lat)
    azi = azimuth(lon0, lat0, lon, lat)
    x = np.cos(np.radians(azi)) * dst
    y = np.sin(np.radians(azi)) * dst
    return x, y
"""

def getcoo(lon, lat):
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:3857')
    xp, yp = transform(inProj, outProj, x, y)
    return xp, yp


def get_rates_density(model, sidx, mmint=0.0, mmaxt=11.0, trt='.*', area=None):
    """
    :parameter model:
        A list of openquake source point instances
    :parameter model:
        A rtree spatial index
    :returns:
        A (key, value) dictionary, where key is the source ID and value
        corresponds to density of the rate of occurrence [eqks/(yr*km2)]
    """
    dens = {}
    if area is None:
        area = get_cell_sizes(model, sidx)
    for cnt, src in enumerate(model):
        if (re.search(trt, src.tectonic_region_type)):
            trates = get_rates_within_m_range(src.mfd, mmint, mmaxt)
            mmin, mmax = src.mfd.get_min_max_mag()
            dens[src.source_id] = trates / area[src.source_id]
            if abs(np.ceil(cnt/10000)*10000)-cnt < 1:
                print (cnt, len(model))
    return dens, area


def find_cell_size(pointl, psrc):
    """
    :parameter pointl:
    :parameter psrc:
    """
    mdx = 1e10
    mdy = 1e10
    dlo = 1e10
    dla = 1e10

    for point in pointl:
        # Check the points are at the same depth
        if abs(psrc.hypocenter_distribution.data[0][1] -
               point.hypocenter_distribution.data[0][1] < 0.5):

            if (np.sign(psrc.location.longitude) !=
                np.sign(point.location.longitude)):
                dltlo = (abs(np.sign(psrc.location.longitude) * 180 -
                             psrc.location.longitude) +
                         abs(np.sign(point.location.longitude) * 180 -
                             point.location.longitude))
            else:
                dltlo = abs(psrc.location.longitude - point.location.longitude)
            dltla = abs(psrc.location.latitude - point.location.latitude)

            # print (dltlo, dltla, dlo, dla)
            if dltlo > 0.025 and dltlo < dlo and dltla < 0.025:
                dx = geodetic_distance(psrc.location.longitude,
                                       psrc.location.latitude,
                                       point.location.longitude,
                                       point.location.latitude)
                mdx = dx if dx < mdx else mdx
                dlo = dltlo

            elif dltla > 0.025 and dltla < dla and dltlo < 0.025:
                dy = geodetic_distance(psrc.location.longitude,
                                       psrc.location.latitude,
                                       point.location.longitude,
                                       point.location.latitude)
                mdy = dy if dy < mdy else mdy
                dla = dltla

    csize = mdx * mdy
    if csize > 2000:
        print (mdx, mdy)
        raise ValueError('Cell size larger than threshold')

    return csize


def get_cell_sizes(pointl, sidx):
    """
    :parameter pointl:
        A list of openquake source point instances
    :parameter model:
        A rtree spatial index
    :return:
        A (key, value) dictionary, where key is the source ID and value
        corresponds to the area of the source.
    """

    logging.info('Computing cells size:', end="")

    lon0 = pointl[0].location.longitude
    lat0 = pointl[0].location.latitude
    csize = {}

    for cnt, psrc in enumerate(pointl):

        lon = psrc.location.longitude
        lat = psrc.location.latitude
        dep = psrc.location.depth

        xc, yc = get_coo(lon, lat, lon0, lat0)

        # Selection area (N.B. units are in km)
        dlt = 50.
        xl = xc - dlt
        yl = yc - dlt
        zl = -1
        xu = xc + dlt
        yu = yc + dlt
        zu = 100

        # Select the grid points at closest distance
        pidx = list(sidx.intersection((xl, yl, zl, xu, yu, zu)))
        # Create the list of points
        tmps = [pointl[iii] for iii in pidx]
        # Find the cell size
        csize[psrc.source_id] = find_cell_size(tmps, psrc)

        if abs(np.ceil(cnt/5000)*5000)-cnt < 1:
            print (cnt, len(pointl))

    return csize


def gcs(points, sidx):
    """
    """
    coo = [(p.location.longitude, p.location.latitude) for p in points]
    coo = np.array(coo)
    lla = sorted(set(coo[:,1]))

    lon0 = points[0].location.longitude
    lat0 = points[0].location.latitude
    csize = {}

    minlo = np.min(coo[:,1])
    maxlo = np.min(coo[:,1])

    # If the model is acreoss the IDL we swap the min and max longitudes
    if np.sign(minlo) != np.sign(maxlo):
        tmp = minlo
        minlo = maxlo
        maxlo = tmp

    for tla in lla:
        xl, yl = get_coo(minlo, tla-0.01, lon0, lat0)
        xu, yu = get_coo(maxlo, tla+0.01, lon0, lat0)
        print (xl, xu)
        pidx = list(sidx.intersection((xl, yl, -10, xu, yu, +500)))
        print ('len:', len(pidx))
