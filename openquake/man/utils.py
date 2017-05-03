import pyproj
import shapely.ops as ops
from functools import partial


def get_area(geom):
    """
    Compute the area of a shapely polygon with lon, lat coordinates.
    See http://tinyurl.com/h35nde4

    :parameter geom:

    :return:
        The area of the polygon in km2
    """
    geom_aea = ops.transform(partial(pyproj.transform,
                                     pyproj.Proj(init='EPSG:4326'),
                                     pyproj.Proj(proj='aea',
                                                 lat1=geom.bounds[1],
                                                 lat2=geom.bounds[3])),
                             geom)
    return geom_aea.area/1e6
