

from openquake.man.single.fault import fault_surface_distance

from openquake.hazardlib.source import SimpleFaultSource
from openquake.hazardlib.mfd import EvenlyDiscretizedMFD
from openquake.hazardlib.scalerel import WC1994
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.geo import Line

class TestFaultSurfaceDistance(unittest.TestCase):

    def setUp(self):
        """
        """
        #
        # Simple fault source
        mfd = EvenlyDiscretizedMFD(6.0, 0.1, [1.])
        msr = WC1994()
        tom = PoissonTOM(1.0)
        trace = Line([Point(10, 45.), Point(10., 45.2)])
        sfs = SimpleFaultSource(source_id='1',
                                name='1',
                                tectonic_region_type='none',
                                mfd=mfd,
                                rupture_mesh_spacing = 5.,
                                magnitude_scaling_relationship=msr,
                                rupture_aspect_ratio=1.0,
                                temporal_occurrence_model=tom,
                                upper_seismogenic_depth=0.,
                                lower_seismogenic_depth=10.,
                                fault_trace=trace,
                                dip=90,
                                rake=90)
        self.srcs = [sfs]
        #
        #
        mfd = EvenlyDiscretizedMFD(6.0, 0.1, [1.])
        msr = WC1994()
        tom = PoissonTOM(1.0)
        trace = Line([Point(10.2, 45.), Point(10.2, 45.2)])
        sfs = SimpleFaultSource(source_id='1',
                                name='1',
                                tectonic_region_type='none',
                                mfd=mfd,
                                rupture_mesh_spacing = 5.,
                                magnitude_scaling_relationship=msr,
                                rupture_aspect_ratio=1.0,
                                temporal_occurrence_model=tom,
                                upper_seismogenic_depth=0.,
                                lower_seismogenic_depth=10.,
                                fault_trace=trace,
                                dip=90,
                                rake=90)
        self.srcsa.append(sfs)

    def test01(self):
        """
        Results computed using the tools available at
        https://geographiclib.sourceforge.io/
        """

        dst = fault_surface_distance(self.srcs), 5.0)
        expected = np.array([43736430.1, 43736430.1, 43736430.1,
                             43774201.0, 43774201.0, 43774201.0,
                             43811937.7, 43811937.7, 43811937.7])/1e6
        computed = get_cell_areas(self.points, self.sidx)

