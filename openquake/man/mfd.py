import numpy

def get_rates_within_m_range(mfd, mmint=0.0, mmaxt=11.0):
    """
    :parameter mfd:
    :parameter mmint:
    :parameter mmaxt:
    """
    rtes = numpy.array(mfd.get_annual_occurrence_rates())
    idx = numpy.nonzero((rtes[:, 0] > mmint) & (rtes[:, 0] < mmaxt))
    return sum(rtes[idx[0], 1])
