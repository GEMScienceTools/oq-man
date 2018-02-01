# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2014-2018 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

url = "https://github.com/GEMScienceTools/oq-mbtk"

README = """
Python and OpenQuake-based Toolkit for the analysis of Seismic Source
Models
Copyright (C) 2017-2018 GEM Foundation
"""

setup(
    name='openquake.man',
    version='0.1.0',
    description=README,
    url=url,
    packages=find_packages(exclude=['tests', 'tests.*']),
    # Minimal requirements, for a complete list see requirements-*.txt
    # matplotlib is brought by the openquake engine
    install_requires=[
        'openquake.engine',
        'pyproj',
    ],
    python_requires='>=3',
    author='GEM Foundation',
    author_email='hazard@globalquakemodel.org',
    maintainer='GEM Foundation',
    maintainer_email='hazard@globalquakemodel.org',
    classifiers=(
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ),
    namespace_packages=['openquake'],
    keywords="seismic hazard",
    license="AGPL3",
    platforms=["any"],
    package_data={"openquake.man": [
        "README.md", "LICENSE"]},
    include_package_data=True,
    zip_safe=False,
)
