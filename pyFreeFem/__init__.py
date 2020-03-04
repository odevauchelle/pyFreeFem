
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Olivier Devauchelle
#
# Many thanks to Eric Lajeunesse, Anaïs Abramian, Valentin Jules & Hugo Chauvet.
#
# When used in a scientific publication, please cite:
#
# Boltzmann Distribution of Sediment Transport, A. Abramian, O. Devauchelle,
# G. Seizilles, E. Lajeunesse, Physical Review Letters, 123, 014501, 2019

"""
pyFreeFem
"""

__author__ = "Olivier Devauchelle"
__copyright__ = "Copyright 2020"
__license__ = "GPL"
__version__ = "0.2"

__all__ = ['TriMesh', 'FreeFemIO', 'FreeFemTools.FreeFemStatics', 'edpScript', 'pyFreeFem.functions']

from pyFreeFem.FreeFemTools.FreeFemStatics import *
from pyFreeFem.TriMesh import *
from pyFreeFem.FreeFemIO import *
from pyFreeFem.edpScript import *
from pyFreeFem.functions import *
