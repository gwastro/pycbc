# Copyright (C) 2017  Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This modules contains information about the announced LIGO/Virgo
compact binary mergers
"""

# For the time being all quantities are the 1-d median value
# FIXME with posteriors when available and we can just post-process that
data = {}

# O1 events, source https://arxiv.org/pdf/1606.04856v3.pdf
event = "GW150914"
data[event] = e = {}
e['time'] = 1126259462.4

# Maybe support multiple data types later and add selection to the accessing class
e['frames'] = {"H1":'https://losc.ligo.org/s/events/GW150914/H-H1_LOSC_4_V2-1126259446-32.gwf',
               "L1":'https://losc.ligo.org/s/events/GW150914/L-L1_LOSC_4_V2-1126259446-32.gwf',
              }
e["median1d"] = d = {}
d["mass1"] = (36.2, -3.8, +5.2)
d["mass2"] = (29.1, -4.4, +3.7)
d["mchirp"] = (28.1, -1.5, +1.8)
d["mtotal"] = (65.3, -3.4, +4.1)
d["chi_eff"] = (-.06, -.14, +.14)
d["redshift"] = (0.09, -.04, +.03)


event = "GW151226"
data[event] = e = {}
e['time'] = 1135136350.65
e['frames'] = {"H1":'https://losc.ligo.org/s/events/GW151226/H-H1_LOSC_4_V2-1135136334-32.gwf',
               "L1":'https://losc.ligo.org/s/events/GW151226/L-L1_LOSC_4_V2-1135136334-32.gwf',
              }
e["median1d"] = d = {}
d["mass1"] = (14.2, -3.7, +8.3)
d["mass2"] = (7.5, -2.3, +2.3)
d["mchirp"] = (8.9, -.3, +.3)
d["mtotal"] = (21.8, -1.7, +5.9)
d["chi_eff"] = (-.21, -.10, +.20)
d["redshift"] = (0.09, -.04, +.03)

event = "LVT151012"
data[event] = e = {}
e['time'] = 1128678900.44
e['frames'] = {"H1":'https://losc.ligo.org/s/events/LVT151012/H-H1_LOSC_4_V2-1128678884-32.gwf',
               "L1":'https://losc.ligo.org/s/events/LVT151012/L-L1_LOSC_4_V2-1128678884-32.gwf',
              }
e["median1d"] = d = {}
d["mass1"] = (23, -6., +18.)
d["mass2"] = (13, -5., +4.)
d["mchirp"] = (15.1, -1.1, +1.4)
d["mtotal"] = (37, -4., +13.)
d["chi_eff"] = (0.0, -0.2, +0.3)
d["redshift"] = (0.20, -0.09, +0.09)

# O2 Events
#https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.118.221101
event = "GW170104"
data[event] = e = {}
e['time'] = 1167559936.6
e['frames'] = {"H1":'https://losc.ligo.org/s/events/GW170104/H-H1_LOSC_4_V1-1167559920-32.gwf',
               "L1":'https://losc.ligo.org/s/events/GW170104/L-L1_LOSC_4_V1-1167559920-32.gwf',
              }
e["median1d"] = d = {}
d["mass1"] = (31.2, -6.0, +8.4)
d["mass2"] = (19.4, -5.9, +5.3)
d["mchirp"] = (21.1, -2.7, +2.4)
d["mtotal"] = (50.7, -5.0, +5.9)
d["chi_eff"] = (-0.12, -0.30, +0.21)
d["redshift"] = (0.18, -0.07, +0.08)

#https://dcc.ligo.org/LIGO-P170814/public/main
event = "GW170814"
data[event] = e = {}
e['time'] = 1186741861.53
e['frames'] = {"H1":"https://dcc.ligo.org/public/0146/P1700341/001/H-H1_LOSC_CLN_4_V1-1186741845-32.gwf",
               "L1":"https://dcc.ligo.org/public/0146/P1700341/001/L-L1_LOSC_CLN_4_V1-1186741845-32.gwf",
               "V1":"https://dcc.ligo.org/public/0146/P1700341/001/V-V1_LOSC_CLN_4_V1-1186741845-32.gwf",
              }
e["median1d"] = d = {}
d["mass1"] = (30.5, -3.0, +5.7)
d["mass2"] = (25.3, -4.2, +2.8)
d["mchirp"] = (24.1, -1.1, +1.4)
d["mtotal"] = (55.9, -2.7, +3.4)
d["chi_eff"] = (0.06, -.12, +.12)
d["redshift"] = (.11, -.04, +0.03)
