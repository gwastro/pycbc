# Copyright (C) 2012  Alex Nitz, Andrew Miller
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
'''
This module tests basic operations including when the starting arrays are on different schemes.
The checkCPU and checkScheme functions must be declared in the child class, because they require
different checks for Arrays, TimeSeries or FrequencySeries.
'''

class array_base(object):
    def setAnswers(self):
        # All the answers are stored here to make it easier to read in the actual tests.
        # Again, it makes a difference whether they are complex or real valued, so there
        # are four sets of possible answers, depending on the dtypes. The fourth value in
        # the list (True or False) lets the checker know whether to use Equal or AlmostEqual
        if self.kind == 'real' and self.okind == 'real':
        
            self.mul = [50, 24, 6, True]
            self.mul_s = [25, 15, 5, True]
                        
            self.add = [15, 11, 7, True]
            self.add_s = [10, 8, 6, True]
                        
            #self.div = [.5, 3./8., 1./6.]
            self.div = [.5, 0.375, .16666666666666666667, False]
            #self.div_s = [1., 3./5., 1./5.]
            self.div_s = [1., 0.6, 0.2, False]
                        
            #self.rdiv = [2., 8./3., 6.]
            self.rdiv = [2., 2.66666666666666666667, 6., False]
            #self.rdiv_s = [1., 5./3., 5.]
            self.rdiv_s = [1., 1.66666666666666666667, 5., False]
            
            self.sub = [-5, -5, -5, True]
            self.sub_s = [0, -2, -4, True]
            
            self.rsub = [5, 5, 5, True]
            self.rsub_s = [0, 2, 4, True]
            
            self.pow1 = [25., 9., 1., False]
            #self.pow2 = [pow(5,-1.5), pow(3,-1.5), pow(1,-1.5)]
            self.pow2 = [0.08944271909999158786, 0.19245008972987525484, 1., False]
            
            self.abs = [5, 3, 1, True]
            
            self.real = [5,3,1, True]
            self.imag = [0, 0, 0, True]
            self.conj = [5, 3, 1, True]
            
            self.sum = 9
            
            self.dot = 80
                        
        if self.kind =='real' and self.okind == 'complex':
            
            self.mul = [50+30j, 24+12j, 6+2j, True]
            self.mul_s = [25+10j, 15+6j, 5+2j, True]
            
            self.add = [15+6j, 11+4j, 7+2j, True]
            self.add_s = [10+2j, 8+2j, 6+2j, True]
            
            #self.div = [25./68.-15.j/68., 3./10.-3.j/20., 3./20.-1.j/20.] 
            self.div = [0.36764705882352941176-0.22058823529411764706j, 0.3-0.15j, 0.15-0.05j, False] 
            #self.div_s = [25./29.-10.j/29., 15./29.-6.j/29., 5./29.-2.j/29.]
            self.div_s = [0.86206896551724137931-0.34482758620689655172j,
                          0.51724137931034482759-0.20689655172413793103j,
                          0.17241379310344827586-0.06896551724137931034j, False]
            
            #self.rdiv = [2.+6.j/5., 8./3.+4.j/3, 6.+2.j]
            self.rdiv = [2.+1.2j, 2.66666666666666666667+1.33333333333333333333j, 6.+2.j, False]
            #self.rdiv_s = [1.+2.j/5., 5./3.+2.j/3., 5.+2.j]
            self.rdiv_s = [1.+0.4j, 1.66666666666666666667+0.666666666666666666667j, 5.+2.j, False]
            
            self.sub = [-5-6j, -5-4j, -5-2j, True]
            self.sub_s = [0-2j, -2-2j, -4-2j, True]
            
            self.rsub = [5+6j, 5+4j, 5+2j, True]
            self.rsub_s = [0+2j, 2+2j, 4+2j, True]
            
            self.pow1 = [25., 9., 1., False]
            #self.pow2 = [pow(5,-1.5), pow(3,-1.5), pow(1,-1.5)]
            self.pow2 = [0.08944271909999158786, 0.19245008972987525484, 1., False]
            
            self.abs = [5, 3, 1, True]
            
            self.real = [5,3,1, True]
            self.imag = [0, 0, 0, True]
            self.conj = [5, 3, 1, True]
            
            self.sum = 9
            
            self.dot = 80+44j
            
        if self.kind == 'complex' and self.okind == 'real':
            
            self.mul = [50+10j, 24+24j, 6+30j, True]
            self.mul_s = [25+5j, 15+15j, 5+25j, True]
            
            self.add = [15+1j, 11+3j, 7+5j, True]
            self.add_s = [10+1j, 8+3j, 6+5j, True]
            
            #self.div = [1./2.+1.j/10., 3./8.+3.j/8., 1./6.+5.j/6.]
            self.div = [0.5+0.1j, 0.375+0.375j, 0.16666666666666666667+0.83333333333333333333j, False]
            #self.div_s = [1.+1.j/5., 3./5.+3.j/5., 1./5.+1.j]
            self.div_s = [1.+0.2j, 0.6+0.6j, 0.2+1.j, False]
            
            #self.rdiv = [25./13.-5.j/13., 4./3.-4.j/3., 3./13.-15.j/13.]
            self.rdiv = [1.92307692307692307692-0.38461538461538461538j,
                         1.33333333333333333333-1.33333333333333333333j,
                         0.23076923076923076923-1.15384615384615384615j, False]
            #self.rdiv_s = [25./26.-5.j/26., 5./6.-5.j/6., 5./26.-25.j/26.]
            self.rdiv_s = [0.96153846153846153846-0.19230769230769230769j,
                           0.83333333333333333333-0.83333333333333333333j,
                           0.19230769230769230769-0.96153846153846153846j, False]
            
            self.sub = [-5+1j, -5+3j, -5+5j, True]
            self.sub_s = [0+1j, -2+3j, -4+5j, True]
            
            self.rsub = [5-1j, 5-3j, 5-5j, True]
            self.rsub_s = [0-1j, 2-3j, 4-5j, True]
            
            self.pow1 = [24.+10.j, 0.+18.j, -24.+10.j, False]
            #self.pow2 = [pow(5+1j,-1.5), pow(3+3j,-1.5), pow(1+5j,-1.5)]
            self.pow2 = [0.08307064054041229214-0.0253416052125975132j,
                         0.04379104225017853491-0.1057209281108342370j,
                        -0.04082059235165559671-0.0766590341356157206j, False]
            
            #self.abs = [pow(26,.5), 3*pow(2,.5), pow(26,.5)]
            self.abs = [5.09901951359278483003,
                        4.24264068711928514641,
                        5.09901951359278483003, False]
            
            self.real = [5,3,1, True]
            self.imag = [1, 3, 5, True]
            self.conj = [5-1j, 3-3j, 1-5j, True]
            
            self.sum = 9+9j
            
            self.dot = 80+64j
            
        if self.kind =='complex' and self.okind =='complex':
            
            self.mul = [44+40j, 12+36j, -4+32j, True]
            self.mul_s = [23+15j, 9+21j, -5+27j, True]
            
            self.add = [15+7j, 11+7j, 7+7j, True]
            self.add_s = [10+3j, 8+5j, 6+7j, True]
            
            #self.div = [7./17.-5.j/34., 9./20.+3.j/20., 2./5.+7.j/10.]
            self.div = [0.41176470588235294118-0.14705882352941176471j, 0.45+0.15j, 0.4+0.7j, False]
            #self.div_s = [27./29.-5.j/29., 21./29.+9.j/29., 15./29.+23.j/29.]
            self.div_s = [0.93103448275862068966-0.17241379310344827586j,
                          0.72413793103448275862+0.31034482758620689655j,
                          0.51724137931034482759+0.79310344827586206897j, False]
            
            #self.rdiv = [28./13.+10.j/13., 2.-2.j/3., 8./13.-14.j/13.]
            self.rdiv = [2.15384615384615384615+0.76923076923076923077j,
                         2.                    -0.66666666666666666667j,
                         0.61538461538461538462-1.07692307692307692308j, False]
            #self.rdiv_s = [27./26.+5.j/26., 7./6.-1.j/2., 15./26.-23.j/26]             
            self.rdiv_s = [1.03846153846153846154+0.19230769230769230769j,
                           1.16666666666666666667-0.5j,
                           0.57692307692307692308-0.88461538461538461538j, False]
            
            self.sub = [-5-5j, -5-1j, -5+3j, True]
            self.sub_s = [0-1j, -2+1j, -4+3j, True]
            
            self.rsub = [5+5j, 5+1j, 5-3j, True]
            self.rsub_s = [0+1j, 2-1j, 4-3j, True]
            
            self.pow1 = [24.+10.j, 0.+18.j, -24.+10.j, False]
            #self.pow2 = [pow(5+1j,-1.5), pow(3+3j,-1.5), pow(1+5j,-1.5)]
            self.pow2 = [0.08307064054041229214-0.0253416052125975132j,
                         0.04379104225017853491-0.1057209281108342370j,
                        -0.04082059235165559671-0.0766590341356157206j, False]
            
            #self.abs = [pow(26,.5), 3*pow(2,.5), pow(26,.5)]
            self.abs = [5.09901951359278483003,
                        4.24264068711928514641,
                        5.09901951359278483003, False]
            
            self.real = [5,3,1, True]
            self.imag = [1, 3, 5, True]
            self.conj = [5-1j, 3-3j, 1-5j, True]
            
            self.sum = 9+9j
            
            self.dot = 52+108j
        self.min = 1
        self.max = 5
            
    def test_mul(self):
        with self.context:
            # CPU with CPU
            c = self.a1 * self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.mul)
            # CPU with Other
            c = self.a2 * self.b1
            self.checkScheme(self.a2, self.b1, self.s, c, self.mul)
            # Other with Other
            c = self.a1 * self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.mul)
            # Other with CPU
            c = self.a1 * self.b2
            self.checkScheme(self.a1, self.b2, self.s, c, self.mul)
            # CPU with scalar
            c = self.a3 * self.s
            self.checkScheme(self.a3, self.b1, self.s, c, self.mul_s)
            # GPU with scalar
            c = self.a1 * self.s
            self.checkScheme(self.a1, self.b1, self.s, c, self.mul_s)
            
            self.assertRaises(TypeError, self.a1.__mul__, self.bad)
            self.assertRaises(ValueError, self.a1.__mul__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1 * self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.mul)
        # CPU with Other
        c = self.a1 * self.b2
        self.checkCPU(self.a1, self.b2, self.s, c, self.mul)
        # CPU with CPU
        c = self.a1 * self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.mul)
        # Other with CPU
        c = self.a2 * self.b1
        self.checkCPU(self.a2, self.b1, self.s, c, self.mul)
        # Other with scalar
        c = self.a3 * self.s
        self.checkCPU(self.a3, self.b1, self.s, c, self.mul_s)
        # CPU with scalar
        c = self.a1 * self.s
        self.checkCPU(self.a1, self.b1, self.s, c, self.mul_s)
        
    def test_rmul(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__rmul__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.mul)
            # CPU with Other
            c = self.a2.__rmul__(self.b1)
            self.checkScheme(self.a2, self.b1, self.s, c, self.mul)
            # Other with Other
            c = self.a1.__rmul__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.mul)
            # Other with CPU
            c = self.a1.__rmul__(self.b2)
            self.checkScheme(self.a1, self.b2, self.s, c, self.mul)
            # CPU with scalar
            c = self.s * self.a3
            self.checkScheme(self.a3, self.b1, self.s, c, self.mul_s)
            # GPU with scalar
            c = self.s * self.a1
            self.checkScheme(self.a1, self.b1, self.s, c, self.mul_s)
            
            self.assertRaises(TypeError, self.a1.__rmul__, self.bad)
            self.assertRaises(ValueError, self.a1.__rmul__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1.__rmul__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.mul)
        # CPU with Other
        c = self.a1.__rmul__(self.b2)
        self.checkCPU(self.a1, self.b2, self.s, c, self.mul)
        # CPU with CPU
        c = self.a1.__rmul__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.mul)
        # Other with CPU
        c = self.a2.__rmul__(self.b1)
        self.checkCPU(self.a2, self.b1, self.s, c, self.mul)
        # Other with scalar
        c = self.s * self.a3
        self.checkCPU(self.a3, self.b1, self.s, c, self.mul_s)
        # CPU with scalar
        c = self.s * self.a1
        self.checkCPU(self.a1, self.b1, self.s, c, self.mul_s)
                
    def test_imul(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 *= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c1, self.mul)
                # CPU with Other
                c2 *= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c2, self.mul)
                # Other with Other
                c4 *= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c4, self.mul)
                # Other with CPU
                c5 *= self.b2
                self.checkScheme(self.a1, self.b2, self.s, c5, self.mul)
                # CPU with scalar
                c3 *= self.s
                self.checkScheme(self.a1, self.b1, self.s, c3, self.mul_s)
                # GPU with scalar
                c6 *= self.s
                self.checkScheme(self.a1, self.b1, self.s, c6, self.mul_s)
                
                self.assertRaises(TypeError, self.a1.__imul__, self.bad)
                self.assertRaises(ValueError, self.a1.__imul__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Other Array and going back to the CPU
            # Other with Other
            c1 *= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c1, self.mul)
            # CPU with Other
            c4 *= self.b2
            self.checkCPU(self.a1, self.b2, self.s, c4, self.mul)
            # CPU with CPU
            c5 *= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c5, self.mul)
            # Other with CPU
            c2 *= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c2, self.mul)
            # Other with scalar
            c3 *= self.s
            self.checkCPU(self.a1, self.b1, self.s, c3, self.mul_s)
            # CPU with scalar
            c6 *= self.s
            self.checkCPU(self.a1, self.b1, self.s, c6, self.mul_s)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__imul__,self.s)
                self.assertRaises(TypeError, self.a1.__imul__,self.b1)
            
    def test_add(self):
        with self.context:
            # CPU with CPU
            c = self.a1 + self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.add)
            # CPU with Other
            c = self.a2 + self.b1
            self.checkScheme(self.a2, self.b1, self.s, c, self.add)
            # Other with Other
            c = self.a1 + self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.add)
            # Other with CPU
            c = self.a1 + self.b2
            self.checkScheme(self.a1, self.b2, self.s, c, self.add)
            # CPU with scalar
            c = self.a3 + self.s
            self.checkScheme(self.a3, self.b1, self.s, c, self.add_s)
            # GPU with scalar
            c = self.a1 + self.s
            self.checkScheme(self.a1, self.b1, self.s, c, self.add_s)
            
            self.assertRaises(TypeError, self.a1.__add__, self.bad)
            self.assertRaises(ValueError, self.a1.__add__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1 + self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.add)
        # CPU with Other
        c = self.a1 + self.b2
        self.checkCPU(self.a1, self.b2, self.s, c, self.add)
        # CPU with CPU
        c = self.a1 + self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.add)
        # Other with CPU
        c = self.a2 + self.b1
        self.checkCPU(self.a2, self.b1, self.s, c, self.add)
        # Other with scalar
        c = self.a3 + self.s
        self.checkCPU(self.a3, self.b1, self.s, c, self.add_s)
        # CPU with scalar
        c = self.a1 + self.s
        self.checkCPU(self.a1, self.b1, self.s, c, self.add_s)
        
    def test_radd(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__radd__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.add)
            # CPU with Other
            c = self.a2.__radd__(self.b1)
            self.checkScheme(self.a2, self.b1, self.s, c, self.add)
            # Other with Other
            c = self.a1.__radd__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.add)
            # Other with CPU
            c = self.a1.__radd__(self.b2)
            self.checkScheme(self.a1, self.b2, self.s, c, self.add)
            # CPU with scalar
            c = self.s + self.a3
            self.checkScheme(self.a3, self.b1, self.s, c, self.add_s)
            # GPU with scalar
            c = self.s + self.a1
            self.checkScheme(self.a1, self.b1, self.s, c, self.add_s)
            
            self.assertRaises(TypeError, self.a1.__radd__, self.bad)
            self.assertRaises(ValueError, self.a1.__radd__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1.__radd__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.add)
        # CPU with Other
        c = self.a1.__radd__(self.b2)
        self.checkCPU(self.a1, self.b2, self.s, c, self.add)
        # CPU with CPU
        c = self.a1.__radd__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.add)
        # Other with CPU
        c = self.a2.__radd__(self.b1)
        self.checkCPU(self.a2, self.b1, self.s, c, self.add)
        # Other with scalar
        c = self.s + self.a3
        self.checkCPU(self.a3, self.b1, self.s, c, self.add_s)
        # CPU with scalar
        c = self.s + self.a1
        self.checkCPU(self.a1, self.b1, self.s, c, self.add_s)
        
    def test_iadd(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 += self.b1
                self.checkScheme(self.a1, self.b1, self.s, c1, self.add)
                # CPU with Other
                c2 += self.b1
                self.checkScheme(self.a1, self.b1, self.s, c2, self.add)
                # Other with Other
                c4 += self.b1
                self.checkScheme(self.a1, self.b1, self.s, c4, self.add)
                # Other with CPU
                c5 += self.b2
                self.checkScheme(self.a1, self.b2, self.s, c5, self.add)
                # CPU with scalar
                c3 += self.s
                self.checkScheme(self.a1, self.b1, self.s, c3, self.add_s)
                # GPU with scalar
                c6 += self.s
                self.checkScheme(self.a1, self.b1, self.s, c6, self.add_s)
                
                self.assertRaises(TypeError, self.a1.__iadd__, self.bad)
                self.assertRaises(ValueError, self.a1.__iadd__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Other Array and going back to the CPU
            # Other with Other
            c1 += self.b1
            self.checkCPU(self.a1, self.b1, self.s, c1, self.add)
            # CPU with Other
            c4 += self.b2
            self.checkCPU(self.a1, self.b2, self.s, c4, self.add)
            # CPU with CPU
            c5 += self.b1
            self.checkCPU(self.a1, self.b1, self.s, c5, self.add)
            # Other with CPU
            c2 += self.b1
            self.checkCPU(self.a1, self.b1, self.s, c2, self.add)
            # Other with scalar
            c3 += self.s
            self.checkCPU(self.a1, self.b1, self.s, c3, self.add_s)
            # CPU with scalar
            c6 += self.s
            self.checkCPU(self.a1, self.b1, self.s, c6, self.add_s)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__iadd__,self.s)
                self.assertRaises(TypeError, self.a1.__iadd__,self.b1)
    
    def test_div(self):
        with self.context:
            # CPU with CPU
            c = self.a1 / self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.div)
            # CPU with Other
            c = self.a2 / self.b1
            self.checkScheme(self.a2, self.b1, self.s, c, self.div)
            # Other with Other
            c = self.a1 / self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.div)
            # Other with CPU
            c = self.a1 / self.b2
            self.checkScheme(self.a1, self.b2, self.s, c, self.div)
            # CPU with scalar
            c = self.a3 / self.s
            self.checkScheme(self.a3, self.b1, self.s, c, self.div_s)
            # GPU with scalar
            c = self.a1 / self.s
            self.checkScheme(self.a1, self.b1, self.s, c, self.div_s)
            
            self.assertRaises(TypeError, self.a1.__div__, self.bad)
            self.assertRaises(ValueError, self.a1.__div__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1 / self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.div)
        # CPU with Other
        c = self.a1 / self.b2
        self.checkCPU(self.a1, self.b2, self.s, c, self.div)
        # CPU with CPU
        c = self.a1 / self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.div)
        # Other with CPU
        c = self.a2 / self.b1
        self.checkCPU(self.a2, self.b1, self.s, c, self.div)
        # Other with scalar
        c = self.a3 / self.s
        self.checkCPU(self.a3, self.b1, self.s, c, self.div_s)
        # CPU with scalar
        c = self.a1 / self.s
        self.checkCPU(self.a1, self.b1, self.s, c, self.div_s)
        
    def test_rdiv(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__rdiv__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.rdiv)
            # CPU with Other
            c = self.a2.__rdiv__(self.b1)
            self.checkScheme(self.a2, self.b1, self.s, c, self.rdiv)
            # Other with Other
            c = self.a1.__rdiv__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.rdiv)
            # Other with CPU
            c = self.a1.__rdiv__(self.b2)
            self.checkScheme(self.a1, self.b2, self.s, c, self.rdiv)
            # CPU with scalar
            c = self.s / self.a3
            self.checkScheme(self.a3, self.b1, self.s, c, self.rdiv_s)
            # GPU with scalar
            c = self.s / self.a1
            self.checkScheme(self.a1, self.b1, self.s, c, self.rdiv_s)
            
            self.assertRaises(TypeError, self.a1.__rdiv__, self.bad)
            self.assertRaises(ValueError, self.a1.__rdiv__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1.__rdiv__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.rdiv)
        # CPU with Other
        c = self.a1.__rdiv__(self.b2)
        self.checkCPU(self.a1, self.b2, self.s, c, self.rdiv)
        # CPU with CPU
        c = self.a1.__rdiv__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.rdiv)
        # Other with CPU
        c = self.a2.__rdiv__(self.b1)
        self.checkCPU(self.a2, self.b1, self.s, c, self.rdiv)
        # Other with scalar
        c = self.s / self.a3
        self.checkCPU(self.a3, self.b1, self.s, c, self.rdiv_s)
        # CPU with scalar
        c = self.s / self.a1
        self.checkCPU(self.a1, self.b1, self.s, c, self.rdiv_s)
        
    def test_idiv(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 /= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c1, self.div)
                # CPU with Other
                c2 /= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c2, self.div)
                # Other with Other
                c4 /= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c4, self.div)
                # Other with CPU
                c5 /= self.b2
                self.checkScheme(self.a1, self.b2, self.s, c5, self.div)
                # CPU with scalar
                c3 /= self.s
                self.checkScheme(self.a1, self.b1, self.s, c3, self.div_s)
                # GPU with scalar
                c6 /= self.s
                self.checkScheme(self.a1, self.b1, self.s, c6, self.div_s)
                
                self.assertRaises(TypeError, self.a1.__idiv__, self.bad)
                self.assertRaises(ValueError, self.a1.__idiv__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Other Array and going back to the CPU
            # Other with Other
            c1 /= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c1, self.div)
            # CPU with Other
            c4 /= self.b2
            self.checkCPU(self.a1, self.b2, self.s, c4, self.div)
            # CPU with CPU
            c5 /= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c5, self.div)
            # Other with CPU
            c2 /= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c2, self.div)
            # Other with scalar
            c3 /= self.s
            self.checkCPU(self.a1, self.b1, self.s, c3, self.div_s)
            # CPU with scalar
            c6 /= self.s
            self.checkCPU(self.a1, self.b1, self.s, c6, self.div_s)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__idiv__,self.s)
                self.assertRaises(TypeError, self.a1.__idiv__,self.b1)
            
    def test_sub(self):
        with self.context:
            # CPU with CPU
            c = self.a1 - self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.sub)
            # CPU with Other
            c = self.a2 - self.b1
            self.checkScheme(self.a2, self.b1, self.s, c, self.sub)
            # Other with Other
            c = self.a1 - self.b1
            self.checkScheme(self.a1, self.b1, self.s, c, self.sub)
            # Other with CPU
            c = self.a1 - self.b2
            self.checkScheme(self.a1, self.b2, self.s, c, self.sub)
            # CPU with scalar
            c = self.a3 - self.s
            self.checkScheme(self.a3, self.b1, self.s, c, self.sub_s)
            # GPU with scalar
            c = self.a1 - self.s
            self.checkScheme(self.a1, self.b1, self.s, c, self.sub_s)
            
            self.assertRaises(TypeError, self.a1.__sub__, self.bad)
            self.assertRaises(ValueError, self.a1.__sub__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1 - self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.sub)
        # CPU with Other
        c = self.a1 - self.b2
        self.checkCPU(self.a1, self.b2, self.s, c, self.sub)
        # CPU with CPU
        c = self.a1 - self.b1
        self.checkCPU(self.a1, self.b1, self.s, c, self.sub)
        # Other with CPU
        c = self.a2 - self.b1
        self.checkCPU(self.a2, self.b1, self.s, c, self.sub)
        # Other with scalar
        c = self.a3 - self.s
        self.checkCPU(self.a3, self.b1, self.s, c, self.sub_s)
        # CPU with scalar
        c = self.a1 - self.s
        self.checkCPU(self.a1, self.b1, self.s, c, self.sub_s)
        
    def test_rsub(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__rsub__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.rsub)
            # CPU with Other
            c = self.a2.__rsub__(self.b1)
            self.checkScheme(self.a2, self.b1, self.s, c, self.rsub)
            # Other with Other
            c = self.a1.__rsub__(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.rsub)
            # Other with CPU
            c = self.a1.__rsub__(self.b2)
            self.checkScheme(self.a1, self.b2, self.s, c, self.rsub)
            # CPU with scalar
            c = self.s - self.a3
            self.checkScheme(self.a3, self.b1, self.s, c, self.rsub_s)
            # GPU with scalar
            c = self.s - self.a1
            self.checkScheme(self.a1, self.b1, self.s, c, self.rsub_s)
            
            self.assertRaises(TypeError, self.a1.__rsub__, self.bad)
            self.assertRaises(ValueError, self.a1.__rsub__, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1.__rsub__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.rsub)
        # CPU with Other
        c = self.a1.__rsub__(self.b2)
        self.checkCPU(self.a1, self.b2, self.s, c, self.rsub)
        # CPU with CPU
        c = self.a1.__rsub__(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.rsub)
        # Other with CPU
        c = self.a2.__rsub__(self.b1)
        self.checkCPU(self.a2, self.b1, self.s, c, self.rsub)
        # Other with scalar
        c = self.s - self.a3
        self.checkCPU(self.a3, self.b1, self.s, c, self.rsub_s)
        # CPU with scalar
        c = self.s - self.a1
        self.checkCPU(self.a1, self.b1, self.s, c, self.rsub_s)
        
    def test_isub(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 -= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c1, self.sub)
                # CPU with Other
                c2 -= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c2, self.sub)
                # Other with Other
                c4 -= self.b1
                self.checkScheme(self.a1, self.b1, self.s, c4, self.sub)
                # Other with CPU
                c5 -= self.b2
                self.checkScheme(self.a1, self.b2, self.s, c5, self.sub)
                # CPU with scalar
                c3 -= self.s
                self.checkScheme(self.a1, self.b1, self.s, c3, self.sub_s)
                # GPU with scalar
                c6 -= self.s
                self.checkScheme(self.a1, self.b1, self.s, c6, self.sub_s)
                
                self.assertRaises(TypeError, self.a1.__isub__, self.bad)
                self.assertRaises(ValueError, self.a1.__isub__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Other Array and going back to the CPU
            # Other with Other
            c1 -= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c1, self.sub)
            # CPU with Other
            c4 -= self.b2
            self.checkCPU(self.a1, self.b2, self.s, c4, self.sub)
            # CPU with CPU
            c5 -= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c5, self.sub)
            # Other with CPU
            c2 -= self.b1
            self.checkCPU(self.a1, self.b1, self.s, c2, self.sub)
            # Other with scalar
            c3 -= self.s
            self.checkCPU(self.a1, self.b1, self.s, c3, self.sub_s)
            # CPU with scalar
            c6 -= self.s
            self.checkCPU(self.a1, self.b1, self.s, c6, self.sub_s)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__isub__,self.s)
                self.assertRaises(TypeError, self.a1.__isub__,self.b1)
        
    def test_pow(self):
        with self.context:
            self.b1 *= 1
            # From CPU
            c1 = self.a1 ** 2
            c2 = self.a2 ** -1.5
            
            self.checkScheme(self.a1, self.b1, self.s, c1, self.pow1)
            self.checkScheme(self.a2, self.b1, self.s, c2, self.pow2)
        # From Other
        c1 = self.a1 ** 2
        c2 = self.a2 ** -1.5
        
        self.checkCPU(self.a1, self.b2, self.s, c1, self.pow1)
        self.checkCPU(self.a2, self.b2, self.s, c2, self.pow2)
        
    def test_abs(self):
        # We want to check that absolute value behaves correctly no matter
        # what quadrant it's in. First we will check with cpu arrays
        t1 = self.a1 * 1
        t2 = self.a1 * -1
        t3 = self.a1 * 1j
        t4 = self.a1 * -1j
        with self.context:
            self.b1 *= 1
            c1 = abs(t1)
            c2 = abs(t2)
            c3 = abs(t3)
            c4 = abs(t4)
            
            self.checkScheme(t1, self.b1, self.s, c1, self.abs)
            self.checkScheme(t1, self.b1, self.s, c2, self.abs)
            self.checkScheme(t1, self.b1, self.s, c3, self.abs)
            self.checkScheme(t1, self.b1, self.s, c4, self.abs)
            # Now coming from the current scheme
            c1 = abs(t1)
            c2 = abs(t2)
            c3 = abs(t3)
            c4 = abs(t4)
            
            self.checkScheme(t1, self.b1, self.s, c1, self.abs)
            self.checkScheme(t1, self.b1, self.s, c2, self.abs)
            self.checkScheme(t1, self.b1, self.s, c3, self.abs)
            self.checkScheme(t1, self.b1, self.s, c4, self.abs)
        # Now taking abs of the current scheme on the CPU
        c1 = abs(t1)
        c2 = abs(t2)
        c3 = abs(t3)
        c4 = abs(t4)
        
        self.checkCPU(t1, self.b2, self.s, c1, self.abs)
        self.checkCPU(t1, self.b2, self.s, c2, self.abs)
        self.checkCPU(t1, self.b2, self.s, c3, self.abs)
        self.checkCPU(t1, self.b2, self.s, c4, self.abs)
        #And finally, from the CPU to the CPU
        c1 = abs(t1)
        c2 = abs(t2)
        c3 = abs(t3)
        c4 = abs(t4)
        
        self.checkCPU(t1, self.b2, self.s, c1, self.abs)
        self.checkCPU(t1, self.b2, self.s, c2, self.abs)
        self.checkCPU(t1, self.b2, self.s, c3, self.abs)
        self.checkCPU(t1, self.b2, self.s, c4, self.abs)
            
        
    def test_real(self):
        with self.context:
            self.b1 *= 1
            # From CPU
            c = self.a1.real()
            self.checkScheme(self.a1, self.b1, self.s, c, self.real)
            
            # From Other
            c = self.a1.real()
            self.checkScheme(self.a1, self.b1, self.s, c, self.real)
        # Now on the CPU, from Other
        c = self.a1.real()
        self.checkCPU(self.a1, self.b2, self.s, c, self.real)
        # And finally CPU on the CPU
        c = self.a1.real()
        self.checkCPU(self.a1, self.b2, self.s, c, self.real)
            

        
    def test_imag(self):
        with self.context:
            self.b1 *= 1
            # From CPU
            c = self.a1.imag()
            self.checkScheme(self.a1, self.b1, self.s, c, self.imag)
            
            # From Other
            c = self.a1.imag()
            self.checkScheme(self.a1, self.b1, self.s, c, self.imag)
        # Now on the CPU, from Other
        c = self.a1.imag()
        self.checkCPU(self.a1, self.b2, self.s, c, self.imag)
        # And finally CPU on the CPU
        c = self.a1.imag()
        self.checkCPU(self.a1, self.b2, self.s, c, self.imag)
        
    def test_conj(self):
        with self.context:
            self.b1 *= 1
            # From CPU
            c = self.a1.conj()
            self.checkScheme(self.a1, self.b1, self.s, c, self.conj)
            
            # From Other
            c = self.a1.conj()
            self.checkScheme(self.a1, self.b1, self.s, c, self.conj)
        # Now on the CPU, from Other
        c = self.a1.conj()
        self.checkCPU(self.a1, self.b2, self.s, c, self.conj)
        # And finally CPU on the CPU
        c = self.a1.conj()
        self.checkCPU(self.a1, self.b2, self.s, c, self.conj)
            
    def test_sum(self):
        with self.context:
            self.b1 *= 1
            # From CPU
            c = self.a1.sum()
            self.checkScheme(self.a1, self.b1, self.s, c, self.sum)
            
            # From Other
            c = self.a1.sum()
            self.checkScheme(self.a1, self.b1, self.s, c, self.sum)
        # Now on the CPU, from Other
        c = self.a1.sum()
        self.checkCPU(self.a1, self.b2, self.s, c, self.sum)
        # And finally CPU on the CPU
        c = self.a1.sum()
        self.checkCPU(self.a1, self.b2, self.s, c, self.sum)
            
    def test_dot(self):
        with self.context:
            # CPU with CPU
            c = self.a1.dot(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.dot)
            # CPU with Other
            c = self.a2.dot(self.b1)
            self.checkScheme(self.a2, self.b1, self.s, c, self.dot)
            # Other with Other
            c = self.a1.dot(self.b1)
            self.checkScheme(self.a1, self.b1, self.s, c, self.dot)
            # Other with CPU
            c = self.a1.dot(self.b2)
            self.checkScheme(self.a1, self.b2, self.s, c, self.dot)
            
            self.assertRaises(TypeError, self.a1.dot, self.bad)
            self.assertRaises(ValueError, self.a1.dot, self.bad2)
            
        # Now taking Other Array and going back to the CPU
        # Other with Other
        c = self.a1.dot(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.dot)
        # CPU with Other
        c = self.a1.dot(self.b2)
        self.checkCPU(self.a1, self.b2, self.s, c, self.dot)
        # CPU with CPU
        c = self.a1.dot(self.b1)
        self.checkCPU(self.a1, self.b1, self.s, c, self.dot)
        # Other with CPU
        c = self.a2.dot(self.b1)
        self.checkCPU(self.a2, self.b1, self.s, c, self.dot)
    
    def test_max(self):
        if self.kind == 'real':
            with self.context:
                self.b1 *= 1
                # From CPU
                c = self.a1.max()
                self.checkScheme(self.a1, self.b1, self.s, c, self.max)
                # From Other
                c = self.a1.max()
                self.checkScheme(self.a1, self.b1, self.s, c, self.max)
            # From Other
            c = self.a1.max()
            self.checkCPU(self.a1, self.b2, self.s, c, self.max)
            # From CPU
            c = self.a1.max()
            self.checkCPU(self.a1, self.b2, self.s, c, self.max)

    def test_min(self):
        if self.kind == 'real':
            with self.context:
                self.b1 *= 1
                # From CPU
                c = self.a1.min()
                self.checkScheme(self.a1, self.b1, self.s, c, self.min)
                # From Other
                c = self.a1.min()
                self.checkScheme(self.a1, self.b1, self.s, c, self.min)
            # From Other
            c = self.a1.min()
            self.checkCPU(self.a1, self.b2, self.s, c, self.min)
            # From CPU
            c = self.a1.min()
            self.checkCPU(self.a1, self.b2, self.s, c, self.min)
                
