#!/usr/bin/env python
import _YUVcomp
import unittest
import random
import math
import os

class VectorTestCase(unittest.TestCase):
    def test_mode_0(self):
    	filename = "AKIYO_352x288_10.yuv"
    	filename2 = "AKIYO_352x288_10_cp.yuv"
    	_YUVcomp.YUV_compress(filename,0,40)
    	a = os.path.getsize(filename)
    	b = os.path.getsize(filename2)
    	self.assertLess(b,a)
    def test_mode_2(self):
        filename = "AKIYO_352x288_10.yuv"
        filename2 = "AKIYO_352x288_10_cp.yuv"
        _YUVcomp.YUV_compress(filename,2,45)
        a = os.path.getsize(filename)
        b = os.path.getsize(filename2)
        self.assertLess(b,a)
    def test_mode_4(self):
        filename = "AKIYO_352x288_10.yuv"
        filename2 = "AKIYO_352x288_10_cp.yuv"
        _YUVcomp.YUV_compress(filename,4,40)
        a = os.path.getsize(filename)
        b = os.path.getsize(filename2)
        self.assertLess(b,a)
    def test_decompress_size(self):
        filename = "AKIYO_352x288_10.yuv"
        filename1 = "AKIYO_352x288_10_cp.yuv"
        filename2 = "AKIYO_352x288_10_cp_de.yuv"
        _YUVcomp.YUV_decompress(filename1)
        a = os.path.getsize(filename)
        b = os.path.getsize(filename2)
        self.assertEqual(b,a)
if __name__ == '__main__':
    unittest.main()
