all:YUV_comp.cpp
	g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` YUV_comp.cpp -o _YUVcomp`python3-config --extension-suffix`
test:
	pytest
run:
	pytest
clean:
	rm -f *.so; rm -rf __pycache__; rm -rf AKIYO_352x288_10_cp.yuv; rm -rf AKIYO_352x288_10_cp_de.yuv
