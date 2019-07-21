# Evolving Cellular Automata rules for image segmentation

This program shows how to evolve rules for a Cellular Automata capable of image segmentation.
Rules are evolved using Genetic Programming.

Each evolved program is a mathematical formula which contains pixel colors (from a given neighbourhood) as variables and 9 operators.
The evolved program is applied repeatedly for each pixel. 
At the end (when **max_num_CA_iterations_with_no_improvements** iterations with no improvements have been made) the obtained image is compared with the mask. The fitness is equal to the distance between those images.

Please note that evolution can take a lot of time (hours to days depending on the parameters).
The user may set the **number of threads** ( _params.num_threads_ ) and the **sub_population size** ( _params.sub_population_size_ ) in order to increase or reduce the computational effort.

## How to compile: 

Add [src/gp_ca.cpp](src/gp_ca.cpp) and [src/qdbmp.cpp](src/qdbmp.cpp) to your C++ console project. 
With _Visual C++_ you may want to add _CRT_SECURE_NO_WARNINGS to compiler settings in order to avoid warnings related to _printf_.

## Other libraries utilized

[QDBMP](http://qdbmp.sourceforge.net) is used for reading images from bmp files.

## Images

Test images are taken from [Berkeley Segmentation Dataset](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/BSDS300/html/dataset/images.html)

## Author: [Mihai Oltean](https://mihaioltean.github.io)

mihai.oltean@gmail.com

## Acknowledgement 

The author acknowledges of having discussions, on this topic, with _Anca Andreica_ during the development of this program.

## License: MIT