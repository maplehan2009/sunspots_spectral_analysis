# Pre-requirement of python:
Install numpy, scipy and matplotlib

# Execution of the code
Change directory to the root folder where find the main.py. Then open your terminal / console, execute this command: python main.py

Or if you use the ipython, just run the main.py in the ipython. The results would be better printed.

# Other files
The other files are support functions like autocovariance, fourier transform, periodogram etc... More details can be found in the explanation part in these functions' python files. Just one thing bizarre: when I used the FFT function provided by numpy.fft, the size of the return array equals that of the input data, this means if the input data has just one element, the FFT will return also only one element in turn. (In fact in the point that lambda = 0). I am not satisfied with FFT of numpy so I wrote the fourier transform myself which is sft.py (slow fourier transform) here. 

# Reference of Time Series Analysis
The main reference is the book Time Series Analysis witten by Mr. François Roueff. In this project, we investigated the spectral analysis of the famous sunspot number series. More details can be found in the description_of_project.pdf.

Voilà, enjoy our code,
Jingtao and Peng
At Ecole Polytechnique, Palaiseau
