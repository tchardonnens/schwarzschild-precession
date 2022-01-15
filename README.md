# schwarzschild-precession
Simulate the trajectory of the star S2 in orbit around the black hole at the center of our galaxy and display in real time the value of the shift between the frequency of the received electromagnetic wave and the emitted wave (distinguishing between the Einstein effect, the Doppler effect of the first and second order).

Resources to create this simulation : 

- Nobel Price 2020, Detection of the Schwarzschild precession in the orbit of the star S2 near the Galactic centre massive black hole : https://arxiv.org/pdf/2004.07187.pdf

- Motion of test particles in the gravitational field of a central fixed mass M : https://en.wikipedia.org/wiki/Schwarzschild_geodesics

<h2>How to make it run?</h2>

The main file is a jupyter notebook. I created a conda env in order to get to work faster ;-)

Install Anaconda:

https://docs.anaconda.com/anaconda/install/windows/
</br>
For Apple Silicon I use Miniconda:
<code>brew install miniforge</code>

Jupyter Notebook: 

https://jupyter.org/install

Install Latex: 

For Windows: https://miktex.org/download
</br>
For MacOS: https://tug.org/mactex/

Commands:
- create the conda env (you have to be in the folder where the .yml file is): 
<code>conda env create -f precession-env.yml</code>

- activate the env to see '(precession2)' before your terminal line:
<code> conda activate precession2 </code>

- spin up the Jupyter Notebook:
<code>jupyter notebook</code>

Don't forget to change the kernel on your Jupyter Notebook to precession2.
https://queirozf.com/entries/jupyter-kernels-how-to-add-change-remove


