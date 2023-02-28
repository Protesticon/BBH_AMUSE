# BBH_AMUSE
Xueyang Hu  
Yi Kang  
Renhao Liang  

### How to run this package?

1. download all files in the same folder;

2. Running **Generate_BBH_seeds.py** to create soft BBHs that orbit around a SMBH.  
It also launches a Nbody simulation that consists of only SMBH and generated BBHs.  
The hardness variations of each BBH would be saved in Figures/Reference/.  
Note that the BBHs here are expected not to interfere with each other in our simulations, hence, if there are significant fluctuations in hardness variations, please consider running Generate_BBH_seeds.py again until all hardnesses ascend smoothly.
The final BBH Particles would be saved in **BBH0.txt**;

3. Then we can initialize the Nbody and Hydro-Gravity Bridge simulations of BBHs in galactic nuclei by **Nbody.py** and **Bridge.py** respectively.  
In **Nbody.py**, we can alternate the range of other particles' masses and particle number.  
In **Bridge.py**, we can alternate the total mass of gas disk and gas particle number, the gas density/distribution follows (Sch\¨odel, R., Eckart, A., Alexander, T., et al. 2007, A&A, 469, 125). And remember to change the result file's name as you like.  

Warning: the simulations of Nbody and Bridge would be VERY time-consuming, and easily getting interrupted, we recommend to copy these codes in Jupyter Notebook before running.

4. More analysis can be achieved using the output files, (e.g., **analysis.py**).

