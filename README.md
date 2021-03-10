## PHYS5394-alejandro

This repository contains classwork for *Statistical Methods*, PHYS 5394.

# Description

* **lab2** -1/27/2021
  - gensgsig.m     - Generate a sine-gaussian signal
  - testgensgsig.m - Test the sine-gaussian (sg) signal and plot it

* **lab4** -2/11/2021
  - gensgsig.m     - No changes from previous lab
  - testgensgsig.m - Sample sg signal at 5x and 0.5x maximum frequency using new function genPlotSig() and plot.

* **lab5** -2/11/2021
  - gensgsig.m     - No changes
  - testgensgsig.m - Perform FFT of sg signal and plot periodogram.
  - filtering.m    - Generate 3 sinusoids at different frequencies, add them, filter each, then perform FFT and compare plots of original vs filtered signal periodograms.

* **lab6** -2/19/2021
  - gensgsig.m     - No changes
  - testgensgsig.m - Added spectrogram computation and plot. Window length and overlap chosen at 0.08 and 0.07 s to better display frequency spread.

* **lab7** -2/23/2021
  - customrand.m     - Generate a pseudo-random number between two open boundaries from uniform distribution.
  - customrandn.m    - Gemerate a pseudo-random number by custom standard deviation and mean from normal distribution.
  - testcustomprng.m - Generate the two probability distribution functions and plot the distribution and histogram.

* **lab8** -3/10/2021
  - whiteningDemo.m  - Reads testData.txt, estimates PSD of signal free segment, designs whitening filter and applies it to entire data. Plots before/after data and spectrograms.

Contact alejandro.reyes01@utrgv.edu
