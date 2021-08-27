# decon
Demonstrate and test TPC signal deconvolution.

The script [run] calls the Root scripr [run.C] which does fast simulation
(convolution) and reconstruction (deconvolution) of one track/event
for each of the three wire views. For help:
<pre>
> ./run -h
</pre>
The following command simulates and reconstructs 100 events with run number (and random seed) 101. The simulation uses 2D convolution with no noise and the reconstruction is 1D DFT deconvolution.
<pre>
> ./run 0.0 conv2d:4:10:0:0 deco1d:dft 101 100
</pre>
The output of the job including plots is copied to ./wsave and the command makeIndex (e.g. in ~dladams/bin) is used to construct html index files to facilitate viewin in a browser. Results are shown for one channel in each plane: 100 for u, 1100 for v and 2100 for collection. Plots include
1. All ROIS for each plane.
2. DFT power before and after deconvolution for each plane.
3. Overlaid, u, v and collection calibrated spectra for each event.
4. Distribution of reconstruction charge (one entry per event) for each plane.

Note that for the example here, the u and v reconstructed charge is much too large because the the 1D deconvolution does not properly account for contributions from neighboring channels.
