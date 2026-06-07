To simulate an electrode in SIMION, a potential array (PA) must be created, either graphically or by means of coding a .GEM file

A .GEM file is the preferred method of creating electrodes because the graphical method is limited to planar, rectangular, or circular symmetry

A .GEM file gives full control over the location and shape of electrodes

After opening SIMION, click "New PA" and "From .GEM File"

Click "Save" and rename the file with extension .pa# which creates a separate numbered .pa file for each electrode programmed in the .GEM file

Click "Refine"

Click "Fast Adjust" to assign voltages to each PA

Then click "View/Open Ion Workbench" and the model is ready to be simulated.

In the .fly2 file, there is a variable at the top of the file that allows you to change the beam energy. It doesn't show up in SIMION when you try to edit the .fly2 file as text so you have to open the .fly2 file in something like Notepad to change the beam energy.
