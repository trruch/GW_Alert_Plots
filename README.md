# GW_Alert_Plots
I'll format this better later on, for now just want to explain how to run the scripts:

### Skymap_Alert_Plot.py
This is meant to be run as soon as the alert is sent out and skymap made available. Eventually, I'll talk to Andre/Nora and get this integrated into their alerts that get sent to Slack so that the script is run automatically and the plot sent to everyone in the channel. To run from the terminal:
1. ``` $ python Skymap_Alert_Plot.py ```
2. Input the url or path to the skymap when prompted
    - ``` Skymap Url (or local path): https://gracedb.ligo.org/api/superevents/S230615az/files/bayestar.fits.gz ```
    - Make sure the url is NOT the multi-order skymap! 
4. Input the Event Name when prompted
    - ``` Event Name: S230615az ```

A png will be saved in the same directory which shows the skymap, 50% and 90% contour regions, name of the event, location of and distance to the max probability pixel.

### Skymap_Hexes_Plot.py
This is meant to be run after the strategy code as finished and an observing json file has been created.
To run from the terminal:
1. ``` $ python Skymap_Hexes_Plot.py ```
2. Input the url or path to the skymap when prompted
    - ``` Skymap Url (or local path): https://gracedb.ligo.org/api/superevents/S230615az/files/bayestar.fits.gz ```
    - Make sure the url is NOT the multi-order skymap! 
4. Input the Event Name when prompted
    - ``` Event Name: S230615az ```
5. Input the path to the json file when prompted
    - ```Json path: ./S230518h_test.json```

This will produce 2 images, saved as png's:
- The same skymap as before, this time overlayed with the DECam hexes (pointing of the telescope)
- A plot detailing how much of the total probability is covered with each hex 


Need to update this once I get the airmass plot up and running.
