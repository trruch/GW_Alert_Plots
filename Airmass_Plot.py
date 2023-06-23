# Orig: Daniel Warshofsky 5-18-2023  (U.Minn)
# Modified: M. Gill 5-19-2023  (U.Minn and FNAL)
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroplan import Observer, FixedTarget
from astroplan import (AltitudeConstraint, AirmassConstraint,AtNightConstraint)
from astroplan import is_observable
from astropy.time import Time
from astroplan.plots import plot_airmass
import matplotlib.pyplot as plt
import sys
​
if __name__=='__main__':
    gw_cand_list=sys.argv[1]
    data=pd.read_csv(gw_cand_list)
    print("Num galaxies in region = ",len(data))
​
    Turbo = Observer(longitude=-93.180333*u.deg, latitude=44.990659*u.deg,
                  elevation=242*u.m, name="Turbo",timezone='US/Central')
​
    CTIO = Observer(longitude=--70.80*u.deg, latitude=-30.17*u.deg,
                  elevation=3000*u.m, name="CTIO",timezone='America/Santiago')
​
    tscope = Turbo ; tscope_str="Turbo"
#    tscope = CTIO ; tscope_str="CTIO"
​
    targets=[FixedTarget(coord=SkyCoord(ra=row[1]['ra']*u.deg,dec=row[1]['dec']*u.deg),name=row[1]['objname']) for row in data.iterrows()]
    constraints = [AltitudeConstraint(10*u.deg, 90*u.deg), AtNightConstraint.twilight_civil()]
​
    time_range=(Time.now(),Time.now()+1*u.day)
    ever_observable = is_observable(constraints, tscope, targets, time_range=time_range)
    data["Observable"]=ever_observable
    observable_data=data[data["Observable"]==True]
    observable_data.to_csv(f'Observable_{gw_cand_list}')
​
    print("Num galaxies in region that are observable = ",len(observable_data))
        
    observable_targets=[FixedTarget(coord=SkyCoord(ra=row[1]['ra']*u.deg,dec=row[1]['dec']*u.deg),name=row[1]['objname']) for row in observable_data.iterrows()]
​
    print("Num galaxies in region of observable targets = ",len(observable_targets))
        
    fig,ax=plt.subplots(figsize=(12,8))
    ax.set_title(f'tscope_str {gw_cand_list.replace(".csv","")}')
    plot_airmass(observable_targets,ax=ax,
            observer=tscope,
            time=tscope.midnight(Time.now()).to_datetime(timezone=tscope.timezone),
            use_local_tz=True,
            brightness_shading=True,
            max_region=3,min_region=1.5)
    plt.legend(loc='best')
    plt.savefig(f'Airmass_plot_{gw_cand_list.replace(".csv",".png")}')
