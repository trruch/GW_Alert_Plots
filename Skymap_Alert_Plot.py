# Orig by T.Ruch
# Changes by MSS Gill
# 7-15-2023

#

#import os
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import ligo.skymap.plot
from matplotlib import pyplot as plt
import healpy as hp
import numpy as np
import ligo.skymap
import pandas as pd
from astroplan import Observer, FixedTarget
from astroplan import (AltitudeConstraint, AirmassConstraint,AtNightConstraint)
from astroplan import is_observable
from astropy.time import Time
from astroplan.plots import plot_airmass
import sys
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun
from astropy.coordinates import get_body
import datetime
import ephem


def make_alert_skymap(map_path):
    skymap = hp.read_map(map_path, field=None)
    prob, distmu, distsigma, distnorm = skymap
    
    npix = len(prob)
    nside = hp.npix2nside(npix)


    maxprobcoord_tup = hp.pix2ang(nside, np.argmax(prob))
    maxprobcoord = [0, 0]
    maxprobcoord[1] = np.rad2deg(0.5*np.pi-maxprobcoord_tup[0])
    maxprobcoord[0] = np.rad2deg(maxprobcoord_tup[1])
    

    def find_area(prob_array, contours, nside = nside):
        areas = []
        for contour in contours:
            sortedprob = np.sort(prob_array)
            probsum = 0
            probcutoff = 1
            npix = 0
            while probsum<contour:
                probsum = probsum+sortedprob[-1]
                probcutoff = sortedprob[-1]
                sortedprob = sortedprob[:-1]
                npix = npix+1

            area = npix * hp.nside2pixarea(nside, degrees=True)
            areas.append(area)
        return ([int(round(i,0)) for i in areas])
    
    
    def ci(level, sorted_prob=np.flip(np.sort(prob))):
        csum = 0
        c = 0
        index = 0
        while csum < level:
            csum += sorted_prob[index]
            c = sorted_prob[index]
            index += 1
        return c

    c90 = ci(.5)
    c50 = ci(.9)
    levels = [c50, c90]
        
    area50, area90 = find_area(prob, contours = [.5, .9])
    
    maxprob_ra = round(maxprobcoord[0],2)
    maxprob_dec = round(maxprobcoord[1],2)
    
    maxprob_dist = int(round(distmu[np.argmax(prob)],0))
    maxprob_distsigma = int(round(distsigma[np.argmax(prob)],0))
    
    return (area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels)

def airmass(event_name,target_coords):
    ''' Target coords is a list of tuples containing ra, dec, and the name of the target'''

    CTIO = Observer(longitude=-70.80*u.deg, latitude=-30.17*u.deg,
                  elevation=3000*u.m, name="CTIO",timezone='America/Santiago')

    tscope = CTIO ; tscope_str = 'CTIO'

    if len(target_coords) > 1:
        targets=[FixedTarget(coord=SkyCoord(ra=coords[0]*u.deg,dec=coords[1]*u.deg),name=coords[2]) for coords in target_coords]
    else:
        targets=FixedTarget(coord=SkyCoord(ra = target_coords[0][0]*u.deg,dec=target_coords[0][1]*u.deg),name='Max Prob Coord')

    constraints = [AltitudeConstraint(10*u.deg, 90*u.deg), AtNightConstraint.twilight_civil()]

    time_range=(Time.now(),Time.now()+1*u.day)
    ever_observable = is_observable(constraints, tscope, targets, time_range=time_range)
    print("Observable?'",ever_observable)
    
    

    fig,ax=plt.subplots(figsize=(12,8))
    ax.set_title('Key Target Airmass Plot')
    plot_airmass(targets,ax=ax,
            observer=tscope,
            time=tscope.midnight(Time.now()).to_datetime(timezone=tscope.timezone),
            use_local_tz=True,
            brightness_shading=True,
            max_region=3,min_region=1.5)
    plt.legend(loc='best')
    plt.savefig(event_name+'_Airmass',dpi=300, bbox_inches = "tight")


def moon(event_name, todays_date):
    date = datetime.date.today()
    m = ephem.Moon(date)
    phase = round(m.moon_phase, 2)
    
    plt.style.use(astropy_mpl_style)
    quantity_support()
    
    
    
    CTIO = EarthLocation(lat=-30.17*u.deg, lon=-70.80*u.deg, height=3000*u.m)
    utcoffset = -4*u.hour  # Eastern Daylight Time
    
    
    
    #midnight = Time('2012-7-13 00:00:00') - utcoffset
    mytime = todays_date + ' 00:00:00'
    midnight = Time(mytime) - utcoffset
    
    
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=CTIO)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
    
    
    
    moon_July12_to_13 = get_body("moon", times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)
    
    
    plt.plot(delta_midnight, moonaltazs_July12_to_13.alt, color='blue', ls='--', label='Moon')
    
    plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -0*u.deg, color='0.9', zorder=0)
    plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -6*u.deg, color='0.8', zorder=0)
    plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -12*u.deg, color='0.7', zorder=0)
    plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -18*u.deg, color='0.6', zorder=0)
    
    plt.text(-11.5, 75, 'Moon phase = {}'.format(phase),bbox=dict(facecolor='red', alpha=0.5))
    plt.title(date)
    plt.legend(loc='upper left')
    plt.xlim(-12*u.hour, 12*u.hour)
    plt.xticks((np.arange(13)*2-12)*u.hour)
    plt.ylim(0*u.deg, 90*u.deg)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    plt.savefig(event_name+'_Moon',dpi=300, bbox_inches = "tight")




if __name__ == "__main__":

    
    print(' \n  ********** Output will be three PNG file saved to disk of: ',
          '\n - The current moon phase, ',
          '\n - A skymap of the probability contours of the localization region of the event ',
          '\n - The airmass to the highest probability pixel of the event as a target, as a function of time \n\n',
         ' You can press return for the next 3 prompts and defaults will be used. \n')

    url = input('Skymap Url (or local path): ')
    name = input(' For plots, event name: ')
    date = input('Todays date (Ex: 2023-6-12): ')

    if url=="":
        url="https://gracedb.ligo.org/api/superevents/S230615az/files/bayestar.fits.gz"
        print("\n Defaulting to event = ",url)

    if name=="":
        name="S230615az"
        print("\n Defaulting to name =  ",name)
        
    if len(date)==0:
        date = datetime.date.today().strftime('%Y-%-m-%-d')
        print("\n Defaulting to today's date =  ",date)

    moon(name, date)
    
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels = make_alert_skymap(url)



    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame

    airmass(name, [(maxprob_ra, maxprob_dec)])

    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.annotate('Event Name: {}'.format(name) + '\n'
                 + r'50% Area: {} deg$^2$'.format(area50) + '\n'
                 + r'90% Area: {} deg$^2$'.format(area90) + '\n'
                 + r'Max Prob Coordinates: ({},{})'.format(maxprob_ra, maxprob_dec) + '\n'
                 + r'Max Prob Distance: {}$\pm${} Mpc'.format(maxprob_dist, maxprob_distsigma)
                 ,(0.9,0.8))
    plt.box(False)
    plt.xticks([])
    plt.yticks([])

    ax = plt.axes(projection='astro hours mollweide')

    ax_inset = plt.axes(
        [0.9, 0.2, 0.2, 0.2],
        projection='astro zoom',
        center=center,
        radius=10*u.deg)

    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(True)
        ax_inset.coords[key].set_ticks_visible(True)
   
    ax.grid()
    ax_inset.grid()
    ax.mark_inset_axes(ax_inset)
    ax.connect_inset_axes(ax_inset, 'upper right')
    ax.connect_inset_axes(ax_inset, 'lower left')
    ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()

    ax.imshow_hpx(url, cmap='cylon')
    cs = ax.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])
    ct = ax_inset.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])

    ax_inset.imshow_hpx(url, cmap='cylon')
    ax_inset.plot(
        maxprob_ra, maxprob_dec,
        transform=ax_inset.get_transform('world'),
        marker=ligo.skymap.plot.reticle(),
        markersize=30,
        markeredgewidth=3)
    ax.plot(
        maxprob_ra, maxprob_dec,
        transform=ax.get_transform('world'),
        marker=ligo.skymap.plot.reticle(inner=0),
        markersize=10,
        markeredgewidth=3)
    
    plt.savefig(name+'_initial_skymap',dpi=300, bbox_inches = "tight")
   
