#import os
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import ligo.skymap.plot
from matplotlib import pyplot as plt
import healpy as hp
import numpy as np
import ligo.skymap
import json
import pandas as pd
from astroplan import Observer, FixedTarget
from astroplan import (AltitudeConstraint, AirmassConstraint,AtNightConstraint)
from astroplan import is_observable
from astropy.time import Time
from astroplan.plots import plot_airmass
import sys


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
    
    return (area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob)

#######################
def ra_dec2theta_phi(ra,dec):
    theta = 0.5 * np.pi - np.pi*dec/180
    phi = np.deg2rad(ra)
    return theta, phi

def get_prob_from_observing_json(json_data, prob_array):
    hex_number = []
    total_prob = []
    for i in range(len(data)):
        ra = data[i]['RA']
        dec = data[i]['dec']
        theta_hex, phi_hex = ra_dec2theta_phi(ra,dec)
        vec = hp.ang2vec(theta_hex, phi_hex)
        decam_hex_disc = hp.query_disc(nside, vec, radius=np.radians(0.9772), nest = True)
        prob_covered = np.sum(prob[decam_hex_disc])
        hex_number.append(i)
        total_prob.append(prob_covered + np.sum(total_prob))

    prob_percent = [i*100 for i in total_prob]
    return hex_number, prob_percent


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
    
    

    fig,ax=plt.subplots(figsize=(12,8))
    ax.set_title('Key Target Airmass Plot')
    plot_airmass(targets,ax=ax,
            observer=tscope,
            time=tscope.midnight(Time.now()).to_datetime(timezone=tscope.timezone),
            use_local_tz=True,
            brightness_shading=True,
            max_region=3,min_region=1.5)
    plt.legend(loc='best')
    plt.savefig(event_name+'_Airmass_Hexes',dpi=300, bbox_inches = "tight")
#######################

if __name__ == "__main__":
  
    
    url = input('Skymap Url (or local path): ')
    name = input('Event Name: ')
    jsonloc = input('Json path: ')
    
    f = open(jsonloc)
    data = json.load(f)

    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)

    airmass_input = [(data[i]['RA'], data[i]['dec'], 'Hex_'+str(i)) for i in range(len(data))]
    airmass(name, airmass_input)

    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame


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

    hex_centers = []
    for i in range(len(data)):
        ra = data[i]['RA']
        dec = data[i]['dec']
        centers = SkyCoord(ra, dec, unit="deg")
        hex_centers.append(centers)
        
        
    ax = plt.axes(
        #[0.05, 0.05, 0.9, 0.9],
        projection='astro hours mollweide')
    
    ax_inset = plt.axes(
        [0.9, 0.2, 0.2, 0.2],
        projection='astro zoom',
        center=hex_centers[0],
        radius=25*u.deg)
    
    ax_hexes = []
    for i in range(len(data)):
        axis = plt.axes(projection='astro degrees zoom')
        ra = data[i]['RA']
        dec = data[i]['dec']
        centers = SkyCoord(ra, dec, unit="deg")
        ax_hexes.append(axis)
        
        
        ax.mark_inset_circle(axis, centers, '.9772 deg', edgecolor = 'green')
        ax_inset.mark_inset_circle(axis, centers, '.9772 deg', edgecolor = 'green')
        axis.remove()


        
        

   
    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(True)
        ax_inset.coords[key].set_ticks_visible(True)
        
        
    ax.grid()
    ax_inset.grid()
    ax.mark_inset_axes(ax_inset)
    ax.connect_inset_axes(ax_inset, 'upper left')
    ax.connect_inset_axes(ax_inset, 'lower left')
    #ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()
    #ax_inset.compass(0.9, 0.1, 0.2)

    ax.imshow_hpx(url, cmap='cylon')
    cs = ax.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])
    ct = ax_inset.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])

    ax_inset.imshow_hpx(url, cmap='cylon')
    plt.savefig(name + '_w_Hexes',dpi=300, bbox_inches = "tight")
   
    
    plt.figure()
    plt.bar(get_prob_from_observing_json(data, prob)[0], get_prob_from_observing_json(data, prob)[1], width = 1, edgecolor = 'k')
    plt.xlabel('Hex #')
    plt.ylabel(r'Cumulative $\%$ of Probability Covered')
    plt.title(name)
    plt.savefig(name+'_Cumulative_Hex_Prob',dpi=300, bbox_inches = "tight")
    
