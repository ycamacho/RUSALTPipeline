#!/usr/bin/env python

##### the following UNIX command must be run before the python script is executed so that
##### the correct installation of pyraf (which has pysalt) is used
##### source /usr/local/astro64/login.csh

'''
The primary data structures used in this pipeline are Python dictionaries:
essentially lists whose elements are "keyword:value" pairs. These dictionaries
(for flats, flat-fielded science images, extracted arc spectra, and so on)
are generally indexed by the header 'GR-ANGLE' keyword value. This makes it 
easy to associate, say, the normalized flats which have some particular GR-ANGLE value
with their corresponding science, standard star, and arc images which have the same
GR-ANGLE value.
'''

import pyfits
import sys
import os
import shutil
import numpy as np
import ds9

from glob import glob
from pyraf import iraf

# importing iraf packages
from iraf import images
from iraf import imutil
from iraf import immatch
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import onedspec
from iraf import twodspec
from iraf import apextract
from iraf import longslit
from iraf import pysalt
from iraf import saltspec

# defining global variables (these shall be the dictionaries)
flats = {}
arcs = {}
standards = {}
sciences = {}
combflats = {}
normflats = {}
auxflatsciences = {}
auxflatarcs = {}
auxflatstandards = {}
flatsciences = {}
flatarcs = {}
flatstandards = {}
wavesols = {}
wavearcs = {}
wavesciences = {}
wavestandards = {}
auxbackgroundsciences = {}
backgroundsciences = {}
auxextractedsciences = {}
auxextractedstandards = {}
auxextractedarcs_std = {}
auxextractedarcs_sci = {}
extractedsciences = {}
extractedstandards = {}
extractedarcs_std = {}
extractedarcs_sci = {}
auxdispsciences = {}
auxdispstandards = {}
dispsciences = {}
dispstandards = {}
stdfiles = {}
sensfiles = {}
auxfluxsciences = {}
fluxsciences = {}

# user menu is at the end of this pipeline

# This function copies all the files in the current directory into a new subdirectory, and changes to that subdirectory.
def clone():
	# this sets the source and destination directories -- assuming all fits files are in current directory
	cwd = os.getcwd()
	dest = cwd+'/pipeline'
	# this begins the cloning process
	shutil.copytree(cwd,dest)
	# this changes the directory to the dest
	os.chdir(dest)
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Cloned and added files from '+cwd+' to '+dest
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

# This function calls all other functions for a full reduction of the data.
def fullreduce():
	# use the global dictionaries for operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# sort through all the fits files and create the initial four dictionaries: flats, arcs, standards, sciences
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Sorting through all the images using: dictionaries()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	dictionaries()
	# combine those flats that have the same GR-ANGLE header keyword value
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Combining flats with same GR-ANGLE using: combineflats()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	combineflats()
	# normalize all the combined flats
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Normalizing each combined flat using: normalizeflats()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	normalizeflats()
	# flat-field the science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flat-fielding each science image with its associated normalized flat using: flattensciences()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	flattensciences()
	# flat-field the arc images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flat-fielding each arc image with its associated normalized flat using: flattenarcs()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	flattenarcs()
	# flat-field the standard star images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flat-fielding each standard star image with its associated normalized flat using: flattenstandards()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	flattenstandards()
	# use pysalt.specidentify to produce wavelength solution files for each flat-fielded arc image
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Identifying the emission lines in each flat-fielded arc image using: specidentifyarcs()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	specidentifyarcs()
	# use wavelength solution files to wavelength-calibrate the 2D arc images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying wavelength solution to each 2D flat-fielded arc image using: wavecalarc()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	wavecalarc()
	# use wavelength solution files to wavelength-calibrate the 2D science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying wavelength solution to each 2D flat-fielded science image using: wavecalsci()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	wavecalsci()
	# use wavelength solution files to wavelength-calibrate the 2D standard star images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying wavelength solution to each 2D flat-fielded standard star image using: wavecalstd()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	wavecalstd()
	# use backgroundsciences() to background-subtract the 2D science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Subtracting background from each 2D wavelength-calibrated science image using: subtractbackground()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	subtractbackground()
	# extract the wavelength-corrected science spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracted each wavelength-calibrated science spectrum using: extractsciences()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	extractsciences()
	# extract the wavelength-corrected standard star spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracting each wavelength-calibrated standard star spectrum using: extractstandards()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	extractstandards()
	# extract the wavelength-corrected arc spectra (twice -- for science and for standard star images)
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracting each wavelength-calibrated arc spectrum twice (for science and for standard star) using: extractarcs()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	extractarcs()
	# run identify on wavelength-corrected, extracted arcs to reduce errors further: just read in linelist, improve fit
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Identifying emission lines in the wavelength-calibrated arc spectra to reduce errors further using: identifyarcs'
	print 'Simply read in the linelist, and improve the fit by deleting outliers'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	identifyarcs()
	# run reidentify for remaining standard-star-associated wavelength-corrected, extracted arcs
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Automatically identifying the emission lines in the wavelength-corrected arcs associated with the standard stars using: reidentifyarcs'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	reidentifyarcs()
	# apply lower-error wavelength solution to wavelength-corrected, extracted science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying new lower-error wavelength solution to each science spectrum using: dispcorsci()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	dispcorsci()
	# apply lower-error wavelength solution to wavelength-corrected, extracted standard star images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying new lower-error wavelength solution to each standard star spectrum using: dispcorstd()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	dispcorstd()
	# flux-calibrate each extracted, dispersion-corrected science image
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flux-calibrating the science spectra using: fluxcal()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	fluxcal()
	# combine all individual fully-calibrated spectra (at each GR-ANGLE) into one final spectrum
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Combining each individual fully-calibrated science spectrum into one final spectrum using: finalize()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	finalize()
	

# This function sorts the fits files in the current directory by their type: flat, arc, standard, science.
# This function should ideally be run first so the pipeline knows what it is working with.
def dictionaries():
	# This returns a list of all the fits filenames in the directory.
	images = glob('m*.fits')
	# This says to use the 4 global dictionaries for flat, arc, standard, and science image names indexed by gr-angle.
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This creates a list of the "default" standard stars in the relevant pysalt subdirectory:
	# /usr/local/astro64/iraf/extern/pysalt/data/standards/spectroscopic/
	cwd = os.getcwd()
	os.chdir('/usr/local/astro64/iraf/extern/pysalt/data/standards/spectroscopic/')
	possiblestandards = glob('m*.dat')
	os.chdir(cwd)
	# this gets rid of the 'm' in the beginning and '.dat' at the end so comparisons can be made with 'OBJECT'
	for n,s in enumerate(possiblestandards):
		snew = s.split('.')
		possiblestandards[n] = snew[0][1:].lower()
	# This loops through all the images and sorts them into one of the 4 dictionaries from above based on type and gr-angle.
	for img in images:
		# if necessary, doing header error corrections
		hdulist = pyfits.open(img,mode='update')
		header = hdulist[0].header
		try:
			if header['MASKTYP'] != 'LONGSLIT':
				header.update('MASKTYP','LONGSLIT')
		except:
			header.update('MASKTYP','LONGSLIT')
		hdulist.flush()
		hdulist.close()
		# sorting into dictionaries
		hdr = pyfits.getheader(img,0)
		try:
			angle = str('%.3f'%hdr['GR-ANGLE'])
			# this truncates the gr-angle so it only has 2 decimal digits (string format)
			# angle = "%.2f"%angle
			a = angle.split('.')
			a[1] = a[1][:2]
			angle = '.'.join(a)
			obj = hdr['OBJECT']
		except:
			continue
		if obj=='FLAT':
			temp = flats.get(angle,[])
			temp.append(img)
			flats[angle] = temp
			print img+' sorted as flat with angle: '+angle
		elif obj=='ARC':
			temp = arcs.get(angle,[])
			temp.append(img)
			arcs[angle] = temp
			print img+' sorted as arc with angle: '+angle
		# check if obj is in list of possible SALT standard stars (lowercase comparisons to ignore case)
		elif obj.lower() in possiblestandards:
			temp = standards.get(angle,[])
			temp.append(img)
			standards[angle] = temp
			print img+' sorted as standard star with angle: '+angle
		# if obj is neither flat nor arc nor in possiblestandards, print obj and ask user to classify manually:
		# as standard (short name) or science (long, weird name maybe with 'RU' in it -- obvious)
		else:
			print 'Could not automatically classify '+img+' as a FLAT, ARC, or SALT standard star'
			print 'The object name is: '+obj
			print 'Please enter 0 if this name resembles an abbreviation for a standard star.'
			print 'Please enter 1 if this name is rather long and resembles a science image.'
			while True:
				answer = raw_input("Enter 0 or 1 (see instructions above): ")
				if answer == '0' or answer == '1':
					break
				else:
					print "Invalid input. You must enter either 0 or 1."
			if answer == 0:
				temp = standards.get(angle,[])
				temp.append(img)
				standards[angle] = temp
				print img+' sorted as standard star with angle: '+angle
			elif answer == 1:
				temp = sciences.get(angle,[])
				temp.append(img)
				sciences[angle] = temp
				print img+' sorted as science with angle: '+angle
	# printing new dictionaries
	print 'The original flat image filenames are:'
	print flats
	print 'The original science image filenames are:'
	print sciences
	print 'The original standard star image filenames are:'
	print standards
	print 'The original arc image filenames are:'
	print arcs

# This function runs flatcombine on flats with the same gr-angle (and thus same exptime)
# default input is the dictionary flats created in the function dictionaries()
# argument dictionary should have gr-angle as keywords, and lists of image file names as the values of those keywords
def combineflats(allflats=flats):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This calls dictionaries() if there does not exist the dictionary flats.
	try:
		len(allflats)
	except:
		dictionaries()
		allflats=flats
	# This will run flatcombine on each list of flats in the dictionary flats indexed by gr-angle.
	# Making an assumption about the flats in each list: same gr-angle => same exptime.
	angles = flats.keys()
	for angle in angles:
		flatlist = allflats[angle]
		# If flatlist doesn't have multiple flats in it, can't "combine", so skip to next flatlist
		if len(flatlist)==0 or len(flatlist)==1:
			continue
		# Loop through the headers of each flat in the list to find the average EXPTIME and AIRMASS
		avgEXPTIME = 0
		avgAIRMASS = 0
		for flat in flatlist:
			print flat
			hdr = pyfits.getheader(flat,0)
			avgEXPTIME = avgEXPTIME + hdr['EXPTIME']
			avgAIRMASS = avgAIRMASS + hdr['AIRMASS']
			hdulist = pyfits.open(flat,mode='update')
			hdrold = hdulist[0].header.copy()
			hdulist.close()
		avgEXPTIME = avgEXPTIME/len(flatlist)
		avgAIRMASS = avgAIRMASS/len(flatlist)
		# output name
		name = 'flt'+str(angle)+'cmb.fits'
		# run iraf.noao.imred.ccdred.flatcombine on this list of flats with appropriate parameters.
		# this calls the run_flatcombine(flatstocombine,combflatname,gratingangle) function
		run_flatcombine(flatstocombine=flatlist,combflatname=name)
		# This adds the newly created combined flat name to the combflats dictionary indexed by gr-angle
		combflats[angle] = name
		# This copies over the header of ONE (the last) image in flatlist
		# and changes the exptime and airmass to the averages above.
		# This is an efficient way to get the header info (date-obs, etc.) pysalt might need.
		# A similar thing is done for the science, arc, and standard images.
		hdulist = pyfits.open(name,mode='update')
		hdulist[0].header = hdrold
		hdr = hdulist[0].header
		hdr.update('EXPTIME',avgEXPTIME)
		hdr.update('AIRMASS',avgAIRMASS)
		hdulist.flush()
		hdulist.close()
		print 'output: '+name
	# printing new dictionaries
	print 'The combined flat image filenames are:'
	print combflats

# This function normalizes the combined flats in the combflats dictionary.
def normalizeflats(combinedflats=combflats):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This calls combineflats() if there does not exist the dictionary combflats.
	try:
		len(combinedflats)
	except:
		combineflats()
		combinedflats=combflats
	# This runs response on each combined flat and adds the normalized flat to the normflats dictionary indexed by gr-angle.
	angles = combinedflats.keys()
	for angle in angles:
		combflat = combinedflats[angle]
		hdulist = pyfits.open(combflat,mode='update')
		hdr = hdulist[0].header
		hdr.update('DISPAXIS',1)
		hdulist.flush()
		hdrold = hdulist[0].header.copy()
		hdulist.close()
		# output name
		name = 'flt'+str(angle)+'nrm.fits'
		# run iraf.noao.twodspec.longslit.response on this combined flat.
		# calls the run_response function below
		run_response(combinedflat=combflat,normflatname=name,gratingangle=angle)
		# Adds the name of the newly created normalized flat to the normflats dictionary indexed by gr-angle
		normflats[angle] = name
		# This copies over the old header
		hdulist = pyfits.open(name,mode='update')
		hdulist[0].header = hdrold
		hdulist.flush()
		hdulist.close()
		print 'output: '+name
	print 'The normalized flat image filenames are:'
	print normflats

# This function divides each science image by the normalized flat with the same gr-angle.
def flattensciences(scienceimages=sciences,normalizedflats=normflats):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This prints an error message and quits function if either of input dictionaries don't exist
	try:
		len(normalizedflats),len(scienceimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# This runs imarith to divide a science image by the normalized flat with the same gr-angle.
	angles = scienceimages.keys()
	for angle in angles:
		sciencelist = scienceimages[angle] # this returns the list of all science images with GR-ANGLE=angle
		if len(sciencelist) == 1:
			scienceimages[angle]=sciencelist[0] # sets scienceimages[angle] equal to the filename of the sole science image
		else:
			avgAIRMASS = 0 # will find the average AIRMASS
			suffix = ''
			for img in sciencelist: # create suffix sequence of sci image numbers and avgAIRMASS
				suffixlist = img.split('.')
				suffix = suffix+suffixlist[0][15:]
				hdr = pyfits.getheader(img,0)
				avgAIRMASS = avgAIRMASS + abs(hdr['AIRMASS'])
			avgAIRMASS = avgAIRMASS/len(sciencelist) # dividing sum of airmasses by number of images
			combsciname = 'sci'+str(angle)+'cmb'+suffix+'.fits'
			auxname = 'sci'+str(angle)+'cmb'+suffix+'AUX.fits'
			imutil.imcopy(input=sciencelist[0]+'[0]',output=combsciname+'[append]')
			imutil.imcopy(input=sciencelist[0]+'[SCI]',output=combsciname+'[append]')
			for i,v in enumerate(sciencelist): # need to append [1] to each filename since imcombine will need that
				sciencelist[i] = v+'[1]'
			sciseq = ','.join(sciencelist)
			run_imcombine(imagestocombine=sciseq,combimgname=auxname)
			hdulist = pyfits.open(auxname)
			hdr = hdulist[0].header
			datnew = hdulist[0].data.copy()
			avgEXPTIME = hdr['EXPTIME']
			hdulist.close()
			hdunew = pyfits.open(combsciname,mode='update')
			hdr0 = hdunew[0].header
			hdr1 = hdunew[1].header
			hdunew[1].data = datnew
			hdr0['EXPTIME'] = avgEXPTIME
			hdr1['EXPTIME'] = avgEXPTIME
			hdr0['AIRMASS'] = avgAIRMASS
			hdr1['AIRMASS'] = avgAIRMASS
			hdunew.flush()
			hdunew.close()
			scienceimages[angle] = combsciname # replaces the list of sciences with the combined science
		science = scienceimages[angle]
		suffixlist = science.split('.')
		suffix = suffixlist[0][15:]
		# this stores the old headers (extension 0 and 1)
		hdulist = pyfits.open(science,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		normflat = normalizedflats.get(angle,-1)
		# if no normflat with same gr-angle found, continue onto next image to be flat-fielded
		if normflat == -1:
			print 'no normalized flat found for '+science+' with angle: '+angle
			print 'skipping flat-fielding of '+science
			continue
		# output names
		namesci = 'sci'+str(angle)+'flt'+suffix+'.fits'
		nameaux = 'sci'+str(angle)+'flt'+suffix+'AUX'+'.fits'
		# run iraf.images.imutil.imarith to simply divide science by normflat
		run_imarith(dividend=science,divisor=normflat,quotient=nameaux)
		auxflatsciences[angle] = nameaux
		# This copies the original file to the final file (to preserve file structure), and replaces data
		imutil.imcopy(input=science+'[0]',output=namesci+'[append]')
		imutil.imcopy(input=science+'[SCI]',output=namesci+'[append]')
		hdulist = pyfits.open(nameaux)
		datnew = hdulist[0].data.copy()
		hdulist.close()
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datnew
		hdunew.flush()
		hdunew.close()
		# This adds the newly flat-fielded science image's name to the dictionary flatsciences.
		flatsciences[angle] = namesci
		# updating user about successful flat-fielding
		print 'output: '+namesci
	print 'The flat-fielded science image filenames are:'
	print flatsciences

# This function divides each arc image by the normalized flat with the same gr-angle.
def flattenarcs(arcimages=arcs,normalizedflats=normflats):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This prints an error message and quits the function if one of the inputted dictionary doesn't exist
	try:
		len(normalizedflats),len(arcimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed'
		return
	# This runs imarith to divide an arc image by the normalized flat with the same gr-angle.
	angles = arcimages.keys()
	for angle in angles:
		arclist = arcimages[angle] # this returns the list of all arc images with GR-ANGLE=angle
		if len(arclist) == 1:
			arcimages[angle]=arclist[0] # sets arcimages[angle] equal to the filename of the sole arc image
		elif len(arclist) == 2:
			for img in arclist: # pick arc with LAMPID!='None'; I set this manually in case of SALT error
				hdr = pyfits.getheader(img,0)
				lampid = hdr['LAMPID'].lower()
				if lampid != 'none':
					arcimages[angle]=img
					break
		else:
			avgAIRMASS = 0 # will find the average AIRMASS
			suffix = ''
			for img in arclist: # create suffix sequence of arc image numbers and avgAIRMASS
				suffixlist = img.split('.')
				suffix = suffix+suffixlist[0][15:]
				hdr = pyfits.getheader(img,0)
				avgAIRMASS = avgAIRMASS + abs(hdr['AIRMASS'])
			avgAIRMASS = avgAIRMASS/len(arclist) # dividing sum of airmasses by number of images
			combarcname = 'arc'+str(angle)+'cmb'+suffix+'.fits'
			auxname = 'arc'+str(angle)+'cmb'+suffix+'AUX.fits'
			imutil.imcopy(input=arclist[0]+'[0]',output=combarcname+'[append]')
			imutil.imcopy(input=arclist[0]+'[SCI]',output=combarcname+'[append]')
			for i,v in enumerate(arclist): # need to append [1] to each filename since imcombine will need that
				arclist[i] = v+'[1]'
			arcseq = ','.join(arclist)
			run_imcombine(imagestocombine=arcseq,combimgname=auxname)
			hdulist = pyfits.open(auxname)
			hdr = hdulist[0].header
			datnew = hdulist[0].data.copy()
			avgEXPTIME = hdr['EXPTIME']
			hdulist.close()
			hdunew = pyfits.open(combarcname,mode='update')
			hdr0 = hdunew[0].header
			hdr1 = hdunew[1].header
			hdunew[1].data = datnew
			hdr0['EXPTIME'] = avgEXPTIME
			hdr1['EXPTIME'] = avgEXPTIME
			hdr0['AIRMASS'] = avgAIRMASS
			hdr1['AIRMASS'] = avgAIRMASS
			hdunew.flush()
			hdunew.close()
			arcimages[angle] = combarcname # replaces the list of arcs with the combined arc
		arc = arcimages[angle]
		# this copies (and stores) the old (original) arc image's header
		suffixlist = arc.split('.')
		suffix = suffixlist[0][15:]
		hdulist = pyfits.open(arc,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		normflat = normalizedflats.get(angle,-1)
		# if no normflat with same gr-angle found, continue onto next image to be flat-fielded
		if normflat == -1:
			print 'no normalized flat found for '+arc+' with angle: '+angle
			continue
		# output name
		namearc = 'arc'+str(angle)+'flt'+suffix+'.fits'
		nameaux = 'arc'+str(angle)+'flt'+suffix+'AUX'+'.fits'
		# run iraf.images.imutil.imarith to simply divide science by normflat
		run_imarith(dividend=arc,divisor=normflat,quotient=nameaux)
		auxflatarcs[angle] = nameaux
		# This copies the original file to the final file (to preserve file structure), and replaces data
		imutil.imcopy(input=arc+'[0]',output=namearc+'[append]')
		imutil.imcopy(input=arc+'[SCI]',output=namearc+'[append]')
		hdulist = pyfits.open(nameaux)
		datnew = hdulist[0].data.copy()
		hdulist.close()
		hdunew = pyfits.open(namearc,mode='update')
		hdunew[1].data = datnew
		hdunew.flush()
		hdunew.close()
		# This adds the newly flat-fielded arc image's name to the dictionary flatarcs.
		flatarcs[angle] = namearc
		print 'output: '+namearc
	print 'The flat-fielded arc image filenames are:'
	print flatarcs

# This function divides each standard star image by the normalized flat with the same gr-angle.
def flattenstandards(standardimages=standards,normalizedflats=normflats):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This prints an error message and quits the function if an inputted dictionary doesn't exist
	try:
		len(normalizedflats),len(standardimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# This runs imarith to divide a standard star image by the normalized flat with the same gr-angle.
	angles = standardimages.keys()
	for angle in angles:
		standardlist = standardimages[angle] # this returns the list of all standard images with GR-ANGLE=angle
		if len(standardlist) == 1:
			standardimages[angle]=standardlist[0] # sets standardimages[angle] equal to the filename of the sole standard image
		else:
			avgAIRMASS = 0 # will find the average AIRMASS
			suffix = ''
			for img in standardlist: # create suffix sequence of standard image numbers and avgAIRMASS
				suffixlist = img.split('.')
				suffix = suffix+suffixlist[0][15:]
				hdr = pyfits.getheader(img,0)
				avgAIRMASS = avgAIRMASS + abs(hdr['AIRMASS'])
			avgAIRMASS = avgAIRMASS/len(standardlist) # dividing sum of airmasses by number of images
			combstdname = 'std'+str(angle)+'cmb'+suffix+'.fits'
			auxname = 'std'+str(angle)+'cmb'+suffix+'AUX.fits'
			imutil.imcopy(input=standardlist[0]+'[0]',output=combstdname+'[append]')
			imutil.imcopy(input=standardlist[0]+'[SCI]',output=combstdname+'[append]')
			for i,v in enumerate(standardlist): # need to append [1] to each filename since imcombine will need that
				standardlist[i] = v+'[1]'
			standardseq = ','.join(standardlist)
			run_imcombine(imagestocombine=standardseq,combimgname=auxname)
			hdulist = pyfits.open(auxname)
			hdr = hdulist[0].header
			datnew = hdulist[0].data.copy()
			avgEXPTIME = hdr['EXPTIME']
			hdulist.close()
			hdunew = pyfits.open(combstdname,mode='update')
			hdr0 = hdunew[0].header
			hdr1 = hdunew[1].header
			hdunew[1].data = datnew
			hdr0['EXPTIME'] = avgEXPTIME
			hdr1['EXPTIME'] = avgEXPTIME
			hdr0['AIRMASS'] = avgAIRMASS
			hdr1['AIRMASS'] = avgAIRMASS
			hdunew.flush()
			hdunew.close()
			standardimages[angle] = combstdname # replaces the list of standards with the combined standard
		standard = standardimages[angle]
		suffixlist = standard.split('.')
		suffix = suffixlist[0][15:]
		hdulist = pyfits.open(standard,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# choose normflat with same rounded gr-angle (string format)
		normflat = normalizedflats.get(angle,-1)
		# if no normflat with same gr-angle found, continue onto next image to be flat-fielded
		if normflat == -1:
			print 'no normalized flat found for '+standard+' with angle: '+angle
			continue
		# output name
		namestd = 'std'+str(angle)+'flt'+suffix+'.fits'
		nameaux = 'std'+str(angle)+'flt'+suffix+'AUX'+'.fits'
		# run iraf.images.imutil.imarith to simply divide standard by normflat
		run_imarith(dividend=standard,divisor=normflat,quotient=nameaux)
		auxflatstandards[angle] = nameaux
		# This copies the original file to the final file (to preserve file structure), and replaces data
		imutil.imcopy(input=standard+'[0]',output=namestd+'[append]')
		imutil.imcopy(input=standard+'[SCI]',output=namestd+'[append]')
		hdulist = pyfits.open(nameaux)
		datnew = hdulist[0].data.copy()
		hdulist.close()
		hdunew = pyfits.open(namestd,mode='update')
		hdunew[1].data = datnew
		hdunew.flush()
		hdunew.close()
		# This adds the newly flat-fielded standard star image's name to the dictionary flatstandards.
		flatstandards[angle] = namestd
		print 'output: '+namestd
	print 'The flat-fielded standard star image filenames are:'
	print flatstandards


# This function runs saltspec.specidentify to produce wavelength solution files.
# This doesn't yet create 2 separate dictionaries -- 1 for the arcs associated with science images, and 1 for
# the arcs associated with standard star images -- partly because SALT doesn't do that (yet).
def specidentifyarcs(arcimages=flatarcs):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This calls flattenarcs() if the dictionary arcimages doesn't exist
	# Should probably just print an error instead
	try:
		len(arcimages)
	except:
		flattenarcs()
		arcimages = flatarcs
	# this prepares for and runs saltspec.specidentify (from pysalt)
	angles = arcimages.keys()
	for angle in angles:
		flatarc = arcimages[angle]
		suffixlist = flatarc.split('.')
		suffix = suffixlist[1][5:]
		# this gets the arc lamp type so the correct linelist can be used by saltspec.specidentify
		hdulist = pyfits.open(flatarc,mode='update')
		hdr = hdulist[0].header
		lamp = hdr['LAMPID']
		hdulist.close()
		# this sets the arguments for run_specidentify function below
		linelistpath = '/usr/local/astro64/iraf/extern/pysalt/data/linelists/'
		# picking the linelist file based on header LAMPID reference variable lamp
		if lamp == 'Th Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'ThAr.txt'
		elif lamp == 'Xe':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Xe.txt'
		elif lamp == 'Ne': 
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Ne.txt'
		elif lamp == 'Cu Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'CuAr.txt'
		elif lamp == 'Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Ar.txt'
		elif lamp == 'Hg Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'HgAr.txt'
		else:
			print 'Could not find the proper linelist for '+lamp+' lamp.'
			continue
		# this sets the name for the output wavelength solution file which will be produced by saltspec.specidentify
		idfilesci = 'arcid'+str(angle)+'sci'+suffix+'.db'
		# idfilestd = 'arcid'+str(angle)+'_std.db'
		# this calls the run_specidentify function found further below
		run_specidentify(arcimage=flatarc,lamplines=linelistpath,idfile=idfilesci)
		# this adds the newly created wavelength solution file to the wavesols dictionary
		wavesols[angle] = idfilesci
	print 'The specidentify wavelength solution filenames are:'
	print wavesols

# This function calls saltspec.specrectify to apply the wavelength solution to arc images.
def wavecalarc(arcimages=flatarcs,sols=wavesols):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(arcimages)
	except:
		print 'The inputted dictionary of flat-fielded arc images does not exist.'
		return
	try:
		len(sols)
	except:
		print 'The inputted dictionary of wavelength solutions does not exist.'
		return
	# This runs saltspec.specrectify on each flat-fielded arc image.
	angles = arcimages.keys()
	for angle in angles:
		# this stores the filenames of an arc image and its associated wavelength solution
		flatarc = arcimages[angle]
		suffixlist = flatarc.split('.')
		suffix = suffixlist[1][5:]
		sol = sols[angle]
		# this stores the old headers of flatarc to copy over to the new image
		hdulist = pyfits.open(flatarc,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 2-D wavelength-corrected arc image
		namearc = 'arc'+str(angle)+'wav'+suffix+'.fits'
		# this calls the run_specrectify function found further below
		run_specrectify(input=flatarc,output=namearc,idfile=sol)
		# this adds the filename of the new 2-D wavelength-corrected image to proper dictionary
		wavearcs[angle] = namearc
		# this notifies the user that an output was created by printing its name
		print 'output: '+namearc
	print 'The wavelength-calibrated arc image filenames are:'
	print wavearcs
		
		
# This function calls saltspec.specrectify to apply the wavelength solution to science images.
def wavecalsci(scienceimages=flatsciences,sols=wavesols):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary of flat-fielded science images does not exist.'
		return
	try:
		len(sols)
	except:
		print 'The inputted dictionary of wavelength solutions does not exist.'
		return
	# This runs saltspec.specrectify on each flat-fielded science image.
	angles = scienceimages.keys()
	for angle in angles:
		# this stores the filenames of a science image and its associated wavelength solution
		flatscience = scienceimages[angle]
		suffixlist = flatscience.split('.')
		suffix = suffixlist[1][5:]
		sol = sols[angle]
		# this stores the old header of flatscience to copy over to the new image
		hdulist = pyfits.open(flatscience,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 2-D wavelength-corrected science image
		namescience = 'sci'+str(angle)+'wav'+suffix+'.fits'
		# this calls the run_specrectify function found further below
		run_specrectify(input=flatscience,output=namescience,idfile=sol)
		# this adds the filename of the new 2-D wavelength-corrected image to proper dictionary
		wavesciences[angle] = namescience
		# this notifies the user that an output was created by printing its name
		print 'output: '+namescience	
	print 'The wavelength-calibrated science image filenames are:'
	print wavesciences

# This function calls saltspec.specrectify to apply the wavelength solution to standard star images.
def wavecalstd(standardimages=flatstandards,sols=wavesols):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(standardimages)
	except:
		print 'The inputted dictionary of flat-fielded standard star images does not exist.'
		return
	try:
		len(sols)
	except:
		print 'The inputted dictionary of wavelength solutions does not exist.'
		return
	# This runs saltspec.specrectify on each flat-fielded standard star image.
	angles = standardimages.keys()
	for angle in angles:
		# this stores the filenames of a standard star image and its associated wavelength solution
		flatstandard = standardimages[angle]
		suffixlist = flatstandard.split('.')
		suffix = suffixlist[1][5:]
		sol = sols[angle]
		# this stores the old header of flatstandard to copy over to the new image
		hdulist = pyfits.open(flatstandard,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output name for the 2-D wavelength-corrected standard star image
		namestandard = 'std'+str(angle)+'wav'+suffix+'.fits'
		# this calls the run_specrectify function found further below
		run_specrectify(input=flatstandard,output=namestandard,idfile=sol)
		# this adds the name of the new 2-D wavelength-corrected image to proper dictionary
		wavestandards[angle] = namestandard
		# this notifies the user that an output was created by printing its name
		print 'output: '+namestandard
	print 'The wavelength-calibrated standard star image filenames are:'
	print wavestandards

# This function subtracts the background (sky lines essentially) from the two-dimensional science images.
def subtractbackground(scienceimages=wavesciences):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks if the inputted dictionary of images exists; if not, it quits the function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary of wavelength-corrected 2D science images does not exist.'
		return
	# this sets the DISPAXIS=1 in the header and calls the run_background function further below
	angles = scienceimages.keys()
	for angle in angles:
		wavescience = scienceimages[angle]
		suffixlist = wavescience.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(wavescience,mode='update')
		hdr = hdulist[1].header
		hdr.update('DISPAXIS',1)
		hdulist.flush()
		hdulist.close()
		# output filenames
		name = 'sci'+str(angle)+'bkg'+suffix+'.fits'
		auxname = 'sci'+str(angle)+'bkg'+suffix+'AUX'+'.fits'
		# this calls the run_background function (NOW COMPLETELY AUTOMATED)
		run_background(twodimage=wavescience,newimage=auxname)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=wavescience+'[0]',output=name+'[append]')
		imutil.imcopy(input=wavescience+'[SCI]',output=name+'[append]')
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# Add background-subtracted science image filename to dictionary backgroundsciences indexed by gr-angle.
		backgroundsciences[angle] = name
		# notify user of successful background-subtraction
		print 'output: '+name
	print 'The background-subtracted science image filenames are:'
	print backgroundsciences
		

# This function runs apall on the wavelength-corrected, background-subtracted science images to extract the SN spectra.
def extractsciences(scienceimages=backgroundsciences):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks if the inputted dictionary of images exists; if not, it quits the function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary of background-subtracted 2D images does not exist.'
		return	
	# This adds the DISPAXIS keyword to each image's header with value 1 and then runs apall.
	angles = scienceimages.keys()
	for angle in angles:
		bkgscience = scienceimages[angle]
		suffixlist = bkgscience.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(bkgscience,mode='update')
		hdr = hdulist[1].header
		hdr.update('DISPAXIS',1)
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this displays the image in ds9 for the user's convenience
		print "You will now use the PyRAF task APALL to extract the spectrum in the image: "+bkgscience
		print "For your convenience, the image has also been opened in ds9 so you can find:"
		print '1. The central row of the supernova spectrum. (Consult your acquisition images if needed.)'
		print '2. Whether the supernova spectrum is bright or faint (barely distinguishable from the background).'
		print 'The image is already background-subtracted so you do NOT need to subtract the background ('b' key).'
		# opens current image in ds9
		d = ds9.ds9(start=True)
		hdulist = pyfits.open(bkgscience)
		d.set_pyfits(hdulist)
		d.set('zoom to fit')
		d.set('zscale')
		# this calls the run_apall function found further below
		# this asks the user if the spectral line looks very faint (look in ds9)
		# if the user says the line is very faint, then apall is run with a different set of parameters
		while True:
			answer = raw_input("Please enter 0 if the spectral line looks very faint, or 1 if it is clearly visible: ")
			if answer == '0' or answer == '1':
				break
			else:
				print "Invalid input. You must enter either 0 or 1."
		# output filenames
		name = 'sci'+str(angle)+'ext'+suffix+'.fits'
		auxname = 'sci'+str(angle)+'ext'+suffix+'AUX'+'.fits'
		run_apall(twodimage=bkgscience,spectrum=auxname,faint=answer)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=bkgscience+'[0]',output=name+'[append]')
		imutil.imcopy(input=bkgscience+'[SCI]',output=name+'[append]')
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# Add extracted science image filename to dictionary extractedsciences indexed by gr-angle.
		extractedsciences[angle] = name
		# notify user of successful extraction
		print 'output: '+name
	print 'The extracted science image filenames are:'
	print extractedsciences
	

# This function runs apall on the wavelength-corrected standard star images to extract the standard star spectra.
def extractstandards(standardimages=wavestandards):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks if the inputted dictionary of images exists; if not, it quits the function
	try:
		len(standardimages)
	except:
		print 'The inputted dictionary of wavelength-corrected 2D images does not exist.'
		return
	# This adds the DISPAXIS keyword to each image's header with value 1 and then runs apall.
	angles = standardimages.keys()
	for angle in angles:
		wavestandard = standardimages[angle]
		suffixlist = wavestandard.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(wavestandard,mode='update')
		hdr = hdulist[1].header
		hdr.update('DISPAXIS',1)
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# setting answer = 1 (bright spectrum) for run_apall since standards will be extracted non-interactively
		answer = 1
		# output filenames
		name = 'std'+str(angle)+'ext'+suffix+'.fits'
		auxname = 'std'+str(angle)+'ext'+suffix+'AUX'+'.fits'
		run_apall(twodimage=wavestandard,spectrum=auxname,faint=answer)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=wavestandard+'[0]',output=name+'[append]')
		imutil.imcopy(input=wavestandard+'[SCI]',output=name+'[append]')
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# Add extracted standard star image filename to dictionary extractedstandards indexed by gr-angle.
		extractedstandards[angle] = name
		# notify user of successful extraction
		print 'output: '+name
	print 'The extracted standard star image filenames are:'
	print extractedstandards


# This function runs apsum on the wavelength-corrected arc images to extract the arc spectra.
# will be run twice for each arc -- once for science image, and once for standard star image
def extractarcs(arcimages=wavearcs):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks if the inputted dictionary of images exists; if not, quits the function
	try:
		len(arcimages)
	except:
		print 'The inputted dictionary of images does not exist.'
		return
	# This runs apsum twice for each arc image (for science-arcs and standard-arcs)
	angles = arcimages.keys()
	for angle in angles:
		wavearc = arcimages[angle]
		suffixlist = wavearc.split('.')
		suffix = suffixlist[1][5:]
		# this sets 'DISPAXIS'=1 and stores the old header so it can be copied over to the new images
		hdulist = pyfits.open(wavearc,mode='update')
		hdr = hdulist[1].header
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output and reference image names for the standard star-reference arc
		namestd = 'arc'+str(angle)+'ext'+suffix+'_std.fits'
		auxnamestd = 'arc'+str(angle)+'ext'+suffix+'_stdAUX.fits'
		refstd = wavestandards[angle]+'[1]' # use wavestandards as standard star reference
		# this calls the run_apsum function found further below for the STD-extracted arc
		run_apsum(twodimage=wavearc,refimage=refstd,spectrum=auxnamestd)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=wavearc+'[0]',output=namestd+'[append]')
		imutil.imcopy(input=wavearc+'[SCI]',output=namestd+'[append]')
		hdulist = pyfits.open(auxnamestd)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namestd,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# Add STD-extracted arc image to dictionary extractedarcs_std indexed by gr-angle.
		extractedarcs_std[angle] = namestd
		# notify user of successful extraction
		print 'output: '+namestd
		
		# this sets the output and reference image filename for the science-reference arc
		namesci = 'arc'+str(angle)+'ext'+suffix+'_sci.fits'
		auxnamesci = 'arc'+str(angle)+'ext'+suffix+'_sciAUX.fits'
		refsci = backgroundsciences[angle]+'[1]' # use backgroundsciences for science reference
		# this calls the run_apsum function found further below for the SCI-extracted arc
		run_apsum(twodimage=wavearc,refimage=refsci,spectrum=auxnamesci)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=wavearc+'[0]',output=namesci+'[append]')
		imutil.imcopy(input=wavearc+'[SCI]',output=namesci+'[append]')
		hdulist = pyfits.open(auxnamesci)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# Add extracted arc image to dictionary extractedarcs_sci indexed by gr-angle.
		extractedarcs_sci[angle] = namesci
		# notify user of successful extraction
		print 'output: '+namesci
	print 'The extracted arc (for science images) image filenames are:'
	print extractedarcs_sci
	print 'The extracted arc (for standard star images) image filenames are:'
	print extractedarcs_std


# This function runs identify on the wavelength-corrected, SCI-extracted arc images to reduce errors further.
def identifyarcs(arcsci=extractedarcs_sci,scienceimages=extractedsciences):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(arcsci),len(scienceimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# this prepares for and runs identify on each extracted arc (associated with a science image)
	angles = arcsci.keys()
	for angle in angles:
		# this stores the filenames of the extracted arc and its associated extracted science
		extarc = arcsci[angle]
		extsci = scienceimages[angle]
		# this gets the arc lamp type so the correct linelist can be used by identify
		hdulist = pyfits.open(extarc,mode='update')
		hdr = hdulist[0].header
		lamp = hdr['LAMPID']
		hdulist.close()
		# this sets the arguments for run_identify function below
		linelistpath = '/usr/local/astro64/iraf/extern/pysalt/data/linelists/'
		# picking the linelist file based on header LAMPID reference variable lamp
		if lamp == 'Th Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'ThAr.txt'
		elif lamp == 'Xe':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Xe.txt'
		elif lamp == 'Ne': # might instead need to use NeAr.salt
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Ne.txt'
		elif lamp == 'Cu Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'CuAr.txt'
		elif lamp == 'Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Ar.txt'
		elif lamp == 'Hg Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'HgAr.txt'
		else:
			print 'Could not find the proper linelist for '+lamp+' lamp.'
			continue
		# this calls the run_identify function found further below on extarc with data extension [1] ([SCI])
		extarc1 = extarc+'[1]'
		run_identify(arcsci=extarc1,lamplines=linelistpath)
		# this adds the 'REFSPEC1' keyword to the associated extracted science image for dispcor later
		hdulist = pyfits.open(extsci,mode='update')
		hdr = hdulist[0].header
		hdr.update('REFSPEC1',extarc)
		hdulist.flush()
		hdulist.close()
		# notifies user of successful identification 
		print 'emission lines in the file '+extarc+' have been identified'
		print 'and it is now linked to its science image: '+extsci
		

# This function runs reidentify on the remaining STD-extracted arcs.
def reidentifyarcs(arcstd=extractedarcs_std,arcsci=extractedarcs_sci,standardimages=extractedstandards):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(arcstd),len(arcsci),len(standardimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# this prepares for and runs reidentify on each extracted arc (associated with a standard star image)
	angles = arcstd.keys()
	for angle in angles:
		# this stores the filenames of the extracted arc, reference sci-extracted arc, and associated extracted std star
		extarcstd = arcstd[angle] # current arc to be identified
		extarcsci = arcsci[angle] # reference arc already identified (with identify)
		extstd = standardimages[angle] # associated standard star image
		# this gets the arc lamp type so the correct linelist can be used by identify
		hdulist = pyfits.open(extarcstd,mode='update')
		hdr = hdulist[0].header
		lamp = hdr['LAMPID']
		hdulist.close()
		# this sets the arguments for run_reidentify function below
		linelistpath = '/usr/local/astro64/iraf/extern/pysalt/data/linelists/'
		# picking the linelist file based on header LAMPID reference variable lamp
		if lamp == 'Th Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'ThAr.txt'
		elif lamp == 'Xe':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Xe.txt'
		elif lamp == 'Ne': # might instead need to use NeAr.salt
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Ne.txt'
		elif lamp == 'Cu Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'CuAr.txt'
		elif lamp == 'Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'Ar.txt'
		elif lamp == 'Hg Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = linelistpath+'HgAr.txt'
		else:
			print 'Could not find the proper linelist for '+lamp+' lamp.'
			continue
		# this calls the run_reidentify function found further below on extarcstd with data extension [1] ([SCI])
		extarcstd1 = extarcstd+'[1]'
		extarcsci1 = extarcsci+'[1]'
		run_reidentify(input=extarcstd1,ref=extarcsci1,lamplines=linelistpath)
		# this adds the 'REFSPEC1' keyword to the associated extracted standard star image for dispcor later
		hdulist = pyfits.open(extstd,mode='update')
		hdr = hdulist[0].header
		hdr.update('REFSPEC1',extarcstd)
		hdulist.flush()
		hdulist.close()
		# notifies user of successful identification 
		print 'emission lines in the file '+extarcstd+' have been identified'
		print 'and it is now linked to its standard star image: '+extstd
		

# This function runs dispcor on the extracted science images to apply the lower-error wavelength solution.
# REFSPEC1 keyword was already added to images by the functions, identifyarcs and reidentifyarcs
def dispcorsci(scienceimages=extractedsciences):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# This runs onedspec.dispcor on each extracted science image.
	angles = scienceimages.keys()
	for angle in angles:
		# this stores the filenames of a science image
		extscience = scienceimages[angle]
		suffixlist = extscience.split('.')
		suffix = suffixlist[1][5:]
		# this stores the old header of extscience to copy over to the new image
		hdulist = pyfits.open(extscience,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 1-D wavelength-corrected science image
		namescience = 'sci'+str(angle)+'dsp'+suffix+'.fits'
		auxnamescience = 'sci'+str(angle)+'dsp'+suffix+'AUX.fits'
		# this calls the run_dispcor function found further below on extscience with data extension [1] ([SCI])
		extscience1 = extscience+'[1]'
		run_dispcor(original=extscience1,corrected=auxnamescience)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=extscience+'[0]',output=namescience+'[append]')
		imutil.imcopy(input=extscience+'[SCI]',output=namescience+'[append]')
		hdulist = pyfits.open(auxnamescience)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namescience,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# this adds the filename of the new 1-D wavelength-corrected image to proper dictionary
		dispsciences[angle] = namescience
		# this notifies the user that an output was created by printing its name
		print 'output: '+namescience
	print 'The dispersion-corrected science image filenames are:'
	print dispsciences


# This function runs dispcor on the extracted standard star images to apply the lower-error wavelength solution.
# REFSPEC1 keyword was already added to images by the functions, identifyarcs and reidentifyarcs
def dispcorstd(standardimages=extractedstandards):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(standardimages)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# This runs onedspec.dispcor on each extracted standard star image.
	angles = standardimages.keys()
	for angle in angles:
		# this stores the filenames of a science image
		extstandard = standardimages[angle]
		suffixlist = extstandard.split('.')
		suffix = suffixlist[1][5:]
		# this stores the old header of extstandard to copy over to the new image
		hdulist = pyfits.open(extstandard,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 2-D wavelength-corrected standard star image
		namestandard = 'std'+str(angle)+'dsp'+suffix+'.fits'
		auxnamestandard = 'std'+str(angle)+'dsp'+suffix+'AUX.fits'
		# this calls the run_dispcor function found further below on extstandard with data extension [1] ([SCI])
		extstandard1 = extstandard+'[1]'
		run_dispcor(original=extstandard1,corrected=auxnamestandard)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=extstandard+'[0]',output=namestandard+'[append]')
		imutil.imcopy(input=extstandard+'[SCI]',output=namestandard+'[append]')
		hdulist = pyfits.open(auxnamestandard)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namestandard,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# this adds the filename of the new 2-D wavelength-corrected image to proper dictionary
		dispstandards[angle] = namestandard
		# this notifies the user that an output was created by printing its name
		print 'output: '+namestandard
	print 'The dispersion-corrected standard star image filenames are:'
	print dispstandards


# This function flux-calibrates the science images using the standard star images with standard, sensfunc, and calibrate.
def fluxcal(dspsci=dispsciences,dspstd=dispstandards):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global auxextractedsciences
	global auxextractedstandards
	global auxextractedarcs_std
	global auxextractedarcs_sci
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global auxdispsciences
	global auxdispstandards
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global auxfluxsciences
	global fluxsciences
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(dspsci),len(dspstd)
	except:
		print 'One of the inputted dictionaries does not exist: cannot proceed.'
		return
	# This runs standard, sensfunc for each standard star spectrum; and then calibrate on the associated science image
	angles = dspstd.keys()
	for angle in angles:
		# this stores the filename for a science image, and its header
		dspscience = dspsci[angle]
		suffixlist = dspscience.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(dspscience,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this stores the filename of the corresponding standard star image and the class of standard star it is
		dspstandard = dspstd[angle]
		hdulist = pyfits.open(dspstandard,mode='update')
		hdr = hdulist[1].header # identify and reidentify messes around with header 0, so added keys to header 1
		starclass = hdr['OBJECT']
		airmass = hdr['AIRMASS']
		exptime = hdr['EXPTIME']
		hdulist.close()
		# this sets the output names for the flux-calibrated science image, the std file, and the sens file
		namesci = 'sci'+str(angle)+'flx'+suffix+'.fits'
		auxnamesci = 'sci'+str(angle)+'flx'+suffix+'AUX.fits'
		stdfile = 'std'+str(angle)+'flx'+suffix
		sensfile = 'sens'+str(angle)+'flx'+suffix
		# this prints the star class and directions for how to use it during the standard task
		print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
		print 'You shall not pass beyond the next task unless you use the following data:'
		print 'This standard star is of type: '+starclass
		print 'When the next task asks you to input the star_name, first input ?,'
		print 'and then input the star_name that the task has chosen to represent '+starclass
		print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
		# this runs onedspec.standard to produce a unique std file (on dspstandard with data extension [1] ([SCI]))
		dspstandard1 = dspstandard+'[1]'
		run_standard(standardimage=dspstandard1,outputname=stdfile,exp=exptime,air=airmass)
		# this adds the std file to the stdfiles dictionary
		stdfiles[angle] = stdfile
		# this runs onedspec.sensfunc to produce a unique sens file
		run_sensfunc(stddata=stdfile,sensname=sensfile)
		# this adds the sens file to the sensfiles dictionary
		sensfiles[angle] = sensfile
		# this runs onedspec.calibrate to calibrate science image using the just-produced sens file (on dspscience[SCI])
		dspscience1 = dspscience+'[SCI]'
		run_calibrate(scienceimage=dspscience,fluximage=auxnamesci,sensfilename=sensfile)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		imutil.imcopy(input=dspscience+'[0]',output=namesci+'[append]')
		imutil.imcopy(input=dspscience+'[SCI]',output=namesci+'[append]')
		hdulist = pyfits.open(auxnamesci)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdunew.close()
		# this adds the flux-calibrated science image filename to dictionary fluxsciences indexed by gr-angle
		fluxsciences[angle] = namesci
		# updates user on successful flux calibration
		print 'output: '+namesci
	print 'The flux-calibrated science image filenames are:'
	print fluxsciences
	
	
# This function runs scombine to finalize the supernovae spectra into one spectrum.
def finalize(fluxspectra=fluxsciences):
	# use the global dictionaries for all operations
	global flats
	global arcs
	global standards
	global sciences   
	global combflats
	global normflats
	global auxflatsciences
	global auxflatarcs
	global auxflatstandards
	global flatsciences
	global flatarcs
	global flatstandards
	global wavesols
	global wavearcs
	global wavesciences
	global wavestandards
	global auxbackgroundsciences
	global backgroundsciences
	global extractedsciences
	global extractedstandards
	global extractedarcs_std
	global extractedarcs_sci
	global dispsciences
	global dispstandards
	global stdfiles
	global sensfiles
	global fluxsciences
	# This prints an error message and quits function if input dictionary doesn't exist.
	try:
		len(fluxspectra)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# output name
	namesci = 'scifinal.fits'
	# create a new list from each of the values in the dictionary fluxsciences to feed into onedspec.scombine
	spectratocombine = []
	angles = fluxspectra.keys()
	for angle in angles:
		spectratocombine.append(fluxspectra[angle]+'[1]')
	# turn list into sequence of comma-separated strings to feed into scombine
	try:
		if len(spectratocombine)>1:
			scombineseq = ','.join(spectratocombine)
			# run scombine interactively to remove cosmic ray leftover effects
			run_scombine(inputlistfile=scombineseq,finalspectrumname=namesci)
			print 'The finalized and combined science image filename is:'
			print namesci
		else:
			print 'There is only one flux-calibrated spectrum: cannot run scombine.'
	except:
		print 'Error running scombine: did not proceed.'


# This function sets the parameters for the pyraf task ccdred.flatcombine and then runs it.
def run_flatcombine(flatstocombine,combflatname):
	# This clears the current parameter values for flatcombine, and sets general parameters again.
	ccdred.flatcombine.unlearn()
	ccdred.flatcombine.combine='average'
	ccdred.flatcombine.reject='avsigclip'
	ccdred.flatcombine.ccdtype=''
	ccdred.flatcombine.scale='mode'
	ccdred.flatcombine.process='no'
	ccdred.flatcombine.subsets='no'
	
	# appends '[1]' to end of each filename since that's the extension for the data to be combined
	for i,v in enumerate(flatstocombine):
		flatstocombine[i] = v+'[1]'
	
	# turn python list into an actual listfile and feed it into ccdred.flatcombine
	f = open('flatcombinelist','w')
	f.write('\n'.join(flatstocombine))
	f.close()
	
	ccdred.flatcombine(input='@flatcombinelist',output=combflatname)

# This function sets the parameters for the pyraf task longslit.response and then runs it.
def run_response(combinedflat,normflatname,gratingangle):
	# this clears current parameter values for longslit.response and sets general parameters.
	longslit.response.unlearn()
	longslit.response.threshold=5.0
	longslit.response.function='spline3'
	longslit.response.order=19
	longslit.response.low_rej=3.0
	longslit.response.high_rej=3.0
	longslit.response.niterate=3
	longslit.response.grow=2
	longslit.response.interactive='no' # make non-interactive for now
	
	# This runs longslit.response (interactively for now)
	longslit.response(calibration=combinedflat,normalization=combinedflat,response=normflatname)
	

# This function sets the parameters for the pyraf task immatch.imcombine and then runs it.
def run_imcombine(imagestocombine,combimgname):
	# This clears the current parameter values for imcombine, and sets general parameters again.
	immatch.imcombine.unlearn()
	immatch.imcombine.project='No'
	immatch.imcombine.combine='average'
	immatch.imcombine.reject='crreject'
	immatch.imcombine.mclip='yes'
	immatch.imcombine.nkeep=1
	immatch.imcombine.rdnoise=0.0
	immatch.imcombine.gain=1.0
	immatch.imcombine.snoise=0.0
	immatch.imcombine.hsigma=3.0
	immatch.imcombine.blank=0.0
	immatch.imcombine.expname='EXPTIME' # average EXPTIME added to output image header
	immatch.imcombine.imcmb='$I,AIRMASS' # copy airmass to one of the IMCMBnnn keywords in output header
	
	# appends '[1]' to end of each filename since that's the extension for the data to be combined
	# for i,v in enumerate(imagestocombine):
		# imagestocombine[i] = v+'[1]'
	
	# turn python list into an actual listfile and feed it into imutil.imcombine
	# f = open('imcombinelist','w')
	# f.write('\n'.join(imagestocombine))
	# f.close()
	
	immatch.imcombine(input=imagestocombine,output=combimgname)
	

# This function sets the parameters for the pyraf task imutil.imarith and then runs it.
def run_imarith(dividend,divisor,quotient):
	# This resets the parameters of imutil.imarith.
	imutil.imarith.unlearn()
	imutil.imarith.divzero=0.0
	# appends [1] to end of dividend image since that's where data lives
	dividend = dividend+'[1]'

	# This runs imutil.imarith
	imutil.imarith(operand1=dividend,op='/',operand2=divisor,result=quotient)	
	
	
# This function sets the parameters for the pysalt task saltspec.specidentify and then runs it.
def run_specidentify(arcimage,lamplines,idfile):
	# this resets the parameters of saltspec.specidentify
	saltspec.specidentify.unlearn()
	saltspec.specidentify.guesstype='rss'
	saltspec.specidentify.automethod='Matchlines'
	saltspec.specidentify.function='legendre'
	saltspec.specidentify.order=3
	saltspec.specidentify.rstep=100
	saltspec.specidentify.rstart='middlerow'
	saltspec.specidentify.mdiff=5
	saltspec.specidentify.inter='yes' # interactive
	saltspec.specidentify.startext=1 # important because SALT data is in ext 1 (name: SCI)
	saltspec.specidentify.clobber='yes'
	saltspec.specidentify.logfile='pysalt.log'
	saltspec.specidentify.verbose='yes'
	
	# this runs saltspec.specidentify
	saltspec.specidentify(images=arcimage,linelist=lamplines,outfile=idfile)


# This function sets the parameters for saltspec.specrectify and runs the task.	
# It outputs 2-D wavelength-corrected images (science, arc, and standard star).
def run_specrectify(input,output,idfile):
	# this resets the parameters of saltspec.specrectify
	saltspec.specrectify.unlearn()
	saltspec.specrectify.outpref=''
	saltspec.specrectify.caltype='line'
	saltspec.specrectify.function='polynomial'
	saltspec.specrectify.order=5
	saltspec.specrectify.inttype='interp'
	saltspec.specrectify.clobber='yes'
	saltspec.specrectify.verbose='yes'
	
	# this runs saltspec.specrectify
	saltspec.specrectify(images=input,outimages=output,solfile=idfile)

	
# This function sets the parameters for the pyraf task twodspec.longslit.background and then runs it.
def run_background(twodimage,newimage):
	# this resets the parameters of longslit.background
	longslit.background.unlearn()
	longslit.background.axis='2'
	longslit.background.interactive='no'
	longslit.background.naverage='-100'
	longslit.background.function='legendre'
	longslit.background.order=2
	longslit.background.low_reject=1.0
	longslit.background.high_reject=1.0
	longslit.background.niterate=10
	longslit.background.grow=1.0
	# since there are 2 extensions (0 and 1), need to specify data operation extension: 1
	twodimage = twodimage+'[1]'
	# This runs longslit.background
	longslit.background(input=twodimage,output=newimage)
	
# This function sets the parameters for the pyraf task apextract.apall and then runs it.
def run_apall(twodimage,spectrum,faint):
	# this resets the parameters of apextract.apall.
	# with as many appropriate parameters set as possible:
	apextract.apall.unlearn()
	apextract.apall.format='multispec'
	apextract.apall.review='no'
	apextract.apall.nsum=50 # set nsum equal to user-inputted colsum
	apextract.apall.line='INDEF' # set line equal to user-inputted col
	apextract.apall.lower=-3.5
	apextract.apall.upper=3.5
	apextract.apall.find='yes'
	apextract.apall.recenter='yes'
	apextract.apall.resize='yes'
	apextract.apall.edit='yes'
	apextract.apall.trace='yes'
	apextract.apall.fittrace='yes'
	apextract.apall.extract='yes'
	apextract.apall.extras='yes'
	apextract.apall.nfind=1 # setting this explicitly in call below so no query (for standard star)
	apextract.apall.b_function='legendre' # b. stuff = background, only used for standards
	apextract.apall.b_order=1
	apextract.apall.b_low_reject=1.0 
	apextract.apall.b_high_reject=1.0
	apextract.apall.b_niterate=1
	apextract.apall.b_grow=1
	apextract.apall.clean='yes'
	apextract.apall.weights='none'
	apextract.apall.t_nsum=25 # increased from 15 to 25 (summing/tracing over more lines)
	apextract.apall.t_nlost=200
	apextract.apall.t_low_reject=1.0 # since there are so many data points (over 1000 generally)
	apextract.apall.t_high_reject=1.0 
	apextract.apall.t_step=3 # decreased from 5 to 3 => more intervals (more points)
	apextract.apall.t_niterate=2 # newly added
	apextract.apall.t_grow=2 # newly added
	apextract.apall.t_function='legendre'
	apextract.apall.t_order=2
	apextract.apall.lsigma=2.0
	apextract.apall.usigma=2.0
	apextract.apall.background='none'
	apextract.apall.readnoise=3 # verify whether rdnoise is nearly the same for all images (SALT CCD)
	apextract.apall.gain=1 # verify whether gain is nearly the same for all images (SALT CCD)
	# interactive apall: yes for science; no for standards
	apextract.apall.interactive='yes'
	# changes some parameters again if user indicated that the line looks really faint
	if faint==0:
		apextract.apall.nsum=-3000
		apextract.apall.b_naverage=-5
		apextract.apall.b_niterate = 2
		apextract.apall.t_nsum=50
		apextract.apall.t_step=50
		apextract.apall.t_nlost=50
	first = twodimage[:3]
	if first == 'sci':
		apextract.apall.interactive='yes'
		print 'running apall interactively on science image: '+twodimage
	elif first=='std':
		apextract.apall.interactive='no'
		apextract.apall.background='fit' # since standards aren't background-subtracted; may not be necessary though
		apextract.apall.weights='variance' # for background
		print 'running apall non-interactively on standard star image: '+twodimage
	# since there are 2 extensions (0 and 1), need to specify extraction extension: 1
	twodimage = twodimage+'[1]'
	# This runs apextract.apall
	apextract.apall(input=twodimage,output=spectrum,nfind=1,trace='yes',fittrace='yes',recenter='yes',resize='yes',edit='yes',extract='yes')


# This function sets the parameters for the pyraf task apextract.apsum and then runs it.  
def run_apsum(twodimage,refimage,spectrum):
	# this resets the parameters of apextract.apsum
	apextract.apsum.unlearn()
	apextract.apsum.interactive='no' # non-interactive
	apextract.apsum.review='no'
	apextract.apsum.background='none'
	apextract.apsum.format='multispec'
	# apextract.apsum.clean='yes'
	# apextract.apsum.weights='variance'
	apextract.apsum.nsum=50
	apextract.apsum.lsigma=2.0
	apextract.apsum.usigma=2.0
	# need to specify extraction extension (1, not 0) since there are 2
	twodimage = twodimage+'[1]'
	# This runs apextract.apsum
	apextract.apsum(input=twodimage,output=spectrum,references=refimage)
	
	
# This function sets the parameters for and runs identify
def run_identify(arcsci,lamplines):
	# this resets the parameters of onedspec.identify
	onedspec.identify.unlearn()
	onedspec.identify.nsum=10
	onedspec.identify.match=-3.0
	onedspec.identify.maxfeatures=50
	onedspec.identify.zwidth=100.0
	onedspec.identify.ftype='emission'
	onedspec.identify.threshold=0.0
	onedspec.identify.function='spline3'
	onedspec.identify.order=1
	onedspec.identify.niterate=0
	onedspec.identify.low_reject=3.0
	onedspec.identify.high_reject=3.0
	onedspec.identify.grow=0.0
	onedspec.identify.autowrite='No'
	onedspec.identify.database='database'
	onedspec.identify.section='middle line'
	
	# this runs longslit.identify
	onedspec.identify(images=arcsci,coordlist=lamplines)
	

# This function sets the parameters for and runs the pyraf task reidentify.
def run_reidentify(input,ref,lamplines):
	# this resets the parameters of onedspec.identify
	onedspec.reidentify.unlearn()
	onedspec.reidentify.interactive='no'
	onedspec.reidentify.newaps='yes'
	onedspec.reidentify.override='no'
	onedspec.reidentify.refit='yes'
	onedspec.reidentify.trace='yes'
	onedspec.reidentify.step=10
	onedspec.reidentify.nsum=10
	onedspec.reidentify.shift=0.
	onedspec.reidentify.search=0.0
	onedspec.reidentify.nlost=0
	onedspec.reidentify.threshold=0.0
	onedspec.reidentify.addfeatures='no'
	onedspec.reidentify.match=-3.0
	onedspec.reidentify.maxfeatures=50
	onedspec.reidentify.verbose='yes'
	# this runs longslit.reidentify non-interactively
	onedspec.reidentify(reference=ref,images=input,coordlist=lamplines)


# This function sets the parameters for the pyraf task onedspec.dispcor and then runs it.
def run_dispcor(original,corrected):
	# this resets the parameters of onedspec.dispcor
	onedspec.dispcor.unlearn()
	onedspec.dispcor.verbose='yes'
	
	# This runs onedspec.dispcor (REFSPEC1 keyword addition was done by identifyarcs and reidentifyarcs)
	onedspec.dispcor(input=original,output=corrected)
	

# This function sets the parameters for the pyraf task onedspec.standard and then runs it.
# The user has already been alerted by fluxcal function above about the star_name.
def run_standard(standardimage,outputname,exp,air):
	# this resets the parameters of onedspec.standard
	onedspec.standard.unlearn()
	onedspec.standard.caldir='noaolib$/onedstds/ctiocal/'
	onedspec.standard.extinction=''
	# onedspec.standard.apertures='1' # only want to interactively choose bandpasses in aperture 1
	onedspec.standard.interact='yes' # interactive defining of bandpasses, need at least 15 bandpasses 
	onedspec.standard.airmass=air
	onedspec.standard.exptime=exp
	onedspec.standard.answer='yes' # user will only manually enter in star_name (will press ? when asked)
	# magnitude of standard star (apparent/absolute) not given to us in SALT headers
	
	# this runs onedspec.standard resulting in an std file named outputname 
	onedspec.standard(input=standardimage,output=outputname)
	

# This function sets the parameters for the pyraf task onedspec.sensfunc and then runs it.
def run_sensfunc(stddata,sensname):
	# this resets the parameters of onedspec.sensfunc
	onedspec.sensfunc.unlearn()
	onedspec.sensfunc.apertures='1' # use only the first aperture (object data)
	onedspec.sensfunc.function='spline3'
	onedspec.sensfunc.order=6
	onedspec.sensfunc.extinction=''
	onedspec.sensfunc.ignoreaps='yes' # create only one sensfunc file
	onedspec.sensfunc.interactive='yes' # interactive fitting 
	onedspec.sensfunc.graphs='sri'
	onedspec.sensfunc.answer='yes' # user shouldn't have to input anything else, just check fit and quit
	
	# runs onedspec.sensfunc resulting in a file named sensname
	onedspec.sensfunc(standards=stddata,sensitivity=sensname)
	

# This function sets the parameters for the pyraf task onedspec.calibrate and then runs it.
def run_calibrate(scienceimage,fluximage,sensfilename):
	# this resets the parameters of onedspec.calibrate
	onedspec.calibrate.unlearn()
	onedspec.calibrate.extinct='no'
	onedspec.calibrate.extinction=''
	onedspec.calibrate.sensitivity=sensfilename # the sensitivity function for flux calibration
	onedspec.calibrate.ignoreaps='yes' # look for and use sens*, not sens*0001, sens*0002, etc.
	# This stores some header information from the science image
	scihdr = pyfits.getheader(scienceimage,1)
	exptime = scihdr['EXPTIME']
	airmass = abs(scihdr['AIRMASS'])
	# this continues to set the necessary parameters of onedspec.standard
	onedspec.calibrate.airmass=airmass
	onedspec.calibrate.exptime=exptime
	
	# runs onedspec.calibrate resulting in a flux-calibrated spectrum with name fluximage (on [SCI] extension)
	scienceimage1 = scienceimage+'[SCI]'
	onedspec.calibrate(input=scienceimage1,output=fluximage)
	

# This function sets the parameters for the pyraf task onedspec.scombine and then runs it.
def run_scombine(inputlistfile,finalspectrumname):
	# this resets the parameters of scombine
	onedspec.scombine.unlearn()
	onedspec.scombine.group='apertures'
	onedspec.scombine.combine='average'
	onedspec.scombine.reject='minmax'
	onedspec.scombine.lthreshold=4.235E-18 #sn2012fr value was 4.235E-18
	onedspec.scombine.lsigma = 3.0
	onedspec.scombine.hsigma = 3.0
	onedspec.scombine.mclip = 'yes' # doesn't matter, not used for reject=minmax
	onedspec.scombine.blank=0.0
	
	# running onedspec.scombine resulting in the final wave- and flux-calibrated spectrum named finalspectrumname
	onedspec.scombine(input=inputlistfile,output=finalspectrumname)
	

# This section prompts the user for inputs asking them what task(s) they want to do.
while True:
	print "Welcome to the Rutgers Supernova Reduction Pipeline for SALT data!"
	print "Please make sure you are in the directory with all of your fits images."
	print "Here are the options for reducing your data."
	print "0. Copy the current directory's files into a new subdirectory."
	print "1. Full reduction of the spectral data."
	print "2. Quit."
	while True:
		answer = raw_input("Please enter a number from the menu: ")
		if answer == '0' or answer == '1' or answer == '2':
			break
		else:
			print "Invalid input. You must enter a number from the menu."
	if answer=='0':
		clone()
	elif answer=='1':
		fullreduce()
	elif answer=='2':
		sys.exit("Thanks for using this pipeline!")
	else:
		print 'You must pick an option from the menu.'
		continue 