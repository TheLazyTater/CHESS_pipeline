#!/usr/bin/env python2.7
# coding: utf-8
#
#
#

import sys
import math
import pickle
import datetime
import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
import astropy.io.fits as fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import matplotlib.text as text
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_gtk3cairo import FigureCanvasGTK3Cairo as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar
from gi.repository import Gtk, GObject, Gdk


class MyWindow(Gtk.Window):

	def __init__(self):
		Gtk.Window.__init__(self, title="Echelle Reduction GUI")

		self.set_default_size(1000, 800)
		self.figure = Figure(figsize=(5,7), dpi=100)
		self.plot_1D = self.figure.add_subplot(212)
		self.plot_2D = self.figure.add_subplot(232)
		self.plot_PHD = self.figure.add_subplot(231)
		self.plot_orders = self.figure.add_subplot(233)
		self.plot_orders.tick_params(axis='both', labelsize=6)
		self.plot_orders.set_title("Orders")
		self.plot_2D.tick_params(axis='both', labelsize=7)
		self.plot_2D.set_title("2D Raw Data")
		self.plot_1D.set_title("1D Extracted Data")
		self.plot_1D.set_xlabel('pixels')
		self.plot_1D.set_ylabel('intensity')
		self.plot_1D.tick_params(axis='both', labelsize=7)
		self.plot_PHD.set_title('Pulse Height Data')
		self.plot_PHD.tick_params(axis='both', labelsize=7)

		self.canvas = FigureCanvas(self.figure)

		self.menubar = Gtk.MenuBar()

		self.menubar_file = Gtk.MenuItem("File")
		self.filemenu = Gtk.Menu()
		self.menubar_file.set_submenu(self.filemenu)
		self.menubar.append(self.menubar_file)

		self.filemenu_open = Gtk.MenuItem("Open")
		self.filemenu_open.connect('activate', self.open_file_dialog)
		self.filemenu.append(self.filemenu_open)

		self.filemenu_exit = Gtk.MenuItem('Exit')
		self.filemenu_exit.connect('activate', Gtk.main_quit)
		self.filemenu.append(self.filemenu_exit)

		self.menubox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
		self.menubox.pack_start(self.menubar, False, False, 0)

		self.toolbar = NavigationToolbar(self.canvas, self)
		self.main_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
		self.add(self.main_box)

		self.statusbar = Gtk.Statusbar()
		self.context_id = self.statusbar.get_context_id("stat bar example")
		self.statusbar.push(0, 'Please Open 2D Fits Data File')

		self.hbutton_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
		self.count_rate_button = Gtk.ToggleButton(label='Raw count rate')
		self.count_rate_button.connect("toggled", self.on_count_rate_button_clicked, self.context_id)
		self.filter_phd_button = Gtk.Button('Filter PHD')
		self.filter_phd_button.connect("clicked", self.on_filter_phd_button_clicked, self.context_id)
		self.gauss_fit_button = Gtk.ToggleButton(label='Fit 1D Gauss')
		self.gauss_fit_button.set_active(False)
		self.gauss_fit_button.connect("toggled", self.on_gauss_fit_button_clicked, self.context_id)
		self.remove_orders_button = Gtk.ToggleButton(label='Remove Orders')
		self.remove_orders_button.connect("toggled", self.on_remove_orders_button_clicked, self.context_id)
		self.add_orders_button = Gtk.ToggleButton(label='Add Orders')
		self.add_orders_button.connect("toggled", self.on_add_orders_button_clicked, self.context_id)
		self.buttonbg = Gtk.Button('Airglow')
		self.buttonbg.connect('clicked', self.buttonbg_clicked, self.context_id)

		self.hbutton_box.pack_start(self.count_rate_button, True, True, 0)
		self.hbutton_box.pack_start(self.filter_phd_button, True, True, 0)
		self.hbutton_box.pack_start(self.gauss_fit_button, True, True, 0)
		self.hbutton_box.pack_start(self.remove_orders_button, True, True, 0)
		self.hbutton_box.pack_start(self.add_orders_button, True, True, 0)
		self.hbutton_box.pack_start(self.buttonbg, True, True, 0)

		self.main_box.pack_start(self.menubox, False, False, 0)
		self.main_box.pack_start(self.toolbar, False, False, 0)
		self.main_box.pack_start(self.canvas, True, True, 0)
		self.main_box.pack_start(self.hbutton_box, False, False, 0)
		self.main_box.pack_start(self.statusbar, False, False, 0)

		#~ self.filename = './2014-03-14-185937.fits'
		#~ self.open_file(self.filename)

	def open_file_dialog(self, widget):
		dialog = Gtk.FileChooserDialog(
			"Please Choose a File",
			self,
			Gtk.FileChooserAction.OPEN,
			(	Gtk.STOCK_CANCEL,
				Gtk.ResponseType.CANCEL,
				Gtk.STOCK_OPEN,
				Gtk.ResponseType.OK ) )
		dialog.set_default_response(Gtk.ResponseType.OK)
		filter = Gtk.FileFilter()
		filter.set_name('fits Files')
		filter.add_mime_type('fits')
		filter.add_pattern('*.fits')
		dialog.add_filter(filter)

		response = dialog.run()

		if response == Gtk.ResponseType.OK :
			self.filename = dialog.get_filename()
			self.statusbar.push(0, 'Opened File: ' + self.filename)
			dialog.destroy()
			self.open_file(self.filename)
		elif response == Gtk.ResponseType.CANCEL:
			dialog.destroy()


	def open_file(self, filename):
		hdulist = fits.open(filename)
		self.photon_list = hdulist[1].data
		hdulist.close()

		self.science(self.photon_list)


	def science(self, photon_list, min_phd=0, max_phd=255):
		PHD = np.array(photon_list['PHD'])
		photon_list_filtered = photon_list[(PHD >= min_phd) & (PHD <= max_phd)]
		image, xedges, yedges = np.histogram2d(photon_list_filtered['X'], photon_list_filtered['Y'], bins=2048, range=[[0,8192],[0,8192]])
		self.update_2D_plot(image)

	# Collapse 2D data into 1D to view orders as peaks
		peaks = []
		for i in range(0, len(image[0])):
			# sums ith row
			peaks.append(np.sum(image[i,:]))

	# Smooth peaks by convolving with a boxcar
		boxcar = np.zeros(300)
		boxcar[100:199] = 1.0
		smooth = np.convolve(peaks, boxcar, mode='same')
		self.peaks_smoothed = peaks / smooth

	# using scipy.signal.find_peaks_cwt() to find centers of orders.
	# this required scipy version .11.0 or greater
		self.orders = signal.find_peaks_cwt(self.peaks_smoothed, np.arange(2, 15))

	# Plot 1D data with lines through the centers of the
	# orders to double check how the peak finder did
		self.update_orders_plot(self.peaks_smoothed, self.orders)

	#
		self.update_PHD_plot(PHD, min_phd, max_phd)

		self.dragbox = []

	### extraction of orders ###
	# find the widths of the orders (difference between peaks)
		peak_widths = []
		for i in range(1, len(self.orders)):
			peak_widths.append(self.orders[i] - self.orders[i - 1])

	# Double last entry to make peak_widths the right length
		peak_widths.append(peak_widths[-1])


	# Add up spectrum lines from 2D plot where y coord is an order +/-FWHM
	# Remember, peak_widths[0] is orders[1]-orders[0].
		FWHM = 0.63 # full width half max
		spectrum = []
		for i in range(len(self.orders)):
			peak = int(self.orders[i])
			#~ width = int(peak_widths[i]*FWHM)
			width = 1
			for j in range(0, width):
				# Find chord length of circular image
				r = len(image[0])/2
				start = int(r - (peak * (2*r - peak))**0.5)
				end = int(2 * r - start)
				for k in range(start, end):
					spectrum.append(image[peak-width/2+j][k])

		self.update_1D_plot(spectrum)

		self.statusbar.push(0, 'Done! ' + self.filename)


	def update_2D_plot(self, image):
		self.plot_2D.cla()
		self.plot_2D.tick_params(axis='both', labelsize=7)
		self.plot_2D.set_title("2D Raw Data")
		max = np.amax(image)/2
		self.plot_2D.imshow(image, norm=LogNorm(vmin=1, vmax=max/2), origin='lower')
		self.canvas.draw()


	def update_orders_plot(self, peaks, orders):
		self.plot_orders.cla()
		self.plot_orders_line, = self.plot_orders.plot([],[])
		self.plot_orders.tick_params(axis='both', labelsize=6)
		self.plot_orders.set_title("Orders")

		self.plot_orders_line.set_xdata(peaks)
		self.plot_orders_line.set_ydata(np.arange(len(peaks)))
		self.plot_orders.hlines(orders, 0, max(peaks), color='purple', label='centers')
		self.plot_orders.relim()
		self.plot_orders.autoscale_view(True, True, True)
		self.canvas.draw()


	# Seems a max of 18980 x values are supported by Cairo.  Since
	# we have more than the max we have to reduce the spectrum
	# resolution to fit.
	# Maybe do this dynamically based on xlim()?
	# xmin, xmax = xlim()   # return the current xlim
	def update_1D_plot(self, spectrum):
		self.plot_1D.cla()
		self.plot_1D_line, = self.plot_1D.plot([],[])
		self.plot_1D.set_title("1D Extracted Data")
		self.plot_1D.set_xlabel('pixels')
		self.plot_1D.set_ylabel('intensity')
		self.plot_1D.tick_params(axis='both', labelsize=7)

		MAX = 18980 #16384 # 2^14
		scale = int(len(spectrum)/MAX) + 1

		print len(spectrum), scale, len(spectrum)/scale

		# When we get wavelength calibration, change this line to
		#~ x = np.linspace(minWL, maxWL, num = len(spectrum)/scale)
		x = np.arange(0, len(spectrum), scale)
		spectrum = [np.sum(spectrum[i:i+scale]) for i in x]

		self.plot_1D_line.set_xdata(x)
		self.plot_1D_line.set_ydata(spectrum)
		self.plot_1D.relim()
		self.plot_1D.autoscale_view(True, True, True)
		self.plot_1D.set_xlabel('pixels')# (x' + str(scale) + ')')
		self.canvas.draw()


	def update_PHD_plot(self, PHD, min_phd, max_phd):
		self.plot_PHD.cla()
		self.plot_PHD.set_title('Pulse Height Data')
		self.plot_PHD.tick_params(axis='both', labelsize=7)
		self.plot_PHD.axvline(x=min_phd, color='purple')
		self.plot_PHD.axvline(x=max_phd, color='purple')
		self.plot_PHD.hist(PHD, bins=256, range=[0,255], histtype='stepfilled')
		self.plot_PHD.relim()
		self.plot_PHD.autoscale_view(True, True, True)
		self.canvas.draw()


## airglow button ##
	def buttonbg_clicked(self, widget, data):
		dialog = Gtk.FileChooserDialog("Please Choose Airglow File", self,
			Gtk.FileChooserAction.OPEN,
			(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
			 Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
		dialog.set_default_response(Gtk.ResponseType.OK)
		filter = Gtk.FileFilter()
		filter.set_name('fits Files')
		filter.add_mime_type('fits')
		filter.add_pattern('*.fits')
		dialog.add_filter(filter)

		response = dialog.run()

		if response == Gtk.ResponseType.OK:
			self.airglow_filename = dialog.get_filename()
			#~ global _AirglowFile
			#~ _AirglowFile = fname
			#print "open file" + fname
			#self.statusbar.push(0,'Opened File:' + fname)
			dialog.destroy()
			self.open_airglow_file(self.airglow_filename)
		elif response == Gtk.ResponseType.CANCEL:
			dialog.destroy()

# opening a‏iglow fits file
	def open_airglow_file(self, airglow_filename):
		#this opens up the fits file
		hdulist = fits.open(airglow_filename)

	## this line will need to be used for real data
		targname = hdulist[0].header['targname']
		#targname = self.targname

	## but for now
		#~ targname = 'HD128627J'

	#  this picks out the actual data from fits file, and turns it into numpy array
		self.glowdata = hdulist[0].data
		hdulist.close()
	# simple subtraction of aurglow image from science image
		scidata = self.scidata - self.glowdata
		self.science(scidata, targname)

## gauss fitting button
	def on_gauss_fit_button_clicked(self, widget, data):

		self.statusbar.push(data,'Ready to fit.  Click on both sides of the emission feature you wish to fit')
		self.xdata = []
		def onclick(event):
			if self.gauss_fit_button.get_active():
				self.xdata.append(event.xdata)
				self.statusbar.push(data, 'one more click...')
				if len(self.xdata) == 2:
					self.statusbar.push(data, 'Ready to fit.  Click on both sides of the emission feature you wish to fit')
					xdata = self.xdata
					self.gauss_fit(xdata)


		#  mouse click event on 1d
		cid = self.canvas.mpl_connect('button_press_event', onclick)
		if [self.gauss_fit_button.get_active()] == [False]:
			self.statusbar.push(0, 'Opened File:' + File)


### gauss fitting ###
	def gauss_fit(self, xdata):

		x = list(self.x)
		xg = []
		xg.append(int(xdata[0]))
		xg.append(int(xdata[1]))
		xg1 = min(xg)
		xg2 = max(xg)
		xgauss = x[xg1:xg2]
		ygauss = self.odo[xg1:xg2]
		right = ygauss[len(xgauss)-4:len(xgauss)]
		left = ygauss[0:4]

	# background subtraction
		averight = sum(right) / len(right)
		aveleft = sum(left) / len(left)
		bg_y = [averight, aveleft]
		rightx = xgauss[len(xgauss)-4:len(xgauss)]
		leftx = xgauss[0:4]
		averightx = sum(rightx) / len(rightx)
		aveleftx = sum(leftx) / len(leftx)
		bg_x = [averightx, aveleftx]

		m,b = np.polyfit(bg_x, bg_y, 1)
		slopex = [i * m for i in xgauss]
		bg_fit = slopex + b



		# making a model gauss
		def gauss(xgauss, MAX, mu, sigma):
			return MAX * np.exp(-(xgauss - mu)**2 / (2.0 * sigma**2))

		avex = sum(xgauss) / len(xgauss)
		guess = [1.0, avex, 1.0]

		# plugging in model to matplotlibs curve_fit()
		coeff, var_matrix = curve_fit(gauss, xgauss, ygauss-bg_fit, p0=guess)
		fit = gauss(xgauss, *coeff)
		sigma = coeff[2]
		FWHM = sigma * 2 * np.sqrt(2 * np.log(2))
		FWHM = round(FWHM, 2)
		fitplot = plt.plot(xgauss, ygauss, color='k')
		plt.plot(xgauss, fit + bg_fit, color='b', linewidth=1.5)
		xpos = xgauss[0] + 0.01 * (coeff[1] - xgauss[0])
		strFWHM = str(FWHM)
		plt.text(xpos, 0.9 * max(ygauss), 'FWHM = ' + strFWHM + '', color='purple', fontweight='bold')

		center = str(round(coeff[1], 2))
		plt.text(xpos, 0.95 * max(ygauss), 'Center = ' + center + '', color='green', fontweight='bold')
		plt.plot(xgauss, bg_fit, 'r--')

		plt.show()
		self.xdata = []


	### count rate button
	def on_count_rate_button_clicked(self, widget, data):
		if self.count_rate_button.get_active():
		   self.statusbar.push(data,'Use zoom feature in navigation bar to select count rate region')
		else:
		   self.statusbar.push(0, 'Opened File:' + self.filename)

		dragbox = []

		def onclick(event):
			if self.count_rate_button.get_active():
				dragbox.append(event.xdata)
				dragbox.append(event.ydata)


		def offclick(event):
			if self.count_rate_button.get_active():
				dragbox.append(event.xdata)
				dragbox.append(event.ydata)
				self.count_rate(dragbox, data)


		self.canvas.mpl_connect('button_press_event', onclick)
		self.canvas.mpl_connect('button_release_event', offclick)


	def count_rate(self, dragbox, data):
		# fake exposure time in seconds
		datafake = './chesstest.fits'
		hdu = fits.open(datafake)
		exptime = hdu[0].header['EXPOSURE']
		dragbox = [int(x) for x in dragbox]
		cntbox = self.scidata[ dragbox[0]:dragbox[2], 1024-dragbox[1]:1024-dragbox[3] ]
		totpix = np.size(cntbox)
		cntrate = np.sum(cntbox) / exptime
		totpix = str(totpix)
		cntrate = str(cntrate)
		self.statusbar.push(data, 'count rate in box = ' + cntrate + ' cnt/sec,    pixels in box = ' + totpix + '')
		return cntrate


	def on_filter_phd_button_clicked(self, widget, data):
		self.phd_window = Gtk.MessageDialog(image = None)
		self.phd_window.set_size_request(500, 100)
		self.phd_window.move(400, 300)
		#self.phd_window.connect("delet_event",lambda w,e:)
		#~ self.phd_window.connect("delete-event", self.phd_window.destroy)

		mainbox = self.phd_window.get_content_area()
		thebox = Gtk.HBox(False, 0)

		label = Gtk.Label("Keep PHD between")
		label2 = Gtk.Label('and')

		label.show()
		label2.show()

		self.ok_button = Gtk.Button('Okay')
		self.ok_button.connect('clicked', self.phd_entry_button_clicked)
		self.min_phd_entry = Gtk.Entry()
		self.min_phd_entry.set_activates_default(True)
		self.max_phd_entry = Gtk.Entry()
		self.max_phd_entry.set_activates_default(True)

		self.min_phd_entry.show()
		self.max_phd_entry.show()
		self.ok_button.show()

		thebox.pack_start(label,False,False,0)
		thebox.pack_start(self.min_phd_entry,False,False,0)

		thebox.pack_start(label2,False,False,0)
		thebox.pack_start(self.max_phd_entry,False,False,0)
		mainbox.pack_start(thebox,True,True,0)
		mainbox.pack_start(self.ok_button,True,False,0)

		mainbox.show()
		thebox.show()
		self.phd_window.show()

	def phd_entry_button_clicked(self, widget):
		min_phd = int(self.min_phd_entry.get_text())
		max_phd = int(self.max_phd_entry.get_text())
		self.phd_window.destroy()
		self.statusbar.push(0, 'Filtering by: ' + str(min_phd) + ' < PHD < ' + str(max_phd))
		self.science(self.photon_list, min_phd, max_phd)


	def on_add_orders_button_clicked(self, widget, data):
		if [self.add_orders_button.get_active()] == [True]:
			self.statusbar.push(data, 'Click where to add an order.')

			def onclick_peak(event):
				if self.add_orders_button.get_active():
					self.add_order(event.ydata)

			self.cid_add = self.canvas.mpl_connect('button_press_event', onclick_peak)

		else:
			self.canvas.mpl_disconnect(self.cid_add)
			del self.cid_add
			self.statusbar.push(0, 'Opened File:' + self.filename)

	def add_order(self, ydata):
		order = ydata

		print 'Number of original orders', len(self.orders)

		self.orders.insert(0, order)
		#~ self.orders.sort() # Don't need to sort yet
		self.update_orders_plot(self.peaks_smoothed, self.orders)

		print 'Corrected number of orders', len(self.orders)


	def on_remove_orders_button_clicked(self, widget, data):
		if [self.remove_orders_button.get_active()] == [True]:
			self.statusbar.push(data, 'Click on the order to remove.')

			def onclick_order(event):
				if self.remove_orders_button.get_active():
					self.remove_order(event.ydata)

			self.cid_remove = self.canvas.mpl_connect('button_press_event', onclick_order)

		else:
			self.canvas.mpl_disconnect(self.cid_remove)
			del self.cid_remove
			self.statusbar.push(0, 'Opened File:' + self.filename)


	def remove_order(self, ydata):
		order = min(self.orders, key=lambda x:abs(x - ydata))

		print 'Number of original orders', len(self.orders)

		self.orders.remove(order)
		self.update_orders_plot(self.peaks_smoothed, self.orders)

		print 'Corrected number of orders', len(self.orders)


if __name__ == '__main__':
	win = MyWindow()
	win.connect("delete-event", Gtk.main_quit)
	win.show_all()
	Gtk.main()
