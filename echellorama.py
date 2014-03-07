#!/usr/bin/env python2
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

		self.plot_1D_line, = self.plot_1D.plot([],[])
		#~ self.plot_2D_line, = self.plot_2D.plot([],[])
		#~ self.plot_PHD_line, = self.plot_PHD.plot([],[])
		self.plot_orders_line, = self.plot_orders.plot([],[])

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
		self.buttonbg = Gtk.Button('Airglow')
		self.buttonbg.connect('clicked', self.buttonbg_clicked, self.context_id)

		self.hbutton_box.pack_start(self.count_rate_button, True, True, 0)
		self.hbutton_box.pack_start(self.filter_phd_button, True, True, 0)
		self.hbutton_box.pack_start(self.gauss_fit_button, True, True, 0)
		self.hbutton_box.pack_start(self.remove_orders_button, True, True, 0)
		self.hbutton_box.pack_start(self.buttonbg, True, True, 0)

		self.main_box.pack_start(self.menubox, False, False, 0)
		self.main_box.pack_start(self.toolbar, False, False, 0)
		self.main_box.pack_start(self.canvas, True, True, 0)
		self.main_box.pack_start(self.hbutton_box, False, False, 0)
		self.main_box.pack_start(self.statusbar, False, False, 0)

		#~ self.open_file('./2014-03-05-230944.fits')

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
			self.statusbar.push(0, 'Opened File:' + self.filename)
			dialog.destroy()
			self.open_file(self.filename)
		elif response == Gtk.ResponseType.CANCEL:
			dialog.destroy()


	def open_file(self, filename):
		"""
		Opens FITS file and sends off the data to be updated
		FITS format:
		hdulist[0].data is a 2D array of photon counts per pixel
		"""
		hdulist = fits.open(filename)
		#~ self.targname = hdulist[0].header['targname']
		#~ self.scidata = hdulist['sci', 1].data
		self.targname = 'test'
		self.image = hdulist[0].data
		self.photon_list = hdulist[1].data
		hdulist.close()

		self.science(self.image, self.photon_list)


	def science(self, image, photon_list):
		self.update_2D_plot(image)

	### smashing 2D data into 1D to view orders as peaks ####
	# filling in y with the sums of the rows of scidata
		peaks = []
		for i in range(0, len(image[0])):
			# sums ith column
			peaks.append(np.sum(image[i,:]))

	# making an x axis with same dimensions as y
		x = np.linspace(0, len(peaks), num = len(peaks))

	# the orders seem to blend around lyman alpha so i have to select
	# the biggest chunk I could use to get a delta y between the
	# orders. i wil have to make this selectable on the GUI some how
		chunk = [0, 500]

	#cutting out the chunk of data that i selected
		self.xchunk = x[-chunk[1]:]
		self.ychunk = peaks[-chunk[1]:]

	# using scipy.signal.find_peaks_cwt() to find centers of orders.
	# this required scipy version .11.0 or greater
		peak_index_list = signal.find_peaks_cwt(self.ychunk, np.arange(3, 15))

	# plotting chunk of 1D data with self.lines through the centers of the
	# orders to double check how the peak finder did
		self.lines = self.xchunk[peak_index_list]
		self.update_orders_plot(self.ychunk, self.xchunk, self.lines)

	### fake PDH stuff ### (fake data for now)
		PHDfake = './chesstest.fits'
		hdu = fits.open(PHDfake)
		PHD = hdu[1].data['PHD']
		self.update_PHDplot(PHD)
		hdu.close()

		self.dragbox = []

	### extraction of orders ###
	# find w, the widths of the orders (difference between peaks)
		peak_widths = []
		for i in range(1, len(peak_index_list)):
			peak_widths.append(peak_index_list[i] - peak_index_list[i - 1])

	# i have to add an extra w at the end of the array to make it the right
	# size i (hopefully this is kosher)
		peak_widths.append(max(peak_widths) - 4)

	# making arrays of 1s and 0s and extracting the
	# 1d orders by matrix multiplication
	#def extraction(self, peak_index_list, x, scidata):
		zeros = np.zeros( (len(x), 1) )
		#~ index = range( 0, len(self.w) )
		#~ reindex = index[::-1]
		global oneDorders
		oneDorders = {}
		for i in range(len(self.w), 0, -1):
			zeros1 = np.copy(zeros)
			zeros1[len(x)-(np.sum(peak_widths[(i):18])) : len(x) - np.sum(peak_widths[(i+1):18])] = 1
			twoD = image * zeros1

			# making 2d orders in to 1d orders
			Y = []
			for j in range(0, len(x)):
				t = np.sum(twoD[:,j])
				Y.append(t)
			# placing 1d orders in dictionary called oneDorders
			oneDorders[str(i)] = Y

		# Send plotting info to update_1D_plot for GUI (for now using just one order until cross correlation is added to script)
		#~ x = np.linspace(0, len(scidata[0]), num=len(scidata[0]))
		self.odo = oneDorders['16']
		self.update_1D_plot(self.odo, x)
		self.save_pickle(oneDorders)


	def update_2D_plot(self, image):
		self.plot_2D.imshow(image, vmin=0, vmax=255, origin='lower')
		self.canvas.draw()

	def update_orders_plot(self, ychunk, xchunk, lines):
		self.plot_orders_line.set_xdata(ychunk)
		self.plot_orders_line.set_ydata(xchunk)
		self.plot_orders.hlines(lines, 0, 2000, color='purple', label='centers')
		self.plot_orders.autoscale_view(False, True, True)
		self.canvas.draw()

	def update_1D_plot(self, odo, x):
		## if you dont want to new airglow subtracted data to over plot but to replot, uncomment this next line
		#~ self.plot_1D.cla()
		#~ self.plot_1D.set_title("1D Extracted Data")
		#~ self.plot_1D.set_xlabel('pixels')
		#~ self.plot_1D.set_ylabel('intensity')
		#~ self.plt = self.plot_1D.plot(x, self.odo)
		self.plot_1D_line.set_xdata(x)
		self.plot_1D_line.set_ydata(odo)
		self.plot_1D.autoscale_view(True, True, True)
		self.canvas.draw()

	def update_PHDplot(self, PHD):
		## if you dont want to new airglow subtracted data to over plot but to replot, uncomment this next line
		self.plot_PHD.cla()
		self.plot_PHD.set_title('Pulse Height Data')
		self.plot_PHD.hist(PHD, bins=80, histtype='stepfilled')
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


### count rate #####
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

#### phd filter button ##
	def on_filter_phd_button_clicked(self, widget, data):
		self.phd_window = Gtk.MessageDialog(image = None)
		self.phd_window.set_size_request(500, 100)
		self.phd_window.move(400, 300)
		#self.phd_window.connect("delet_event",lambda w,e:)

		mainbox = self.phd_window.get_content_area()
		self.phd_window.add(mainbox)
		thebox = Gtk.HBox(False, 0)

		label = Gtk.Label("Discard PHD between")
		label2 = Gtk.Label('and')

		label.show()
		label2.show()

		self.okbutton = Gtk.Button('Okay')
		self.okbutton.connect('clicked', self.phd_entry_button)

		self.entry = Gtk.Entry()
		self.entry.set_activates_default(True)

		self.entry2 = Gtk.Entry()
		self.entry2.set_activates_default(True)
		self.entry.show()
		self.entry2.show()
		self.okbutton.show()

		thebox.pack_start(label,False,False,0)
		thebox.pack_start(self.entry,False,False,0)

		thebox.pack_start(label2,False,False,0)
		thebox.pack_start(self.entry2,False,False,0)
		mainbox.pack_start(thebox,True,True,0)
		mainbox.pack_start(self.okbutton,True,False,0)


		mainbox.show()
		thebox.show()
		self.phd_window.show()

	def phd_entry_button(self, widget):
		minphd = self.entry.get_text()
		maxphd = self.entry2.get_text()
		phdfilt = [minphd,maxphd]
		self.phd_window.destroy()
		self.filter_phd(phdfilt)

### phd filter function ##
	def filter_phd(self, phdfilt):
		phdfilt = [int(x) for x in phdfilt]

		fakedata = './chesstest.fits'
		hdu = fits.open(fakedata)
		PHD = hdu[1].data['PHD']
		PHD = np.array(PHD)
		data  = hdu[1].data
		newdata = data[(PHD > phdfilt[0]) & (PHD < phdfilt[1])]

		plt.subplot(221)
		oldplot = plt.plot(data['X'], data['Y'], linestyle='', marker='.')

		plt.subplot(222)
		newplot = plt.plot(newdata['X'], newdata['Y'], linestyle='', marker='.')

		plt.show()

### mouse click on remove orders ###
	def on_remove_orders_button_clicked(self, widget, data):
		self.statusbar.push(data, 'Click on the order to exclude.')

		# y data values of clicked locations (not pixels!)
		#~ orders = []

		def onclick_order(event):
			if self.remove_orders_button.get_active():
				#~ orders.append(event.ydata)
				#~ self.remove_orders(orders, self.lines)
				self.remove_orders(event.ydata, self.lines)

		self.canvas.mpl_connect('button_press_event', onclick_order)

		if [self.remove_orders_button.get_active()] == [False]:
			self.statusbar.push(0, 'Opened File:' + self.filename)

### removing orders
	def remove_orders(self, ydata, lines):
		#~ for order in orders:
			#~ # Find closest order to clicked location
			#~ bad = min(lines, key=lambda x:abs(x - order))
			#~ bad_orders.append(bad)
			#~ bad_orders_index.append(np.where(lines == bad))

		order = min(lines, key=lambda x:abs(x - ydata))
		#~ newlines = [line for line in lines if line == order]
		#~ newlines = filter(lambda x:x == order, lines)
		#~ print newlines
		#~ self.update_orders_plot(self.ychunk, self.xchunk, newlines)

		print 'Number of original orders', len(oneDorders)

		for key, value in oneDorders.items():
			#~ print value
			if value == order:
				del oneDorders[key]

		self.update_orders_plot(self.ychunk, self.xchunk, oneDorders['16'])

		print 'Corrected number of orders', len(oneDorders)

		self.save_pickle(oneDorders)
		#~ self.new_dictionary(bad_orders_index)


	def new_dictionary(self, bad_orders_index):

		maxint = len(oneDorders)
		for i in range(0, len(bad_orders_index)):
			num = int(bad_orders_index[i][0])

			if oneDorders.get(str(num), None):
				oneDorders.pop(str(num))
			for key in oneDorders.keys():
				k = int(key)
				if k > num:
					oneDorders[str(k - 1)] = oneDorders.get(key)
				if k == maxint - 1:
					oneDorders.pop(key)

		print 'Corrected number of orders', len(oneDorders)

		self.save_pickle(oneDorders)
		#for key in oneDorders.keys(): print key


	def save_pickle(self, oneDorders):
		order_dict = oneDorders
		now = datetime.datetime.now()
		date = now.strftime("%m_%d_%Y")
		targname = str(self.targname)
		pickle.dump(order_dict, open('' + targname + '_1D_' + date + '.p', 'wb'))

if __name__ == '__main__':
	win = MyWindow()
	win.connect("delete-event", Gtk.main_quit)
	win.show_all()
	Gtk.main()
