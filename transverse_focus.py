import sys
#sys.path.append('/home/mvarvera/HiPACE++/PositronPWFA/analysis/')
sys.path.append('/home/mvarvera/src/hipace/tools/')

import numpy as np
import read_insitu_diagnostics as diag
#import defs
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import constants


# This version uses `OpenPMDTimeSeries` and is a little faster because it doesn't calculate a bunch of stuff upon initialization the way `defs` does.

from openpmd_viewer import OpenPMDTimeSeries

def transverse_slice(timeseries, insitu_path, sigma_x, name, ext='png', z_choice = None):
	ExmBy,info = timeseries.get_field(field='ExmBy_lev1', iteration=0)
	witnessInsitu = diag.read_file(insitu_path + 'reduced_witness.0000.txt')
	
	zAx = diag.z_axis(witnessInsitu)
	beam = zAx[(witnessInsitu['sum(w)']!=0)[0]]

	if z_choice:
		choice = np.argmin(np.abs(zAx - z_choice))

	centroid = np.argmin(np.abs(zAx - -12.89))
	head = np.argmin(np.abs(zAx - -11.9))
	tail = np.argmin(np.abs(zAx - -13.42))
	
	plt.figure(figsize = (10, 6))
	plt.plot(info.x, ExmBy[head], 'tab:blue', ls = '-', label = f'Head $k_p\\xi = {zAx[head]:.2f}$')
	plt.plot(info.x, ExmBy[centroid], 'purple', ls = '-', label = f'Centroid $k_p\\xi = {zAx[centroid]:.2f}$')
	plt.plot(info.x, ExmBy[tail], 'tab:red', ls = '-', label = f'Tail $k_p\\xi = {zAx[tail]:.2f}$')
#	plt.plot(info.x, ExmBy[head], 'tab:green', ls = '-', label = f'Head 2')
#	plt.plot(info.x, ExmBy[centroid], 'k', ls = '-', label = f'Centroid 2')
#	plt.plot(info.x, ExmBy[tail], 'b', ls = '-', label = f'Tail 2')
	# plt.plot(info.x, ExmBy[choice], 'm', ls = '-', label = f'$k_p\\xi = {zAx[choice]:.2f}$')

	#ax = plt.gca()
	#ax.axhline(0, c='gray', ls='--', lw=1)
	#ax.axvline(0, c='gray', ls='--', lw=1)
	plt.grid(lw=0.4)

	x = np.linspace(-5*sigma_x, 5*sigma_x, 100)
	plt.plot(x, np.exp( -1/2 * np.power(x, 2) / sigma_x**2 ) - 2., 'gray', label='Radial distribution')

	plt.ylabel('$(E_x-B_y)/E_0$', fontsize=15)
	plt.xlabel('$k_px$', fontsize=15)
	plt.ylim(-2., 2.)
	plt.xlim(info.xmin, info.xmax)
	#plt.xlim(-5*sigma_x, 5*sigma_x)
	plt.legend(loc = 'best', fontsize=13)
	plt.savefig(f'{name}.{ext}', dpi = 300, bbox_inches = 'tight')
	plt.close()

###
def transverse_slices(timeseries, insitu_path, sigma_x, name, ext='png'):
	fs = 26

	ExmBy,info = timeseries.get_field(field='ExmBy_lev1', iteration=0)
	witnessInsitu = diag.read_file(insitu_path + 'reduced_witness.0000.txt')
	
	zAx = diag.z_axis(witnessInsitu)
	beam = zAx[(witnessInsitu['sum(w)']!=0)[0]]

	head = np.argmin(np.abs(zAx - max(beam)))
	tail = np.argmin(np.abs(zAx - min(beam)))
	#head = np.argmin(np.abs(zAx - -11.9))
	#tail = np.argmin(np.abs(zAx - -13.42))
	#print(f'Head idx = {head}\t Tail idx = {tail}')

	cmap = plt.get_cmap('jet')
	colors = cmap(np.linspace(0,1, head-tail+1))

	fig = plt.figure(figsize=(10, 6))

	for field_idx,idx in zip(range(tail, head+1), range(head-tail+1)):
		plt.plot(info.x, ExmBy[field_idx], ls = '-', color=colors[idx], label='_nolegend_')

	sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=zAx[head], vmax=zAx[tail]))
	sm.set_array([])  # Required for ScalarMappable
	ax = plt.gca()
	cbar = plt.colorbar(sm, ax=ax, pad=.01)
	cbar.set_label(r'$k_p\xi$', fontsize=fs)
	cbar.ax.set_title(r'$\mathbf{Head}$', fontsize=20)
	cbar.ax.set_xlabel(r'$\mathbf{Tail}$', fontsize=20)

	#plt.grid(lw=0.4)
	ax.grid(True, lw=0.4, zorder=1)

	ylim = -3.8 # -5 # -3.8 # -2.95 # -2.9 # -2.25
	
	x = np.linspace(info.xmin, info.xmax, 100)
	#x = np.linspace(-5*sigma_x, 5*sigma_x, 100)
	plt.plot(x, np.exp( -1/2 * np.power(x, 2) / sigma_x**2 ) + ylim, c='k', ls='--', lw=0.75, label=r'$\mathrm{Radial\ distribution}$')

	#plt.axhline(y=0, color='k', ls='-', lw=0.35, label='_')
	
	#plt.vlines(x=-sigma_x, ymin=0.03, ymax=1.3, color='k', ls='-.', lw=2, label=r'$k_p x = -\sigma_{x,p}$')
	plt.vlines(x=-sigma_x, ymin=0.03, ymax=1.2, color='k', ls='-.', lw=2, label=r'$k_p x = -\sigma_{x,p}$')
	#plt.vlines(x=-sigma_x, ymin=0.03, ymax=1.2, color='k', ls='-.', lw=2, label=r'$\mathrm{Inset\ slice}$')
	#plt.vlines(x=-sigma_x, ymin=0.17, ymax=1.18, color='k', ls='-.', lw=2, label=r'$\mathrm{Inset\ slice}$')
	#plt.vlines(x=-sigma_x, ymin=0.17, ymax=1.5, color='k', ls='-.', lw=2, label=r'$\mathrm{Inset\ slice}$')
	#plt.vlines(x=-sigma_x, ymin=0.17, ymax=1.12, color='k', ls='-.', lw=2, label=r'$\mathrm{Inset\ slice}$')

	plt.ylabel('$(E_x-B_y)/E_0$', fontsize=fs)
	plt.xlabel('$k_px$', fontsize=fs)
	#plt.ylim(-2., 2.)
	plt.ylim(bottom=ylim, top=None) 
	plt.xlim(info.xmin, info.xmax)
	#plt.xlim(-5*sigma_x, 5*sigma_x)
	plt.legend(loc = 'best')
	
	# These are in unitless percentages of the figure size. (0,0 is bottom left)
	#left, bottom, width, height = [0.21, 0.3, 0.25, 0.25]
	left, bottom, width, height = [0.21, 0.34, 0.25, 0.25] # [0.18, 0.35, 0.25, 0.25] # LUMI_front
	#left, bottom, width, height = [0.205, 0.305, 0.25, 0.25] # [0.185, 0.33, 0.25, 0.25] # [0.185, 0.31, 0.25, 0.25] # [0.18, 0.35, 0.25, 0.25] # [0.18, 0.26, 0.25, 0.25]
	ax2 = fig.add_axes([left, bottom, width, height], zorder=2)

	sigma_x_idx = np.argmin(np.abs(info.x + sigma_x))
	ax2.scatter(zAx[tail:head+1], [ ExmBy[i][sigma_x_idx] for i in range(tail, head+1) ], c=colors, edgecolors='k', s=6, lw=0.1)
	#plt.text(0.03, 0.98, f'$k_p x = -\sigma_x$', ha='left', va='top', transform = ax2.transAxes, fontsize=13)
	#plt.text(0.05, 0.95, r'$k_p x = -\sigma_{x,p}$', ha='left', va='top', transform = ax2.transAxes, fontsize=20)
	#plt.text(0.05, 0.95, f'$k_p x = {info.x[sigma_x_idx]:.2f}$', ha='left', va='top', transform = ax2.transAxes)
	ax2.set_ylabel(r'$(E_x-B_y)/E_0$', fontsize=20)
	ax2.set_xlabel(r'$k_p\xi$', fontsize=20)
	ax2.set_yticks([round(ExmBy[tail][sigma_x_idx],1), 0.5, 1.0])
	#ax2.set_ylim(bottom=round(ExmBy[tail][sigma_x_idx],1), top=None)
	#ax2.yaxis.set_major_formatter('${x:.2g}$')
	#ax2.set_yticks([i/2 for i in range(5)])
	#ax2.grid(lw=0.4)

	plt.savefig(f'{name}.{ext}', dpi = 500, bbox_inches = 'tight')

	plt.close()

def view_focus(timeseries, insitu_path, X_SLICE, sigma_x, name, ext='png'):
	'''
	Parameters
        ----------
	X_SLICE : list of floats
		X-coordinates to slice across (normalized to plasma skin depth). For reference, current sims have sigma_x = .03 k_p
	'''
	s = 1.05 # ylim scale from field extrema
	zmin = -14.5
	zmax = -11.
	cmap = plt.get_cmap('jet')
	colors = cmap(np.linspace(0,1,len(X_SLICE)))
	xmin,xmax= (-5.1, 1.)
	
	# calculate relative x coordinate for slicing and truncate at domain boundaries (assumes x domain goes from {-|xmax|, ..., |xmax|})
	_,info_pre = timeseries.get_field(field='ExmBy_lev1', iteration=0, slice_across='z', slice_relative_position=0)
	
	plt.close()
	fig, axs = plt.subplots(1, 1, figsize=(10, 6))
	
	min_field_bound,max_field_bound = (0,0)
	for idx,x_slice in enumerate(X_SLICE):
		x_rel = np.sign(x_slice) * min(abs(x_slice/info_pre.xmax), 1.)

		field,info = timeseries.get_field(field='ExmBy_lev1', iteration=0, slice_across='x', slice_relative_position=x_rel)

		zmask = np.logical_and(zmin <= info.z, info.z <= zmax)
		field = field[zmask]
		info.z = info.z[zmask] 

		min_field_bound,max_field_bound = ( min(min_field_bound, s * min(field)), max(max_field_bound,s * max(field)) )
		
		axs.plot(info.z, field, color=colors[idx], alpha=.5)
	
	p_beam_profile = diag.per_slice_charge(diag.read_file(insitu_path + 'reduced_witness.0000.txt'))[0] * constants.c # initial on-axis longitudinal current profile of e+ beam	
	scaled_p_beam_profile = 2.5e-7*p_beam_profile[zmask] + xmin # min_field_bound
	axs.plot(info.z, scaled_p_beam_profile, 'k', alpha = .5, linewidth=.5, label='_exclude')
	
	# Annotation placement
	max_idx = np.argmax(scaled_p_beam_profile)  # Index of maximum value
	annotation_x = info.z[max_idx]  # Corresponding z-value
	annotation_y = 0.95*xmin  # *min_field_bound  # scaled_p_beam_profile[max_idx]*2  # Maximum value in the profile
	axs.text(annotation_x, annotation_y, r'On-axis $j_z$', fontsize=14, color='k', ha='center', va='center', alpha=0.6)
	
	axs.set_xlim(zmin, zmax)
	axs.set_xlabel(r'$k_p\xi$', fontsize = 18)

	axs.set_ylim(xmin, xmax)
	#axs.set_ylim(min_field_bound, max_field_bound)
	axs.set_ylabel(r'$(E_x - B_y)/E_0$', fontsize = 18)
	
	# Add colorbar
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min(X_SLICE)/sigma_x, vmax=max(X_SLICE)/sigma_x))
	sm.set_array([])  # Required for ScalarMappable
	cbar = plt.colorbar(sm, ax=axs, pad=.01)
	cbar.set_label(r'$k_p (x/\sigma_x)$', fontsize=18)
	
	plt.savefig(f'{name}.{ext}', dpi = 300, bbox_inches = 'tight')
	plt.close()

if __name__ == '__main__':

	n0 = 5e16 # plasma density
	iteration = 0
	sigma_x = 0.03
#	x = [ 0., sigma_x/4, sigma_x/2, sigma_x, 2*sigma_x ] # desired x-coordinate (normalized to k_p^{-1}) to slice across for plot of focusing field
#	step = sigma_x/32
#	x = np.arange(0., 5*sigma_x + step, step)
	#file_num = 8
	#for file_num in range(8,15):
	
	dirname = f'LUMI_eta16'

	p  = f'/data/data4/mvarvera/hdf5/{dirname}/' # path to field data
	ip = f'/data/data4/mvarvera/insitu/{dirname}/' # path to insitu diagnostics

	ts = OpenPMDTimeSeries(path_to_dir=p, check_all_files=True)
	#view_focus(timeseries=ts, insitu_path=ip, X_SLICE=x, sigma_x=sigma_x,  name=f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_focus')
	#transverse_slice(timeseries=ts, insitu_path=ip, sigma_x=sigma_x,  name=f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_transverse_focus')
	transverse_slices(timeseries=ts, insitu_path=ip, sigma_x=sigma_x,  name=f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_transverse_focus')

