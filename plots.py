import sys
sys.path.append('/home/mvarvera/HiPACE++/PositronPWFA/analysis/')
sys.path.append('/home/mvarvera/src/hipace/tools/')

import read_insitu_diagnostics as diag
import defs
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import constants

def ExmBy_plot_v2(data_class, name, ext='png'):
	s = 2.
	col = 'm'
	pd = 0.025
	prfl = 'k'

	plt.close()
	fig, axs = plt.subplots(1, 1, sharex='col', figsize=(6, 6)) # , gridspec_kw={'width_ratios': [1.75, 1]})
	# plt.subplots_adjust(hspace=0.075)

	im = axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.ExmBy.T, cmap = 'RdBu', vmin = -s, vmax = s) # level 0
	axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.ExmBy_lev1.T, cmap = 'RdBu', vmin = -s, vmax = s) # level 1
	axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.jz_beam.T, cmap = 'seismic_rT', vmin = -7.5e1, vmax = 7.5e1) # level 0

	axs.plot(data.info.z, 4.e-8 * data.profile[data.iteration] - 4, prfl, alpha = .5, linewidth=.5)

	fs = 26
	
	#divider = make_axes_locatable(axs)
	#cax = divider.append_axes("top", size = "5%", pad = "4%")
	#cb = plt.colorbar(im, cax = cax, extend='both', orientation='horizontal', location='top')
	cb = plt.colorbar(im, extend='both', orientation='horizontal', pad=pd, location='top')
	cb.set_label(r'$(E_x - B_y)/E_0 $', fontsize = fs, labelpad=8)
	cb.ax.xaxis.set_ticks_position('top')
	cb.ax.xaxis.set_label_position('top')
	
	ax2 = axs.twinx()
	ax2.plot(data_class.info_lev1.z, data_class.Ez_lev1, color = col, linewidth=.5) # level 1
	ax2.set_ylim(-.8, .8)
	ax2.set_ylabel(r'$E_z/E_0$',  labelpad = 0.5, fontsize = fs, color = col)
	ax2.spines["right"].set_color(col)
	# ax2.spines["left"].set_visible(False)
	ax2.tick_params(axis='y', colors=col)
	#divider2 = make_axes_locatable(axs)
	
	#cax2 = divider2.append_axes("top", size = "5%", pad = pd)
	#divider3 = make_axes_locatable(ax2)
	#cax3 = divider3.append_axes("right", size = "4%", pad = pd)
	#cax3.remove()
	
	
	


	axs.set_xlim(data_class.info.zmin, data_class.info.zmax)
	#axs.set_xticks([i for i in range(int(data_class.info.zmin), int(data_class.info.zmax))])
	axs.set_xlabel(r'$k_p\xi$', fontsize = fs)

	axs.set_xlim(-15,6)
	axs.set_ylim(-4, 4)
	#axs.set_ylim(-6, 6)
	axs.set_ylabel(r'$k_px$', fontsize = fs)

	plt.savefig(f'{name}_v2.{ext}', dpi = 500, bbox_inches = 'tight')
	plt.close()


def ExmBy_plot(data_class, name, ext='png'):
	s = 2.
	col = 'm'
	pd = 1.1 # .9
	prfl = 'k'

	plt.close()
	fig, axs = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(10, 6)) # , gridspec_kw={'width_ratios': [1.75, 1]})
	# plt.subplots_adjust(hspace=0.075)

	im = axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.ExmBy.T, cmap = 'RdBu', vmin = -s, vmax = s) # level 0
	axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.ExmBy_lev1.T, cmap = 'RdBu', vmin = -s, vmax = s) # level 1
	axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.jz_beam.T, cmap = 'seismic_rT', vmin = -7.5e1, vmax = 7.5e1) # level 0
	
#	xd,wd,zd = data.ts.get_particle(var_list=['x','w','z'], iteration=data.iteration, species="drive")	
#	xp,wp,zp = data.ts.get_particle(var_list=['x','w','z'], iteration=data.iteration, species="witness")	
#	xr,wr,zr = data.ts.get_particle(var_list=['x','w','z'], iteration=data.iteration, species="recovery")
#	
#	data.customCMAP('seismic')
#	vlim = 10
#	axs.hist2d(zd[wd>0], xd[wd>0], bins=2000, cmap='seismic_rT', vmin=-vlim, vmax=vlim)
#	axs.hist2d(zp[wp>0], xp[wp>0], bins=2000, cmap='seismic_T', vmin=-vlim, vmax=vlim)
#	axs.hist2d(zr[wr>0], xr[wr>0], bins=2000, cmap='seismic_rT', vmin=-vlim, vmax=vlim)

	#axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.jz_beam_lev1.T, cmap = 'seismic_rT', vmin = -7.5e1, vmax = 7.5e1) # level 1

	axs.plot(data.info.z, 4.e-8 * data.profile[data.iteration] - 4, prfl, alpha = .5, linewidth=.5)
	#axs.plot(data_class.info.z, 1e-7 * data_class.profile[data_class.iteration] - 6, 'k', alpha = .5, linewidth=.5)

	# axs.vlines(-11., data_class.info.xmin, data_class.info.xmax, color = 'k', linestyle = '--', alpha = .5)
	fs = 26

	ax2 = axs.twinx()
	# ax2.hlines(.4, data_class.info.zmin, data_class.info.zmax, color = col, linestyle = '--', alpha = .5)
	# ax2.plot(data_class.info.z, data_class.Ez, color = col) # level 0
	ax2.plot(data_class.info_lev1.z, data_class.Ez_lev1, color = col, linewidth=.5) # level 1
	ax2.set_ylim(-.8, .8)
	ax2.set_ylabel(r'$E_z/E_0$',  labelpad = 1, fontsize = fs, color = col)
	ax2.spines["right"].set_color(col)
	# ax2.spines["left"].set_visible(False)
	ax2.tick_params(axis='y', colors=col)
	divider2 = make_axes_locatable(axs)
	cax2 = divider2.append_axes("right", size = "4%", pad = pd)
	divider3 = make_axes_locatable(ax2)
	cax3 = divider3.append_axes("right", size = "4%", pad = pd)
	cax3.remove()
	cb2 = plt.colorbar(im, cax = cax2, extend='both')
	cb2.set_label(r'$(E_x - B_y)/E_0 $', fontsize = fs)

	axs.set_xlim(data_class.info.zmin, data_class.info.zmax)
	#axs.set_xticks([i for i in range(int(data_class.info.zmin), int(data_class.info.zmax))])
	axs.set_xlabel(r'$k_p\xi$', fontsize = fs)

	axs.set_xlim(-15,6)
	axs.set_ylim(-4, 4)
	#axs.set_ylim(-6, 6)
	axs.set_ylabel(r'$k_px$', fontsize = fs)

	plt.savefig(f'{name}.{ext}', dpi = 500, bbox_inches = 'tight')
	plt.close()

def Ez_plot(data_class, name, ext='png'):
	s = 1.
	col = 'm'
	pd = .1 # .5
	prfl = 'k'

	fig, axs = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(10, 6)) # , gridspec_kw={'width_ratios': [1.75, 1]})
	# plt.subplots_adjust(hspace=0.075)

	im = axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.ts.get_field(field = 'Ez_lev0', iteration = 0)[0].T, cmap = 'RdBu', vmin = -s, vmax = s) # level 0
	axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.ts.get_field(field = 'Ez_lev1', iteration = 0)[0].T, cmap = 'RdBu', vmin = -s, vmax = s) # level 1
	axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.jz_beam.T, cmap = 'seismic_rT', vmin = -7.5e1, vmax = 7.5e1) # level 0
	#axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.jz_beam_lev1.T, cmap = 'seismic_rT', vmin = -5e1, vmax = 5e1) # level 1

	divider2 = make_axes_locatable(axs)
	cax2 = divider2.append_axes("right", size = "4%", pad = pd)
	# divider3 = make_axes_locatable(ax2)
	# cax3 = divider3.append_axes("right", size = "4%", pad = pd)
	# cax3.remove()
	cb2 = plt.colorbar(im, cax = cax2, extend='both')
	cb2.set_label(r'$E_z/E_0 $', fontsize = 18)


	axs.set_xlim(data_class.info.zmin, data_class.info.zmax)
	axs.set_xlabel(r'$k_p\xi$', fontsize = 18)

	axs.set_ylim(-6, 6)
	axs.set_ylabel(r'$k_px$', fontsize = 18)

	plt.savefig(f'{name}.{ext}', dpi = 300, bbox_inches = 'tight')
	plt.close()

def rho_plot(data_class, name, ext='png'):
	s = 6.5e-2
	col = 'm'
	pd = .1
	prfl = 'k'

	fig, axs = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(10, 6)) # , gridspec_kw={'width_ratios': [1.75, 1]})
	# plt.subplots_adjust(hspace=0.075)

	im = axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.ts.get_field(field = 'rho_lev0', iteration = 0)[0].T*(constants.e * data_class.n0), cmap = 'RdBu', vmin = -s, vmax = s) # level 0
	axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.ts.get_field(field = 'rho_lev1', iteration = 0)[0].T*(constants.e * data_class.n0), cmap = 'RdBu', vmin = -s, vmax = s) # level 1
	# axs.pcolormesh(data_class.info.z, data_class.info.x, data_class.jz_beam.T * IA, cmap = 'RdBuT', vmin = -1e15, vmax = 1e15) # level 0
	# axs.pcolormesh(data_class.info_lev1.z, data_class.info_lev1.x, data_class.jz_beam_lev1.T, cmap = 'seismic_rT', vmin = -1e1, vmax = 1e1) # level 1

	divider2 = make_axes_locatable(axs)
	cax2 = divider2.append_axes("right", size = "4%", pad = pd)
	# divider3 = make_axes_locatable(ax2)
	# cax3 = divider3.append_axes("right", size = "4%", pad = pd)
	# cax3.remove()
	cb2 = plt.colorbar(im, cax = cax2, extend='both')
	cb2.set_label(r'$\rho/(e \cdot n_0) $', fontsize = 18)


	axs.set_xlim(data_class.info.zmin, data_class.info.zmax)
	axs.set_xlabel(r'$k_p\xi$', fontsize = 18)


	axs.set_ylim(-6, 6)
	axs.set_ylabel(r'$k_px$', fontsize = 18)

	plt.savefig(f'{name}.{ext}', dpi = 300, bbox_inches = 'tight')
	plt.close()

if __name__ == '__main__':
	
	step = False
	r = True
	n = True
	
	#for file_num in range(8,15):
	#	dirname = f'LUMI_bigDr_gP_{file_num}'
	#	p = f'/data/data4/mvarvera/hdf5/{dirname}/'
	#	ip = f'/data/data4/mvarvera/insitu/{dirname}/'
	
	dirname = f'front_eta22_13'
	p = f'/data/data4/mvarvera/hdf5/{dirname}/'
	ip = f'/data/data4/mvarvera/insitu/{dirname}/'


	data = defs.Functions(path = p, insitu_path = ip, n0 = 5e16, iteration = 0, normalized = n, recovery = r, mesh_refinement = True, src_path='/home/mvarvera/src/hipace/')

	data.customCMAP()
	
#	print(f'###########################\n\t{file_num}\t\n###########################')
	print(f'Drive {diag.emittance_x(data.driveInsitu["average"])[0]*data.kp_inv:.3e}\t Witness {diag.emittance_x(data.witnessInsitu["average"])[0]*data.kp_inv:.3e}\tRecovery {diag.emittance_x(data.recoveryInsitu["average"])[0]*data.kp_inv:.3e}')
	print(f'Drive {data.charge(q = diag.total_charge(data.driveInsitu))[0] * 1e9:.3f} nC\t Witness {data.charge(q = diag.total_charge(data.witnessInsitu))[0] * 1e9:.3f} nC\tRecovery {data.charge(q = diag.total_charge(data.recoveryInsitu))[0] * 1e9:.3f} nC')
	print(f'Efficiency: {data.quickEfficiency(0):.2f} %')
#	print(f'kp_inv: {data.kp_inv}')
	#print(data.diag.position_std(data.recoveryInsitu['average'])[0])
	#print(data.kp_inv * data.diag.position_std(data.recoveryInsitu['average'])[0])
	#qfactor = constants.e * (1e6 * data.n0) * data.kp_inv**3
	#print(f'Charge conversion factor: {qfactor}')
	#print(f'3 nC --> {3e-9/qfactor}')
	
	ExmBy_plot_v2(data,f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_emb')
	#ExmBy_plot(data,f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_emb')
	#Ez_plot(data,f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_ez')
	#rho_plot(data,f'/home/mvarvera/HiPACE++/Stability/images/{dirname}_rho')
