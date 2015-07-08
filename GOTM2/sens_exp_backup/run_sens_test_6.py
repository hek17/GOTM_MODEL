#!/usr/bin/env python

"""
 Purpose:
 
   To run a series of sensitivity tests which access the impact of chlorophyll concentration
   on light attentuation in water and to quantify the subsequent change in heat flux.
 
 To run the FORTRAN code ~dkn/Documents/projects/Ocean_heat_flux/code/new_version/dwcpn.f90
 
 Procedure:
    1) Runs BIRD/PENGUIN I model to create I(z) profiles
    2) Reads output intensity file to extract the I(z) profile for each of the 12 time slots
    3) Calls optimisation script optLightZ to fit the attenuation coefficients B, K1, K2
    4) Write input files for GOTM run (including optimized extinction parameters)
    5) Runs GOTM

 Inputs:

 Outputs: 
 
    Output files are saved in the following directory structure:
 
    |-- 01_intensity_profiles
    |-- 02_optimisation_output
    |-- 03_GOTM_input_ext
    `-- 04_GOTM_output
    
    Directory contents:
    
    * 01_intensity_profiles:   Profiles of intensity with zepth in the water column, as output 
                               by BIRD/PENGUIN.
    * 02_optimisation_output:  Output coefficients from I profile optimisation.
    
    * 03_GOTM_input_ext:       Modified extinction profiles (which are then symlinked from the 
                               GOTM case directory).
    * 04_GOTM_output:          Output netcdf files from the GOTM run.
    
 Variable definitions:
    z   : Depth profile 
    I   : Intensity profile corresponding to depths in z
    B   : Fraction of high absorption (red) PAR extinction
    K1  : High absorption (red) attenuation coefficient
    K2  : Low absorption (blue) attenuation coefficient

 Example run commands:
 
    === Command line ===
    run_sens_test_5.py

    === iPython ===
    ipython 
    %load_ext autoreload
    %autoreload 2
    %run run_sens_test_5.py
    
Written by: Diane Knappett, Plymouth Marine Laboratory
Date: 22/09/2014
"""
import numpy as np
import os
import pdb
import pylab as pl
import subprocess

from scipy.optimize import fmin_slsqp, curve_fit, leastsq

def run_penguin(B_0, input_file):
   """ 
   Write input file for BIRD/PENGUIN, then run the main fortran90 script
   """
   
   # Define inputs
   bottom_z = '  110.0'   # Bottom depth (f7.1)
   lat      = ' 00.00'    # Latitude     (f6.2)
   lon      = ' -30.00'   # Longitude    (f7.2)
   day      = '001'       # Day of year  (i3)
   alphaB   = '0.0914'    # Biological properties (not required)    (f6.4)
   P_mB     = '  2.76'    # Biological properties (not required)    (f6.2)
   z_m      = '  10.00'   # Biological properties (not required)    (f7.2)
   h        = '  0.00'    #                                         (f6.2)
   sigma    = '   1.00'   # Chlorophyll profile                     (f7.2)
   cloud    = '  0.00'    # Cloud cover (0 for sensitivity test)    (f6.2)
   yelsub   = '0.30'      # Standart yellow substance concentration (f4.2)

   # Write input file
   f = open(input_file, 'w')   
   f.write('bottom_z Latitude Longitude Day alphaB  P_mB  z_m B_0  h  sigma  Cloud(%) Yelsub\n')
   f.write('(f7.1,x,f6.2,x,f7.2,x,i3,x,f6.4,x,f6.2,x,f7.2,x,f8.4,x,f6.2,x,f7.2,x,f6.2,x,f4.2)\n')
   f.write(bottom_z +' '+ lat +' '+ lon +' '+ day +' '+ alphaB +' '+ P_mB +' '+ z_m +' '+ B_0 +' '+ h +' '+ sigma +' '+ cloud +' '+ yelsub +'\n')
   f.close()

   # Run fortran script
   bashCommand = fortran_script +' '+ input_file +' > '+ log_file
   os.system(bashCommand)

   output_file = input_file.replace('.dat', '.out')

   return output_file


def readIz_multi(file_Iz):
   d = open(file_Iz, 'r')     
   vals = d.readlines()
   d.close()  

   intensities = dict()
   for line in vals:
      if line.startswith('Time'):
         currenttime = line.split()[1]
         intensities[currenttime] = []
      elif line.strip() == "":
         continue
      else:
         thisitem = line.strip().split()
         intensities[currenttime].append(thisitem)
    # end if
   
   return intensities  


def funcI(p):
   B, K1, K2 = p
   print p
      
   err =  I/I[0]-(B*np.exp(-K1*z) + (1.-B)*np.exp(-K2*z))
   return np.sqrt(np.mean(err**2))


def plotIz_fit(z, I, I_fit, B, K1, K2, save_file=False):
   fig = pl.figure(figsize=(6,10))
   ax = fig.add_subplot(111)
   ax.scatter(I, -z, edgecolors=None)
   ax.plot(I_fit, -z, lw=2, c='r')
   ax.set_ylim(np.min(-z), 0)
   ax.set_xlim(0,np.max(I)+1)
   ax.set_xlabel(r'I(z) $\rm[W/m^2]$', fontsize=18)
   ax.set_ylabel(r'z $\rm[m]$', fontsize=18)
   
   fig.text(0.15,0.33,r'$I(z) = I_0(B\,exp(-K_1\,z) + (1-B)\,exp(-K_2\,z))$' % B, fontsize=16)
   fig.text(0.65,0.3,'B  = %5.2f' % B, fontsize=16)
   fig.text(0.65,0.27,'K1 = %5.2f' % K1, fontsize=16)
   fig.text(0.65,0.24,'K2 = %5.2f' % K2, fontsize=16)
   
   if save_file != False:
      fig.savefig(save_file, dpi=300)
   
   
def plotIz(z, I0, jerlov_param, our_param):
   
   A, g1, g2 = jerlov_param
   B, K1, K2 = our_param
   
   I_jer = I0*(A*np.exp(-z/g1)+(1.-A)*np.exp(-z/g2))
   
   I_new = I0*(A*np.exp(-z/g1)+(1-A)*(B*np.exp(-K1*z) + (1-B)*np.exp(-K2*z)))
   
   fig = pl.figure(figsize=(6,10))
   ax = fig.add_subplot(111)
   ax.plot(I_jer, -z, lw=2, c='r')
   ax.plot(I_new, -z, lw=2, c='k')


def write_opt_Iz(optI_file):
   f = open(optI_file, 'w')   
   f.write('Optimised parameters (B, K1, K2): \n')
   f.write('\n')
   f.write("%2s: %14.12f" % ('B', B) +'\n')
   f.write("%2s: %14.12f" % ('K1', K1) +'\n')
   f.write("%2s: %14.12f" % ('K2', K2) +'\n')
   f.write('\n')
   f.write("%5s %13s %13s" % ('z', 'I', 'I_fit') + '\n')
   #f.write('z    I              I_fit \n')
   for depth in range(0,len(z),1):
      f.write( "%5.1f %13.9f %13.9f" % (z[depth], I[depth], I_fit[depth]) + '\n')
      #f.write("{:.2f}".format(3.1415926)
      #f.write(str(z[depth]) +'  '+ str(I[depth]) +'  '+ str(I_fit[depth]) + '\n')
   f.close()


def write_gotm_files(gotm_file, form, year, month, day, time, line):
      
   f = open(gotm_file, 'w') 
   
   for y in year:
      for m in month:
         for d in day:
            for t in time:
               date = (y, m, d, t)
               f.write(form % (date + line))
               f.write('\n')
   f.close()       


def rewrite_obs_nml(orig_nml, obs_nml, ext_method, ext_file):
   
   # Read existing obs.nml file
   d = open(orig_nml, 'r')     
   info = d.readlines()
   d.close()  

   data = info
   for i, line in enumerate(info):
      if line.startswith('&extinct'):
         data[i+1]= "   extinct_method = "+ ext_method +", \n"
         data[i+2]= "   extinct_file = '"+ ext_file +"', \n"
      else:
         continue
   # end if
      
   # Write modified obs.nml file
   w = open(obs_nml, 'w')   
   for line in info:
      w.write(line)
   w.close()


def rewrite_gotm_run(orig_nml, gotm_run, yy, mm, dd, tt, outname):
   
   # Read existing obs.nml file
   d = open(orig_nml, 'r')     
   info = d.readlines()
   d.close()  
  
   lyy = len(yy)
   lmm = len(mm)
   ldd = len(dd)
   ltt = len(tt)

   data = info
   for i, line in enumerate(info):
      if line.startswith('&output'):
         data[i+3]= "   out_fn = '"+ outname +"',\n"
      elif line.startswith('&time'):
         data[i+3]= "   start = '%4i-%02d-%02d %02d:00:00', \n" % (yy[0],mm[0],dd[0],tt[0])
         data[i+4]= "   stop = '%4i-%02d-%02d %02d:00:00', \n" % (yy[len(yy)-1],mm[len(mm)-1],dd[len(dd)-1],tt[len(tt)-1])
      else:
         continue
   # end if
   
   # Write modified gotmrun.nml file
   w = open(gotm_run, 'w')   
   for line in info:
      w.write(line)
   w.close()


def run_command(args):
   p        = subprocess.Popen(args)
   out, err = p.communicate()
   
   
# ====== MAIN PROGRAM =====

# Define directories
main_dir       = '/users/rsg/dkn/Documents/projects/Ocean_heat_flux'
sens_dir       = main_dir + '/sens_test'
fortran_script = main_dir + '/code/new_version/x86_64-linux/dwcpn'
log_file       = main_dir + '/code/python_scripts/processing_logs/bird-penguin_log.txt'

gotm_dir       = '/users/rsg/dkn/scratch_local/GOTM/developers_version_mod'
case_dir       = gotm_dir + '/gotm-cases/sens_test'

# Define chlorophyll array
chl_arr  = ['    0.01', '    0.10', '    1.00', '    10.0'] 
name_arr = ['0-01','0-10','1-00','10-0']


for i in range(0,len(chl_arr),1):

   #==============
   # BIRD/PENGUIN
   #==============

   input_file =  sens_dir + '/01_intensity_profiles/sensitivity_chl_'+ name_arr[i] +'.dat'

   # Run BIRD/PENGUIN to create I(z) profiles
   Iz_file = run_penguin(chl_arr[i], input_file)
  
   # Read output I file
   intensities = readIz_multi(Iz_file)
   
   for t, time in enumerate(intensities):
      
      z_str = np.array([ list[:][0] for list in intensities[time]])
      I_str = np.array([ list[:][1] for list in intensities[time]])
      
      z     = z_str.astype(np.float)
      I     = I_str.astype(np.float)
      
      # Skip to next iteration if intensity array contains no values
      if len(I) == 0 : continue

      #==============
      # Optimisation
      #==============

      # Optimisation output file
      optI_file = sens_dir + '/02_optimisation_output/optI_file-chl_'+str.lstrip(chl_arr[i])+'-time_'+str(time)[0:4]+'.txt'

      # Define I0
      I0 = I[0]
      
      # Run optimisation                 
      popt= fmin_slsqp(funcI, x0=(0.1, 0.2, 0.2), bounds=((0.,1.),(-1,10),(-1,10)))
      B, K1, K2 = popt
      
      I_fit = I0*(B*np.exp(-K1*z) + (1.-B)*np.exp(-K2*z))
      
      # Plot output parameters
      png_out  = os.path.basename(Iz_file.replace('.out', '.png'))
      png_path = sens_dir + '/02_optimisation_output/' + png_out
      #save_file = Iz_file.replace('.out', '.png')
      plotIz_fit(z, I, I_fit, B, K1, K2, save_file = png_path)
      
      # Write optimized parameters to ASCII file
      write_opt_Iz(optI_file)

      #========================
      # Write GOTM input files
      #========================

      # Define run period
      yy  = [1998] 
      mm  = [1]
      dd  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] 
      tt  = [0, 6, 12, 18] 
      
      
      # METEOROLOGICAL FILE
      #---------------------
      # meteo.dat columns:
      #   date: yyyy-mm-dd hh:mm:ss
      #   x-component of wind (10 m) in m s-1
      #   y-component of wind (10 m) in m s-1
      #   air pressure (2 m) in hectopascal
      #   dry air temperature (2 m) in Celsius
      #   rel. hum. in % or wet bulb temp. in C or dew point temp. in C
      #   cloud cover in 1/10

      meteo_file   = sens_dir + '/03_GOTM_input_ext/meteo/meteo_file.dat'
      meteo_format = "%4i-%02d-%02d %02d:00:00  %5.2f  %5.2f  %7.2f  %5.2f %5.2f  %5.2f"
      meteo_var    = (1.00, 5.00, 1013.0, 20.0, 80.0, 0.0)
      
      write_gotm_files(meteo_file, meteo_format, yy, mm, dd, tt, meteo_var)
      
      # Create symbolic link in case directory
      args = ["ln","-s", meteo_file, case_dir + "/meteo_file.dat"]
      run_command(args)
      
      #   EXTINCTION FILE
      #---------------------   
      # extinction.dat columns:
      #   A:  observed light extinction: non-visible fraction
      #   g1: observed light extinction: e-folding depth of non-visible fraction
      #   B:  observed light extinction: red fraction
      #   k1: observed light extinction: e-folding depth of red fraction
      #   k2: observed light extinction: e-folding depth of green-blue fraction
      
      ext_name     = 'extinction-chl_'+str.lstrip(chl_arr[i])+'-time_'+str(time)[0:4]+'.txt'
      ext_file     = sens_dir + '/03_GOTM_input_ext/extinction/' + ext_name
      ext_format   = "%4i/%02d/%02d %02d:00:00   %5.3f   %5.3f   %5.3f   %5.3f   %5.3f"
      A            = 0.62
      g1           = 0.60
      ext_var      = (A, g1, B, K1, K2)
      
      write_gotm_files(ext_file, ext_format, yy, mm, dd, tt, ext_var)   
      
      # Create symbolic link in case directory
      args = ["ln","-s", ext_file, case_dir +'/'+ ext_name]
      run_command(args)
      
      # MODIFY '.nml' FILES 
      #---------------------
      
      # Specify extinction method (0= user define, 3=default)
      ext_method = '0'
      
      # GOTM RUN FILE
      # Attributes unique file name to GOTM output
      orig_nml = case_dir + '/orig_nml/gotmrun.nml' 
      gotm_run = sens_dir + '/03_GOTM_input_ext/gotmrun/gotmrun_'+str.lstrip(chl_arr[i])+'-time_'+str(time)[0:4]+'.nml'
      
      if ext_method == '0':
         gotm_outfile = 'sens_test_extinct_case_'+ ext_method +'-chl_'+str.lstrip(chl_arr[i])+'-time_'+str(time)[0:4]
      else:
         gotm_outfile = 'sens_test_extinct_case_' + ext_method
            
      rewrite_gotm_run(orig_nml, gotm_run, yy, mm, dd, tt, gotm_outfile)  
      
      # Create symbolic link in case directory
      args = ["cp", gotm_run,  case_dir + '/gotmrun.nml']
      #args = ["ln","-s", gotm_run,  case_dir + '/gotmrun.nml']
      run_command(args)
      
      
      # OBS.NML 
      # Write GOTM obs.nml
      orig_nml = case_dir + '/orig_nml/obs.nml'
      obs_nml  = sens_dir + '/03_GOTM_input_ext/obs/obs_'+str.lstrip(chl_arr[i])+'-time_'+str(time)[0:4]+'.nml'
      
      rewrite_obs_nml(orig_nml, obs_nml, ext_method, case_dir +'/'+ ext_name)
      
      # Create symbolic link in case directory
      args = ["ln","-s", obs_nml,  case_dir + '/obs.nml']
      run_command(args)
      
      
      #========================
      #        Run GOTM
      #========================
      
      # cd to GOTM directory
      bashCommand = 'cd ' + case_dir
      os.system(bashCommand)
      
      # Run GOTM
      # Run run_megs.sh to generate the l2 file
      args = [case_dir + "/gotm"]
      run_command(args)
      
      # Move output file to sens test directory
      gotm_out = sens_dir + '/04_GOTM_output/' + gotm_outfile
            
      args = ["mv", case_dir +"/"+ gotm_outfile, gotm_out]
      run_command(args)     
      
      #========================
      #    Remove symlinks
      #========================     
     
      # Remove symbolic links
      extf   = '/'+ ext_name
      rm_arr = ["/meteo_file.dat", extf, "/gotmrun.nml", "/obs.nml" ]
      
      for trash in rm_arr:
         print 'Removing file: ', case_dir + trash
         args = ["rm", case_dir + trash]
         run_command(args)   
      
      #pdb.set_trace()
         