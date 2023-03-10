#!/usr/bin/env python
# coding: utf-8


# In[11]:


import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import hvsrpy
from hvsrpy import utils


# ## Time Domain Settings
# ---

# In[12]:


# Input file name (may be a relative or full path).
# file_name = "UT.STN11.A2_C50.miniseed"
file_name = "UT.STN11.A2_C150.miniseed"
file_name = "SSHV/M1P_N2.mseed"
    #file_name = 'SSHV/UT.STN11.A2_C50.miniseed'
    #import obspy
    #stream=obspy.read(file_name)
    #stream2=stream[0:3]
    #stream2.write('M1B_1.mseed',format='mseed')
    # file_name = "UT.STN12.A2_C50.miniseed"
    # file_name = "UT.STN12.A2_C150.miniseed"

    # Window length in seconds. In general low frequency peaks require longer window lengths.
    # See the SESAME guidelines for specific window length recommendations.

def do_hvsr(file_name):
    windowlength = 60

    # Boolean to control whether Butterworth filter is applied.
    # Geopsy does not apply a bandpass filter.
    filter_bool = False
    # Low-cut frequency for bandpass filter.
    filter_flow = 0.1
    # High-cut frequency for bandpass filter.
    filter_fhigh = 50
    # Filter order.
    filter_order = 5

    # Width of cosine taper {0. - 1.}. Geopsy default of 0.05 is equal to 0.1 -> 0.1 is recommended
    width = 0.1


    # ## Frequency Domain Settings
    # ---

    # In[13]:


    # Konno and Ohmachi smoothing constant. 40 is recommended.
    bandwidth = 40

    # Minimum frequency after resampling
    resample_fmin = 0.4
    # Maximum frequency after resampling
    resample_fmax = 20
    # Number of frequencies after resampling
    resample_fnum = 40
    # Type of resampling {'log', 'linear'}
    resample_type = 'log'

    # Upper and lower frequency limits to restrict peak selection. To use the entire range use `None`.
    peak_f_lower = None
    peak_f_upper = None


    # ## HVSR Settings
    # ---

    # In[14]:


    # Method for combining horizontal components {"squared-average", "geometric-mean", "single-azimuth"}.
    # Geopsy's default is "squared-average" -> "geometric-mean" is recommended.
    method = "geometric-mean"
    # If method="single-azimuth", set azimuth in degree clock-wise from north. If method!="single-azimuth", value is ignored.
    azimuth = 0

    # Boolean to control whether frequency domain rejection proposed by Cox et al. (2020) is applied.
    # Geopsy does not offer this functionality.
    rejection_bool = True
    # Number of standard deviations to consider during rejection. Smaller values will reject more windows -> 2 is recommended.
    n = 2
    # Maximum number of iterations to perform for rejection -> 50 is recommended.
    max_iterations = 50

    # Distribution of f0 {"lognormal", "normal"}. Geopsy default "normal" -> "lognormal" is recommended.
    distribution_f0 = "lognormal"
    # Distribution of mean curve {"lognormal", "normal"}. Geopsy default "lognormal" -> "lognormal" is recommended.
    distribution_mc = "lognormal"


    # ## Plot Settings
    # ---

    # In[15]:


    # Manually set the ylimits of the HVSR figures. Default is None so limits will be set automatically.
    ymin, ymax = 0, 10


    # ## Perform Calculation
    # ---

    # In[16]:


    fig = plt.figure(figsize=(6,6), dpi=150)
    gs = fig.add_gridspec(nrows=6,ncols=6)

    ax0 = fig.add_subplot(gs[0:2, 0:3])
    ax1 = fig.add_subplot(gs[2:4, 0:3])
    ax2 = fig.add_subplot(gs[4:6, 0:3])

    if rejection_bool:
        ax3 = fig.add_subplot(gs[0:3, 3:6])
        ax4 = fig.add_subplot(gs[3:6, 3:6])
    else:
        ax3 = fig.add_subplot(gs[0:3, 3:6])
        ax4 = False

    start = time.time()
    sensor = hvsrpy.Sensor3c.from_mseed(file_name)
    bp_filter = {"flag":filter_bool, "flow":filter_flow, "fhigh":filter_fhigh, "order":filter_order}
    resampling = {"minf":resample_fmin, "maxf":resample_fmax, "nf":resample_fnum, "res_type":resample_type}
    hv = sensor.hv(windowlength, bp_filter, width, bandwidth, resampling, method, f_low=peak_f_lower, f_high=peak_f_upper, azimuth=azimuth)
    end = time.time()
    print(f"Elapsed Time: {str(end-start)[0:4]} seconds")

    individual_width = 0.3
    median_width = 1.3
    for ax, title in zip([ax3, ax4], ["Before Rejection", "After Rejection"]):
        # Rejected Windows
        if title=="After Rejection":
            if len(hv.rejected_window_indices):
                label = "Rejected"
                for amp in hv.amp[hv.rejected_window_indices]:
                    ax.plot(hv.frq, amp, color='#00ffff', linewidth=individual_width, zorder=2, label=label)
                    label=None

        # Accepted Windows
        label="Accepted"
        for amp in hv.amp[hv.valid_window_indices]:
            ax.plot(hv.frq, amp, color='#888888', linewidth=individual_width,
                    label = label if title=="Before Rejection" else "")
            label=None

        # Window Peaks
        ax.plot(hv.peak_frq, hv.peak_amp, linestyle="", zorder=2,
                marker='o', markersize=2.5, markerfacecolor="#ffffff", markeredgewidth=0.5, markeredgecolor='k',
                label="" if title=="Before Rejection" and rejection_bool else r"$f_{0,i}$")

        # Peak Mean Curve
        print('wow',hv.mc_peak_frq(distribution_mc),hv.mc_peak_amp(distribution_mc))
        ax.plot(hv.mc_peak_frq(distribution_mc), hv.mc_peak_amp(distribution_mc), linestyle="", zorder=4,
                marker='D', markersize=4, markerfacecolor='#66ff33', markeredgewidth=1, markeredgecolor='k',
                label = "" if title=="Before Rejection" and rejection_bool else r"$f_{0,mc}$")

        # Mean Curve
        label = r"$LM_{curve}$" if distribution_mc=="lognormal" else "Mean"
        ax.plot(hv.frq, hv.mean_curve(distribution_mc), color='k', linewidth=median_width,
                label="" if title=="Before Rejection" and rejection_bool else label)

        # Mean +/- Curve
        label = r"$LM_{curve}$"+" ± 1 STD" if distribution_mc=="lognormal" else "Mean ± 1 STD"
        ax.plot(hv.frq, hv.nstd_curve(-1, distribution_mc),
                color='k', linestyle='--', linewidth=median_width, zorder=3,
                label = "" if title=="Before Rejection" and rejection_bool else label)
        ax.plot(hv.frq, hv.nstd_curve(+1, distribution_mc),
                color='k', linestyle='--', linewidth=median_width, zorder=3)

        # f0 +/- STD
        if ymin is not None and ymax is not None:
            ax.set_ylim((ymin, ymax))
        label = r"$LM_{f0}$"+" ± 1 STD" if distribution_f0=="lognormal" else "Mean f0 ± 1 STD"
        _ymin, _ymax = ax.get_ylim()
        ax.plot([hv.mean_f0_frq(distribution_f0)]*2, [_ymin, _ymax], linestyle="-.", color="#000000")
        ax.fill([hv.nstd_f0_frq(-1, distribution_f0)]*2 + [hv.nstd_f0_frq(+1, distribution_f0)]*2, [_ymin, _ymax, _ymax, _ymin],
                color = "#ff8080",
                label="" if title=="Before Rejection" and rejection_bool else label)
        ax.set_ylim((_ymin, _ymax))

        ax.set_xscale('log')
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("HVSR Amplitude")
        if rejection_bool:
            if title=="Before Rejection":
                print("\nStatistics before rejection:")
                hv.print_stats(distribution_f0)
                c_iter = hv.reject_windows(n, max_iterations=max_iterations,
                                           distribution_f0=distribution_f0, distribution_mc=distribution_mc)
            elif title=="After Rejection":
                fig.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.51, 0), columnspacing=2)

                print("\nAnalysis summary:")
                display(pd.DataFrame(columns=[""], index=["Window length", "No. of windows", "Number of iterations to convergence", "No. of rejected windows"],
                        data=[f"{windowlength}s", str(sensor.ns.nseries), f"{c_iter} of {max_iterations} allowed", str(sum(hv.rejected_window_indices))]))
                print("\nStatistics after rejection:")
                hv.print_stats(distribution_f0)
        else:
            display(pd.DataFrame(columns=[""], index=["Window length", "No. of windows"],
                             data=[f"{windowlength}s", str(sensor.ns.nseries)]))
            hv.print_stats(distribution_f0)
            fig.legend(loc="upper center", bbox_to_anchor=(0.77, 0.4))
            break
        ax.set_title(title)

    norm_factor = sensor.normalization_factor
    for ax, timerecord, name in zip([ax0,ax1,ax2], [sensor.ns, sensor.ew, sensor.vt], ["NS", "EW", "VT"]):
        ctime = timerecord.time
        amp = timerecord.amp/norm_factor
        ax.plot(ctime.T, amp.T, linewidth=0.2, color='#888888')
        ax.set_title(f"Time Records ({name})")
        ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xlim(0, windowlength*timerecord.nseries)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Normalized Amplitude')
        ax.plot(ctime[hv.rejected_window_indices].T, amp[hv.rejected_window_indices].T, linewidth=0.2, color="cyan")

    if rejection_bool:
        axs = [ax0, ax3, ax1, ax4, ax2]
    else:
        axs = [ax0, ax3, ax1, ax2]

    for ax, letter in zip(axs, list("abcde")):
        ax.text(0.97, 0.97, f"({letter})", ha="right", va="top", transform=ax.transAxes, fontsize=12)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)


    fig.tight_layout(h_pad=1, w_pad=2, rect=(0,0.08,1,1))
    plt.show()
    plt.close()

    # ## SESAME (2004) Reliability and Clarity Criteria
    # ---
    #
    # Please note that this functionality is still in the beta stage of development.
    # So please review the classification carefully before incorporating the results into your work.
    # If you believe you have found an error or if you wish to provide other feedback related to this functionality it is both welcomed and appreciated.

    # In[17]:


    reliability = utils.sesame_reliability(hv.meta["Window Length"], len(hv.valid_window_indices), hv.frq, hv.mean_curve(), hv.std_curve(), search_limits=(peak_f_lower, peak_f_upper), verbose=1)
    clarity = utils.sesame_clarity(hv.frq, hv.mean_curve(), hv.std_curve(), hv.std_f0_frq(distribution="normal"), search_limits=(peak_f_lower, peak_f_upper), verbose=1)


    # ## Save Figure to File
    # ---

    # In[18]:
    figure_name_out = file_name[:-6]+'.png'
    #figure_name_out=file_name.split('/')[0].split('.')[0]
    #figure_name_out=file_name.split('/')[0] + '/SSHV/Figures/'+file_name.split('/')[2].split('.')[0]+'.png'        figure_name_out=file_name[:-6]+'.png'
    print('figname',figure_name_out)
    print('flex',figure_name_out)
    fig.savefig(figure_name_out, dpi=300, bbox_inches='tight')
    #plt.close()
    print("Figure saved successfully!")


    # ## Save Results to Text File
    # ---
    # In[19]:


    #file_name_out = "example_output_hvsrpy.hv"
   #file_name_out=file_name.split('/')[0] + '/SSHV/Results/' + file_name.split('.')[0].split('/')[2] + '.hv'
    file_name_out=file_name[:-6]+'.hv'
    print('lal',file_name.split('.')[0])
    hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="hvsrpy")
    #print("Results saved successfully!")   
    print("Results saved successfully!")
    return hv,reliability,clarity

    # if np.sum(reliability) > 2 and np.sum(clarity) > 4:
    #     figure_name_out = "example_hvsr_figure.png"
    #     #figure_name_out=file_name.split('/')[0].split('.')[0]
    #     #figure_name_out=file_name.split('/')[0] + '/SSHV/Figures/'+file_name.split('/')[2].split('.')[0]+'.png'
    #     figure_name_out=file_name[:-6]+'.png'
    #     print('figname',figure_name_out)
    #     print('flex',figure_name_out)
    #     fig.savefig(figure_name_out, dpi=300, bbox_inches='tight')
    #     #plt.close()
    #     print("Figure saved successfully!")


    #     # ## Save Results to Text File
    #     # ---

    #     # In[19]:


    #     #file_name_out = "example_output_hvsrpy.hv"
    #     #file_name_out=file_name.split('/')[0] + '/SSHV/Results/' + file_name.split('.')[0].split('/')[2] + '.hv'
    #     file_name_out=file_name[:-6]+'.png'
    #     print('lal',file_name.split('.')[0])
    #     hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="hvsrpy")
    #     print("Results saved successfully!")
    # else:
    #     print('do nothing')
        # figure_name_out = "example_hvsr_figure.png"
        # # figure_name_out=file_name.split('/')[0].split('.')[0]
        # figure_name_out = file_name.split('/')[0] + '/SSHV/Figures/' + file_name.split('/')[2].split('.')[0] + '.png'
        # print('flex', figure_name_out)
        # fig.savefig(figure_name_out, dpi=300, bbox_inches='tight')
        # # plt.close()
        # print("Figure saved successfully!")
        #
        # # ## Save Results to Text File
        # # ---
        #
        # # In[19]:
        #
        # # file_name_out = "example_output_hvsrpy.hv"
        # file_name_out = file_name.split('/')[0] + '/SSHV/Results/' + file_name.split('.')[0].split('/')[2] + '.hv'
        # print('lal', file_name.split('.')[0])
        # hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="hvsrpy")
        # print("Results saved successfully!")
        # #print('dont save file, hv curve is not reliable')

# ## Save Results to Geopsy-Style Text File
# ---

# In[20]:


# file_name_out = "example_output_geopsy.hv"
#
# hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="geopsy")
# print("Results saved successfully!")


# In[ ]:




