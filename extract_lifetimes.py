#exec(open("extract_lifetimes.py").read())
import ROOT
import numpy as np
import sys

particle_type = 'proton' #proton or alpha
particle_E = 2218
particle_dE = 500#particle_E*0.01 #gate for deviation from center of energy peak


file_path = r'F:\e21010\pxct\Adams_5days_9000pps.root'
file = ROOT.TFile.Open(file_path, "READ")
tree = file.Get("tree")

#show PID plots
pid_canvas = ROOT.TCanvas("pid_plot", "Particle Identification")
tree.Draw("MSD12_e:MSD26_e>>(1000,0,8000,500,0,4000)", "1==1", "colz")
pid_canvas.SetLogz()

proton_selector = "MSD12_e>121 && MSD12_e <853"
alpha_selector = "MSD12_e>1300"

pid_selector = (proton_selector if particle_type == 'proton' else alpha_selector) + '&& (MSD12_e + MSD26_e) > %f && (MSD12_e+MSD26_e < %f)'%(particle_E-particle_dE, particle_E + particle_dE)
selected_pid_canvas = ROOT.TCanvas('selected_pid', 'selected particles')
num_selected_particles = tree.Draw('MSD12_e:MSD26_e>>(1000,0,8000,500,0,4000)', pid_selector , 'colz')
selected_pid_canvas.SetLogz()

#note: these numbers are still hard coded some places, so be careful if changing.
min_xray_energy = 2
max_xray_energy = 12



'''
Transition | Range for mu 
NiKL | 7459-7478
NiKM | 8261-8339
CuKL | 8027-8049
CuKM | 8903-8978
ZnKL | 8615-8640
ZnKM | 9569-9652
'''

#set up PDF for gated xray spectrum which we'll fit with unbinned maximum likilihood
normalized_pdf = ROOT.TF1("normalized_pdf", "[19]*(1+[0]*x + gaus(1) + gaus(4) + gaus(7) + gaus(10) + gaus(13) + gaus(16))/(1*10 + 0.5*[0]*(12*12 - 2*2) + sqrt(2*pi)*([1]*[3] + [4]*[6]+[7]*[9]+[10]*[12]+[13]*[15] + [16]*[18]))", min_xray_energy, max_xray_energy);
normalized_pdf.SetParName(19, "normalization")
normalized_pdf.SetParLimits(19,1,1)
normalized_pdf.SetParameter(19,1)


normalized_pdf.SetNpx(1000)
normalized_pdf.SetParName(0, "m")
normalized_pdf.SetParLimits(0,-0.1,0.1)
for i in range(6):
    normalized_pdf.SetParName(i*3+1, 'a%d'%i)
    normalized_pdf.SetParName(i*3+2, 'mu%d'%i)
    normalized_pdf.SetParName(i*3+3, 'sigma%d'%i)


#constrain central values
n = 0
normalized_pdf.SetParLimits(2, 7.3, 7.6)
normalized_pdf.SetParLimits(5, 8.1,8.5)
normalized_pdf.SetParLimits(8, 7.9, 8.1)
normalized_pdf.SetParLimits(11, 8.8, 9.1)
normalized_pdf.SetParLimits(14, 8.5, 8.7)
normalized_pdf.SetParLimits(17, 9.4, 9.7)
#magnitude
lower_b = 0.01
upper_b = 1000
for n in range(6):
    normalized_pdf.SetParLimits(n*3+1, lower_b, upper_b)
#sigma
#(~260 eV FWHM or 110 eV @1 sigma expected from detector resolution, although there will be more than one peak which will broaden the gaussian)
upper_b = 0.5
lower_b = 0.05
for n in range(6):
    normalized_pdf.SetParLimits(3*n+3, lower_b, upper_b)

#set up canvas for drawing proton gated xray spectrum
bins = 200
xray_canvas_ungated = ROOT.TCanvas('all_xrays', 'all xrays')
dssd3e_hist_ungated  = dssd3e_hist2 = ROOT.TH1F("dssd3e_hist_ungated", "ungated xrays", bins, min_xray_energy, max_xray_energy)
num_all_particle_xrays = tree.Draw("LEGe_e>>dssd3e_hist_ungated", proton_selector)
norm = num_all_particle_xrays/bins*(max_xray_energy - min_xray_energy)
while(tree.UnbinnedFit("normalized_pdf", "LEGe_e", 'LEGe_e>=%f && LEGe_e <= %f'%(min_xray_energy, max_xray_energy))):
    #renormalize to draw on top of plot, draw, and then set back to 1
    normalized_pdf.SetParLimits(19, norm, norm)
    normalized_pdf.SetParameter(19, norm)
    normalized_pdf.DrawCopy("LSAME")
    ROOT.gStyle.SetOptFit(11111111)
    normalized_pdf.SetParLimits(19, 1, 1)
    normalized_pdf.SetParameter(19, 1)
    
normalized_pdf.SetParLimits(19, norm, norm)
normalized_pdf.SetParameter(19, norm)
#restrict FWHM and location of peaks in future fits
for i in range(6):
    mu = normalized_pdf.GetParameter(i*3+2)
    normalized_pdf.SetParLimits(i*3+2, mu, mu)
    sigma = normalized_pdf.GetParameter(i*3+3)
    normalized_pdf.SetParLimits(i*3+3, sigma, sigma)
#restrict background slope
m = normalized_pdf.GetParameter(0)
normalized_pdf.SetParLimits(0,m,m)

xray_canvas = ROOT.TCanvas('xrays_on_protons_canvas', 'xray spectrum gated on %s with energy %f +/- %f keV'%(particle_type, particle_E, particle_dE))
dssd3e_hist2 = ROOT.TH1F("dssd3e_hist2", "xrays gated on protons", bins, min_xray_energy, max_xray_energy)
unbinned_fit_proton_selector = '%s && LEGe_e>=%f && LEGe_e <= %f'%(pid_selector, min_xray_energy, max_xray_energy)
num_gated_xrays = tree.Draw("LEGe_e>>dssd3e_hist2", unbinned_fit_proton_selector)
#try doing an unbinned fit
norm = num_gated_xrays/bins*(max_xray_energy - min_xray_energy)
while(tree.UnbinnedFit("normalized_pdf", "LEGe_e", unbinned_fit_proton_selector)):
    #renormalize to draw on top of plot, draw, and then set back to 1
    normalized_pdf.SetParLimits(19, norm, norm)
    normalized_pdf.SetParameter(19, norm)
    normalized_pdf.Draw("SAME")
    ROOT.gStyle.SetOptFit(11111111)
    normalized_pdf.SetParLimits(19, 1, 1)
    normalized_pdf.SetParameter(19, 1)
normalized_pdf.SetParLimits(19, norm, norm)
normalized_pdf.SetParameter(19, norm)
#tree.UnbinnedFit("normalized_pdf", "LEGe_e", unbinned_fit_proton_selector, 'M')

#calculate number of counts in each peak
#TODO: do this including correlations between parameters
index_to_peak_name = ['NiKL', 'NiKM', 'CuKL', 'CuKM', 'ZnKL', 'ZnKM']
peak_name_to_index = {index_to_peak_name[i]:i for i in range(len(index_to_peak_name))}
counts_in_peak = {}#record (num, +/-)
ps = [normalized_pdf.GetParameter(i) for i in range(19)]
overall_norm = 1/(1*10 + 0.5*ps[0]*(12*12 - 2*2) + np.sqrt(2*np.pi)*(ps[1]*ps[3] + ps[4]*ps[6]+ps[7]*ps[9]+ps[10]*ps[12]+ps[13]*ps[15] + ps[16]*ps[18]))
for peak_name in peak_name_to_index:
    i = peak_name_to_index[peak_name]
    magnitude = normalized_pdf.GetParameter(i*3+1)*overall_norm
    magnitude_error = normalized_pdf.GetParError(i*3+1)*overall_norm
    sigma = normalized_pdf.GetParameter(i*3+3)
    sigma_error = normalized_pdf.GetParError(i*3+3)
    counts = num_gated_xrays*sigma*magnitude*(2*np.pi)**0.5
    counts_error = counts*magnitude_error/magnitude#no longer propagating uncertainty in gaussian width, since that parameter is 0
    counts_in_peak[peak_name] = (counts, counts_error)

#calculate x-rays from accidental coincidences
coincidence_window = 200e-9 #s
beam_rate = 9000 #pps
#TODO: this method only works for low enough rates, fix this
coincidence_prob = beam_rate*coincidence_window
#for E in keV, between 
xray_detector_efficiency = lambda E: 7.8/100#(1.14*E+2.59)/100 #TODO: simulate rathe than read off Lijie's paper
#make dictionairy of xrays per background decay
#Zn would be from 60Ga decay, Cu from Zn, and Ni from Cu
#format is (fraction of decays which emit this xray, uncertainty)
#Cu and Ni are from https://www.nndc.bnl.gov/nudat3/dec_searchi.jsp
#Zn is from 0.092% of 60Ga decays which electron capture + radiative yields from EADL, since our simulation currently doesn't include IC for 60Ga decay
#TODO: incorporate uncertainty for 60Ga decay radiative yield
background_xray_yield = {'ZnKL':(0.110022*0.46587, 0), 'ZnKM':(0.110022*(2.0209e-3+1.0634e-2+1.0237e-2), 0), 'CuKL':(0.0318, 0.00327), 'CuKM':(0.00379, 0.0034), 'NiKL':(0.0255, 0.00103), 'NiKM':(0.00307,0.0001166)}
#peak energies below are only used to determine detector efficiency
peak_energies = {'NiKL':7.47, 'NiKM':8.3, 'CuKL':8.03, 'CuKM':8.94, 'ZnKL':8.62, 'ZnKM':9.6}
coincidence_counts = {}
for peak in background_xray_yield:
    expected_counts = num_selected_particles*coincidence_prob*background_xray_yield[peak][0]*xray_detector_efficiency(peak_energies[peak])
    #TODO: uncertainty from uncertainty in radiative yields
    coincidence_counts[peak] = (expected_counts, expected_counts**0.5)
print('peak, measured counts, accidental coincidences, real counts')
actual_counts = {}
for peak in counts_in_peak:
    actual_counts[peak] = (counts_in_peak[peak][0] - coincidence_counts[peak][0], (counts_in_peak[peak][1]**2 + coincidence_counts[peak][1]**2)**0.5)
    print('%s, %f+/-%f, %f+/-%f, %f+/-%f'%(peak, counts_in_peak[peak][0], counts_in_peak[peak][1], \
            coincidence_counts[peak][0], coincidence_counts[peak][1], \
            actual_counts[peak][0], actual_counts[peak][1]))

#calculate width of state
hbar = 6.582199E-16 #eV-s
kshell_lifetime = 4.22e-16 #s
first_peak = 'ZnKL' #xray peak of element before charged particle emmission
second_peak = 'CuKL' if particle_type == 'proton' else 'NiKL' #xray peak of element after charged particle emmission
background_peak = 'CuKL' if particle_type == 'alpha' else 'NiKL'
#radiative transitioin probabilities from EADL
radiative_transition_prob = {'NiKL':.12106+.236419,\
                    'CuKL':0.131119+.255668, \
                    'ZnKL':0.14062+0.273809}

vacancies_filled_before_cpe = [actual_counts[first_peak][0]/radiative_transition_prob[first_peak], actual_counts[first_peak][1]/radiative_transition_prob[first_peak]]
vacancies_filled_before_cpe[1] = (vacancies_filled_before_cpe[1]**2 + vacancies_filled_before_cpe[0])**0.5 #sqrtN uncertainty for number of events which occur during experiment
print('K-shell vacancies filled prior to charged particle emmission: %f+/-%f'%tuple(vacancies_filled_before_cpe))
vacancies_filled_after_cpe = [actual_counts[second_peak][0]/radiative_transition_prob[second_peak], actual_counts[second_peak][1]/radiative_transition_prob[second_peak]]
vacancies_filled_after_cpe[1] = (vacancies_filled_after_cpe[1]**2 + vacancies_filled_after_cpe[0])**0.5 #sqrtN uncertainty for number of events which occur during experiment
print('K-shell vacancies filled after charged particle emmission: %f+/-%f'%tuple(vacancies_filled_after_cpe))
lifetime = (vacancies_filled_before_cpe[0]/vacancies_filled_after_cpe[0]*kshell_lifetime, \
            vacancies_filled_before_cpe[0]/vacancies_filled_after_cpe[0]*kshell_lifetime*\
                ((vacancies_filled_before_cpe[1]/vacancies_filled_before_cpe[0])**2 + (vacancies_filled_after_cpe[1]/vacancies_filled_after_cpe[0])**2)**0.5)
print('lifetime = %e+/-%e s'%lifetime)

ROOT.gApplication.Run()