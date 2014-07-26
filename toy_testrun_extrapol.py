


from ROOT import TGraphErrors,TCanvas

shift = [0.,10.,50.0,80.0,100.0]
#testrun
mean = [-670.57, -671.80, -678.53, -683.57, -686.50]
mean_err = [0.7, 0.7, 0.71, 0.72, 0.72]
sigma = [21.48, 21.54, 21.91,22.12,22.13]
sigma_err = [0.71, 0.69,0.71, 0.70, 0.71]

#full
mean_full = [0.08, -0.68, -3.07, -4.89, -6.18]
mean_err_full = [0.13, 0.13, 0.14, 0.14, 0.15]
sigma_full = [4.08, 4.00, 4.12, 4.27, 4.41]
sigma_err_full = [0.13, 0.12, 0.13, 0.14, 0.14]


gr_mean_testrun = TGraphErrors(len(shift))
gr_sigma_testrun = TGraphErrors(len(shift))
gr_rel_sigma_testrun = TGraphErrors(len(shift))

gr_mean_full = TGraphErrors(len(shift))
gr_sigma_full = TGraphErrors(len(shift))
gr_rel_sigma_full = TGraphErrors(len(shift))

for i in range(len(shift)):
    gr_mean_testrun.SetPoint(i,shift[i],mean[i])
    gr_mean_testrun.SetPointError(i,0.,2.0) #mean_err[i])
    gr_sigma_testrun.SetPoint(i,shift[i],sigma[i])
    gr_sigma_testrun.SetPointError(i,0.,sigma_err[i])
    gr_rel_sigma_testrun.SetPoint(i,shift[i],(sigma[i]-sigma[0])/sigma_err[0])
    gr_rel_sigma_testrun.SetPointError(i,0.,0.)


    gr_mean_full.SetPoint(i,shift[i],mean_full[i])
    gr_mean_full.SetPointError(i,0.,mean_err_full[i])
    gr_sigma_full.SetPoint(i,shift[i],sigma_full[i])
    gr_sigma_full.SetPointError(i,0.,sigma_err_full[i])
    gr_rel_sigma_full.SetPoint(i,shift[i],(sigma_full[i]-sigma_full[0])/sigma_err_full[0])
    gr_rel_sigma_full.SetPointError(i,0.,0.)

    print gr_sigma_full.GetErrorY(i)






c_testrun_sigma = TCanvas('c_testrun_sigma','c_testrun_sigma',10,10,700,500)
gr_sigma_testrun.SetTitle('; L1 top misalignment in y (#mum);#sigma(z_{vtx}) (mm)')
gr_sigma_testrun.SetMarkerStyle(20)
gr_sigma_testrun.SetMarkerStyle(20)
gr_sigma_testrun.Draw('ALP')

c_testrun_mean = TCanvas('c_testrun_mean','c_testrun_mean',10,10,700,500)
gr_mean_testrun.SetTitle('; L1 top misalignment in y (#mum); <z_{vtx}> (mm)')
gr_mean_testrun.SetMarkerStyle(20)
gr_mean_testrun.SetMarkerStyle(20)
gr_mean_testrun.Draw('ALP')

#c_testrun_rel_sigma = TCanvas('c_testrun_rel_sigma','c_testrun_rel_sigma',10,10,700,500)
#gr_rel_sigma_testrun.SetMarkerStyle(20)
#gr_rel_sigma_testrun.Draw('ALP')

c_full_sigma = TCanvas('c_full_sigma','c_full_sigma',10,10,700,500)
gr_sigma_full.SetTitle('; L1 top misalignment in y (#mum);#sigma(z_{vtx}) (mm)')
gr_sigma_full.SetMarkerStyle(20)
gr_sigma_full.Draw('ALP')


c_full_mean = TCanvas('c_full_mean','c_full_mean',10,10,700,500)
gr_mean_full.SetTitle('; L1 top misalignment in y (#mum); <z_{vtx}> (mm)')
gr_mean_full.SetMarkerStyle(20)
gr_mean_full.Draw('ALP')

#c_full_rel_sigma = TCanvas('c_full_rel_sigma','c_full_rel_sigma',10,10,700,500)
#gr_rel_sigma_full.SetMarkerStyle(20)
#gr_rel_sigma_full.Draw('ALP')


c_rel_sigma = TCanvas('c_rel_sigma','c_rel_sigma',10,10,700,500)
gr_rel_sigma_full_cp = gr_rel_sigma_full.Clone('gr_rel_sigma_full_cp')
gr_rel_sigma_full_cp.SetTitle('Relative z_{vtx} width difference; L1 top misalignment in y (#mum);(#sigma-#sigma_{ideal})/#delta(#sigma_{ideal})')
gr_rel_sigma_full_cp.SetLineColor(2)
gr_rel_sigma_full_cp.SetMarkerColor(2)
gr_rel_sigma_full_cp.SetMarkerStyle(20)
gr_rel_sigma_full_cp.Draw('ALP')
gr_rel_sigma_testrun.SetMarkerStyle(20)
gr_rel_sigma_testrun.Draw('LP,same')


ans = raw_input('pause')



