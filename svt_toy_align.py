#!/usr/bin/python


import os,sys

from math import sqrt,log,tan,atan,cos
import utils
from ROOT import TLegend,TCanvas,SetOwnership,TH2F,TH1F,Double,TFile,TGraph,TGraphErrors
import toy_mult_scat as ms
import argparse


def getArgs():
  parser = argparse.ArgumentParser(description='Run toy fitter')
  parser.add_argument('--targetz',required=True,type=float, help='Target position w.r.t. nominal 0')
  parser.add_argument('--debug','-d',action='store_true',help='Debug output.')
  parser.add_argument('--useinvmass',action='store_true',help='Use inv mass.')
  parser.add_argument('--name',default='tmp',help='Name to add to results')
  parser.add_argument('--l1sigma',type=float, default=-1.0, help='hit error L1')
  parser.add_argument('--l2sigma',type=float, default=-1.0, help='hit error L2')
  parser.add_argument('--l3sigma',type=float, default=-1.0, help='hit error L3')
  parser.add_argument('--beamspotwidth',type=float, default=0.5, help='BS width')
  parser.add_argument('--pspectrum', default='p_spectrum.root', help='momentum spectrum file')
  parser.add_argument('--geometry', default='default', help='noL1, full,default')
  parser.add_argument('--l1misaligny',type=float, default=0.0, help='L1 misalign in y')
  parser.add_argument('--l2misaligny',type=float, default=0.0, help='L2 misalign in y')
  parser.add_argument('--l3misaligny',type=float, default=0.0, help='L3 misalign in y')
  #parser.add_argument('--events','-n',type=int,default=-1,help='Max events to process.')
  #parser.add_argument('--half',required=True,help='Top or bottom half tracks')
  #parser.add_argument('--mc','-m',action='store_true',help='Simulation input')
  #parser.add_argument('--save','-s',action='store_true',help='Save output')
  print sys.argv
  args = parser.parse_args();
  print args
  return args


if __name__=='__main__':
    print 'toy align'

    args = getArgs()
    debug = args.debug
    useInvM = args.useinvmass
    
    target_z = args.targetz
    hit_l1_sigma = args.l1sigma
    hit_l2_sigma = args.l2sigma
    hit_l3_sigma = args.l3sigma

    # simulate beam spot (assume Gaussian in y: wrong but try it for now)
    beamspot_sigma_y = args.beamspotwidth
    pspectrum_filename = args.pspectrum
    mom_multiplier = 1.0
    geometry = args.geometry
    
    utils.misalignments[1] = args.l1misaligny
    utils.misalignments[2] = args.l2misaligny
    utils.misalignments[3] = args.l3misaligny

    
    y_beamspot = 0.
    if geometry=='noL1':
        print 'noL1 geom chosen'
        l1_z = 200.
        l1_y = 35.8-20.0 
        nom_angle = atan((l1_y-y_beamspot)/(l1_z-target_z))
        l2_z = 300.
        hit_l3_z = 500.
    elif geometry=='default':
        print 'default geom chosen'
        l1_z = 100.
        l1_y = 36.9-20.0 
        nom_angle = atan((l1_y-y_beamspot)/(l1_z-target_z))
        l2_z = 200.
        l2_y = y_beamspot+tan(nom_angle)*(l2_z-target_z)
        l3_z = 300.
    elif geometry=='full':
        print 'full geom chosen'
        l1_z = 100.
        l1_y = 1.5 
        nom_angle = atan((l1_y-y_beamspot)/(l1_z-target_z))
        l2_z = 200.
        l2_y = y_beamspot+tan(nom_angle)*(l2_z-target_z)
        l3_z = 300.
    else:
        print 'Not a valid geometry: %s' % geometry
        sys.exit(1)
    
    h_p = None
    if pspectrum_filename!='':
        f = TFile(pspectrum_filename,"READ")
        if 'invM' in pspectrum_filename:
            useInvM = True
            h_p = f.Get('h_invM')
        else:
            h_p = f.Get("htmp")
        #f.Close()
    
    # initialize some variables
    plotDir = 'plots'
    

    h_zvtx = TH1F('h_zvtx','z position at crossing;z (mm);Entries',100,target_z-50.-abs(target_z*0.5),target_z+50.+abs(target_z*0.5))
    #h_zvtx = TH1F('h_zvtx','z position at crossing;z (mm);Entries',100,-950,-550)
    h_yvtx = TH1F('h_yvtx','y position at crossing;y (mm);Entries',100,-3.,3.)
    h_zposAtY0 = TH1F('h_zposAtY0','z position at y=0;z (mm);Entries',100,target_z-10.-abs(target_z*0.75),target_z+10.+abs(target_z*0.75))
    h_ypos_target = TH1F('h_ypos_target','y position at target;y (mm);Entries',100,-4.,4.)
    h_ypos_target_p = TH2F('h_ypos_target_p','Width of y position at target;Track momentum (GeV);#sigma(y) (mm)',20,0.,3.0*mom_multiplier,50,-2.5,2.5)
    h_ypos_target_slope = TH2F('h_ypos_target_slope','Y position at target;Track slope (rad);y (mm)',50,0.01,0.05,50,-10.,10.)
    h_zposAtY0_2 = TH1F('h_zposAtY0_2','z position at y=0 (trk2);z (mm);Entries',100,target_z-10.-abs(target_z*0.75),target_z+10.+abs(target_z*0.75))
    h_ypos2_target = TH1F('h_ypos2_target','y position at target (trk2);y (mm);Entries',100,-4.,4.)
    h_ypos2_target_p = TH2F('h_ypos2_target_p','Width of y position at target (trk2);Track momentum (GeV);#sigma(y) (mm)',20,0.,3.0*mom_multiplier,50,-2.5,2.5)
    h_ypos2_target_slope = TH2F('h_ypos2_target_slope','Y position at target;Track slope (rad);y (mm)',50,0.01,0.05,50,-10.,10.)
    h_angle = TH1F('h_angle','Angle of tracks;Angle [mrad];Entries',100,atan(nom_angle)*1.0e3-30,atan(nom_angle)*1.0e3+30)
    h_angle_2 = TH1F('h_angle_2','Angle of tracks (trk2);Angle [mrad];Entries',100,-1*(atan(nom_angle)*1.0e3)-30,-1.0*(atan(nom_angle)*1.0e3)+30)
    h_y1 = TH1F('h_y1','h_y1',50,-1.0*(atan(nom_angle)*(l1_z-target_z))-5.,atan(nom_angle)*(l1_z-target_z)+5.)
    h_y2 = TH1F('h_y2','h_y2',50,-1.0*(atan(nom_angle)*(l2_z-target_z))-5.,atan(nom_angle)*(l2_z-target_z)+5.)
    h_y3 = TH1F('h_y3','h_y3',50,-1.0*(atan(nom_angle)*(l3_z-target_z))-5.,atan(nom_angle)*(l3_z-target_z)+5.)
    h_invMass = TH1F('h_invMass','Invariant mass;Invariant mass (GeV);Entries',50,0.,0.3)
    h_invMass_truth = TH1F('h_invMass_truth','Invariant mass;Invariant mass (GeV);Entries',50,0.,0.3)
    h_chi2 = TH1F('h_chi2','h_chi2',50,0.,5)
    h_chi2ndf = TH1F('h_chi2ndf','h_chi2ndf',50,0.,5)
    h_chi2_2 = TH1F('h_chi2_2','h_chi2',50,0.,5)
    h_chi2ndf_2 = TH1F('h_chi2ndf_2','h_chi2ndf',50,0.,5)
    h_hit_sigma = TH2F('h_hit_sigma',';z (mm);Hit uncertainty (incl. MS) (mm)',9,-50,850,50,0,0.5);
    #h_tmp = TH2F('h_tmp',';z (mm);Hit position y (mm)',10,target_z-100.,l3_z+50.,10,-23,23.)
    h_tmp = TH2F('h_tmp',';Hit z (mm);Hit y (mm)',10,target_z-100.,l3_z+50.,10,-(y_beamspot+tan(nom_angle)*(l3_z-target_z))*1.5,(y_beamspot+tan(nom_angle)*(l3_z-target_z))*1.5)
    h_trk_p = TH1F('h_trk_p',';Track momentum (GeV);Entries',100,0.,3.*mom_multiplier)
    h_tmp.SetStats(False)
    ct = TCanvas('ct','ct',10,10,700,500)
    h_tmp.Draw()
    #fit_vec = []
    gr_theta_p = TGraph()

    fits = []
    grbs = []
    for i in range(5000):
        # create a basic track at nominal angle from a beamspot position
        y_beamspot = utils.rnd.Gaus(0.,beamspot_sigma_y)
        if h_p!=None:
            if useInvM:
                invariant_mass = h_p.GetRandom()
                print 'M ',invariant_mass,' nom_angle ', nom_angle
                p = sqrt(invariant_mass*invariant_mass/(2.0*(1-cos(nom_angle*2.0))))
            else:
                p = h_p.GetRandom()
            p = p*mom_multiplier
            h_trk_p.Fill(p)
            p = p*1.0e3 #MeV
            E = sqrt(p*p+ms.me*ms.me)
            beta = p/E
            x = ms.t/ms.X0
            theta = ms.MSAngle(x,ms.z,beta,p)
            hit_l1_sigma = ms.single_hit_res
            hit_l2_sigma = sqrt((ms.single_hit_res*ms.single_hit_res) + (l2_z-l1_z)*(l2_z-l1_z)*theta*theta)            
            hit_l3_sigma = sqrt((ms.single_hit_res*ms.single_hit_res) + (l3_z-l1_z)*(l3_z-l1_z)*theta*theta + (l3_z-l2_z)*(l3_z-l2_z)*theta*theta)
            gr_theta_p.SetPoint(i,p,theta)
        h_hit_sigma.Fill(l1_z,hit_l1_sigma)
        h_hit_sigma.Fill(l2_z,hit_l2_sigma)
        h_hit_sigma.Fill(l3_z,hit_l3_sigma)

        #place the hits w.r.t. to L1 hit at innermost strip in v3 detector
        # randomly choose angle of particle
        rndm_angle = 0. # utils.rnd.Rndm()*0.02
        angle = atan((l1_y-y_beamspot)/(l1_z-target_z)) + rndm_angle
	angle_2 = angle+utils.rnd.Gaus(0.,theta)
	angle_3 = angle_2+utils.rnd.Gaus(0.,theta)
	hit_l1_y = y_beamspot+tan(angle)*(l1_z-target_z)
	hit_l2_y = hit_l1_y+tan(angle_2)*(l2_z-l1_z)
	hit_l3_y = hit_l2_y+tan(angle_3)*(l3_z-l2_z)
	base_track = utils.Track()
	base_track.addHit(l1_z,hit_l1_y,hit_l1_sigma,ms.single_hit_res)
	base_track.addHit(l2_z,hit_l2_y,hit_l2_sigma,ms.single_hit_res)
	base_track.addHit(l3_z,hit_l3_y,hit_l3_sigma,ms.single_hit_res)

        #hit_l2_y = y_beamspot+tan(angle)*l2_z
        #hit_l3_y = y_beamspot+tan(angle)*l3_z
        #base_track = utils.Track()
        #base_track.addHit(l1_z,y_beamspot+(l1_z-target_z)*atan(angle),hit_l1_sigma)
        #base_track.addHit(l2_z,y_beamspot+(l2_z-target_z)*atan(angle),hit_l2_sigma)
        #base_track.addHit(l3_z,y_beamspot+(l3_z-target_z)*atan(angle),hit_l3_sigma)
        # get the new track
        track = base_track.getWiggledTrack(True)
        fit = utils.Fit(track,'pol1')
        fit.SetRange(target_z-800.,l3_z+800.)

        if debug: print 'i=%d angle %f' % (i,angle)

        angle = atan((-1.0*l1_y-y_beamspot)/(l1_z-target_z)) - rndm_angle
	angle_2 = angle+utils.rnd.Gaus(0.,theta)
	angle_3 = angle_2+utils.rnd.Gaus(0.,theta)
	hit_l1_y = y_beamspot+tan(angle)*(l1_z-target_z)
	hit_l2_y = hit_l1_y+tan(angle_2)*(l2_z-l1_z)
	hit_l3_y = hit_l2_y+tan(angle_3)*(l3_z-l2_z)
	base_track = utils.Track()
	base_track.addHit(l1_z,hit_l1_y,hit_l1_sigma,ms.single_hit_res)
	base_track.addHit(l2_z,hit_l2_y,hit_l2_sigma,ms.single_hit_res)
	base_track.addHit(l3_z,hit_l3_y,hit_l3_sigma,ms.single_hit_res)

        track2 = base_track.getWiggledTrack(False)
        fit2 = utils.Fit(track2,'pol1')
        #fit2.mirrorX(target_z)
        fit2.SetRange(target_z-800.,l3_z+800.)
        if i<5 and debug: 
            fit.printProp()
            fit2.printProp()
        h_zposAtY0.Fill(fit.getZAtY(0.))
        h_ypos_target.Fill(fit.getYAtZ(target_z))
        h_ypos_target_p.Fill(p*1.0e-3,fit.getYAtZ(target_z))
        target_z_displ = 0.0 #mm
        y_tmp = -1.
        if fit.getSlope() > 0:
          y_tmp = fit.getYAtZ(target_z + target_z_displ)
        else:
          y_tmp = -1. * fit.getYAtZ(target_z + target_z_displ)          
        h_ypos_target_slope.Fill(abs(fit.getSlope()),y_tmp)
        h_zposAtY0_2.Fill(fit2.getZAtY(0.))
        h_ypos2_target.Fill(fit2.getYAtZ(target_z))
        h_ypos2_target_p.Fill(p*1.0e-3,fit2.getYAtZ(target_z))
        y_tmp = -1.
        if fit2.getSlope() > 0:
          y_tmp = fit2.getYAtZ(target_z + target_z_displ)
        else:
          y_tmp = -1. * fit2.getYAtZ(target_z + target_z_displ)          
        h_ypos2_target_slope.Fill(abs(fit2.getSlope()),y_tmp)
        zvtx = utils.getPol1Crossing(fit,fit2)
        yvtx = fit.getYAtZ(zvtx)
        h_zvtx.Fill(zvtx)
        h_yvtx.Fill(yvtx)
        h_angle.Fill(fit.getSlope()*1.0e3)
        h_angle_2.Fill(fit2.getSlope()*1.0e3)
        h_y1.Fill(fit.track.hits[0].y)
        h_y2.Fill(fit.track.hits[1].y)
        h_y3.Fill(fit.track.hits[2].y)
        h_y1.Fill(fit2.track.hits[0].y)
        h_y2.Fill(fit2.track.hits[1].y)
        h_y3.Fill(fit2.track.hits[2].y)
        h_chi2.Fill(fit.getChi2())
        h_chi2ndf.Fill(fit.getChi2()/fit.getNDF())
        h_chi2_2.Fill(fit2.getChi2())
        h_chi2ndf_2.Fill(fit2.getChi2()/fit2.getNDF())
        
        opening_angle = fit.getSlope()-fit2.getSlope()
        M = 1.0e-3*sqrt(2.0*p*p*(1-cos(opening_angle)))
        if debug: print 'M ',M
        h_invMass.Fill(M)
        if useInvM:
            M_truth = invariant_mass
        else:
            M_truth = 1.0e-3*sqrt(2.0*p*p*(1-cos((nom_angle+rndm_angle)*2)))
        h_invMass_truth.Fill(M_truth)

        # plot a few of the fits
        ct.cd()
        if i<10:
            if i==1:
                grb = track.getGraph()
                grb.SetMarkerStyle(19+i)
                grb.SetMarkerSize(1.5)
                grb.Draw('same,P')            
                grbs.append(grb)
                grb2 = track2.getGraph()
                grb2.SetMarkerStyle(19+i)
                grb2.SetMarkerSize(1.5)
                grb2.Draw('same,P')            
                grbs.append(grb2)
            fnc = fit.getFunction()
            fnc.Draw("same,L")
            fnc2 = fit2.getFunction()
            fnc2.Draw("same,L")
        fits.append(fit)
        fits.append(fit2)
                
    
    if h_p!=None:
        saveStr = '%s_bsy%.1f_geom%s_momdistr'%(args.name,beamspot_sigma_y,geometry)
    else: 
        saveStr = '%s_bsy%.1f_%.0f_%.0f_%.0f_%.0f'%(args.name,beamspot_sigma_y,target_z,hit_l1_sigma*1.0e3,hit_l2_sigma*1.0e3,hit_l3_sigma*1.0e3)

    ct.SaveAs('%s/fits_3hits_%s.png'%(plotDir,saveStr),'png')

    #ans = raw_input('pause')

    #sys.exit(0)

    cp = TCanvas('cp','cp',10,10,700,500)
    h_trk_p.Draw()
    cp.SaveAs('%s/trk_mom_distr_3hits_%s.png'%(plotDir,saveStr),'png')

    chs = TCanvas('chs','chs',10,10,700,500)
    chs.Divide(2,2)
    chs.cd(1)
    h_hit_sigma.SetTitle('Hit uncertainty (incl. MS);Hit position z (mm); Hit uncertainty (mm)')
    h_hit_sigma.Draw("colz")
    chs.cd(2)
    #h_hit_sigma.ProjectionY('pr_l1',2,2).Draw()
    h_sig_l2 = h_hit_sigma.ProjectionY('pr_l2',h_hit_sigma.GetXaxis().FindBin(l2_z),h_hit_sigma.GetXaxis().FindBin(l2_z))
    mean_l2_sigma = h_sig_l2.GetMean()
    h_sig_l2.Draw()
    h_sig_l3 = h_hit_sigma.ProjectionY('pr_l3',h_hit_sigma.GetXaxis().FindBin(l3_z),h_hit_sigma.GetXaxis().FindBin(l3_z))
    h_sig_l3.SetLineColor(2)
    mean_l3_sigma = h_sig_l3.GetMean()
    h_sig_l3.Draw('same')
    utils.myText(0.6,0.6,'<#sigma>_{L2}=%.3fmm'%mean_l2_sigma,0.05,1)
    utils.myText(0.6,0.53,'<#sigma>_{L3}=%.3fmm'%mean_l3_sigma,0.05,2)
    chs.cd(3)
    gr_theta_p.SetTitle('Multiple scattering angle;Track momentum (MeV);#theta_{0}')
    gr_theta_p.Draw('AXP')
    
    chs.SaveAs('%s/ms_distr_3hits_%s.png'%(plotDir,saveStr),'png')

    c2 = TCanvas('c2','c2',10,10,700,500)
    c2.Divide(2,2)
    c2.cd(1)
    h_y1.Draw()
    c2.cd(2)
    h_y2.Draw()
    c2.cd(3)
    h_y3.Draw()


    c3 = TCanvas('c3','c3',10,10,700,700)
    c3.Divide(1,2)
    c3.cd(1)
    h_angle.Fit('gaus')
    h_angle.SetStats(False)
    h_angle.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2f'%(h_angle.GetFunction('gaus').GetParameter(2)),0.05,1)
    c3.cd(2)
    h_angle_2.Fit('gaus')
    h_angle_2.SetLineColor(3)
    h_angle_2.Draw('')
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2f'%(h_angle.GetFunction('gaus').GetParameter(2)),0.05,1)

    c3.SaveAs('%s/angle_distr_3hits_%s.png'%(plotDir,saveStr),'png')

    c1 = TCanvas('c1','c1',10,10,700,500)
    c1.Divide(2,2)
    c1.cd(1)
    h_zposAtY0.Fit('gaus')
    h_zposAtY0.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2fmm'%(h_zposAtY0.GetFunction('gaus').GetParameter(2)),0.05,1)
    c1.cd(2)
    h_zposAtY0_2.Fit('gaus')
    h_zposAtY0_2.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2fmm'%(h_zposAtY0_2.GetFunction('gaus').GetParameter(2)),0.05,1)
    c1.cd(3)
    h_ypos_target.Fit('gaus')
    h_ypos_target.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2fmm'%(h_ypos_target.GetFunction('gaus').GetParameter(2)),0.05,1)
    c1.cd(4)
    h_ypos2_target.Fit('gaus')
    h_ypos2_target.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2fmm'%(h_ypos2_target.GetFunction('gaus').GetParameter(2)),0.05,1)

    c1.SaveAs('%s/pos_3hits_%s.png'%(plotDir,saveStr),'png')


    c123 = TCanvas('c123','c123',10,10,700,500)
    c123.Divide(2,2)
    c123.cd(1)
    h_ypos_target_slope.Draw("colz")
    c123.cd(2)
    h_ypos2_target_slope.Draw("colz")
    c123.cd(3)
    h_ypos_target_slope.ProfileX().Draw("");
    c123.cd(4)
    h_ypos2_target_slope.ProfileX().Draw("")

    c123.SaveAs('%s/pos_y_slope_3hits_%s.png'%(plotDir,saveStr),'png')
    


    c111 = TCanvas('c111','c111',10,10,700,500)
    c111.Divide(2,2)
    c111.cd(1)
    h_ypos_target_p.Draw('colz')
    #h_ypos_target_p.ProfileY().Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit #sigma=%.2fmm'%(h_ypos_target.GetFunction('gaus').GetParameter(2)),0.05,1)
    c111.cd(2)
    h_ypos2_target_p.Draw('colz')

    c111.SaveAs('%s/sigma_pos_vs_p_3hits_%s.png'%(plotDir,saveStr),'png')
    
    c1111 = TCanvas('c1111','c1111',10,10,700,500)
    #h_ysigmapos_prj_p []
    gr_ysigmapos_prj_p  = TGraphErrors()
    ib=0
    for b in range(h_ypos_target_p.GetNbinsX()):
      x_b = h_ypos_target_p.GetXaxis().GetBinCenter(b+1)
      h_syp_b = h_ypos_target_p.ProjectionY('%s_prjy'%h_ypos_target_p.GetName(),b+1,b+1) 
      print 'b %d x_b %f N %f' % (b,x_b,h_syp_b.Integral(-1,9999))
      if h_syp_b.Integral(-1,9999) > 30:
        h_syp_b.Fit('gaus')
        #gr_ysigmapos_prj_p.SetPoint(b,x_b,h_syp_b.GetRMS())
        fg = h_syp_b.GetFunction('gaus')
        gr_ysigmapos_prj_p.SetPoint(ib,x_b,fg.GetParameter(2))
        gr_ysigmapos_prj_p.SetPointError(ib,0.,fg.GetParError(2))
        ib+=1
    gr_ysigmapos_prj_p.SetMarkerStyle(20)
    gr_ysigmapos_prj_p.SetMarkerSize(1.2)
    #gr_ysigmapos_prj_p.SetTitle(';Track momentum (GeV);#sigma(YOCA) (#mum)')
    #h_tmp_prj_p = TH2F('h_tmp_prj_p',';Track momentum (GeV);#sigma(YOCA) (#mum)',10,0,2,10,0,460)
    h_tmp_prj_p = TH2F('h_tmp_prj_p',';Track momentum (GeV);#sigma(YOCA) (mm)',10,0,2,10,0,0.46) #1.5)
    h_tmp_prj_p.SetStats(False)
    h_tmp_prj_p.Draw('')
    gr_ysigmapos_prj_p.Draw('P,same')
    
    c1111.SaveAs('%s/sigma_fit_pos_vs_p_3hits_%s.png'%(plotDir,saveStr),'png')


    #ans = raw_input('pause')

    c11 = TCanvas('c11','c11',10,10,1000,500)
    c11.Divide(2,1)
    c11.cd(1)
    h_zvtx.Fit('gaus')
    h_zvtx.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit <m>=%.2f#pm%.2fmm'%(h_zvtx.GetFunction('gaus').GetParameter(1),h_zvtx.GetFunction('gaus').GetParError(1)),0.05,1)
    utils.myText(0.2,0.47,'Fit #sigma=%.2f#pm%.2fmm'%(h_zvtx.GetFunction('gaus').GetParameter(2),h_zvtx.GetFunction('gaus').GetParError(2)),0.05,1)
    c11.cd(2)
    h_yvtx.Fit('gaus')
    h_yvtx.Draw()
    utils.myText(0.2,0.81,'Target z=%.2fmm'%target_z,0.05,1)
    utils.myText(0.2,0.74,'#sigma_{L1}=%.0f#mum'%(hit_l1_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.67,'#sigma_{L2}=%.0f#mum'%(mean_l2_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.6,'#sigma_{L3}=%.0f#mum'%(mean_l3_sigma*1.0e3),0.05,1)
    utils.myText(0.2,0.53,'Fit <m>=%.2f#pm%.2fmm'%(h_yvtx.GetFunction('gaus').GetParameter(1),h_yvtx.GetFunction('gaus').GetParError(1)),0.05,1)
    utils.myText(0.2,0.47,'Fit #sigma=%.2f#pm%.2fmm'%(h_yvtx.GetFunction('gaus').GetParameter(2),h_yvtx.GetFunction('gaus').GetParError(2)),0.05,1)

    c11.SaveAs('%s/vtxpos_3hits_%s.png'%(plotDir,saveStr),'png')

    c4 = TCanvas('c4','c4',10,10,700,500)
    h_chi2.Draw()
    h_chi2_2.SetLineColor(3)
    h_chi2_2.Draw('same')
    utils.myText(0.2,0.81,'<#Chi^{2}>=%.2f (top)'%h_chi2.GetMean(),0.05,1)
    utils.myText(0.2,0.73,'<#Chi^{2}>=%.2f (bot)'%h_chi2_2.GetMean(),0.05,1)
    c4.SaveAs('%s/chi2_%s.png'%(plotDir,saveStr),'png')

    c_invm = TCanvas('c_invm','c_invm',10,10,700,500)
    h_invMass_truth.Draw()
    h_invMass.SetLineColor(2)
    h_invMass.SetMarkerColor(2)
    h_invMass.SetMarkerStyle(20)
    h_invMass.SetMarkerSize(1.2)
    h_invMass.Draw('E,same')
    leg_inv = TLegend(0.6,0.6,0.9,0.9)
    leg_inv.SetFillColor(0)
    leg_inv.AddEntry(h_invMass_truth,'Generated','L')
    leg_inv.AddEntry(h_invMass,'Fitted','P')
    leg_inv.Draw()
    c_invm.SaveAs('%s/invM_3hits_%s.png'%(plotDir,saveStr),'png')

    ans = raw_input('pause')
