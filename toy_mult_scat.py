#!/usr/bin/python

from math import sqrt,log
from ROOT import TCanvas,TGraph
from utils import myText

def MSAngle(x,z,beta,p):
    theta = 13.6/(beta*p)*z*sqrt(x)*(1.0+0.038*log(x))
    return theta

me = 0.511 #MeV/c2
X0 = 9.370 #cm for Si
z = 1 #electrons
t = 2.0*320.0e-4 #cm
single_hit_res = 0.006 #mm

if __name__=='__main__':
    print 'Do MS calc'
    
    
    gr_th_p = TGraph()
    
    for i in range(1,20):
        p = float(i)
        p = p/10.0*1.0e3 #MeV
        E = sqrt(p*p+me*me)
        beta = p/E
        x = t/X0
        theta = MSAngle(x,z,beta,p)
        gr_th_p.SetPoint(i-1,p,theta)
    
    c1 = TCanvas('c1','c1',10,10,700,500)
    c1.Divide(2,1)
    c1.cd(1)
    gr_th_p.SetTitle('Multiple scattering angle;Momentum (MeV);#theta_0')
    gr_th_p.Draw('AXL')
    myText(0.2,0.81,'t=%.0f#mum,X_{0}=%.2fcm'%(t*1.0e4,X0),0.05,1)


    layers = {1:100.,2:200.,3:300.,4:500.,5:700.}
    mserr = {}
    for l,d in layers.iteritems():
        mserr[d] = 0.
    
    p = 1200.0 #MeV
    E = sqrt(p*p+me*me)
    beta = p/E
    t = 320.0e-4 #cm
    x = t/X0
    theta = MSAngle(x,z,beta,p)
    

    ilayer=0
    for l,d in layers.iteritems():
        print 'Layer %d' % l
        # uncorrelated: sum over previous hit uncertainties
        for ll,dd in layers.iteritems():
            if dd<d: 
                dl = d-dd #distance to the current layer
                delta = theta*dl
                mserr[d] = mserr[d] + delta*delta
                print 'Adding error (%f from theta=%f,dl=%f) from layer %d to hit error ' % (delta,theta,dl,ll)
    
    gr_mserr_layer = TGraph()
    i=0
    for d in mserr.keys():
        mserr[d] = sqrt(mserr[d])
        gr_mserr_layer.SetPoint(i,d,mserr[d])
        i=i+1

    c1.cd(2)
    gr_mserr_layer.SetTitle('Multiple scattering uncertainty;Layer position;Uncertainty (mm)')
    gr_mserr_layer.SetMarkerStyle(20)
    gr_mserr_layer.Draw('AXP')
    myText(0.2,0.81,'p=%.0fMeV'%p,0.05,1)

    c1.SaveAs('svt_toy_ms.png','png')
    ans = raw_input('pause')
