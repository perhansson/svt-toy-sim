from ROOT import TGraphErrors,TLatex,TRandom
import math


rnd = TRandom()
misalignments = {0:0.,1:0.,2:0.}


class Hit:
    def __init__(self):
        self.y = 0.0
        self.z = 0.0
        self.sigma_y = 0.0
	self.hit_res = 0.0
    def __init__(self,z,y,sigma_y,hit_res):
        self.y = y
        self.z = z
        self.sigma_y = sigma_y
	self.hit_res = hit_res
    def getWiggledHit(self,misalign_y):
        if misalign_y != 0.0:
            print 'hit at %.3f: applying misalign %.3f' % (self.z,misalign_y)
        newhit = Hit(self.z,self.y+misalign_y+rnd.Gaus(0.,self.hit_res),self.sigma_y,self.hit_res)
        return newhit
    def getTotError(self):
        return math.sqrt(self.hit_res*self.hit_res+self.sigma_y*self.sigma_y)


class Track:
    def __init__(self):
        self.hits = []
    def printProp(self):
        print 'Track with hits at: '
        for h in self.hits: 
            print '(%f,%f+-%f)' % (h.z,h.y,h.sigma_y)
    def addHit(self,z,y,sigma_y,hit_res):
        newhit = Hit(z,y,sigma_y,hit_res)
        exists = False
        for h in self.hits:
            if h.y==newhit.y or h.z==newhit.z:
                print 'Error: hit exists'
                sys.exit(1)
                exists=True
        if not exists: 
            self.hits.append(newhit)
    def getWiggledTrack(self,applyMisAlignment):
        trk = Track()
        for layer in range(len(self.hits)):
            h = self.hits[layer]
            if applyMisAlignment:
                misalign = misalignments[layer+1]
            else:
                misalign=0.
            newhit = h.getWiggledHit(misalign)
            trk.addHit(newhit.z,newhit.y,newhit.sigma_y,newhit.hit_res)
        return trk
    def getGraph(self):
        gr = TGraphErrors()
        ihit=0
        for h in self.hits:
            gr.SetPoint(ihit,h.z,h.y)
            gr.SetPointError(ihit,0.,h.getTotError())
            ihit=ihit+1
        return gr

def fitTrack(track,fnctype):
    #print 'fit track'
    #use graph for fit
    gr = TGraphErrors()
    ihit=0
    for h in track.hits:
        gr.SetPoint(ihit,h.z,h.y)
        gr.SetPointError(ihit,0.,h.getTotError())
        ihit=ihit+1
        #print ' set hit %d at %f,%f' % (ihit,h.z,h.y)
    #print ' %d hits' % gr.GetN()
    gr.Fit(fnctype,'Q0')
    return gr

class Fit:
    def __init__(self,track,fnctype):
        self.track = track
        self.fnctype = fnctype
        self.gr = fitTrack(self.track,self.fnctype)
    def printProp(self):
        f = self.gr.GetFunction(self.fnctype)
        #print 'Track fit: k=%f m=%f ' % (f.GetParameter(0),f.GetParameter(1))
        self.track.printProp()
    def getSlope(self):
        return self.gr.GetFunction(self.fnctype).GetParameter(1)
    def getYAtZ(self,z):
        return self.gr.GetFunction(self.fnctype).Eval(z)
    def getZAtY(self,y):
        return self.gr.GetFunction(self.fnctype).GetX(y)
    def SetRange(self,xmin,xmax):
        self.gr.GetFunction(self.fnctype).SetRange(xmin,xmax)
    def getChi2(self):
        return self.gr.GetFunction(self.fnctype).GetChisquare()
    def getNDF(self):
        return self.gr.GetFunction(self.fnctype).GetNDF()
    def mirrorX(self,start_pos):
        slope = self.gr.GetFunction(self.fnctype).GetParameter(1)
        m = self.gr.GetFunction(self.fnctype).GetParameter(0)
        #m_new = m - 2.0*slope*(-start_pos)
        m_new = -1.0*m
        self.getFunction().SetParameter(0,m_new)
        self.getFunction().SetParameter(1,slope*-1)
    def getFunction(self):
        return self.gr.GetFunction(self.fnctype)

def getSimpleSigmaZ(lt,l2,s1,h1,h2):
    sig_z = ((lt+l2)*(h1+s1)-h2*lt)/(h2-(h1+s1))
    return sig_z

def myText(x,y,text, tsize,color):
    l = TLatex()
    l.SetTextSize(tsize); 
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x,y,text);

def getPol1Crossing(fit1,fit2):
    if fit1.fnctype!='pol1':
        print 'ERROR: only works for pol1!!!'
        sys.exit(1)
    f1 = fit1.gr.GetFunction(fit1.fnctype)
    f2 = fit2.gr.GetFunction(fit2.fnctype)
    if (f2.GetParameter(0)-f1.GetParameter(0))==0:
        print 'WARNING pol1 crossing had the same slope'
        return 0.
    k1 = f1.GetParameter(1)
    k2 = f2.GetParameter(1)
    m1 = f1.GetParameter(0)
    m2 = f2.GetParameter(0)
    x = (m1-m2)/(k2-k1)
    #print 'k1=%f m1=%f k2=%f m2=%f => x=%f' % (k1,m1,k2,m2,x)
    return x
