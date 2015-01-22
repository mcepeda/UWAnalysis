'''
Look Up Table for plotting info
Author: T.M.Perry UW
'''
import fnmatch as fnm

def ranger(leaf): 
 steps = None
 xmin = None
 xmax = None
 setLogY = False
 if fnm.fnmatch(leaf,'*_eta*'):
   steps = 20
   xmin = -3.
   xmax = 5.
   xtitle = leaf
   xunits = None

 elif fnm.fnmatch(leaf,'*_phi*'):
  # steps = 24
  # xmin = -3.14159
  # xmax = 3.14159
   #steps = 31
   #xmin = -3.45
   #xmax = 5.85
   steps = 35
   xmin = -3.4033
   xmax = 5.7594
   xtitle = leaf
   xunits = None

 elif leaf=='Wpt' or fnm.fnmatch(leaf,'J*_pt'):
   steps = 12
   xmin = 25
   xmax = 205.
   xtitle = leaf
   xunits = 'GeV'

 elif leaf=='muon_pt':
   steps = 50
   xmin = 0.
   xmax = 200.
   xtitle = 'p_{T}^{#mu}'
   xunits = 'GeV'

 elif leaf=='ptEle':
   steps = 50
   xmin = 0.
   xmax = 200.
   xtitle = 'p_{T}^{e}'
   xunits = 'GeV'

 elif leaf=='muon_charge':
   steps = 5
   xmin = -2.
   xmax = 2.
   xtitle = leaf
   xunits = 'e'

 elif fnm.fnmatch(leaf,'nJets*'):
   steps = 8
   xmin = 2.
   xmax = 10.
   xtitle = leaf
   xunits = None
   setLogY = True

 elif fnm.fnmatch(leaf,'mt*'):
   # Andrea
   #steps = 20
   #xmin = 45.
   #xmax = 205.
   steps = 16 
   #steps = 24 
   xmin = 0.
   xmax = 240.
   xtitle = 'Transverse Mass'
   xunits = 'GeV'

 elif  leaf == 'met_pt':
   steps = 24
   xmin = 0.
   xmax = 240.
   xtitle = 'MET'
   xunits = 'GeV'

 elif leaf == 'vertices':
   steps = 40
   xmin = 0.
   xmax = 40.
   xtitle = 'Nr. Vertices'
   xunits = None

 elif leaf == 'DiMuonPt':
   steps = 100
   xmin = 0.
   xmax = 200.
   xtitle = 'p_{T}^{#mu#mu}'
   xunits = 'GeV'

 elif leaf == 'DiMuonMass':
   steps = 100
   xmin = 0.
   xmax = 200.
   xtitle = 'm^{#mu#mu}'
   xunits = 'GeV'

 elif leaf == 'ht':
   steps = 20 
   #steps = 25 
   #steps = 50
   xmin = 50.
   xmax = 550.
   #xmax = 400.
   xtitle = 'h_{T}'
   xunits = 'GeV'

 elif leaf == 'mJJ':
   steps = 25
   xmin = 0.
   xmax = 500.
   xtitle = 'm^{JJ}'
   xunits = 'GeV'

 elif fnm.fnmatch(leaf,'mJ3J4'):
   steps = 15
   xmin = 20.
   xmax = 170.
   xtitle = 'm^{J3J4}'
   xunits = 'GeV'

 elif leaf == 'ptJJ':
   steps = 50
   xmin = 0.
   xmax = 300.
   xtitle = 'p_{T}^{JJ}'
   xunits = 'GeV'

 elif leaf == 'LHEProduct':
   steps = 7
   xmin = 4.
   xmax = 11.
   xtitle = leaf
   xunits = None

 elif fnm.fnmatch(leaf,"J*fullId*Pu"):
   steps = 21
   xmin = 0.
   xmax = 1.1
   xtitle = leaf
   xunits = None
   setLogY = True

 elif fnm.fnmatch(leaf,"J*btag") or fnm.fnmatch(leaf,"*Beta*"):
   steps = 20
   xmin = 0.
   xmax = 1.
   xtitle = leaf
   xunits = None
   setLogY = True

 elif leaf=='J1DR'or leaf=='J2DR':
   steps = 80
   xmin = 0.3
   xmax = 0.7
   xtitle = leaf
   xunits = None
   setLogY = True

 elif fnm.fnmatch(leaf,'J*JetParton'):
   steps = 14
   xmin = -7.
   xmax = 7.
   xtitle = leaf
   xunits = 'flavor'

 elif fnm.fnmatch(leaf,'(J1_mass_SSV*+J2_mass_SSV*)'):
   steps = 45
   xmin = 0.
   xmax = 15.
   xtitle = leaf
   xunits = 'GeV'

 elif fnm.fnmatch(leaf,'J*_pt_SSV*'):
   steps = 25
   xmin = 0.
   xmax = 100.
   xtitle = leaf
   xunits = 'GeV'

 elif fnm.fnmatch(leaf,'J*_mass_*SV*'):
   steps = 32
   xmin = 0.
   xmax = 8.
   xtitle = leaf
   xunits = 'GeV'

 elif leaf == '(((J1_eta-J2_eta)*(J1_eta-J2_eta))+((J1_phi-J2_phi)*(J1_phi-J2_phi)))':
   steps = 20
   xmin = 0.
   xmax = 5.
   xtitle = '#Delta #phi^{2}_{j1j2}+#Delta #eta^{2}_{j1j2}'
   xunits = None

 elif leaf=='dz':
   steps = 40
   xmin = -0.2
   xmax = 0.2
   xtitle = leaf
   xunits = None

 elif fnm.fnmatch(leaf,'J*_nTracksSSV') or fnm.fnmatch(leaf,'(J*_nTracksSSV)'):
   steps = 30
   xmin = 0.
   xmax = 30
   xtitle = leaf
   xunits = None

 elif fnm.fnmatch(leaf,'J*_pt_PV') or fnm.fnmatch(leaf,'J*_mass_PV'):
   steps = 50
   xmin = 0.
   xmax = 400.
   xtitle = leaf
   xunits = 'GeV'

 else :
   print("\n\n")
   print "YOU ARE A NINNY! The leaf's name isn't valid"
   print("\n\n")
   steps = 20
   xmin = 0.
   xmax = 1.
   xtitle = leaf
   xunits = None

# print leaf
# print steps
# print xmin
# print xmax

 return steps,xmin,xmax,xtitle,xunits,setLogY
