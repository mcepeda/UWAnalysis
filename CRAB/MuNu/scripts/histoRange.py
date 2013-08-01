'''
Look Up Table for plotting info
Author: T.M.Perry UW
'''

def ranger(leaf): 
 steps = None
 xmin = None
 xmax = None
 setLogY = False
 if leaf=='muonPt':
   steps = 50
   xmin = 0.
   xmax = 200.
   xtitle = 'p_{T}^{#mu}'
   xunits = 'GeV'
 elif leaf=='Mt':
   steps =50 
   xmin = 0.
   xmax = 200.
   xtitle = 'Transverse Mass'
   xunits = 'GeV'
 elif leaf == 'MET' or leaf == 'met':
   steps = 50
   xmin = 0.
   xmax = 200.
   xtitle = 'MET'
   xunits = 'GeV'
 elif leaf == 'vertices':
   steps = 60
   xmin = 0.
   xmax = 40.
   xtitle = 'Nr. Vertices'
   xunits = 'xx'
 elif leaf == 'DiMuonMass':
   steps = 100
   xmin = 0.
   xmax = 200.
   xtitle = 'm^{#mu#mu}'
   xunits = 'GeV'
 elif leaf == 'muonEta':
   steps = 20
   xmin = -3.
   xmax = 5.
   xtitle = '#eta^{#mu}'
   xunits = 'xx'
 elif leaf == 'muonPhi':
   steps = 20
   xmin = -4.
   xmax = 6.
   xtitle = '#phi^{#mu}'
   xunits = 'xx'
 elif leaf == 'ht':
   steps = 50
   xmin = 50.
   xmax = 400.
   xtitle = 'h_{T}'
   xunits = 'GeV'
 elif leaf == 'highestJetEta' or leaf=='secondJetEta' or leaf == 'thirdJetEta':
   steps = 20
   xmin = -3.
   xmax = 5.
   xtitle = leaf
   xunits = 'xx'
 elif leaf == 'highestJetPhi' or leaf=='secondJetPhi' or leaf=='thirdJetPhi':
   steps = 20
   xmin = -4.
   xmax = 6.
   xtitle = leaf
   xunits = 'xx'
 elif leaf == 'highestJetPt' or leaf=='secondJetPt' or leaf=='thirdJetPt':
   steps = 50
   xmin = 0.
   xmax = 200.
   xtitle = leaf
   xunits = 'GeV'
 elif leaf == 'mJJ':
   steps = 100
   xmin = 1.
   xmax = 300.
   xtitle = 'm^{JJ}'
   xunits = 'GeV'
 elif leaf == 'mJ3J4':
   steps = 100
   xmin = 1.
   xmax = 300.
   xtitle = 'm^{J_{3}J_{4}}'
   xunits = 'GeV'
 elif leaf == 'ptJJ':
   steps = 100
   xmin = 0.
   xmax = 300.
   xtitle = 'p_{T}^{JJ}'
   xunits = 'GeV'
 elif leaf == 'WPt':
   steps = 100
   xmin = 0.
   xmax = 300.
   xtitle = 'p_{T}^{W}'
   xunits = 'GeV'
 elif leaf == 'LHEProduct':
   steps = 7
   xmin = 4.
   xmax = 11.
   xtitle = leaf
   xunits = 'xx'
 elif leaf=='J1CSVbtag'or leaf=='J1CSVMVAbtag'or leaf=='J2CSVbtag'or leaf=='J2CSVMVAbtag':
   steps = 20
   xmin = 0.
   xmax = 1.
   xtitle = leaf
   xunits = 'xx'
   setLogY = True
 elif leaf=='J1DR'or leaf=='J2DR':
   steps = 20
   xmin = 0.
   xmax = 1.
   xtitle = leaf
   xunits = 'xx'
 elif leaf=='J1JetParton' or leaf =='J2JetParton':
   steps = 14
   xmin = -7.
   xmax = 7.
   xtitle = leaf
   xunits = 'flavor'
 elif leaf=='J1SVMassb' or leaf =='J2SVMassb':
   steps = 10
   xmin = 0.
   xmax = 15.
   xtitle = leaf
   xunits = 'GeV'
 elif leaf=='muonCharge':
   steps = 5
   xmin = -2.
   xmax = 2.
   xtitle = leaf
   xunits = 'e'
 elif leaf == '(((highestJetEta-secondJetEta)*(highestJetEta-secondJetEta))+((highestJetPhi-secondJetPhi)*(highestJetPhi-secondJetPhi)))':
   steps = 20
   xmin = 0.
   xmax = 2.
   xtitle = '#Delta #phi^{2}+#Delta #eta^{2}'
   xunits = 'xx'
 else :
   print("\n\n")
   print "YOU ARE A NINNY! The leaf's name isn't valid"
   print("\n\n")
   steps = 20
   xmin = 0.
   xmax = 1.
   xtitle = leaf
   xunits = 'xx'

# print leaf
# print steps
# print xmin
# print xmax

 return steps,xmin,xmax,xtitle,xunits,setLogY
