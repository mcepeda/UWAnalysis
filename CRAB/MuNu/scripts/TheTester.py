'''
Tests if all the variable names are valid
'''

import histoRange as hr 
import TheParameters as p
import TheCuts as ct

lumi,bNr,btype,jNr,njetcut,jetcut,I,F,iso_value,antiIso_value,path,extraName,leafs,drawW,drawZ,drawQCD,drawData,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,eventTreeLocation,extraCut,Ctl_Andrea = p.arams() 

CutsMCn, CutsMCnW, CutsMCi,CutsDatan,CutsDatai,CutsMCnwl,CutsMCiwl,CutsMCnwc,CutsMCiwc,CutsMCnwcc,CutsMCiwcc,CutsMCnwbb,CutsMCiwbb,CutsMCnT,CutsMCiT,CutsMCnTup,CutsMCiTup,CutsMCnTdn,CutsMCiTdn, cutMcNonIso_Bup, cutMcNonIsoW_Bup, cutMcIso_Bup, cutMcWlNonIso_Bup, cutMcWlIso_Bup, cutMcWcNonIso_Bup, cutMcWcIso_Bup, cutMcWccNonIso_Bup, cutMcWccIso_Bup, cutMcWbbNonIso_Bup, cutMcWbbIso_Bup, cutMcNonIso_Bdn, cutMcNonIsoW_Bdn, cutMcIso_Bdn, cutMcWlNonIso_Bdn, cutMcWlIso_Bdn, cutMcWcNonIso_Bdn, cutMcWcIso_Bdn, cutMcWccNonIso_Bdn, cutMcWccIso_Bdn, cutMcWbbNonIso_Bdn, cutMcWbbIso_Bdn = ct.cutmaker(
leafs[0],iso_value,antiIso_value,lumi,bNr,btype,jNr,njetcut,jetcut,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,extraCut,Ctl_Andrea
)


for leaf in leafs:
 steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)
 print(leaf)
 print('%s %s'%(xtitle,xunits))
 print('%s(%s,%s)'%(steps,xmin,xmax))
 print('Log Scale: %s'%setLogY)
 print('')

print("Data Cuts: \n %s\n"%(CutsDatai))
print("MC Cuts: \n %s\n"%(CutsMCi))
print("Bup Cuts: \n %s\n"%(cutMcIso_Bup))
print("Bdn Cuts: \n %s\n"%(cutMcIso_Bdn))
