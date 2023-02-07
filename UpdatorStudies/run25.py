import os
energies = ["100","50","20","10"]
etas = ["25"]


for i,eta in enumerate(etas):
    for en in energies:
        cmd = "cmsRun python/ConfigFile_cfg.py %s %s" %(eta,en)
        os.system(cmd)