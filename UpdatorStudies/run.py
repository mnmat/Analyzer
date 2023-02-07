import os
energies = ["10","50","20","10"]
etas = ["16"]#,"21","25","29"]


for i,eta in enumerate(etas):
    for en in energies:
        cmd = "cmsRun python/ConfigFile_cfg.py %s %s" %(eta,en)
        os.system(cmd)