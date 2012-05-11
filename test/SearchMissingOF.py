#!/usr/bin/env python
import re
import sys
import os
import ConfigParser
import optparse
import subprocess

ListStringdCACHE=[]
#ListStringOut=[]

ListStringZahldCACHE=[]
ListStringZahlOut=[]

ListdCacheMustDelete=[]


usage = "usage: %prog [options]"
optParser = optparse.OptionParser(usage)

optParser.add_option("-u", "--user", help=r"user  (default = %s)" %(""), dest="user", default="")
optParser.add_option("-t", "--tag", help="which tag do you want",dest="tag", default="")
optParser.add_option("-v", "--version", help="Skimmer version",dest="version", default="")
optParser.add_option("-b", "--basedir", help=r"Basedir (default = %s)" %(os.getcwd()),dest="basedir", default=os.getcwd())
(options, args) = optParser.parse_args()
if (options.user=="" or options.tag=="" or options.version==""):
    optParser.error("Can't handle no tag or user or version at the momment  try -h")
    sys.exit()

cmndlst = [ "uberftp","grid-ftp.physik.rwth-aachen.de",r"'active; ls /pnfs/physik.rwth-aachen.de/cms/store/user/%s/output/%s/%s'"%(options.user,options.version,options.tag)]
p = subprocess.Popen( cmndlst,stdout=subprocess.PIPE )
StringdCACHE, err = p.communicate()

ListStringdCACHE = StringdCACHE.split("\n")
ListStringdCACHE=ListStringdCACHE[2:-1]
#print(StringdCACHE)

resDir=os.listdir(options.basedir+r"/CRAB-%s-%s/"%(options.tag,options.version))

for i in resDir:
    if "crab_" in i:
        ListStringOut=os.listdir(options.basedir+r"/CRAB-%s-%s/%s/res/"%(options.tag,options.version,i))


for i in ListStringdCACHE:

	dC = re.findall(r".*out_([^_]*)_.*",i)

	if(len(dC)!=0):
		ListStringZahldCACHE.append(dC[0])

for i in ListStringOut:

	Ou = re.findall(r".*CMSSW_(.*)\.stdo.*",i)
	
	if(len(Ou)!=0):
		ListStringZahlOut.append(Ou[0])


#print(ListStringZahlOut)
#print(ListStringZahldCACHE)

pos=0

print("crab -resubmit "),
for i in ListStringZahldCACHE:

	if((i not in ListStringZahlOut)):
		print(i+","),
		#print(ListStringdCACHE[pos])

		ListdCacheMustDelete.append(ListStringdCACHE[pos].split()[-1])
		#print("rm "+ListdCacheMustDelete[-1])

	pos=pos+1

print("")

for i in ListdCacheMustDelete:

	print(r"uberftp grid-ftp.physik.rwth-aachen.de 'active; rm /pnfs/physik.rwth-aachen.de/cms/store/user/%s/output/%s/%s/"%(options.user,options.version,options.tag)+i+"'")

print("\n")

if (len(ListdCacheMustDelete) == 0):
	print("No root file without stdout file")
else:
	print(r"There are %i root files with no stdout file"%(len(ListdCacheMustDelete)))	
