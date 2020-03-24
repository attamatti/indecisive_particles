#!/usr/bin/env python
import sys
import glob
import matplotlib.pyplot as plt
import itertools
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile_new(f):
    inhead = True
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    count = 0
    labcount = 0
    for i in alldata:
        if '_rln' in i:
            labelsdic[i.split()[0]] = labcount
            labcount +=1
        if inhead == True:
            header.append(i.strip("\n"))
            if '_rln' in i and '#' in i and  '_rln' not in alldata[count+1] and '#' not in alldata[count+1]:
                inhead = False
        elif len(i.split())>=1:
            data.append(i.split())
        count +=1
    
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#

errmsg = 'USAGE: indecisive_parts.py <class3d files dir> <start iteration> <% switches allowed> <particles file>'
### find the star files

try:
    startiter = int(sys.argv[2])
    allowed = float(sys.argv[3])
    if allowed > 1.0:
        allowed = allowed/100.0
    finlabels,finheader,findata = read_starfile_new(sys.argv[4])
    extra = False
    if '--extra' in sys.argv:
            extra = True
except:
    sys.exit(errmsg)
    
all_dat_files  = glob.glob('{0}*_data.star'.format(sys.argv[1]))
datafiles = {}          ##{iteration number: filename}
for i in all_dat_files:
    fname = i.split('_')
    for j in fname:
        if 'it0' in j:
            itno = int(j.replace('it',''))
            datafiles[itno] = i

iters = list(datafiles)
iters.sort()  

def get_classes_for_parts(starfile,partsdic,classes):
    labels,header,data = read_starfile_new(starfile)
    print(starfile)
    for i in data:
        try:
            partsdic[i[labels['_rlnImageName']]].append(int(i[labels['_rlnClassNumber']]))
        except:
            partsdic[i[labels['_rlnImageName']]] = [int(i[labels['_rlnClassNumber']])]
        if i[labels['_rlnClassNumber']] not in classes:
            classes.append(i[labels['_rlnClassNumber']])
    return(partsdic,classes)
    
        

### read intial file and get the parts and classes
partsdic = {}          #{partname:[class iter 1, class iter 2, ..., class iter n]}
classes = []
for i in range(startiter,len(datafiles)):
    partsdic,classes = get_classes_for_parts(datafiles[i],partsdic,classes)
classes.sort()
classes = [int(x) for x in classes]

### add in the _rlnLogLikeliContribution _rlnMaxValueProbDistribution from final iteration
labels,header,data = read_starfile_new(datafiles[len(datafiles)-1])
for i in data:
    partsdic[i[labels['_rlnImageName']]] = [partsdic[i[labels['_rlnImageName']]],i[labels['_rlnLogLikeliContribution']],i[labels['_rlnMaxValueProbDistribution']]]

### track each particles's switches
movescount = {}         ##{particle: # of class changes movements}
partsmoves = {}         ##{particle: # of moves}
for i in partsdic:
    partsmoves[i] = 0
    movescount[i] = 0
    n= 1
    for it in partsdic[i][0]:
        try:
            if it != partsdic[i][0][1]:
                movescount[i] += 1
                partsmoves[i] += 1
        except:
            pass

moves = [movescount[x] for x in movescount]
plt.hist(moves,range(len(datafiles)))
plt.savefig('total_switches.png')
plt.close()

if extra == True:

    #### relationships bewteen classes
    def get_class_relations(classno,partsdic):
        thisclasscounts = []
        for i in partsdic:
            if partsdic[i][0][-1] == classno:
                for sw in partsdic[i][0][:-1]:
                    thisclasscounts.append(sw)
        return(thisclasscounts)
    
    def get_defined_switches(classno,partsdic):
        defined_switches = {}           ### {(from,to):count}
        possible_switches = itertools.permutations(classes,2)
        for i in possible_switches:
            defined_switches[i] = 0
        for i in partsdic:
            if partsdic[i][0][-1] == classno:
                n=1
                for sw in partsdic[i][0]:
                    try:
                        switch = (sw,partsdic[i][0][n])
                        n+=1
                        defined_switches[switch] +=1
                    except:
                        pass
        return(defined_switches)
    
    ### total time in any given class
    total_time = {}        #{class:[locations durning previous iters]}
    for i in classes:
        print('class number {0}'.format(i))
        total_time[i] = get_class_relations(i,partsdic)
        for j in classes:
            timeinclass = total_time[i].count(j)/float(len(total_time[i]))
            print(j,round(timeinclass,3))
    
    #### switches
    defined_switches = {}       ## {classno:{switch:count}}
    for i in classes:
        defined_switches[i] = get_defined_switches(i,partsdic)
    dswitchlist = list(defined_switches)
    dswitchlist.sort()
    for i in dswitchlist:
        print(i,defined_switches[i])
        
    ### analysis of defined switches not in or out of final class
    ### not an escepcially intereting analysis right now
    non_fin_switches = {}           ### {(from,to):count} not counting when a particle moved in or of its final class
    possible_switches = itertools.permutations(classes,2)
    for i in possible_switches:
        non_fin_switches[i] = 0
    for i in defined_switches:
        for j in defined_switches[i]:
            if i not in j:
                non_fin_switches[j] += defined_switches[i][j]
    for i in non_fin_switches:
        print(i,non_fin_switches[i])

#### plot the number of switches vs maxvalprobdist and loglikelicontrib
nswitches_dic = {}      ## {number of switches:[all LLCs],[all MVPDS]}
for i in range(len(datafiles)-startiter):
    nswitches_dic[i] = [[],[]]
for i in partsmoves:
    nswitches_dic[int(partsmoves[i])][0].append(float(partsdic[i][1]))
    nswitches_dic[int(partsmoves[i])][1].append(float(partsdic[i][2]))

f,(ax1,ax2) = plt.subplots(1,2)
for i in range(len(datafiles)-startiter):    
    meanllc = np.mean(nswitches_dic[i][0])
    stdllc = np.std(nswitches_dic[i][0])
    meanmvpd = np.mean(nswitches_dic[i][1])
    stdmvpd = np.std(nswitches_dic[i][1])
    print('')
    print('number of switches',i,'{0}% of iterations'.format(int((i/float(len(datafiles)-startiter))*100.0)))
    print('number of parts',len(nswitches_dic[i][0]))
    print('mean,std LLC',meanllc,stdllc)
    print('mean,std MVPD',meanmvpd,stdmvpd)
    ax1.scatter(i,meanmvpd,s=50)
    ax1.plot([i,i],[meanmvpd+stdmvpd,meanmvpd-stdmvpd])
    ax2.scatter(i,meanllc,s=50)
    ax2.plot([i,i],[meanllc+stdllc,meanllc-stdllc])
plt.savefig('LLC_MVPD.png')

##### filter the particles
count = 0
output = open('decisive_particles.star','w')
for i in finheader:
    output.write('{0}\n'.format(i))
print('\nfiltering particles over the switch threshold')
for i in findata:
    pct = movescount[i[finlabels['_rlnImageName']]]/float(len(datafiles)-startiter)
    if pct <= allowed:
        output.write('{0}\n'.format('   '.join(i)))
        count+=1
print('wrote decisive_particles.star with {0}/{1} particles'.format(count,len(data)))