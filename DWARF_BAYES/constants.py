import math
import csv

NE=31
Estart=100.0
Eincrement=0.15
# NE=8
# Estart=80.0
# Eincrement=2.0
Erange=[]
for i in range(NE):
    if (i==0):
        Erange.append(Estart)
    else:
        Erange.append(Erange[i-1]+Eincrement*Erange[i-1])

print "Erange from constants.py: ", Erange

NM=30
Mstart=115.0
Mincrement=0.15
Mlist=[]
for i in range(NM):
    if (i==0):
        Mlist.append(Mstart)
    else:
        Mlist.append(Mlist[i-1]+Mlist[i-1]*Mincrement)


#(Former HACK) Sadly, this number is necessarily (not true anymore after introducing "big numbers" class -- bvv).
#It is canceled, but if I had to write this software again it wouldn't
#be necessary. It becomes necessary to set this number arbitrarily high
#when calculating lots of data due to the smallness of the P_Npi_S
#probabilities.
#P_Npi_S_factor=100
#P_Npi_S_factor=15
P_Npi_S_factor=10.0

#Slist=[0,1,2,3,4,5,6,7,8,9,10]
#Slist=[10**(-26+(26-19)*n/100) for n in range(101)]
#Slist=[math.pow(10,(-29+(29-21)*(1.0*n)/100)) for n in range(101)]
#Slist=[math.pow(10,(-30+(30-17)*(1.0*n)/100)) for n in range(101)]

#use this:
#Slist=[math.pow(10,(-40+(40-19)*(1.0*n)/100)) for n in range(101)]
NS=100

# TODO: Introduce some name for the lower boundary:
f=lambda n : math.pow(10,(-26+(26-19)*(1.0*n)/NS))
Slist=[f(n) for n in range(NS+1)]
Slist_sum=math.fsum(Slist) #bvv: why would we need this?

#or...
#start=10e-26
#finish=10e-22
#increment=(finish-start)/100
#Slist=[start+i*increment for i in range(100)]

#TODO: Is this right?
Slist_diff=[f(n+1)-f(n) for n in range(NS+1)]
Slist_integration_factors={}
for n in range(NS+1):
    #Slist_integration_factors[f(n)]=Slist_diff[n]/math.fsum(Slist_diff)
    Slist_integration_factors[f(n)]=Slist_diff[n]

def get_event_line(filename, separator):
    header_count = 2
    myreader = csv.reader(open(filename), delimiter=separator, skipinitialspace=True)
    for i in range(0, header_count):
        myreader.next()
    for line in myreader:
        if line:
            yield line


#J-factors... units are GeV^2cm^-5
# All taken from http://veritash.sao.arizona.edu:8081/DM-AsPEN-SWG/150714_032238/DMAspenSlidesForIACTPlenary_v1.pdf
# except for WilmanI.

J = {
    # 10^19.4 => 2.51188643151e19
    "Segue": 2.51e19, \

    # 10^18.8 => 6.3095734448e18
    "Draco": 6.31e18, \

    #10^18.9 => 7.94328234724e18
    "UrsaMinor": 7.94e18, \

    # 10^18.2 => 1.58489319246e18
    "Bootes1": 1.58e18, \

    # It was just here before, most likely written by Ryan.
    "WilmanI": 8.4e18
}

ALEX_J_CONV_PSF_PATHS = {"d": "Alex/intJPSFs.txt", "f": "Alex/intJPSFs_pass5f_chainrow1039.txt"} 
BAYES_SRC_NAME_FROM_ALEX = {"Energy": "Energy", "BooI": "Bootes1", "Seg1": "Segue", "Draco": "Draco", "UMi": "UrsaMinor"}

# ALEX_JS_CONV_PSF =
# {"d": {src: [J1, J2, ... , JN]},
#  "f": {src: [J1, J2, ... , JN]}}
# where JN - number of energy bins in intJPSF.dat
# JN is supposed to be equal to NE.
ALEX_JS_CONV_PSF = {} 
for k, v in ALEX_J_CONV_PSF_PATHS.iteritems():
    with open(v) as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        reader.next()
        header = reader.next()
        header = ["Energy"] + header
        ALEX_JS_CONV_PSF[k] = {}
        for h in header:
            ALEX_JS_CONV_PSF[k][BAYES_SRC_NAME_FROM_ALEX[h]] = []

        for row in reader:
            for h, v in zip(header, row):
                ALEX_JS_CONV_PSF[k][BAYES_SRC_NAME_FROM_ALEX[h]].append(float(v))
    # This is not present in intJPSFs.txt:
    ALEX_JS_CONV_PSF[k]["WilmanI"] = 31 * [J["WilmanI"]]
