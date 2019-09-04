import numpy
from matplotlib import pyplot as pyplot
from matplotlib import animation
import csv
from matplotlib import gridspec
from matplotlib.lines import Line2D
import matplotlib.text as text

# this plots the paper wave plate

ExFieldFSFile = 'Documents/SSProject/Debug/a4mmDATAscrappy/ExPaperRef.csv'
EyFieldFSFile = 'Documents/SSProject/Debug/a4mmDATAscrappy/EyPaperRef.csv'
EzFieldFSFile = 'Documents/SSProject/Debug/a4mmDATAscrappy/EzPaperRef.csv'

ExFieldPAFile = 'Documents/SSProject/Debug/a4mmDATAscrappy/ExPaperSample.csv'
EyFieldPAFile = 'Documents/SSProject/Debug/a4mmDATAscrappy/EyPaperSample.csv'
EzFieldPAFile = 'Documents/SSProject/Debug/a4mmDATAscrappy/EzPaperSample.csv'

materialField = 'Documents/SSProject/Debug/a4mmDATAscrappy/PaperMaterial.csv'

arrowDim = 32
#cutFactor = 13
start = 0
end = -2
hMax=0.6
colorMax=0.5

#maxSteps=100

deltaT = 8e-14
speedC = 299792458.0
deltaS = 60e-6/2

deltaT = deltaS/(speedC*1.8)

domainX = 0.006 # metres
domainY = domainX # metres
domainZ = domainX # metres
timePeriod =5e-11 # seconds

GRIDX = int(domainX/deltaS)-2
GRIDY = int(domainY/deltaS)-2
GRIDZ = int(domainZ/deltaS)-2
TSTEPS = int(timePeriod/deltaT) +1
print GRIDX
print GRIDY
print TSTEPS

cutFactor = GRIDX/arrowDim

def charCsvToNumpy(filename):
    field = numpy.zeros((GRIDX,GRIDY))
    with open(filename, 'rb') as material:
        Reader = csv.reader(material, delimiter=';', quotechar='|')
        rows=[]
        
        for row in Reader:
            thisRow = row[:-1]
            thisRowAsFloat = numpy.zeros(len(thisRow))
            if row!=[]:
                for xi in range(len(thisRow)):
                    #print ('char: %c,,,'%thisRow[xi])
                    if (thisRow[xi]=='v'):
                        thisRowAsFloat[xi]=0
                    else:
                        thisRowAsFloat[xi]=1
                        
                rows.append(thisRowAsFloat)
            else:
                field=numpy.array(rows).astype(numpy.float)
    return field

def csv2numpy(filename):
    field = numpy.zeros((TSTEPS,GRIDX,GRIDY))
    
    with open(filename, 'rb') as material:
        Reader = csv.reader(material, delimiter=';', quotechar='|')
        rows=[]
        tstep=0
        for row in Reader:
            #if tstep>maxSteps:
                #break
            thisRow = row[:-1]
            if row!=[]:
                rows.append(thisRow)
            else:
                field[tstep]=numpy.array(rows).astype(numpy.float)
                tstep+=1
                rows=[]
    return field
    
def cropField(field,cells):
    field=field[:,cells:-cells,cells:-cells]
    
    
# get fields
EyPAField = csv2numpy(EyFieldPAFile)
EzPAField = csv2numpy(EzFieldPAFile)

materialV = charCsvToNumpy(materialField)



imShowFieldPA = numpy.zeros((TSTEPS,GRIDX,GRIDY))
for i in range(TSTEPS):
    imShowFieldPA[i] = EyPAField[i] - materialV*hMax/2 # may need to explicitly for loop




unitScale=1000
domain = [0, domainX*unitScale, 0, domainY*unitScale]

fig = pyplot.figure(figsize =(12,7.2))

pyplot.xlabel('x/mm')
pyplot.ylabel('y/mm')
pyplot.axis(domain)



objects=[]
xvals=numpy.linspace(0,domainX*unitScale,GRIDX)

for tim in xrange(0,TSTEPS,1):
    timeStamp=(tim*deltaT * 1e12)-2.2
    
    im1 = pyplot.imshow(imShowFieldPA[tim,:,:].T, cmap='gist_heat', alpha=1,interpolation = 'nearest',\
    extent=domain,origin='lower',vmin=-hMax,vmax=hMax)

    
    timer = pyplot.text(0.18,3.58,'%.1fps'%timeStamp)
    objects.append([timer,im1])



ani = animation.ArtistAnimation(fig,objects, interval=10,blit=False)
pyplot.show()