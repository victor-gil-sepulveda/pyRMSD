import matplotlib.pyplot
from matplotlib.ticker import FormatStrFormatter
import pylab
import numpy

files =  ["out/cuda5-35.out","out/openmp5-35.out","out/serial5-35.out"]
tags = {"out/cuda5-35.out":"Cuda","out/openmp5-35.out":"OpenMP","out/serial5-35.out":"Serial"}
numberOfStructures = range(5,40,5)

file_times = {}
# Parse files
for filename in files:
	handler = open(filename,"r")
	times = []	
	for line in handler:
		if line[:4]== "With":
			times.append(float(line.split()[4]));
	file_times[tags[filename]] = numpy.array(times)


# Calculate speedups
speedups = {}
for i in range(len(tags)):
	first = tags[tags.keys()[i]]
	for j in range(i+1,len(tags)):
		second = tags[tags.keys()[j]]
		speedups[(first,second)] = file_times[first] /file_times[second]

for key in speedups:
	print key
	print speedups[key]

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)

# Draw the plot
colors = {"Cuda": 'r',"OpenMP": 'b' , "Serial":'g'}
axis = matplotlib.pyplot.subplot(111)
#pylab.gcf().set_dpi(50)
pylab.title ("QCP Calculator Performance")
legend_lines = []
legend_keys = []
for key in tags:
	line, = axis.plot(numberOfStructures, file_times[tags[key]],colors[tags[key]])
	matplotlib.pyplot.setp(line, linewidth = 3)
	legend_lines.append(line)
	legend_keys.append(tags[key])
	line, = axis.plot(numberOfStructures, file_times[tags[key]],colors[tags[key]]+'o')
	matplotlib.pyplot.setp(line, linewidth = 3)
# Axis labels
matplotlib.pyplot.xlabel('Number Of Structures')
matplotlib.pyplot.ylabel('Calculation Time')
matplotlib.pyplot.gca().yaxis.set_major_formatter(FormatStrFormatter('%d s'))
matplotlib.pyplot.gca().xaxis.set_major_formatter(FormatStrFormatter('%d k'))
matplotlib.pyplot.gca().yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))

# Grid
matplotlib.pyplot.grid(True)

# Legend 
matplotlib.pyplot.legend(legend_lines, legend_keys, loc = 2)
# Show it
matplotlib.pyplot.show()
