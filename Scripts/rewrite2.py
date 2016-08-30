import sys

fw = open('Fieldrewrite.dat','w');

name = "PICfield"
exten = ".dat"
filenum = "0";
fullname = name + filenum + exten;
fr = open(fullname,'r');

print "Filename: ", sys.argv[0]
print "Number: ", sys.argv[1]

for line in fr:
        vals = line.split()
        print "Iterator: ", vals[0], vals[1]
