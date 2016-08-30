fw = open('velspace.dat','w');
name = "elecParticles"
exten = ".dat"

for i in range(0,100):
        filenum = str(2*i+1);
        fullname = name + filenum+ exten;
        print "Iterator: ", i
	fr = open(fullname,'r');
	x = fr.readline();
	fw.write(x);
